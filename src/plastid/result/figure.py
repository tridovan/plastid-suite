"""Figures: a pure function of a :class:`Result`, drawn through one theme (Phase 5, D5-01/D5-03).

Three rules shape this module, and each is tested:

* **A figure reads cells, it never recomputes.** Every series is built by a ``series()`` function
  that returns plain data taken from table cells, so a figure and its table cannot disagree —
  a measured failure of the early toolkit in a new costume. ``Figure.source_tables`` records
  what it read.
* **One theme.** Fonts, colours and sizes live in :data:`THEME`, not in each figure, so eleven
  figure types read as one system.
* **Byte-identical across runs.** matplotlib stamps a date and a per-process random salt into its
  vector output; both are pinned here, measured rather than assumed. SVG needs
  the salt *and* the date pinned, PDF has a metadata hook, and the PostScript backend has none —
  its ``%%CreationDate`` is rewritten to a fixed stamp.

matplotlib is an optional dependency (``pip install 'plastid-suite[figures]'``): it is imported
inside the function, so the minimum install stays pure Python and a missing extra is a named
error, never a traceback.
"""

from __future__ import annotations

import io
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

from plastid.core.deps import require_module
from plastid.result.result import Result

__all__ = [
    "FIGURE_FORMATS",
    "Figure",
    "FigureSpec",
    "THEME",
    "Theme",
    "export_figure",
    "missing_tables",
    "themed",
]

# Vector only: SVG for editing, PDF for submission, EPS when a journal still demands it.
FIGURE_FORMATS = ("svg", "pdf", "eps")

_EPS_DATE = re.compile(rb"%%CreationDate: [^\n]*")
_EPS_DATE_FIXED = b"%%CreationDate: Thu Jan  1 00:00:00 1970"

# Everything the theme touches in matplotlib's global state, restored on the way out.
_THEMED_RCPARAMS = (
    "svg.hashsalt",
    "svg.fonttype",
    "pdf.fonttype",
    "font.family",
    "font.size",
    "axes.titlesize",
    "axes.labelsize",
    "xtick.labelsize",
    "ytick.labelsize",
    "axes.edgecolor",
    "axes.labelcolor",
    "text.color",
    "xtick.color",
    "ytick.color",
    "axes.spines.top",
    "axes.spines.right",
    "grid.color",
    "grid.linewidth",
    "lines.linewidth",
    "figure.facecolor",
    "savefig.bbox",
)


@dataclass(frozen=True, slots=True)
class Theme:
    """One place for every font, colour and size a figure uses (D5-01)."""

    font_family: str = "DejaVu Sans"  # bundled with matplotlib — the same face on every machine
    font_size: float = 9.0
    line_width: float = 1.2
    line: str = "#1f3b57"
    accent: str = "#b4472e"
    muted: str = "#8a8f98"
    grid: str = "#d9dde2"
    background: str = "#ffffff"
    # Subtle fills behind the region bands. IRa and IRb share one colour: they are the same
    # sequence, and two shades would read as two different regions.
    region_fills: tuple[tuple[str, str], ...] = (
        ("LSC", "#f4f6f8"),
        ("IRa", "#e3eaf1"),
        ("SSC", "#f8f3ea"),
        ("IRb", "#e3eaf1"),
    )

    def region_fill(self, region: str) -> str:
        for name, colour in self.region_fills:
            if name == region:
                return colour
        return self.background


THEME = Theme()


@dataclass(frozen=True, slots=True)
class Figure:
    """A rendered figure: what it read, and its bytes in each vector format."""

    id: str
    title: str
    source_tables: tuple[str, ...]  # "<module>.<table>" — the D5-03 audit trail
    data: tuple[tuple[str, bytes], ...]

    def __post_init__(self) -> None:
        if not self.id or "." in self.id:
            raise ValueError(f"figure id must be a plain name without dots: {self.id!r}")
        if not self.data:
            raise ValueError(f"figure {self.id!r} carries no payload")
        formats = [fmt for fmt, _ in self.data]
        if len(set(formats)) != len(formats):
            raise ValueError(f"figure {self.id!r} repeats a format: {formats}")
        for fmt, blob in self.data:
            if fmt not in FIGURE_FORMATS:
                raise ValueError(f"figure {self.id!r}: {fmt!r} is not a vector format")
            if not blob:
                raise ValueError(f"figure {self.id!r}: the {fmt} payload is empty")

    @property
    def formats(self) -> tuple[str, ...]:
        return tuple(fmt for fmt, _ in self.data)

    def payload(self, fmt: str) -> bytes:
        for name, blob in self.data:
            if name == fmt:
                return blob
        raise KeyError(f"figure {self.id!r} has no {fmt!r} payload (have: {self.formats})")


def missing_tables(
    requires: Sequence[tuple[str, str]], results: Mapping[str, Result]
) -> tuple[str, ...]:
    """The ``(module, table)`` pairs this run does not hold — named as ``module.table``.

    Shared by :class:`FigureSpec` and the R figures (``result/rfigure.py``): the two catalogues
    must name a missing source identically, or the same run would explain the same absence twice.
    """
    return tuple(
        f"{module}.{table}"
        for module, table in requires
        if module not in results or table not in results[module].tables
    )


@dataclass(frozen=True, slots=True)
class FigureSpec:
    """A figure type: which tables it needs, and how to build it from a run's results.

    ``build`` receives the matplotlib **figure**, not an axes: a figure may be one panel, several,
    or a polar projection, and only the figure can make all three. ``requires`` is keyed by
    **(module, table)** and never by table id alone — ``regions`` exists
    in two modules with different columns, and a figure that resolved the wrong one would draw a
    picture of a different table without saying so.
    """

    id: str
    title: str
    requires: tuple[tuple[str, str], ...]
    build: Callable[[Mapping[str, Result]], Figure]
    formats: tuple[str, ...] = FIGURE_FORMATS  # a page is PDF only; the others carry all three

    def missing(self, results: Mapping[str, Result]) -> tuple[str, ...]:
        """The required tables this run does not hold — named as ``module.table``."""
        return missing_tables(self.requires, results)


@contextmanager
def themed(theme: Theme = THEME) -> Iterator[Any]:
    """Apply the theme to matplotlib's global state and put it back — no hidden globals.

    ``svg.hashsalt`` is what makes the SVG's element ids stable: without it matplotlib salts them
    with a per-process random value and two runs of the same input differ (measured).
    """
    pyplot = require_module("matplotlib.pyplot", extra="figures", purpose="Figure export")
    matplotlib = require_module("matplotlib", extra="figures", purpose="Figure export")
    previous = {key: matplotlib.rcParams[key] for key in _THEMED_RCPARAMS}
    matplotlib.rcParams.update(
        {
            "svg.hashsalt": "plastid",
            # SVG is the *editable* copy, so its labels stay text rather than outlines; PDF is the
            # submission copy, and journals want TrueType (42) rather than matplotlib's Type 3.
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.family": theme.font_family,
            "font.size": theme.font_size,
            "axes.titlesize": theme.font_size + 1,
            "axes.labelsize": theme.font_size,
            "xtick.labelsize": theme.font_size - 1,
            "ytick.labelsize": theme.font_size - 1,
            "axes.edgecolor": theme.muted,
            "axes.labelcolor": theme.line,
            "text.color": theme.line,
            "xtick.color": theme.muted,
            "ytick.color": theme.muted,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": theme.grid,
            "grid.linewidth": 0.6,
            "lines.linewidth": theme.line_width,
            "figure.facecolor": theme.background,
            "savefig.bbox": "tight",
        }
    )
    try:
        yield pyplot
    finally:
        matplotlib.rcParams.update(previous)


def _save(fig: Any, fmt: str) -> bytes:
    """One format, with everything matplotlib would otherwise vary pinned."""
    buffer = io.BytesIO()
    if fmt == "svg":
        fig.savefig(buffer, format="svg", metadata={"Date": None})
        return buffer.getvalue()
    if fmt == "pdf":
        fig.savefig(buffer, format="pdf", metadata={"CreationDate": None})
        return buffer.getvalue()
    if fmt == "eps":
        fig.savefig(buffer, format="eps")
        # The PostScript backend takes no metadata, so the only way to a fixed file is to write
        # the date out of it — matplotlib's own test suite strips the same line.
        return _EPS_DATE.sub(_EPS_DATE_FIXED, buffer.getvalue(), count=1)
    raise ValueError(
        f"unsupported figure format {fmt!r} — figures are vector only: {FIGURE_FORMATS}"
    )


def export_figure(
    *,
    figure_id: str,
    title: str,
    source_tables: Sequence[str],
    build: Callable[[Any], None],
    figsize: tuple[float, float],
    formats: Sequence[str] = FIGURE_FORMATS,
) -> Figure:
    """Draw once, export each format, and return the payloads as data.

    ``build`` is handed the figure and makes its own axes — ``figure.subplots()`` for one panel,
    ``figure.subplots(1, 2)`` for two, ``figure.add_subplot(projection="polar")`` for a circle.
    """
    with themed() as pyplot:
        fig = pyplot.figure(figsize=figsize)
        try:
            build(fig)
            data = tuple((fmt, _save(fig, fmt)) for fmt in formats)
        finally:
            pyplot.close(fig)
    return Figure(id=figure_id, title=title, source_tables=tuple(source_tables), data=data)
