"""Dot plot of one region against a reference record — the picture behind an orientation verdict.

Source: ``reference_audit.region_matches`` — every exact match the audit found between a record's
region and the record the set agrees on, in **both** orientations. The figure reads those rows and
draws them: it matches nothing itself, so the picture and the table cannot disagree (D5-03).

The two orientations are drawn in different colours on purpose. A region that runs the same way as
the subject fills one diagonal; a region whose SSC is inverted fills the other, and that is what
the eye is here to check. The audit reached the same verdict from the record's gene strands
(:func:`plastid.core.reference_audit.ssc_orientation`); when a person has to confirm a flip-flop
once and pin it (D5-08), this is the second thing they look at.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result

_MATCHES = "region_matches"
_COLUMNS = ("region", "query", "subject", "orientation", "query_pos", "subject_pos", "length")

# Forward matches are the ordinary case; inverted ones are the finding. The inverted colour is
# the theme's accent — the same red the toolkit uses for "look here" everywhere else.
_COLOURS = {"+": THEME.line, "-": THEME.accent}
_LABELS = {"+": "same orientation as the subject", "-": "opposite orientation"}


@dataclass(frozen=True, slots=True)
class Dot:
    """One exact-match run, exactly as the table carries it: 1-based positions within the region."""

    subject_pos: int
    query_pos: int
    length: int
    orientation: str


@dataclass(frozen=True, slots=True)
class Panel:
    """One query record against the subject — the dots, and the labels the caption needs."""

    query: str
    region: str
    subject: str
    dots: tuple[Dot, ...]

    @property
    def inverted(self) -> int:
        return sum(1 for dot in self.dots if dot.orientation == "-")

    @property
    def forward(self) -> int:
        return sum(1 for dot in self.dots if dot.orientation == "+")


@dataclass(frozen=True, slots=True)
class DotplotSeries:
    """Everything the figure draws, read from the matches table — the D5-03 audit trail."""

    panels: tuple[Panel, ...]


def series(results: Mapping[str, Result]) -> DotplotSeries:
    """Read ``reference_audit.region_matches`` into one panel per query record."""
    result = results.get("reference_audit")
    if result is None:
        raise ValueError(
            f"dotplot needs the 'reference_audit' module's result; have {sorted(results)}"
        )
    if _MATCHES not in result.tables:
        raise ValueError(f"dotplot needs table reference_audit.{_MATCHES}")
    table = result.tables[_MATCHES]
    absent = [name for name in _COLUMNS if name not in {column.name for column in table.columns}]
    if absent:
        raise ValueError(f"dotplot: {table.id} has no column(s) {absent}")
    panels: list[Panel] = []
    for query in dict.fromkeys(row["query"] for row in table.rows):
        rows = [row for row in table.rows if row["query"] == query]
        panels.append(
            Panel(
                query=query,
                region=rows[0]["region"],
                subject=rows[0]["subject"],
                dots=tuple(
                    Dot(
                        subject_pos=row["subject_pos"],
                        query_pos=row["query_pos"],
                        length=row["length"],
                        orientation=row["orientation"],
                    )
                    for row in rows
                ),
            )
        )
    return DotplotSeries(panels=tuple(panels))


def caption(data: DotplotSeries) -> str:
    """The figure's own title, built from what the rows say — never from a recomputation."""
    first = data.panels[0]
    if len(data.panels) == 1:
        return f"{first.region} dot plot — {first.query} against {first.subject}"
    return f"{first.region} dot plot — {len(data.panels)} records against {first.subject}"


def _panel(axes: Any, panel: Panel) -> None:
    """One scatter, both orientations, with the run lengths carried as dot areas."""
    for orientation, colour in _COLOURS.items():
        dots = [dot for dot in panel.dots if dot.orientation == orientation]
        if not dots:
            continue
        axes.scatter(
            [dot.subject_pos for dot in dots],
            [dot.query_pos for dot in dots],
            s=[max(1.0, min(12.0, dot.length / 20)) for dot in dots],
            color=colour,
            linewidths=0,
            label=f"{_LABELS[orientation]} ({len(dots)})",
        )
    axes.set_xlabel(f"{panel.subject} position (bp, 1-based)")
    axes.set_ylabel(f"{panel.query} position (bp, 1-based)")
    axes.set_title(
        panel.query if len(panel.dots) else panel.query,
        loc="left",
        fontsize=THEME.font_size - 1,
        color=THEME.muted,
        pad=4,
    )
    axes.grid(color=THEME.grid)
    axes.set_axisbelow(True)
    axes.legend(loc="upper right", frameon=False, fontsize=THEME.font_size - 2, markerscale=3)


def _build(data: DotplotSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        panels = data.panels
        if len(panels) == 1:
            axes_list = [figure.subplots()]
        else:
            axes_list = figure.subplots(1, len(panels), squeeze=False)[0]
        figure.suptitle(title, x=0.02, ha="left", fontsize=THEME.font_size, color=THEME.line)
        for axes, panel in zip(axes_list, panels, strict=True):
            _panel(axes, panel)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the table, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="dotplot",
        title=text,
        source_tables=(f"reference_audit.{_MATCHES}",),
        build=_build(data, text),
        figsize=(4.2 * len(data.panels), 3.8),
    )


SPEC = FigureSpec(
    id="dotplot",
    title="Region dot plot against a reference record",
    requires=(("reference_audit", _MATCHES),),
    build=build,
)
