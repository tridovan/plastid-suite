"""The R figures: a second catalogue, drawn by the manuscript's own scripts (M10.3/M10.4).

**Why this is not :class:`~plastid.result.figure.FigureSpec`.** Three reasons, each measured:

* ``Figure.__post_init__`` accepts only the vector formats in :data:`FIGURE_FORMATS`, and a test
  asserts it — the R engine writes PNG and PDF, and a raster format would be a different kind of
  thing rather than a narrower one.
* ``FigureSpec.build`` receives the results mapping alone, so a figure drawn that way would call R
  with **no ``RunContext``** — and provenance (name · version · argv · exit code · logs) is not
  optional in this toolkit.
* One script writes several named images (F5 writes three PNG and three PDF), which ``Figure``'s
  one-payload-per-format shape cannot express.

The bundle says the same thing in its own words: ``03_tich_hop_PTS/figure_spec_mau.py`` — the
scripts cannot meet PTS's byte-identical guarantee, so the honest home is a supplementary export
step beside the figure step, not a slot in ``_SPECS``.

**What an R figure still is.** A figure **reads cells and never recomputes**: every number it hands
R comes from a table cell, through that column's own ``format`` rendering, so nothing is rounded
twice and a figure cannot disagree with the table it read. The pivot lives in
``result/figures/r_*.py`` as a pure ``frame(results) -> RFrame``; the R call is injected by
``cli.py``, because this layer may not import ``external/``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Protocol, runtime_checkable

from plastid.result.figure import missing_tables
from plastid.result.result import Result

__all__ = ["RFrame", "RFigureSkip", "RFigureSpec", "scope_lines"]


class RFigureSkip(Exception):
    """A **declared** skip: the run holds the tables, but not the shape this figure needs.

    Not every absence is a missing table — a run may hold ``substitutions`` in the CDS scope, which
    has no ``substitutions`` table at all. That is a named skip with a reason a user can act on, not
    a crash: a bug is not a skip, and anything other than this exception propagates.
    """


@runtime_checkable
class RFrame(Protocol):
    """What a pivot returns: the files R will read, and the caveats its caption must carry."""

    def inputs(self) -> tuple[tuple[str, str], ...]:
        """``(filename, text)`` for every file the script reads — one, or two for F6."""
        ...

    def caption(self) -> str:
        """The figure's own lines, before the source tables' scope labels are appended."""
        ...


@dataclass(frozen=True, slots=True)
class RFigureSpec:
    """One R-drawn figure: which tables it needs, which script draws it, and how it is pivoted.

    ``script`` is a key into :data:`plastid.external.rscript.R_SCRIPTS`. That binding is checked by
    a test rather than here, because this layer may not import ``external/``.
    """

    id: str
    title: str
    requires: tuple[tuple[str, str], ...]
    script: str
    input_mode: str  # "file" — one --input; "dir" — an --inputdir holding several files
    frame: Callable[[Mapping[str, Result]], RFrame]
    notes: tuple[str, ...] = ()
    # A parameter naming a file the script reads **directly**, rather than one this toolkit pivots
    # out of a table. `R_F9_phylogeny_tree.R` is the one case: the tree is the user's own, built
    # outside this toolkit, so there is no cell to read and nothing to write — only a path to be
    # named, and a digest to record. Its caption says which file drew the figure.
    source_param: str | None = None

    def __post_init__(self) -> None:
        if not self.id or "." in self.id:
            raise ValueError(f"R figure id must be a plain name without dots: {self.id!r}")
        if self.input_mode not in ("file", "dir"):
            raise ValueError(
                f"R figure {self.id!r}: input_mode must be 'file' or 'dir', got {self.input_mode!r}"
            )

    def missing(self, results: Mapping[str, Result]) -> tuple[str, ...]:
        """The required tables this run does not hold — named as ``module.table``."""
        return missing_tables(self.requires, results)

    def scope_lines(self, results: Mapping[str, Result]) -> tuple[str, ...]:
        """``<module>.<table> — scope: <label>`` for every table this figure read.

        The label is the one the **run** stored, so a data label applied by the runner reaches the
        caption too: a figure drawn from a stand-in says so on its own face.
        """
        return scope_lines(self.requires, results)


def scope_lines(
    requires: Sequence[tuple[str, str]], results: Mapping[str, Result]
) -> tuple[str, ...]:
    """The scope label of every ``(module, table)`` a figure read, in the order it declared them."""
    return tuple(
        f"{module}.{table} — scope: {results[module].tables[table].scope.label}"
        for module, table in requires
    )


def as_text(value: Any) -> str:
    """One cell as text: the empty string for a blank, never the word ``None``."""
    return "" if value is None else str(value)
