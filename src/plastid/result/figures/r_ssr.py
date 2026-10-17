"""Paper Figure 5 — SSR distribution, three panels, drawn by ``R_F5_SSR_distribution.R``.

The script reads one row per species with sixteen columns, and pivots them itself into its three
panels: (A) LSC/SSC/IR, (B) the six unit sizes, (C) gene/tRNA/rRNA/intron/IGS. Four PTS tables
feed that row — ``ssr_summary``, ``ssr_by_region``, ``ssr_by_unit`` and the per-SSR ``ssrs`` table
— and all four describe **one genome** (``ssr`` is an arity-1 module), so the figure is one row.

**``tRNA`` and ``rRNA`` are written blank, not zero.** PTS's ``ssrs.context`` is only
``Gene``/``Intron``/``IGS`` (``modules/ssr.py``: a gene occupies its exons, and the rule does not
split gene kinds), so the toolkit does not hold the number the manuscript's two extra columns ask
for. A zero would be a claim that no SSR sits in a tRNA or rRNA — which the toolkit has not
measured. The caption says so.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import csv_text
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

_SUMMARY = "ssr_summary"
_BY_REGION = "ssr_by_region"
_BY_UNIT = "ssr_by_unit"
_SSRS = "ssrs"
_FILENAME = "ssr_summary_for_r.csv"

# The script's own column names, in its own order — it selects them by name.
_UNIT_COLUMNS = (
    ("mono", "Mono"),
    ("di", "Di"),
    ("tri", "Tri"),
    ("tetra", "Tetra"),
    ("penta", "Penta"),
    ("hexa", "Hexa"),
)
_CONTEXT_COLUMNS = ("Gene", "Intron", "IGS")
_UNSOURCED = ("tRNA", "rRNA")


@dataclass(frozen=True, slots=True)
class SsrFrame:
    """One species' row, exactly the sixteen columns the script selects."""

    values: tuple[str, ...]
    unsourced: tuple[str, ...]

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = (
            "Species",
            "Total_SSRs",
            "LSC",
            "SSC",
            "IR",
            "Mono",
            "Di",
            "Tri",
            "Tetra",
            "Penta",
            "Hexa",
            "Gene",
            "tRNA",
            "rRNA",
            "Intron",
            "IGS",
        )
        return ((_FILENAME, csv_text(header, (self.values,))),)

    def caption(self) -> str:
        return (
            "One genome per run (the ssr module takes one): the paper's 34-species panel needs 34 "
            f"runs and this toolkit combines none of them. {', '.join(self.unsourced)} are blank "
            "because PTS records an SSR's context as Gene/Intron/IGS only and does not split gene "
            "kinds — blank is 'not measured', never zero."
        )


def _cell(table, ref: str) -> str:
    return table.format(ref)


def frame(results: Mapping[str, Result]) -> SsrFrame:
    """The four SSR tables folded into the script's one wide row — counts, never recounts."""
    result = results["ssr"]
    summary = result.tables[_SUMMARY]
    by_region = {row["region"]: row for row in result.tables[_BY_REGION].rows}
    by_unit = {row["unit"]: row for row in result.tables[_BY_UNIT].rows}
    contexts: dict[str, int] = dict.fromkeys(_CONTEXT_COLUMNS, 0)
    for row in result.tables[_SSRS].rows:
        if row["context"] in contexts:
            contexts[row["context"]] += 1

    # The one aggregation in this figure: the paper draws a single IR bar, PTS reports its two
    # copies separately. Asserted by a test rather than left implicit.
    ir = sum(by_region[name]["total"] for name in ("IRa", "IRb") if name in by_region)
    values = (
        summary.value("accession"),
        _cell(summary, "total"),
        _cell(result.tables[_BY_REGION], "LSC.total") if "LSC" in by_region else "",
        _cell(result.tables[_BY_REGION], "SSC.total") if "SSC" in by_region else "",
        str(ir),
        *[str(by_unit[name]["n"]) if name in by_unit else "0" for name, _label in _UNIT_COLUMNS],
        str(contexts["Gene"]),
        "",  # tRNA — not held by this toolkit
        "",  # rRNA — not held by this toolkit
        str(contexts["Intron"]),
        str(contexts["IGS"]),
    )
    return SsrFrame(values=values, unsourced=_UNSOURCED)


SPEC = RFigureSpec(
    id="r_ssr",
    title="SSR distribution (paper figure 5, drawn by R)",
    requires=(
        ("ssr", _SUMMARY),
        ("ssr", _BY_REGION),
        ("ssr", _BY_UNIT),
        ("ssr", _SSRS),
    ),
    script="ssr",
    input_mode="file",
    frame=frame,
)
