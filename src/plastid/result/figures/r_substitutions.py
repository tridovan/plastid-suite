"""Paper Figure 7 — substitution types and Ts/Tv, drawn by ``R_F7_pairwise_substitution.R``.

The script reads ``Substitution_Type`` rows against one column per genome: twelve ``A_to_G``-style
rows with ``Category`` ``Substitution``, then one ``Ts/Tv Ratio`` row with ``Category`` ``Ratio``.
It sums the complementary pairs itself, so the twelve **directed** changes are what this pivot
supplies.

**This figure had no source until now.** Measured 17/10/2026: ``core/align.py`` counted only
transitions and transversions, so the toolkit did not hold the six types the manuscript draws.
``substitutions.substitutions_by_type`` is that table, added in the same change as this figure.

The values are each query genome **against the reference** (row 0 of the alignment) — the reference
is not a column, so the caption names it. Ts/Tv comes from whichever per-genome table the run's
scope produced: ``substitutions`` in the whole-alignment scope, ``substitutions_by_genome`` in the
CDS scope. A run holding neither is skipped by name rather than drawn without its ratio panel.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.core.align import SUBSTITUTION_TYPES
from plastid.result.figures.r_common import csv_text
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSkip, RFigureSpec

_SUMMARY = "substitutions_summary"
_BY_TYPE = "substitutions_by_type"
_RATIO_TABLES = ("substitutions", "substitutions_by_genome")
_FILENAME = "substitutions.csv"


@dataclass(frozen=True, slots=True)
class SubstitutionsFrame:
    """The wide table F7 reads: twelve directed types, one ratio row, one column per genome."""

    genomes: tuple[str, ...]
    counts: dict[tuple[str, str], str]  # (genome, type) → the count column's rendering
    ratios: dict[str, str]  # genome → the ts_tv column's rendering
    ratio_table: str
    reference: str

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = ("Substitution_Type", "Category", *self.genomes)
        rows = [
            (
                type_name,
                "Substitution",
                *[self.counts[(genome, type_name)] for genome in self.genomes],
            )
            for type_name in SUBSTITUTION_TYPES
        ]
        rows.append(("Ts/Tv Ratio", "Ratio", *[self.ratios[genome] for genome in self.genomes]))
        return ((_FILENAME, csv_text(header, rows)),)

    def caption(self) -> str:
        return (
            f"each of {len(self.genomes)} genome(s) against the reference {self.reference}; "
            "the twelve directed changes are the rows and the script sums the six complementary "
            f"pairs itself. Ts/Tv is read from substitutions.{self.ratio_table}. "
            "A blank Ts/Tv means no transversion was found — it is undefined, not infinite."
        )


def frame(results: Mapping[str, Result]) -> SubstitutionsFrame:
    """The type table pivoted to the script's shape, with Ts/Tv from the run's own scope."""
    result = results["substitutions"]
    by_type = result.tables[_BY_TYPE]
    summary = result.tables[_SUMMARY]
    ratio_table = next((name for name in _RATIO_TABLES if name in result.tables), None)
    if ratio_table is None:
        raise RFigureSkip(
            f"this run's substitutions scope has neither {' nor '.join(_RATIO_TABLES)}, so the "
            f"Ts/Tv panel has no cells to read — no substitution figure for this run"
        )
    ratio = {
        row["genome"]: result.tables[ratio_table].format(f"{row['genome']}.ts_tv")
        for row in result.tables[ratio_table].rows
    }
    genomes = tuple(dict.fromkeys(row["genome"] for row in by_type.rows))
    counts = {
        (row["genome"], row["substitution_type"]): by_type.format(f"{row['row_id']}.count")
        for row in by_type.rows
    }
    missing = [genome for genome in genomes if genome not in ratio]
    if missing:
        raise RFigureSkip(
            f"substitutions.{ratio_table} has no row for {', '.join(missing)} — the type table and "
            f"the ratio table describe different genomes, so this figure would be half-empty"
        )
    return SubstitutionsFrame(
        genomes=genomes,
        counts=counts,
        ratios=ratio,
        ratio_table=ratio_table,
        reference=summary.value("accession") or summary.value("organism"),
    )


SPEC = RFigureSpec(
    id="r_substitutions",
    title="Substitution types and Ts/Tv (paper figure 7, drawn by R)",
    requires=(("substitutions", _SUMMARY), ("substitutions", _BY_TYPE)),
    script="substitutions",
    input_mode="file",
    frame=frame,
)
