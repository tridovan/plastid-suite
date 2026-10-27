"""Paper Figure 4 — amino-acid composition, drawn by ``R_F4_amino_acid_heatmap.R``.

The script reads a wide table: ``AA_Code`` (one-letter), ``Amino_Acid``, then one column per
species of percentages. PTS's ``amino_acid.amino_acid_composition`` is already one row per residue
(``residue, count, percent, hydropathy``) but holds **one genome**, so the pivot is a column rename
plus the single species column. The script draws **every** species column in its input (M11.7), so
a wider table is a valid figure; this pivot supplies the one genome a run holds. The species column
keeps its case, and ``Amino_Acid`` is written but the script overwrites it from ``AA_Code``; it is
filled in anyway so the file says what it holds.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import AMINO_ACID_NAMES, csv_text
from plastid.result.figures.r_common import species_column as _species_column
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

_COMPOSITION = "amino_acid_composition"
_SUMMARY = "amino_acid"
_COMPOSITION_COLUMNS = ("residue", "count", "percent", "hydropathy")
_SUMMARY_COLUMNS = ("accession", "organism", "n_cds", "n_residues", "gravy")
_SCRIPT = "R_F4_amino_acid_heatmap.R"
_FILENAME = "amino_acid_composition.csv"


@dataclass(frozen=True, slots=True)
class AminoAcidFrame:
    """The wide table F4 reads — one residue per row, one genome's column."""

    species: str
    residues: tuple[tuple[str, str, str], ...]  # (one-letter, full name, percent rendering)
    n_residues: str
    gravy: str

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = ("AA_Code", "Amino_Acid", self.species)
        return ((_FILENAME, csv_text(header, self.residues)),)

    def caption(self) -> str:
        return (
            f"{len(self.residues)} residues of one genome ({self.species}); "
            f"{self.n_residues} residues counted, GRAVY {self.gravy}. "
            "The script draws every species column in its input, in that order or in "
            "--species-order; this run supplies one."
        )


def frame(results: Mapping[str, Result]) -> AminoAcidFrame:
    """The composition table, one column wide — percents are the column's own rendering."""
    result = results["amino_acid"]
    table = result.tables[_COMPOSITION]
    summary = result.tables[_SUMMARY]
    for name in _COMPOSITION_COLUMNS:
        if name not in {column.name for column in table.columns}:
            raise ValueError(
                f"amino_acid.{_COMPOSITION} has no {name!r} column — cannot draw the composition"
            )
    for name in _SUMMARY_COLUMNS:
        if name not in {column.name for column in summary.columns}:
            raise ValueError(f"amino_acid.{_SUMMARY} has no {name!r} column — cannot draw it")

    organism = summary.value("organism")
    species = _species_column(organism)
    return AminoAcidFrame(
        species=species,
        residues=tuple(
            (
                row["residue"],
                AMINO_ACID_NAMES.get(row["residue"], row["residue"]),
                table.format(f"{row['residue']}.percent"),
            )
            for row in table.rows
        ),
        n_residues=summary.format("n_residues"),
        gravy=summary.format("gravy"),
    )


SPEC = RFigureSpec(
    id="r_amino_acid",
    title="Amino-acid composition (paper figure 4, drawn by R)",
    requires=(("amino_acid", _COMPOSITION), ("amino_acid", _SUMMARY)),
    script="amino_acid",
    input_mode="file",
    frame=frame,
)
