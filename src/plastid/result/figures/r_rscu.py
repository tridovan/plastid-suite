"""Paper Figure 3 — RSCU by codon, drawn by ``R_F3_RSCU_codon_usage.R``.

The script reads a **wide** table: one row per sense codon, one column per species. PTS's
``codon_usage.rscu`` is long (``codon, amino_acid, count, rscu``) and covers **one genome**, so the
pivot here is the long→wide step and it produces exactly one species column.

Two things about that shape are load-bearing and neither is obvious from the script:

* the species column is ``Genus_species`` **with the case kept** — the script prints the header it
  was given;
* ``Amino_Acid`` is the **full** name (``Alanine``), which the script maps to a three-letter code
  for the facet labels.

The script draws **every** species column in its input (M11.7), so a wider table — one column per
sample — is a valid figure with no list to intersect; `plastid figures` supplies the one genome a
run holds until the batch path does the joining.

Every RSCU value is the ``rscu`` column's own rendering, so the figure cannot round a cell twice —
and a blank RSCU (undefined for a family with no observations) stays blank rather than becoming a
zero-height bar that reads as "never used".
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import AMINO_ACID_NAMES, csv_text
from plastid.result.figures.r_common import species_column as _species_column
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

_RSCU = "rscu"
_SUMMARY = "codon_usage"
_RSCU_COLUMNS = ("codon", "amino_acid", "count", "rscu")
_SUMMARY_COLUMNS = ("accession", "organism", "n_cds", "enc")
_SCRIPT = "R_F3_RSCU_codon_usage.R"
_FILENAME = "codon_usage_rscu.csv"


@dataclass(frozen=True, slots=True)
class RscuFrame:
    """The wide table F3 reads, plus what its caption has to say about the run it came from."""

    species: str  # the organism as the script's column spelling
    codons: tuple[tuple[str, str, str], ...]  # (full amino-acid name, codon, RSCU rendering)
    n_cds: str
    enc: str

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = ("Amino_Acid", "Codon", self.species)
        return ((_FILENAME, csv_text(header, self.codons)),)

    def caption(self) -> str:
        return (
            f"{len(self.codons)} sense codons of one genome ({self.species}); "
            f"{self.n_cds} CDS, ENC {self.enc}. "
            "The script draws every species column in its input, in that order or in "
            "--species-order; this run supplies one."
        )


def frame(results: Mapping[str, Result]) -> RscuFrame:
    """The ``rscu`` table, one column wide — every cell read, none recomputed."""
    result = results["codon_usage"]
    table = result.tables[_RSCU]
    summary = result.tables[_SUMMARY]
    for name in _RSCU_COLUMNS:
        if name not in {column.name for column in table.columns}:
            raise ValueError(f"codon_usage.{_RSCU} has no {name!r} column — cannot draw RSCU")
    for name in _SUMMARY_COLUMNS:
        if name not in {column.name for column in summary.columns}:
            raise ValueError(f"codon_usage.{_SUMMARY} has no {name!r} column — cannot draw RSCU")

    organism = summary.value("organism")
    species = _species_column(organism)
    return RscuFrame(
        species=species,
        codons=tuple(
            (
                AMINO_ACID_NAMES.get(row["amino_acid"], row["amino_acid"]),
                row["codon"],
                table.format(f"{row['codon']}.rscu"),
            )
            for row in table.rows
        ),
        n_cds=summary.format("n_cds"),
        enc=summary.format("enc"),
    )


SPEC = RFigureSpec(
    id="r_rscu",
    title="RSCU by codon (paper figure 3, drawn by R)",
    requires=(("codon_usage", _RSCU), ("codon_usage", _SUMMARY)),
    script="rscu",
    input_mode="file",
    frame=frame,
)
