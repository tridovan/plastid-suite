"""Paper Figure 6 — nucleotide diversity π, drawn by ``R_F6_nucleotide_diversity_pi.R``.

The script reads **two** files, each ``Region<TAB>Value``: panel A is the coding/RNA features,
panel B the non-coding ones (introns and intergenic spacers). It is therefore the one R figure
that must be given ``--inputdir`` rather than ``--input``: ``duong_dan_vao()`` returns the
``--input`` file for *every* read, so passing one file would make both panels read panel A.

Unlike the other five, this figure comes from an **arity-2** module: ``diversity`` compares a
reference against every query genome in the run, so π is measured over all of them and the figure
genuinely covers the run's whole panel. The features are the reference's own annotation — one row
per reference gene, intron or spacer — which is what the manuscript's panel B draws.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from plastid.result.figures.r_common import tsv_text
from plastid.result.result import Result
from plastid.result.rfigure import RFigureSpec

_SUMMARY = "diversity_summary"
_GENES = "diversity_genes"
_INTRONS = "diversity_introns"
_IGS = "diversity_igs"
_GENES_FILE = "genes_for_R_plot.txt"
_NONCODING_FILE = "noncoding_for_R_plot.txt"


@dataclass(frozen=True, slots=True)
class DiversityFrame:
    """The two files F6 reads: coding features, then non-coding ones."""

    genes: tuple[tuple[str, str], ...]  # (label, π rendering)
    noncoding: tuple[tuple[str, str], ...]
    n_genomes: str
    window_bp: str

    def inputs(self) -> tuple[tuple[str, str], ...]:
        header = ("Region", "Value")
        return (
            (_GENES_FILE, tsv_text(header, self.genes)),
            (_NONCODING_FILE, tsv_text(header, self.noncoding)),
        )

    def caption(self) -> str:
        return (
            f"{len(self.genes)} coding and {len(self.noncoding)} non-coding features of the "
            f"reference annotation, π over {self.n_genomes} aligned genomes "
            f"({self.window_bp} bp windows for the summary's own number). "
            "Panel A is the reference's genes, panel B its introns then its intergenic spacers — "
            "the manuscript's own split. Panel B's labels are the tables' own keys: an intron is "
            "`gene.copy.n`, so the two inverted-repeat copies of one intron are two rows and a "
            "gene with two introns has two, and only the key tells those two cases apart."
        )


def frame(results: Mapping[str, Result]) -> DiversityFrame:
    """The three diversity tables as the script's two ``Region/Value`` files."""
    result = results["diversity"]
    summary = result.tables[_SUMMARY]
    genes = result.tables[_GENES]
    introns = result.tables[_INTRONS]
    spacers = result.tables[_IGS]
    # The labels are the tables' own keys, never a name built from parts. Measured 17/10/2026:
    # a gene-name label collides, and the script's `factor(Region, levels = Region)` then halts
    # with "factor level [6] is duplicated" — the figure does not draw at all. Two different
    # collisions produce it, and neither is visible from the gene name alone: `rps12.1.2` and
    # `rps12.2.2` are the two **IR copies of one** intron, while `clpP.1.1` and `clpP.1.2` are two
    # **different** introns of one gene. Only `intron_id` (``gene.copy.n``) tells them apart.
    return DiversityFrame(
        genes=tuple((row["gene"], genes.format(f"{row['gene']}.pi")) for row in genes.rows),
        noncoding=tuple(
            [(row["intron_id"], introns.format(f"{row['intron_id']}.pi")) for row in introns.rows]
            + [(row["igs_id"], spacers.format(f"{row['igs_id']}.pi")) for row in spacers.rows]
        ),
        n_genomes=summary.format("n_genomes"),
        window_bp=summary.format("window_bp"),
    )


SPEC = RFigureSpec(
    id="r_diversity",
    title="Nucleotide diversity π (paper figure 6, drawn by R)",
    requires=(
        ("diversity", _GENES),
        ("diversity", _INTRONS),
        ("diversity", _IGS),
        ("diversity", _SUMMARY),
    ),
    script="diversity",
    input_mode="dir",
    frame=frame,
)
