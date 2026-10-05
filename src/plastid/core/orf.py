"""ORF checks for protein-coding features (M3.3, P3-01).

The check answers one question per CDS — *does this annotation translate the way it says it
does* — and answers it from the record's own declarations (``/transl_table``, ``/codon_start``),
never from a default that happens to fit most species.

Three rules keep it honest:

* **never rewrite a CDS to make it pass.** A failure is reported; repairing it would erase the
  evidence that something upstream is wrong;
* **refuse rather than guess.** A translation table whose amino acids we do not know, or a
  ``/codon_start`` that is not 1, raises :class:`OrfError` — translating the wrong frame is a
  silent-wrong, and this project has met enough of those;
* **a flagged start is not a failure.** Measured on the golden genome: 5 of 86 CDS do not start
  with ATG, and every one of them is a curated, published annotation.
"""

from __future__ import annotations

from dataclasses import dataclass

from plastid.core.model import Gene, Genome, coding_sequence

__all__ = [
    "CANONICAL_START",
    "FAILURE_STATUSES",
    "NON_CANONICAL_STARTS",
    "ORF_STATUSES",
    "START_CODONS",
    "STOP_CODONS",
    "OrfCheck",
    "OrfError",
    "annotated_not_coding",
    "check_cds",
    "check_genome",
]

STATUS_OK = "ok"
STATUS_ALTERNATIVE_START = "alternative_start"
STATUS_NO_START = "no_start"
STATUS_NO_STOP = "no_stop"
STATUS_INTERNAL_STOP = "internal_stop"
STATUS_NOT_MULTIPLE_OF_THREE = "not_multiple_of_three"
STATUS_PSEUDOGENE = "pseudogene"

# Report order: the problems first, then the two classes that are not failures. Every status is
# listed here, so a table built from it shows a zero rather than hiding the row.
ORF_STATUSES: tuple[str, ...] = (
    STATUS_NOT_MULTIPLE_OF_THREE,
    STATUS_NO_START,
    STATUS_NO_STOP,
    STATUS_INTERNAL_STOP,
    STATUS_ALTERNATIVE_START,
    STATUS_PSEUDOGENE,
    STATUS_OK,
)

# Statuses a reader should act on. ``alternative_start`` is deliberately not one of them.
FAILURE_STATUSES: frozenset[str] = frozenset(
    {STATUS_NOT_MULTIPLE_OF_THREE, STATUS_NO_START, STATUS_NO_STOP, STATUS_INTERNAL_STOP}
)

STOP_CODONS: frozenset[str] = frozenset({"TAA", "TAG", "TGA"})
CANONICAL_START = "ATG"

# Translation tables this toolkit can translate: tables 1 and 11 share every amino-acid
# assignment and differ only in initiation codons (D3-03), so one code table answers both.
# Any other table is refused — see OrfError.
START_CODONS: dict[int, tuple[str, ...]] = {
    1: ("TTG", "CTG", "ATG"),  # NCBI table 1 (standard)
    11: ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG", "ACG"),  # plastid, plus ACG below
}

# ACG is NOT in NCBI table 11's initiator list, yet curated plastid records use it for rpl2
# (monocots) and ndhD (dicots), where C-to-U RNA editing converts it to AUG in the transcript.
# Measured 2026-10-03: 6 ACG starts across the project's panel, 2 of them in the golden genome.
# Classing it "no_start" would be a false alarm on a published record, so it is accepted as a
# start — and reported as non-canonical so nobody reads it as an ordinary one.
NON_CANONICAL_STARTS: frozenset[str] = frozenset({"ACG"})

START_CLASS_CANONICAL = "canonical"
START_CLASS_ALTERNATIVE = "alternative"
START_CLASS_NON_CANONICAL = "non_canonical"
START_CLASS_INVALID = "invalid"


class OrfError(ValueError):
    """An ORF check cannot be made honestly for this feature."""


@dataclass(frozen=True, slots=True)
class OrfCheck:
    """One CDS and what its own sequence says about it. Facts first, ``status`` derived from them.

    ``start_class`` and the codon facts stay visible whichever way ``status`` falls: a CDS that
    fails on its frame may still carry a perfectly good start codon, and a reader checking the
    report should see that rather than a single word.
    """

    gene: Gene
    transl_table: int
    codon_start: int
    length_bp: int
    start_codon: str
    stop_codon: str
    start_class: str  # canonical | alternative | non_canonical | invalid
    internal_stop_positions: tuple[int, ...] = ()  # 1-based codon indices inside the CDS

    @property
    def length_is_multiple_of_three(self) -> bool:
        return self.length_bp % 3 == 0

    @property
    def status(self) -> str:
        """The headline class. A pseudogene is reported as one, never as a failed CDS."""
        if self.gene.pseudo:
            return STATUS_PSEUDOGENE
        if not self.length_is_multiple_of_three:
            return STATUS_NOT_MULTIPLE_OF_THREE
        if self.start_class == START_CLASS_INVALID:
            return STATUS_NO_START
        if self.stop_codon not in STOP_CODONS:
            return STATUS_NO_STOP
        if self.internal_stop_positions:
            return STATUS_INTERNAL_STOP
        if self.start_class in (START_CLASS_ALTERNATIVE, START_CLASS_NON_CANONICAL):
            return STATUS_ALTERNATIVE_START
        return STATUS_OK

    @property
    def is_failure(self) -> bool:
        return self.status in FAILURE_STATUSES


def check_cds(genome: Genome, gene: Gene, *, default_table: int = 1) -> OrfCheck:
    """Check one CDS. Raises :class:`OrfError` when the check cannot be made honestly."""
    if gene.kind != "CDS":
        raise OrfError(f"{gene.name}: kind is {gene.kind!r}, not CDS — no reading frame to check")
    table = gene.transl_table if gene.transl_table is not None else default_table
    if table not in START_CODONS:
        raise OrfError(
            f"{gene.name}: translation table {table} is not one this toolkit can translate "
            f"(supported: {sorted(START_CODONS)}) — refusing to translate a code we do not know"
        )
    if gene.codon_start != 1:
        raise OrfError(
            f"{gene.name}: /codon_start={gene.codon_start} — only 1 is supported; any other "
            f"value shifts the reading frame, and translating it as 1 answers a question the "
            f"record did not ask"
        )

    sequence = coding_sequence(genome, gene)
    complete_codons = len(sequence) // 3
    internal = tuple(
        index + 1
        for index in range(1, complete_codons - 1)
        if sequence[index * 3 : index * 3 + 3] in STOP_CODONS
    )
    return OrfCheck(
        gene=gene,
        transl_table=table,
        codon_start=gene.codon_start,
        length_bp=len(sequence),
        start_codon=sequence[:3],
        stop_codon=sequence[-3:],
        start_class=_start_class(sequence[:3], table),
        internal_stop_positions=internal,
    )


def check_genome(genome: Genome, *, default_table: int = 1) -> tuple[OrfCheck, ...]:
    """Every ``CDS`` feature, both IR copies, pseudogenes included (and marked)."""
    return tuple(
        check_cds(genome, gene, default_table=default_table)
        for gene in genome.genes
        if gene.kind == "CDS"
    )


def annotated_not_coding(genome: Genome) -> tuple[Gene, ...]:
    """``gene`` features no same-named CDS/tRNA/rRNA implements — the truncated-ycf1 shape.

    A ``gene`` feature is a claim; the feature under it is what makes the claim real. Both
    halves of the match matter, and the golden genome is why (measured 2026-10-03):

    * by **name**, because the truncated ycf1 copy shares its name with the real ycf1 CDS
      elsewhere in the record — matching on position alone would report nothing, and a
      ``pseudo`` ycf1 with no CDS of its own is exactly the case (PROJECT.md §9.B-10);
    * by **position**, because the truncated copy's 900 bp overlap the first 21 bp of *ndhF*
      — matching on name alone would excuse it, and matching on "any CDS overlaps" would let a
      neighbouring gene's frame stand in for the missing one.
    """
    implementing = tuple(
        (gene.name, gene) for gene in genome.genes if gene.kind in ("CDS", "tRNA", "rRNA")
    )
    return tuple(
        gene
        for gene in genome.genes
        if gene.kind == "gene"
        and not any(
            name == gene.name
            and any(part.overlaps(other) for part in gene.parts for other in feature.parts)
            for name, feature in implementing
        )
    )


def _start_class(codon: str, table: int) -> str:
    if codon == CANONICAL_START:
        return START_CLASS_CANONICAL
    if codon not in START_CODONS[table]:
        return START_CLASS_INVALID
    if codon in NON_CANONICAL_STARTS:
        return START_CLASS_NON_CANONICAL
    return START_CLASS_ALTERNATIVE
