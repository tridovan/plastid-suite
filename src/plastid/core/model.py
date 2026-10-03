"""Core data model shared by every module.

Coordinate convention: internal coordinates are **0-based half-open** ``[start, end)``
(the same convention BioPython uses). GenBank and every exported table use **1-based
inclusive**; convert only through :meth:`Interval.to_1based` /
:meth:`Interval.from_1based` — no other module converts coordinates.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from plastid.core.sequence import reverse_complement

__all__ = [
    "GENE_KINDS",
    "Gene",
    "GeneCounts",
    "Genome",
    "Interval",
    "cds_features",
    "coding_sequence",
    "count_genes",
]

# The feature kinds counted as "genes" — the fixed rule is 132 = 86 CDS + 38 tRNA + 8 rRNA.
GENE_KINDS: tuple[str, ...] = ("CDS", "tRNA", "rRNA")


@dataclass(frozen=True, slots=True)
class Interval:
    """A 0-based half-open interval ``[start, end)`` on one strand.

    ``partial_start`` / ``partial_end`` say that the *coordinate* is incomplete in the record
    (GenBank ``<`` / ``>``, i.e. BioPython's ``BeforePosition`` / ``AfterPosition``) — they are
    about the coordinate, not about the biological 5'/3' end, because on the minus strand those
    are opposite ends. The submission writer needs them: dropping a partial marker silently
    turns "this gene runs off the contig" into "this gene starts here", which is a different
    claim in a submitted record.
    """

    start: int
    end: int
    strand: int = 1
    partial_start: bool = False
    partial_end: bool = False

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"end ({self.end}) must be >= start ({self.start})")
        if self.strand not in (1, -1):
            raise ValueError(f"strand must be +1 or -1, got {self.strand}")

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlaps(self, other: Interval) -> bool:
        """True when the two half-open intervals share at least one base."""
        return self.start < other.end and other.start < self.end

    def to_1based(self) -> tuple[int, int]:
        """Return GenBank-style inclusive coordinates."""
        return (self.start + 1, self.end)

    @classmethod
    def from_1based(cls, start: int, end: int, strand: int = 1) -> Interval:
        """Build from GenBank-style inclusive coordinates."""
        return cls(start=start - 1, end=end, strand=strand)


@dataclass(frozen=True, slots=True)
class Gene:
    """A feature: one part for simple locations, several for ``join()`` / trans-splicing."""

    name: str
    kind: str
    parts: tuple[Interval, ...]
    product: str | None = None
    pseudo: bool = False
    # the record's own /trans_splicing qualifier — the annotation's word, not our inference
    trans_splicing: bool = False
    # the record's own /transl_table; None means the record declares none (fall back to defaults).
    # Never assume a table: measured 2026-10-03, every CDS in the project's dataset declares 11,
    # and P3-01 refuses a table it cannot translate rather than translating the wrong one.
    transl_table: int | None = None
    # the record's own /codon_start; anything but 1 shifts the reading frame (D3-03)
    codon_start: int = 1
    # the record's own /exception — a statement only the submitter can make ("RNA editing",
    # "trans-splicing", …). Kept because dropping it is not neutral: measured 2026-10-03, a CDS
    # whose start codon is created by RNA editing is rejected by NCBI's validator unless its
    # /exception travels with it (PROJECT.md §8.18).
    exception: str | None = None

    def __post_init__(self) -> None:
        if not self.parts:
            raise ValueError(f"gene {self.name!r} has no location parts")

    @property
    def strand(self) -> int:
        """+1 / -1 when all parts agree; 0 for mixed-strand (trans-spliced) features."""
        strands = {part.strand for part in self.parts}
        return strands.pop() if len(strands) == 1 else 0

    @property
    def span(self) -> tuple[int, int]:
        """Outer bounds, 1-based inclusive."""
        return (
            min(part.start for part in self.parts) + 1,
            max(part.end for part in self.parts),
        )

    @property
    def length(self) -> int:
        """Sum of part lengths — exons only for spliced features."""
        return sum(part.length for part in self.parts)


@dataclass(slots=True)
class Genome:
    """One organelle genome. ``sequence`` is uppercase; ``genes`` keeps every copy (both IRs)."""

    name: str
    sequence: str
    accession: str | None = None
    topology: str = "circular"
    genes: tuple[Gene, ...] = ()
    repeat_regions: tuple[Interval, ...] = ()
    # the record's own source qualifiers (organism · organelle · mol_type · cultivar · …). They
    # travel with the record for the same reason /transl_table does: a submission has to state
    # them, and inventing an organism would be worse than refusing. Measured 2026-10-03 — a
    # table without them fails NCBI's own validator twice (§8.18).
    source: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.sequence)

    def subseq(self, interval: Interval) -> str:
        """Forward-strand slice; circular genomes wrap around the origin."""
        length = len(self.sequence)
        if interval.length > length:
            raise ValueError(f"interval length {interval.length} exceeds genome length {length}")
        if self.topology != "circular":
            return self.sequence[interval.start : interval.end]
        return "".join(self.sequence[i % length] for i in range(interval.start, interval.end))


@dataclass(frozen=True, slots=True)
class GeneCounts:
    """The fixed counting rule: 132 = 86 CDS + 38 tRNA + 8 rRNA.

    ``gene_features`` counts ``gene`` features and is a **diagnostic only** — it must
    never be reported as the gene count (CLAUDE.md invariant).
    """

    cds: int
    trna: int
    rrna: int
    gene_features: int
    unique_names: int
    repeated_names: int

    @property
    def total(self) -> int:
        return self.cds + self.trna + self.rrna


def cds_features(genome: Genome) -> tuple[Gene, ...]:
    """The counted CDS set: every ``CDS`` feature, both IR copies, pseudogenes excluded.

    One source of truth for the gene-based analyses (P2-03 codon usage, P2-04 composition) — the
    golden genome's set is 86, the number ``count_genes`` reports.
    """
    return tuple(gene for gene in genome.genes if gene.kind == "CDS" and not gene.pseudo)


def coding_sequence(genome: Genome, gene: Gene) -> str:
    """A feature's sequence 5'→3': every part read on **its own** strand, in annotation order.

    A trans-spliced CDS may mix strands and GenBank writes a join() in transcript order, so
    neither re-sorting the parts nor reading them all on the gene's strand is safe. Measured on
    both golden rps12 copies (2026-10-03): only this reading yields the canonical ``ATG…TAA``
    ORF — re-sorting gives ``AAT…TAT``, reading everything forward gives ``ATA…ATT``.
    """
    chunks = [
        genome.subseq(part) if part.strand > 0 else reverse_complement(genome.subseq(part))
        for part in gene.parts
    ]
    return "".join(chunks)


def count_genes(genome: Genome) -> GeneCounts:
    """Count genes by the fixed rule, over *all* copies (both IR copies included)."""
    by_kind = Counter(gene.kind for gene in genome.genes if gene.kind in GENE_KINDS)
    gene_features = [gene for gene in genome.genes if gene.kind == "gene"]
    names = Counter(gene.name for gene in gene_features)
    return GeneCounts(
        cds=by_kind.get("CDS", 0),
        trna=by_kind.get("tRNA", 0),
        rrna=by_kind.get("rRNA", 0),
        gene_features=len(gene_features),
        unique_names=len(names),
        repeated_names=sum(1 for count in names.values() if count > 1),
    )
