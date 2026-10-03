"""Core data model shared by every module.

Coordinate convention: internal coordinates are **0-based half-open** ``[start, end)``
(the same convention BioPython uses). GenBank and every exported table use **1-based
inclusive**; convert only through :meth:`Interval.to_1based` /
:meth:`Interval.from_1based` — no other module converts coordinates.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

__all__ = ["GENE_KINDS", "Gene", "GeneCounts", "Genome", "Interval", "count_genes"]

# The feature kinds counted as "genes" — the fixed rule is 132 = 86 CDS + 38 tRNA + 8 rRNA.
GENE_KINDS: tuple[str, ...] = ("CDS", "tRNA", "rRNA")


@dataclass(frozen=True, slots=True)
class Interval:
    """A 0-based half-open interval ``[start, end)`` on one strand."""

    start: int
    end: int
    strand: int = 1

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
