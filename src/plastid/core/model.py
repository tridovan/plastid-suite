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
    "RepeatRegion",
    "cds_features",
    "coding_sequence",
    "count_genes",
    "gene_qualifiers",
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
    # /exception travels with it.
    exception: str | None = None
    # the record's own /translation — the curator's protein, kept **because it is the thing an
    # annotation is transferred from** (D5-08): the search layer aligns a query against
    # these, never against a translation this toolkit made. It is read, never written: /translation
    # is NCBI's to assign on a new record (see ``gene_qualifiers``).
    translation: str | None = None

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

    @property
    def partial(self) -> bool:
        """Any part carries a partial coordinate — a truncated 5' or 3' end.

        The flag is on the part, not the gene, so any part order gives the same answer. Every CDS
        analysis states whether such genes are included, and counts them (M11.17); silently
        merging a truncated gene into a measured value is what this property exists to stop.
        """
        return any(part.partial_start or part.partial_end for part in self.parts)


def gene_qualifiers(gene: Gene) -> dict[str, str | None]:
    """The qualifiers this toolkit writes for a feature, in a stable order.

    Shared by both writers — the five-column table and GenBank — so one feature cannot come out
    differently in the two formats. Two of the choices here are decisions rather than formatting:

    * **``/gene`` is written only when it says something ``/product`` does not.** Measured
      2026-10-03: the stand-in mitogenome ``OX459128.1`` carries no ``/gene`` on any feature, so
      our reader names those features by their product sentence. Writing both would put a
      sentence where a gene *symbol* belongs; writing only the product reads back identically.
    * **``None`` means a bare qualifier** (``/pseudo``, ``/trans_splicing``). Every real record
      writes those with no value, and BioPython writes them bare for ``None`` — an empty string
      comes out as ``=""``, which is a different file.

    ``/locus_tag``, ``/protein_id`` and ``/translation`` are deliberately absent: the first is a
    submitter's assignment this toolkit does not make, and the other two are NCBI's to assign.
    """
    qualifiers: list[tuple[str, str | None]] = []
    if gene.product is None or gene.name != gene.product:
        qualifiers.append(("gene", gene.name))
    if gene.product:
        qualifiers.append(("product", gene.product))
    if gene.exception:
        qualifiers.append(("exception", gene.exception))
    if gene.kind == "CDS":
        if gene.transl_table is not None:
            qualifiers.append(("transl_table", str(gene.transl_table)))
        if gene.codon_start != 1:
            qualifiers.append(("codon_start", str(gene.codon_start)))
    if gene.trans_splicing:
        qualifiers.append(("trans_splicing", None))
    if gene.pseudo:
        qualifiers.append(("pseudo", None))
    return dict(qualifiers)


@dataclass(frozen=True, slots=True)
class RepeatRegion:
    """One ``repeat_region`` feature: one copy, even when its location is a ``join``.

    A GenBank location that crosses the origin is written as ``join(m..L,1..k)``; BioPython
    hands those over as separate parts. Keeping the parts together is what makes the copy whole
    (M11.3): a copy one base short no longer pairs with the other, and a flattened list cannot say
    which piece belongs to which copy. ``parts`` is the one source; everything else is derived.

    ``start``/``end``/``length`` are the outer bounds in the record's frame. For a wrapped copy
    they are not a useful linear interval — use :meth:`arc` with the genome length when the
    circular extent matters.
    """

    parts: tuple[Interval, ...]

    def __post_init__(self) -> None:
        if not self.parts:
            raise ValueError("a repeat_region has no location parts")

    @property
    def length(self) -> int:
        return sum(part.length for part in self.parts)

    @property
    def strand(self) -> int:
        strands = {part.strand for part in self.parts}
        return strands.pop() if len(strands) == 1 else 0

    @property
    def start(self) -> int:
        return min(part.start for part in self.parts)

    @property
    def end(self) -> int:
        return max(part.end for part in self.parts)

    @property
    def partial_start(self) -> bool:
        first = min(self.parts, key=lambda part: part.start)
        return first.partial_start

    @property
    def partial_end(self) -> bool:
        last = max(self.parts, key=lambda part: part.end)
        return last.partial_end

    def arc(self, genome_length: int) -> tuple[int, int]:
        """``(start, length)`` of the single circular arc the parts cover.

        A joined copy is contiguous on the circle: its parts may be written in either order, and
        the gap between them (at most one) is the origin. More than one gap means the location is
        not one repeat — a malformed annotation, refused by the caller rather than silently
        stitched. ``start`` is the first base after the gap, so a wrapped copy has
        ``start + length > genome_length``.
        """
        if genome_length <= 0:
            raise ValueError("genome_length must be positive")
        parts = sorted(self.parts, key=lambda part: (part.start, part.end))
        for earlier, later in zip(parts, parts[1:], strict=False):
            if earlier.end > later.start:
                raise ValueError(
                    f"repeat_region parts overlap: {earlier.start + 1}..{earlier.end} and "
                    f"{later.start + 1}..{later.end}"
                )
            if later.end > genome_length:
                raise ValueError(
                    f"repeat_region part {later.start + 1}..{later.end} is past the sequence "
                    f"end ({genome_length} bp)"
                )
        gaps = [
            (earlier.end, later.start)
            for earlier, later in zip(parts, parts[1:], strict=False)
            if later.start > earlier.end
        ]
        circular_gap = parts[-1].end, parts[0].start + genome_length
        if circular_gap[1] > circular_gap[0]:
            gaps.append(circular_gap)
        if len(gaps) > 1:
            raise ValueError(
                "repeat_region parts do not form one circular arc "
                f"({len(gaps)} gaps: {[(a + 1, b) for a, b in gaps]})"
            )
        if not gaps:
            # The parts tile the whole circle. That is not an IR copy; keep the answer explicit.
            return 0, genome_length
        gap_end = gaps[0][1] % genome_length
        return gap_end, self.length


@dataclass(slots=True)
class Genome:
    """One organelle genome. ``sequence`` is uppercase; ``genes`` keeps every copy (both IRs)."""

    name: str
    sequence: str
    accession: str | None = None
    topology: str = "circular"
    genes: tuple[Gene, ...] = ()
    # One entry per ``repeat_region`` feature. A bare Interval is accepted for callers that build
    # a single-part region and is wrapped here; there is no second, independently mutable form.
    repeat_regions: tuple[RepeatRegion, ...] = ()
    # the record's own source qualifiers (organism · organelle · mol_type · cultivar · …). They
    # travel with the record for the same reason /transl_table does: a submission has to state
    # them, and inventing an organism would be worse than refusing. Measured 2026-10-03 — a
    # table without them fails NCBI's own validator twice.
    source: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Single-part callers (tests, annotation writers) may pass an Interval. Normalizing here
        # keeps one representation in the model instead of two that can drift.
        self.repeat_regions = tuple(
            region if isinstance(region, RepeatRegion) else RepeatRegion((region,))
            for region in self.repeat_regions
        )

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
