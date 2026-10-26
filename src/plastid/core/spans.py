"""Derived spans — introns, trans-splicing junctions and intergenic spacers (P2B-01, D-g).

Two modules need the same boundaries: ``introns`` publishes the spans, ``diversity`` measures π
over them. The derivation lives here so the two cannot disagree about where a boundary is — the
reason D-g fixed the IGS definition once instead of leaving it to each reading module.

Every span is stated in the **genome's own (GenBank) frame** — the coordinates a record's
features are written in — while the region label comes from the normalized frame. That is why an
:class:`~plastid.core.ir.IRResult` is an argument rather than something a caller could forget to
pass: :meth:`~plastid.core.normalize.NormalizedFrame.region_at` converts between the two frames,
and comparing them directly is right only when the record happens to start at the LSC (measured
10/10/2026: the golden's rotation offset is 0, which is exactly how such a mistake survives).

Definitions:

* an **intron** is the gap between two consecutive parts of one CDS, both ends in one region
  (17 on the golden — the number the ``introns`` module publishes, both IR copies included);
* a **trans-splicing junction** is that same gap when it crosses a region boundary: reported,
  never counted as an intron;
* an **intergenic spacer** (IGS) is the gap between two consecutive **exons** of *different*
  called genes lying inside one single-copy region. A gene occupies its exons, not the space
  between them: taking the outer span instead would let a trans-spliced ``rps12`` swallow the
  genes that sit between its LSC exon and its IR exons, and let ``trnK`` swallow ``matK``, which
  is annotated inside its intron — measured on the golden, both are real spacers that the outer
  span silently loses. The IR is excluded (D-g), because its two copies are the same sequence
  annotated twice and the annotation is not always mirrored base for base — measured on the
  golden: 3 of the 18 spacer pairs do not mirror exactly (ycf1's truncated copy carries no CDS;
  two rRNA spacers differ by one base). IR-resident, boundary-crossing and overlapping gaps are
  counted and named, never silently dropped.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from plastid.core.ir import IRResult
from plastid.core.model import GENE_KINDS, Gene, Genome, Interval, cds_features
from plastid.core.normalize import NormalizedFrame, normalize

__all__ = [
    "Gap",
    "IgsSpan",
    "IntronSet",
    "IntronSpan",
    "JunctionSpan",
    "SpacerSet",
    "introns_of",
    "spacers_of",
]


@dataclass(frozen=True, slots=True)
class IntronSpan:
    """One intron: the gap between two consecutive CDS parts, both ends in one region."""

    intron_id: str  # "ndhA.1.1" — gene.copy.number, the id the `introns` table publishes
    gene: str
    start: int  # 0-based half-open, the genome's own (GenBank) frame
    end: int
    phase: int  # bases of the preceding exon(s) modulo three, in transcript order
    region: str  # "LSC" | "IRa" | "SSC" | "IRb", from the normalized frame

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_1based(self) -> tuple[int, int]:
        return (self.start + 1, self.end)


@dataclass(frozen=True, slots=True)
class JunctionSpan:
    """A gap between two CDS parts whose ends sit in **different** regions: trans-splicing.

    The two parts are joined after transcription, so the span between them is not an intron and
    has no phase — the record's own ``/trans_splicing`` qualifier is cross-checked against this
    geometry by the ``introns`` module, which reports a disagreement rather than smoothing it.
    """

    junction_id: str  # numbered with the introns of the same gene copy: "rps12.1.1"
    gene: str
    start: int
    end: int
    regions: tuple[str, str]  # the regions the two parts sit in, in transcript order

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_1based(self) -> tuple[int, int]:
        return (self.start + 1, self.end)


@dataclass(frozen=True, slots=True)
class IgsSpan:
    """One intergenic spacer: the gap between the exons of two neighbouring called genes."""

    igs_id: str  # "psaB-rps12" — the two flanking genes; ".2" when a pair repeats
    left: str  # upstream gene in genome order
    right: str  # downstream gene
    start: int  # 0-based half-open, the genome's own (GenBank) frame
    end: int
    region: str

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_1based(self) -> tuple[int, int]:
        return (self.start + 1, self.end)


@dataclass(frozen=True, slots=True)
class Gap:
    """A gap between two genes that is **not** a spacer: it crosses a boundary, or overlaps."""

    left: str
    right: str
    length: int  # negative when the two genes overlap


@dataclass(frozen=True, slots=True)
class IntronSet:
    """Every intron and trans-splicing junction of one genome, with the region of each."""

    introns: tuple[IntronSpan, ...]  # both IR copies, as the `introns` module publishes them
    junctions: tuple[JunctionSpan, ...]


@dataclass(frozen=True, slots=True)
class SpacerSet:
    """Every published intergenic spacer, plus everything looked at and not published."""

    igs: tuple[IgsSpan, ...]  # gaps inside LSC or SSC, in genomic order — the published set
    in_ir: tuple[IgsSpan, ...]  # gaps inside the IR: excluded (D-g), counted and named
    crossing: tuple[Gap, ...]  # gaps crossing a region boundary — not spacers
    no_spacer: tuple[Gap, ...]  # gene pairs that overlap or abut — no gap to measure


def introns_of(genome: Genome, ir: IRResult) -> IntronSet:
    """Derive every intron and trans-splicing junction from the CDS part boundaries."""
    frame = normalize(genome, ir)
    introns: list[IntronSpan] = []
    junctions: list[JunctionSpan] = []
    copies: Counter[str] = Counter()
    for gene in cds_features(genome):
        copies[gene.name] += 1
        preceding = 0
        pairs = zip(gene.parts, gene.parts[1:], strict=False)  # one shorter by design
        for number, (first, second) in enumerate(pairs, start=1):
            preceding += first.length
            start = min(first.end, second.end)
            end = max(first.start, second.start)
            # gene.copy.n: the IR carries a second copy of several genes, and a repeated table
            # key would make every reference to it ambiguous
            span_id = f"{gene.name}.{copies[gene.name]}.{number}"
            _check_gap(gene, first, second, end - start)
            head, tail = _region_of(frame, start), _region_of(frame, end - 1)
            if head == tail:
                introns.append(
                    IntronSpan(
                        intron_id=span_id,
                        gene=gene.name,
                        start=start,
                        end=end,
                        phase=preceding % 3,
                        region=head,
                    )
                )
            else:
                junctions.append(
                    JunctionSpan(
                        junction_id=span_id,
                        gene=gene.name,
                        start=start,
                        end=end,
                        regions=(head, tail),
                    )
                )
    return IntronSet(introns=tuple(introns), junctions=tuple(junctions))


def spacers_of(genome: Genome, ir: IRResult) -> SpacerSet:
    """Derive every intergenic spacer, and count the gaps that are not spacers.

    Called genes are the counted set — CDS · tRNA · rRNA, pseudogenes excluded — and each gene
    contributes **each of its exons** as an occupied interval. Two exons of the *same* feature
    are never a spacer: the gap between them is that feature's intron (or its trans-splicing
    junction), which :func:`introns_of` owns.
    """
    frame = normalize(genome, ir)
    length = len(genome)
    occupied = sorted(
        (
            (part.start, part.end, gene, index)
            for index, gene in enumerate(genome.genes)
            if gene.kind in GENE_KINDS and not gene.pseudo
            for part in gene.parts
        ),
        key=lambda entry: (entry[0], entry[1], entry[2].name),
    )
    if len(occupied) < 2:
        raise ValueError(
            f"{genome.accession or genome.name}: {len(occupied)} exon(s) — a spacer needs two "
            f"exons, so this record's annotation cannot produce one"
        )

    igs: list[IgsSpan] = []
    in_ir: list[IgsSpan] = []
    crossing: list[Gap] = []
    no_spacer: list[Gap] = []
    seen: Counter[str] = Counter()
    for index, (_start, end, gene, feature) in enumerate(occupied):
        following_start, _following_end, following, following_feature = occupied[
            (index + 1) % len(occupied)
        ]
        if index + 1 == len(occupied):  # the gap that wraps the origin
            following_start += length
        gap_length = following_start - end
        if feature == following_feature:  # an intron or a junction of one gene — not a spacer
            continue
        if gap_length <= 0:
            no_spacer.append(Gap(gene.name, following.name, gap_length))
            continue
        head = _region_of(frame, end % length)
        tail = _region_of(frame, (following_start - 1) % length)
        if head != tail:
            crossing.append(Gap(gene.name, following.name, gap_length))
            continue
        stem = f"{gene.name}-{following.name}"
        seen[stem] += 1
        suffix = "" if seen[stem] == 1 else f".{seen[stem]}"
        span = IgsSpan(
            igs_id=f"{stem}{suffix}",
            left=gene.name,
            right=following.name,
            start=end,
            end=following_start,
            region=head,
        )
        (in_ir if head in ("IRa", "IRb") else igs).append(span)
    return SpacerSet(
        igs=tuple(igs),
        in_ir=tuple(in_ir),
        crossing=tuple(crossing),
        no_spacer=tuple(no_spacer),
    )


def _region_of(frame: NormalizedFrame, position: int) -> str:
    region = frame.region_at(position)
    if region is None:  # a normalized frame covers every position; this is a broken frame
        raise ValueError(f"position {position} is in no region of the normalized frame")
    return region.name


def _check_gap(gene: Gene, first: Interval, second: Interval, gap: int) -> None:
    """A gap between two parts of a CDS is at least one base — otherwise the record is wrong.

    This replaces the algebraic "exons + gaps tile the span" check the ``introns`` module carried
    until 10/10/2026: measured on that date, the identity holds for **any**
    pair of parts, so that check could never fail. What it silently allowed is the case it was
    meant to catch — parts that overlap or nest give a negative "intron", which was published as
    a length. The gap below is a fact about the record, so this check can fail.
    """
    if gap < 1:
        raise ValueError(
            f"{gene.name}: parts {first.to_1based()} and {second.to_1based()} leave a gap of "
            f"{gap} bp — the record's join() overlaps or is out of order, and no intron can be "
            f"that short"
        )
