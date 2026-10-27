"""The reference's own exon structure, placed on the target by aligning its own DNA (P3B-10).

**Why this exists.** A call is transferred from a reference *protein*, and blastx seeds on
amino-acid words: an exon shorter than a word is invisible to it, and so is any boundary the
alignment happens to run straight through. Measured on the corpus (34 records), that
is where most of this annotator's remaining differences came from and none of them were noise:

* the first exons of ``petB`` (6 bp) · ``petD`` (8) · ``rpl16`` (9) · ``rps16`` (40) are missed in
  **every** record — 136 conclusion-level differences, the top of the extent table;
* ``rps12``'s third exon (26 bp) is missed wherever the gene is written with three exons;
* and where an intron sits *inside* one alignment the protein search has no evidence at all about
  where it starts: ``ycf3`` landed 5 and 2 bp away from the published boundaries, ``rpoC1`` 15 and
  12, ``rpl2`` 4 — the local alignment ends where the score stops paying, not where the exon does.

A 6 bp exon **cannot** be found by searching for it: six bases match a 3 kb window by chance about
twice per frame, so a search for the exon would be guessing with extra steps. What can be brought
is the reference's **own structure** — the exon boundaries its record already states — and its own
**DNA**, which is what places them.

**How.** For every internal boundary of a reference gene, this module builds a *junction*: the
exon before, the whole intron, and the exon after, written 5'→3' in the orientation **those two
exons share**. The target is searched with that (``external/blast.py``, ``run_blastn_trace``), and
because the hit carries its alignment, each boundary maps to an exact position on the target — gaps
included, which coordinate arithmetic alone would get wrong by the intron's own indels. Per pair,
never per gene: a *trans*-spliced gene's exons may lie on both strands of an inverted repeat
(``NC_072562.1`` writes ``rps12`` as ``-1, +1, +1``), and the pair across the trans-splice has no
junction while the pair that carries the 26 bp exon still does.

A junction whose intron diverged in length can come back as **two** alignments, one exon each
(measured 27/10/2026: ``MW201730``'s ``rps16``, 766 bp + 421 bp). Each exon is still mapped end to
end by an alignment of the junction's own DNA, so the pair is placed under the same rules as one.

**What it refuses to do.** Place anything an alignment did not put there: a boundary is only moved
when the junction aligned *and* both ends of the exon mapped, and a hit is only used when the
exons it places actually land on the call — which is what stops an inverted-repeat gene from being
moved onto the other copy. Nothing is inferred from an intron's length or from a reading frame.

**Conventions.** Everything here is **1-based inclusive**, the GenBank and BLAST convention, in
genome coordinates. Offsets *within* a junction are positions in its own sequence, also 1-based.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from io import StringIO
from typing import Protocol

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plastid.core.model import Gene, Genome, Interval
from plastid.core.sequence import reverse_complement

__all__ = [
    "MAX_INTRON",
    "Alignment",
    "Junction",
    "Placement",
    "Transfer",
    "junctions",
    "place",
    "render_junctions",
]

# The longest intron a junction is built for. Measured 10/10/2026 over every multi-exon CDS of the
# pinned set: the longest real intron is **2,815 bp** (``rpoC1``), and the one thing far above it
# is ``rps12``'s **69,056 bp** gap between its LSC exon and its IR exons — which is not an intron
# at all but the distance across the repeat between two transcripts that are spliced in *trans*.
# A junction that long is one alignment across the whole inverted repeat, and it would answer the
# wrong question; skipping it costs nothing, because the boundary it would place is already
# inside a protein alignment (measured: every record's ``rps12`` exon 1 ends exactly where the
# published record ends it).
MAX_INTRON = 10_000


class Alignment(Protocol):
    """What the transfer needs of one hit: a subject position's place on the target, and its score.

    Stated as a protocol so this module depends on the *shape* of an alignment rather than on the
    program that produced it — the layering rule keeps ``core`` free of ``external``.
    """

    bitscore: float

    def map_subject(self, position: int) -> int | None:
        """Where a 1-based position of the junction sits on the target, or ``None`` if unaligned."""
        ...


@dataclass(frozen=True, slots=True)
class Junction:
    """One internal boundary of one reference gene, as DNA a target can be searched with.

    ``sequence`` is the reference's own exon *k*, its whole intron, and its exon *k+1*, written
    5'→3' in the gene's orientation, so an alignment of it answers one question: where are these
    two exons on the target, and the intron between them. ``index`` is exon *k*'s position in the
    gene's exon list, which is how a placement is put back in transcript order.
    """

    key: str
    sequence: str
    left: tuple[int, int]  # exon k's offsets within `sequence`, 1-based inclusive
    right: tuple[int, int]  # exon k+1's offsets, likewise
    boundary: int  # the offset of the intron's first base
    index: int  # exon k's index in the gene's exon list
    exons: int  # how many exons the gene has
    intron: int  # the intron's length in the reference
    first: bool  # exon k is the gene's first exon
    last: bool  # exon k+1 is the gene's last exon


@dataclass(frozen=True, slots=True)
class Transfer:
    """One call's own alignment, as the transfer needs to know it — genome coordinates, 1-based."""

    symbol: str
    strand: int
    spans: tuple[tuple[int, int], ...]  # the call's exons, transcript order
    partial_start: bool
    partial_end: bool


@dataclass(frozen=True, slots=True)
class Placement:
    """The reference's own structure as it sits on the target, and what had to be said about it.

    An end is *established* when the placed run contains the reference's own first (or last) exon:
    the gene then starts (ends) where an alignment of the reference's own start (stop) codon put
    it, which is what the caller's partial flags mean. A run that stops short of it leaves the
    target's gene partial there — the structure says the gene continues, and this module did not
    find where.
    """

    spans: tuple[tuple[int, int], ...]  # 1-based inclusive, transcript order
    partial_start: bool
    partial_end: bool
    notes: tuple[str, ...]


def junctions(genome: Genome, gene: Gene, prefix: str) -> tuple[Junction, ...]:
    """Every internal boundary of one reference gene whose two exons sit on one strand.

    ``prefix`` is the first two fields of the search id and the start of the third; each junction
    appends its own index, so a hit comes back naming the gene *and* the boundary it is about.

    **Per pair, never per gene.** A record may spell one gene's exons on more than one strand —
    measured 10/10/2026, ``NC_072562.1`` and ``NC_072563.1`` write ``rps12`` as ``-1, +1, +1``,
    which is what a *trans*-spliced gene looks like when its halves come from opposite strands of
    an inverted repeat. Two exons of one pair are still colinear, and their junction is still the
    DNA between them; only the pair that spans the trans-splice is skipped, because the DNA
    between *those* exons is 69 kb of repeat and not an intron.

    A boundary is also skipped when the gap is longer than :data:`MAX_INTRON`, or when the two
    exons are not in gene order at all — searching either would place a boundary from an alignment
    of something that is not the intron.
    """
    out: list[Junction] = []
    for index in range(len(gene.parts) - 1):
        pair = _pair(genome, gene.parts[index], gene.parts[index + 1])
        if pair is None:
            continue
        sequence, left, right, intron = pair
        out.append(
            Junction(
                key=f"{prefix}{index}",
                sequence=sequence,
                left=left,
                right=right,
                boundary=left[1] + 1,
                index=index,
                exons=len(gene.parts),
                intron=intron,
                first=index == 0,
                last=index + 1 == len(gene.parts) - 1,
            )
        )
    return tuple(out)


def _pair(
    genome: Genome, left: Interval, right: Interval
) -> tuple[str, tuple[int, int], tuple[int, int], int] | None:
    """Exon k + the intron + exon k+1, read on the strand the two exons share, or ``None``.

    Offsets are into the junction's own sequence, 1-based inclusive, and the two exons must come
    in gene order: anything else is not an intron this module may place a boundary from.
    """
    if left.strand != right.strand or left.strand == 0:
        return None
    if left.strand > 0:
        if right.start < left.end:
            return None
        sequence = genome.sequence[left.start : right.end]
        first = (1, left.end - left.start)
        second = (right.start - left.start + 1, right.end - left.start)
    else:
        if right.end > left.start:
            return None
        sequence = reverse_complement(genome.sequence[right.start : left.end])
        high = left.end
        first = (high - left.end + 1, high - left.start)
        second = (high - right.end + 1, high - right.start)
    intron = second[0] - first[1] - 1
    if intron < 0 or intron > MAX_INTRON:
        return None
    return sequence, first, second, intron


def render_junctions(junctions: Iterable[Junction]) -> str:
    """The junctions as the FASTA the search runs against — BioPython, as every writer here.

    Each record's id is the junction's own :attr:`Junction.key`, which is what a hit comes back
    as; written through the same three-field rule the reference headers follow, so a key that
    survives BLAST's truncation is the key that was written.
    """
    handle = StringIO()
    SeqIO.write(
        [
            SeqRecord(Seq(junction.sequence), id=junction.key, description="")
            for junction in junctions
        ],
        handle,
        "fasta",
    )
    return handle.getvalue()


def place(
    transfer: Transfer,
    junctions: Sequence[Junction],
    hits: Mapping[str, Iterable[Alignment]],
) -> Placement | None:
    """The reference's structure on the target, or ``None`` when there is nothing to transfer.

    Every boundary is taken from an alignment of the reference's own DNA; a boundary whose
    junction did not align, or whose exon did not map end to end, keeps the coordinate the protein
    alignment gave it. The placement is anchored to the call: a junction is only used when the
    exons it maps land on the call's own exons, so a hit to a gene's other inverted-repeat copy —
    or to a paralog — is rejected rather than believed.
    """
    if not junctions:
        return None
    spans: dict[int, tuple[int, int]] = {}
    notes: list[str] = []
    for junction in junctions:
        chosen = _best(junction, hits.get(junction.key, ()), transfer)
        if chosen is None:
            notes.append(
                f"{transfer.symbol}: the reference's own intron {junction.intron} bp long did not "
                f"align to the target — the protein alignment's own boundary stands, and the "
                f"published record may differ from it by however far that alignment ran on"
            )
            continue
        left, right = chosen
        spans[junction.index] = left
        spans[junction.index + 1] = right
    if not spans:
        return None
    anchors = _anchors(transfer, spans)
    run = _run(transfer, spans)
    if run is None:
        return None
    exons = junctions[0].exons
    # ``index`` is the reference exon the span is, or ``None`` for one of the call's own exons
    # that no junction reached — a boundary whose reference exon is across ``rps12``'s 69 kb
    # trans-splice, or a junction that simply did not align. Those are kept exactly as the
    # protein alignment called them: dropping one would delete a real exon from the record.
    ordered = _interleave(transfer, spans, anchors, run)
    # The call's own ends are kept where it established them: its initiator and stop tests are
    # evidence about *this* genome, and the placement is evidence about where the reference's
    # structure sits. They only disagree when the alignment ran short, which is when the
    # placement is the better of the two.
    index, span = ordered[0]
    if index == 0 and anchors.get(0) == 0 and not transfer.partial_start:
        ordered[0] = (0, _with_end(span, transfer.spans[0], transfer.strand, five_prime=True))
    index, span = ordered[-1]
    ends_the_gene = index == exons - 1 and anchors.get(index) == len(transfer.spans) - 1
    if ends_the_gene and not transfer.partial_end:
        ordered[-1] = (
            index,
            _with_end(span, transfer.spans[-1], transfer.strand, five_prime=False),
        )
    placed = sum(end - start + 1 for _index, (start, end) in ordered)
    was = sum(end - start + 1 for start, end in transfer.spans)
    notes.append(
        f"{transfer.symbol}: the reference's own exon structure was placed on the target — "
        f"{len(ordered)} exon(s), {placed} bp against the protein alignment's {was} bp; the "
        f"boundaries it moved come from an alignment of the reference's own introns, gaps included"
    )
    return Placement(
        spans=tuple(span for _index, span in ordered),
        partial_start=transfer.partial_start if ordered[0][0] is None else ordered[0][0] > 0,
        partial_end=transfer.partial_end if ordered[-1][0] is None else ordered[-1][0] < exons - 1,
        notes=tuple(notes),
    )


def _interleave(
    transfer: Transfer,
    spans: Mapping[int, tuple[int, int]],
    anchors: Mapping[int, int],
    run: Sequence[int],
) -> list[tuple[int | None, tuple[int, int]]]:
    """The call's own exons and the placed ones, in transcript order.

    The call's exons are the spine, because *their* order is the one the protein alignment
    established; a placed exon is written where the reference says it belongs among them, which is
    what puts ``petB``'s 6 bp first exon in front of the alignment's exon and leaves ``rps12``'s
    first exon exactly where the caller put it.
    """
    matched = {part: index for index, part in anchors.items()}
    out: list[tuple[int | None, tuple[int, int]]] = []
    cursor = run[0] - 1
    for part, span in enumerate(transfer.spans):
        index = matched.get(part)
        if index is None:
            out.append((None, span))
            continue
        out.extend((step, spans[step]) for step in range(cursor + 1, index + 1) if step in spans)
        cursor = index
    out.extend((step, spans[step]) for step in range(cursor + 1, run[-1] + 1) if step in spans)
    return out


def _best(
    junction: Junction, hits: Iterable[Alignment], transfer: Transfer
) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """The target spans of the junction's two exons — what answers *this* call's question.

    One alignment carrying both exons is the ordinary shape and is preferred. The search also
    splits a junction into two alignments where the intron's own length diverges — measured
    27/10/2026, ``MW201730``'s ``rps16`` junction came back as 766 bp + 421 bp, one exon each —
    and each exon is still mapped end to end by an alignment of the junction's own DNA. That pair
    is accepted under the same rules as one hit: the two exons must sit in gene order and at least
    one of them must land on the call.
    """
    hits = list(hits)
    whole: tuple[int, float, tuple[int, int], tuple[int, int]] | None = None
    lefts: list[tuple[float, tuple[int, int]]] = []
    rights: list[tuple[float, tuple[int, int]]] = []
    for hit in hits:
        left = _mapped(hit, junction.left)
        right = _mapped(hit, junction.right)
        if left is not None:
            lefts.append((hit.bitscore, left))
        if right is not None:
            rights.append((hit.bitscore, right))
        if left is None or right is None or not _ordered(transfer.strand, left, right):
            continue
        overlap = sum(_anchored(span, transfer.spans) for span in (left, right))
        if not overlap:
            continue
        if whole is None or (overlap, hit.bitscore) > (whole[0], whole[1]):
            whole = (overlap, hit.bitscore, left, right)
    if whole is not None:
        return whole[2], whole[3]
    pair: tuple[int, float, tuple[int, int], tuple[int, int]] | None = None
    for left_score, left in lefts:
        for right_score, right in rights:
            if not _ordered(transfer.strand, left, right):
                continue
            overlap = sum(_anchored(span, transfer.spans) for span in (left, right))
            if not overlap:
                continue
            if pair is None or (overlap, left_score + right_score) > (pair[0], pair[1]):
                pair = (overlap, left_score + right_score, left, right)
    if pair is None:
        return None
    return pair[2], pair[3]


def _mapped(hit: Alignment, offsets: tuple[int, int]) -> tuple[int, int] | None:
    """The target span the given junction offsets sit on, or ``None`` if either end is unaligned."""
    first, second = hit.map_subject(offsets[0]), hit.map_subject(offsets[1])
    if first is None or second is None:
        return None
    return (min(first, second), max(first, second))


def _ordered(strand: int, left: tuple[int, int], right: tuple[int, int]) -> bool:
    """Do the two exons sit in gene order on the target, an intron apart and not overlapping?"""
    return right[0] > left[1] if strand > 0 else right[1] < left[0]


def _overlap(span: tuple[int, int], part: tuple[int, int]) -> int:
    return max(0, min(span[1], part[1]) - max(span[0], part[0]) + 1)


def _anchored(span: tuple[int, int], parts: Sequence[tuple[int, int]]) -> int:
    """How much of the shorter of the two a placed exon and a call exon share.

    The same rule the corpus measures a match by (half the shorter exon, never by envelope): a
    junction that aligns to a *homologous* place shares a few bases with the call by accident, and
    a shared 26 bases of a 249 bp exon is not two annotations of one exon. Measured 10/10/2026:
    without this, ``NC_072562.1``'s mixed-strand ``rps12`` junction anchored on 26 shared bases and
    would have moved the golden's call 700 bp.
    """
    shorter = span[1] - span[0] + 1
    return max(
        (
            overlap
            for part in parts
            if (overlap := _overlap(span, part)) * 2 >= min(shorter, part[1] - part[0] + 1)
        ),
        default=0,
    )


def _anchors(transfer: Transfer, spans: Mapping[int, tuple[int, int]]) -> dict[int, int]:
    """``reference exon index -> the call's own exon it lands on``, for the exons that do.

    A placed exon that lands on no call exon is a **new** exon (``petB``'s 6 bp one) and is not an
    anchor — which is what keeps it out of the run's own extent, since the run is grown around the
    exons the call and the reference agree about.
    """
    out: dict[int, int] = {}
    for index, span in spans.items():
        best = max(
            range(len(transfer.spans)),
            key=lambda part: _overlap(span, transfer.spans[part]),
            default=None,
        )
        if best is not None and _anchored(span, [transfer.spans[best]]):
            out[index] = best
    return out


def _run(transfer: Transfer, spans: Mapping[int, tuple[int, int]]) -> list[int] | None:
    """The exons to write: the contiguous stretch of placed exons the call's own exons sit in.

    A junction that failed leaves a hole; the exons beyond it are placed but disconnected from the
    call, and writing them would put a second, unattached piece of gene in the record.
    """
    anchored = _anchors(transfer, spans)
    if not anchored:
        return None
    first, last = min(anchored), max(anchored)
    while first - 1 in spans and _ordered(transfer.strand, spans[first - 1], spans[first]):
        first -= 1
    while last + 1 in spans and _ordered(transfer.strand, spans[last], spans[last + 1]):
        last += 1
    return list(range(first, last + 1))


def _with_end(
    placed: tuple[int, int], called: tuple[int, int], strand: int, *, five_prime: bool
) -> tuple[int, int]:
    """``placed`` with its 5' (or 3') coordinate taken from the call's own exon."""
    lo, hi = placed
    if (strand > 0) == five_prime:
        lo = called[0]
    else:
        hi = called[1]
    return (min(lo, hi), max(lo, hi))
