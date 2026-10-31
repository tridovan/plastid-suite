"""Reading an assembly without trusting the assembler: circular identity and read support.

An assembler's "complete" says that a path closed through its graph. It does not say that the
path is right, and a closed path over a wrong join is the failure this file exists to expose. Two
independent readings are made here, both pure and both free of the assembler's own code:

* **identity with a reference**, as a circle: any origin, either strand. A plastome has no
  privileged first base, so two records of one molecule can differ in every position of a linear
  comparison while being the same sequence;
* **read support**, from the reads mapped back onto the assembled sequence. One definition covers
  every position of the circle, so the four inverted-repeat junctions and the origin are lookups
  in the same profile rather than special cases, and a plastome with no inverted repeat is read by
  the same rule.

Boundaries. A boundary ``b`` lies between bases ``b - 1`` and ``b``, 0-based, on a circle of length
``L``: ``b`` runs over ``[0, L)`` and boundary 0 is the join of the last base to the first. Reads
are mapped to ``S + S[:w]`` so a read can cross the origin; the extended coordinates are folded
back with ``b % L``, and each read is counted once because only primary alignments are read.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

from Bio.Align import PairwiseAligner

from plastid.core.model import Interval
from plastid.core.sequence import reverse_complement

__all__ = [
    "ASSEMBLY_STATUSES",
    "STATUS_COMPLETE",
    "STATUS_FRAGMENTED",
    "STATUS_GAPPED",
    "STATUS_NO_RESULT",
    "SUPPORT_NOT_MEASURED",
    "SUPPORT_SUPPORTED",
    "SUPPORT_UNDETERMINED",
    "AlignedRead",
    "AssembledPath",
    "AssemblyError",
    "BoundaryCall",
    "ReferenceComparison",
    "SequenceDifference",
    "SupportProfile",
    "SupportRule",
    "WeakRegion",
    "build_profile",
    "call_boundary",
    "compare_to_reference",
    "ir_junctions",
    "orient_to",
    "parse_sam",
    "same_circular",
    "weak_regions",
]

# What the assembler claims. Only "complete" is a closed sequence with every base known, and it is
# a claim about a path, which the read support below tests.
STATUS_COMPLETE = "complete"
STATUS_GAPPED = "gapped"  # a closed path with ambiguous bases
STATUS_FRAGMENTED = "fragmented"  # no closed path: separate scaffolds
STATUS_NO_RESULT = "no_result"  # no plastome sequence at all, which is not a measured zero
ASSEMBLY_STATUSES = (STATUS_COMPLETE, STATUS_GAPPED, STATUS_FRAGMENTED, STATUS_NO_RESULT)

SUPPORT_SUPPORTED = "supported"
SUPPORT_UNDETERMINED = "undetermined"
SUPPORT_NOT_MEASURED = "not measured"

_SKIPPED_FLAGS = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary
# A read does not vouch for what it jumps over. An aligner bridges a stretch of the target that the
# molecule lacks with a deletion, so a read over a wrong join can look as if it crossed it, in bp.
_LONG_INDEL = 10
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_REFERENCE_CONSUMING = frozenset("=XDN")

# Anchors for orienting a sequence. A k-mer present once on one strand places that strand; 31 bp
# occurs by chance far less than once in a genome of 10^5 bases.
_ANCHOR_KMER = 31
_ANCHOR_STEP = 2_000

# Walking two near-identical sequences: after a difference, the next k-mer both share within this
# window is where they run together again, in bp. The window bounds the largest indel that can be
# placed; a difference wider than it leaves the comparison unaligned rather than guessed.
_RESYNC_KMER = 31
_RESYNC_WINDOW = 6_000
_MAX_SEGMENT_CELLS = 4_000_000


class AssemblyError(Exception):
    """An assembly's input or output cannot be read as this step's contract requires."""


@dataclass(frozen=True, slots=True)
class AssembledPath:
    """One way through the assembly graph: one sequence when closed, several contigs when not.

    A plastome with an inverted repeat has two closed paths, because flip-flop recombination gives
    two isomers and reads shorter than the repeat cannot tell them apart. Both are kept.
    """

    name: str
    contigs: tuple[str, ...]
    circular: bool

    @property
    def length(self) -> int:
        """Total bases over the contigs."""
        return sum(len(contig) for contig in self.contigs)


def same_circular(first: str, second: str) -> bool:
    """Whether two sequences are one circular molecule: any origin, either strand."""
    if len(first) != len(second):
        return False
    doubled = second + second
    return first in doubled or reverse_complement(first) in doubled


def _anchors(reference: str, kmer: int) -> list[tuple[int, str]]:
    """Reference k-mers of known bases, at a fixed step: ``(position, k-mer)``, 0-based."""
    picked: list[tuple[int, str]] = []
    for position in range(0, max(1, len(reference) - kmer), _ANCHOR_STEP):
        word = reference[position : position + kmer]
        if len(word) == kmer and set(word) <= set("ACGT"):
            picked.append((position, word))
    return picked


def orient_to(query: str, reference: str, *, kmer: int = _ANCHOR_KMER) -> str | None:
    """``query`` on ``reference``'s strand and origin; None when no anchor places it.

    The strand is the majority of the single-copy anchors, because an isomer of the same plastome
    has its small single-copy region on the other strand and those anchors vote against the rest;
    the origin is the first anchor that agrees with that majority. An anchor inside the inverted
    repeat occurs twice and is not used. The result is good to within the indels before that
    anchor, which an end-to-end alignment absorbs.
    """
    if len(query) <= kmer or len(reference) <= kmer:
        return None
    forward = query + query[: kmer - 1]
    backward_query = reverse_complement(query)
    backward = backward_query + backward_query[: kmer - 1]
    forward_hits: list[tuple[int, int]] = []
    backward_hits: list[tuple[int, int]] = []
    for position, word in _anchors(reference, kmer):
        in_forward, in_backward = forward.count(word), backward.count(word)
        if in_forward == 1 and in_backward == 0:
            forward_hits.append((position, forward.index(word)))
        elif in_backward == 1 and in_forward == 0:
            backward_hits.append((position, backward.index(word)))
    if not forward_hits and not backward_hits:
        return None
    if len(forward_hits) >= len(backward_hits):
        chosen, hits = query, forward_hits
    else:
        chosen, hits = backward_query, backward_hits
    reference_position, query_position = hits[0]
    start = (query_position - reference_position) % len(chosen)
    return chosen[start:] + chosen[:start]


@dataclass(frozen=True, slots=True)
class SequenceDifference:
    """One place two sequences part, 0-based half-open on the oriented query and the reference."""

    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    mismatches: int
    gap_columns: int


@dataclass(frozen=True, slots=True)
class ReferenceComparison:
    """One assembled sequence against one reference, every aligned column counted once."""

    identical: bool  # the same circular molecule, any origin, either strand
    query_length: int
    reference_length: int
    identity: float | None  # matched / every compared column, gaps included; None if unaligned
    matched: int | None
    mismatches: int | None
    gap_columns: int | None
    origin_gap_columns: int | None  # the part of gap_columns at the reference's own join
    method: str
    differences: tuple[SequenceDifference, ...] = ()

    @property
    def length_delta(self) -> int:
        """Query minus reference, in bp."""
        return self.query_length - self.reference_length


def _common_prefix(first: str, second: str, i: int, j: int) -> int:
    """How far ``first[i:]`` and ``second[j:]`` run equal. Doubling, then bisection on slices."""
    limit = min(len(first) - i, len(second) - j)
    if limit <= 0 or first[i] != second[j]:
        return 0
    low, step = 1, 1024
    while low < limit:
        high = min(low + step, limit)
        if first[i + low : i + high] != second[j + low : j + high]:
            limit = high
            break
        low = high
        step *= 2
    while low < limit:
        middle = (low + limit + 1) // 2
        if first[i + low : i + middle] == second[j + low : j + middle]:
            low = middle
        else:
            limit = middle - 1
    return low


def _next_sync(first: str, second: str, i: int, j: int) -> tuple[int, int] | None:
    """The nearest point past ``(i, j)`` where both sequences share a k-mer, or None.

    The first sequence is stepped base by base and each k-mer searched in the window of the second,
    so a substitution resynchronises at once and an indel on either side within the window is
    found. A k-mer of one base repeated says nothing about position, and is stepped over.
    """
    last = min(len(first) - _RESYNC_KMER, i + _RESYNC_WINDOW)
    for offset in range(i, last + 1):
        word = first[offset : offset + _RESYNC_KMER]
        if len(set(word)) == 1:
            continue
        found = second.find(word, j, j + _RESYNC_WINDOW + _RESYNC_KMER)
        if found >= 0:
            return offset, found
    return None


def _segment_columns(a: str, b: str) -> tuple[int, int, int] | None:
    """``(matched, mismatches, gap columns)`` of a global alignment; None when too large."""
    if not a or not b:
        return 0, 0, len(a) + len(b)
    if len(a) * len(b) > _MAX_SEGMENT_CELLS:
        return None
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    alignment = aligner.align(a, b)[0]
    matched = mismatches = gaps = 0
    for x, y in zip(alignment[0], alignment[1], strict=True):
        if x == "-" or y == "-":
            gaps += 1
        elif x == y:
            matched += 1
        else:
            mismatches += 1
    return matched, mismatches, gaps


def _walk(
    query: str, reference: str
) -> tuple[int, int, int, int, tuple[SequenceDifference, ...]] | None:
    """Walk two oriented sequences together: exact runs, and an alignment at each parting.

    Returns ``(matched, mismatches, gap columns, of which at the origin, differences)``, or None
    when a parting is wider than the window, which is what the other isomer of a plastome looks
    like. Both sequences are circles, so what is left over at either end of the walk is not an
    edge effect: it is sequence one circle has at the other's join, and it counts as a difference.
    """
    i = j = matched = mismatches = gaps = origin_gaps = 0
    differences: list[SequenceDifference] = []
    while i < len(query) and j < len(reference):
        run = _common_prefix(query, reference, i, j)
        matched += run
        i += run
        j += run
        if i >= len(query) or j >= len(reference):
            break
        sync = _next_sync(query, reference, i, j)
        if sync is None:
            return None
        next_i, next_j = sync
        columns = _segment_columns(query[i:next_i], reference[j:next_j])
        if columns is None:
            return None
        matched += columns[0]
        mismatches += columns[1]
        gaps += columns[2]
        if i == 0 and j == 0:
            origin_gaps += columns[2]
        differences.append(
            SequenceDifference(
                query_start=i,
                query_end=next_i,
                reference_start=j,
                reference_end=next_j,
                mismatches=columns[1],
                gap_columns=columns[2],
            )
        )
        i, j = next_i, next_j
    left_over = (len(query) - i) + (len(reference) - j)
    if left_over:
        gaps += left_over
        origin_gaps += left_over
        differences.append(
            SequenceDifference(
                query_start=i,
                query_end=len(query),
                reference_start=j,
                reference_end=len(reference),
                mismatches=0,
                gap_columns=left_over,
            )
        )
    return matched, mismatches, gaps, origin_gaps, tuple(differences)


def compare_to_reference(query: str, reference: str) -> ReferenceComparison:
    """Compare two circular sequences, and refuse to invent a number when they cannot be placed.

    The two must be near identical: the other isomer of a reference has a whole region reversed,
    and is reported as unaligned, which is true, rather than as thousands of mismatches, which is
    not.
    """
    unaligned = {
        "identity": None,
        "matched": None,
        "mismatches": None,
        "gap_columns": None,
        "origin_gap_columns": None,
    }
    if same_circular(query, reference):
        return ReferenceComparison(
            identical=True,
            query_length=len(query),
            reference_length=len(reference),
            identity=1.0,
            matched=len(query),
            mismatches=0,
            gap_columns=0,
            origin_gap_columns=0,
            method="exact",
        )
    oriented = orient_to(query, reference)
    if oriented is None:
        return ReferenceComparison(
            identical=False,
            query_length=len(query),
            reference_length=len(reference),
            method="no anchor",
            **unaligned,
        )
    walked = _walk(oriented, reference)
    if walked is None:
        return ReferenceComparison(
            identical=False,
            query_length=len(query),
            reference_length=len(reference),
            method="unaligned region",
            **unaligned,
        )
    matched, mismatches, gaps, origin_gaps, differences = walked
    compared = matched + mismatches + gaps
    return ReferenceComparison(
        identical=False,
        query_length=len(query),
        reference_length=len(reference),
        identity=matched / compared if compared else None,
        matched=matched,
        mismatches=mismatches,
        gap_columns=gaps,
        origin_gap_columns=origin_gaps,
        method="walked",
        differences=differences,
    )


@dataclass(frozen=True, slots=True)
class AlignedRead:
    """One primary alignment: a 0-based half-open span on the extended sequence, and its CIGAR."""

    start: int
    end: int
    left_clip: int  # soft-clipped read bases before ``start``
    right_clip: int  # soft-clipped read bases after ``end``
    matched: int
    mismatched: int
    inserted: int
    deleted: int
    blocks: tuple[tuple[int, int], ...]  # the aligned (=/X) reference spans, half-open
    segments: tuple[tuple[int, int], ...]  # the span, cut at every long indel


def _parse_cigar(position: int, cigar: str) -> AlignedRead:
    cursor = segment_start = position
    left_clip = right_clip = matched = mismatched = inserted = deleted = 0
    blocks: list[tuple[int, int]] = []
    segments: list[tuple[int, int]] = []
    seen_aligned = False
    for count_text, operation in _CIGAR.findall(cigar):
        count = int(count_text)
        if operation == "M":
            raise AssemblyError(
                f"CIGAR {cigar!r} uses M: this reader needs =/X operations (minimap2 --eqx) to "
                "count mismatches"
            )
        if operation == "S":
            if seen_aligned:
                right_clip = count
            else:
                left_clip = count
        elif operation == "H":
            continue
        elif operation in "=X":
            seen_aligned = True
            if operation == "=":
                matched += count
            else:
                mismatched += count
            blocks.append((cursor, cursor + count))
            cursor += count
        elif operation == "I":
            inserted += count
            if count >= _LONG_INDEL:
                segments.append((segment_start, cursor))
                segment_start = cursor
        elif operation in _REFERENCE_CONSUMING:
            deleted += count
            if count >= _LONG_INDEL:
                segments.append((segment_start, cursor))
                segment_start = cursor + count
            cursor += count
    if not blocks:
        raise AssemblyError(f"CIGAR {cigar!r} has no aligned bases")
    segments.append((segment_start, cursor))
    return AlignedRead(
        start=position,
        end=cursor,
        left_clip=left_clip,
        right_clip=right_clip,
        matched=matched,
        mismatched=mismatched,
        inserted=inserted,
        deleted=deleted,
        blocks=tuple(blocks),
        segments=tuple(segments),
    )


def parse_sam(lines: Iterable[str]) -> Iterator[AlignedRead]:
    """Primary alignments of a SAM stream. Unmapped, secondary and supplementary lines are dropped.

    A read counts once: an inverted-repeat read has two equally good placements, and reading one of
    them spreads the repeat's reads over its two copies, which is what leaves the depth level
    across the circle.
    """
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.split("\t", 6)
        if len(fields) < 6:
            raise AssemblyError(f"not a SAM alignment line: {line[:80]!r}")
        if int(fields[1]) & _SKIPPED_FLAGS or fields[5] == "*":
            continue
        yield _parse_cigar(int(fields[3]) - 1, fields[5])


@dataclass(frozen=True, slots=True)
class SupportProfile:
    """What the mapped reads say at every position of one circular sequence.

    ``spanning[b]`` counts reads with at least ``flank`` aligned bases of reference on both sides
    of boundary ``b``, with no long indel between them. ``clipped[b]`` counts reads whose alignment
    stops there with at least ``min_clip`` read bases left over. ``depth[i]`` counts aligned bases
    at base ``i``. All three are indexed on the circle, ``[0, length)``.
    """

    length: int
    flank: int
    min_clip: int
    spanning: tuple[int, ...]
    clipped: tuple[int, ...]
    depth: tuple[int, ...]
    reads: int
    matched: int
    mismatched: int
    inserted: int
    deleted: int

    @property
    def mean_depth(self) -> float:
        """Mean aligned depth, in reads per base."""
        return sum(self.depth) / self.length

    @property
    def sd_depth(self) -> float:
        """Population standard deviation of the aligned depth, in reads per base."""
        mean = self.mean_depth
        return (sum((value - mean) ** 2 for value in self.depth) / self.length) ** 0.5

    @property
    def error_rate(self) -> float | None:
        """(mismatched + inserted + deleted) bases per matched base; None with no matched base."""
        if not self.matched:
            return None
        return (self.mismatched + self.inserted + self.deleted) / self.matched


def _accumulate(deltas: Sequence[int]) -> list[int]:
    running = 0
    values: list[int] = []
    for delta in deltas:
        running += delta
        values.append(running)
    return values


def _fold(values: Sequence[int], length: int) -> tuple[int, ...]:
    """Extended coordinates back onto the circle: position ``i`` adds to ``i % length``."""
    folded = [0] * length
    for index, value in enumerate(values):
        folded[index % length] += value
    return tuple(folded)


def build_profile(
    alignments: Iterable[AlignedRead],
    *,
    length: int,
    extension: int,
    flank: int,
    min_clip: int,
) -> SupportProfile:
    """Accumulate alignments into the circle's profile, in time linear in reads plus positions."""
    if length < 1:
        raise AssemblyError(f"a profile needs a sequence, got length {length}")
    if flank < 1 or min_clip < 1:
        raise AssemblyError(f"flank and min_clip must be positive, got {flank} and {min_clip}")
    size = length + extension + 2
    spanning = [0] * size
    depth = [0] * size
    clipped = [0] * size
    reads = matched = mismatched = inserted = deleted = 0
    for read in alignments:
        if read.end > length + extension:
            raise AssemblyError(
                f"an alignment ends at {read.end}, past the extended sequence "
                f"({length} + {extension})"
            )
        reads += 1
        matched += read.matched
        mismatched += read.mismatched
        inserted += read.inserted
        deleted += read.deleted
        for start, end in read.segments:
            if start + flank <= end - flank:
                spanning[start + flank] += 1
                spanning[end - flank + 1] -= 1
        for start, end in read.blocks:
            depth[start] += 1
            depth[end] -= 1
        if read.left_clip >= min_clip:
            clipped[read.start] += 1
        if read.right_clip >= min_clip:
            clipped[read.end] += 1
    return SupportProfile(
        length=length,
        flank=flank,
        min_clip=min_clip,
        spanning=_fold(_accumulate(spanning)[: length + extension + 1], length),
        clipped=_fold(clipped[: length + extension + 1], length),
        depth=_fold(_accumulate(depth)[: length + extension], length),
        reads=reads,
        matched=matched,
        mismatched=mismatched,
        inserted=inserted,
        deleted=deleted,
    )


@dataclass(frozen=True, slots=True)
class SupportRule:
    """What counts as support: reads that actually cross a boundary, in reads.

    Only crossing is evidence for a join. Reads that stop at a boundary say little: on correct
    sequence half the reads can stop at a motif the sequencer reads badly, and at a wrong join on an
    inverted repeat they stop on one side only, like the rest. They are reported and do not decide.
    """

    min_spanning: int
    clip_window: int  # bp either side summed for the report, as a clip lands around a break


@dataclass(frozen=True, slots=True)
class BoundaryCall:
    """One named boundary, read from the profile."""

    name: str
    boundary: int  # 0-based, between bases boundary - 1 and boundary
    spanning: int
    clipped: int
    fraction: float | None  # spanning / (spanning + clipped), reported; None when no read speaks
    depth_before: float  # mean depth over the window before the boundary
    depth_after: float
    call: str


def _window_mean(values: Sequence[int], start: int, width: int) -> float:
    return sum(values[(start + offset) % len(values)] for offset in range(width)) / width


def _clips_around(profile: SupportProfile, point: int, window: int) -> int:
    return sum(
        profile.clipped[(point + offset) % profile.length] for offset in range(-window, window + 1)
    )


def call_boundary(
    profile: SupportProfile,
    rule: SupportRule,
    *,
    name: str,
    boundary: int,
    depth_window: int = 100,
) -> BoundaryCall:
    """Supported, or undetermined: absence of evidence is never reported as a verified join."""
    point = boundary % profile.length
    spanning = profile.spanning[point]
    clipped = _clips_around(profile, point, rule.clip_window)
    said = spanning + clipped
    fraction = spanning / said if said else None
    supported = spanning >= rule.min_spanning
    width = min(depth_window, profile.length)
    return BoundaryCall(
        name=name,
        boundary=point,
        spanning=spanning,
        clipped=clipped,
        fraction=fraction,
        depth_before=_window_mean(profile.depth, point - width, width),
        depth_after=_window_mean(profile.depth, point, width),
        call=SUPPORT_SUPPORTED if supported else SUPPORT_UNDETERMINED,
    )


@dataclass(frozen=True, slots=True)
class WeakRegion:
    """A run of boundaries too few reads cross; ``start`` after ``end`` wraps the origin."""

    start: int  # first weak boundary, 0-based
    end: int  # last weak boundary, 0-based, inclusive
    length: int
    lowest: int  # fewest spanning reads inside the run
    most_clipped: int  # most reads stopping within the clip window, inside the run


def weak_regions(profile: SupportProfile, rule: SupportRule) -> tuple[WeakRegion, ...]:
    """Every run of the circle too few reads cross, which catches a join nobody asked about."""
    length = profile.length
    clips = [_clips_around(profile, point, rule.clip_window) for point in range(length)]
    weak = [spanning < rule.min_spanning for spanning in profile.spanning]

    def region(start: int, run: int) -> WeakRegion:
        points = [(start + step) % length for step in range(run)]
        return WeakRegion(
            start=start,
            end=points[-1],
            length=run,
            lowest=min(profile.spanning[point] for point in points),
            most_clipped=max(clips[point] for point in points),
        )

    if all(weak):
        return (region(0, length),)
    if not any(weak):
        return ()
    # Start the walk at a boundary that is not weak, so a run across the origin stays one run.
    origin = weak.index(False)
    regions: list[WeakRegion] = []
    run_start: int | None = None
    for offset in range(1, length + 1):
        point = (origin + offset) % length
        if weak[point]:
            if run_start is None:
                run_start = offset
            continue
        if run_start is not None:
            regions.append(region((origin + run_start) % length, offset - run_start))
            run_start = None
    return tuple(sorted(regions, key=lambda item: item.start))


def ir_junctions(
    first: Interval, second: Interval, *, length: int, frame_offset: int = 0
) -> tuple[tuple[str, int], ...]:
    """The four junctions of a quadripartite circle: ``(name, boundary)``, walking the circle.

    ``first`` and ``second`` are the two copies in ascending order, in a frame rotated by
    ``frame_offset`` from the sequence the boundaries are returned on. The longer single-copy arc
    is the LSC. Named by the convention LSC → IRb → SSC → IRa → LSC: JLB, JSB, JSA, JLA.
    """
    inner = second.start - first.end
    outer = length - second.end + first.start
    if inner > outer:
        # The first copy follows the SSC: SSC → first → LSC → second → SSC.
        pairs = (
            ("JSA", first.start),
            ("JLA", first.end),
            ("JLB", second.start),
            ("JSB", second.end),
        )
    else:
        pairs = (
            ("JLB", first.start),
            ("JSB", first.end),
            ("JSA", second.start),
            ("JLA", second.end),
        )
    return tuple((name, (point + frame_offset) % length) for name, point in pairs)
