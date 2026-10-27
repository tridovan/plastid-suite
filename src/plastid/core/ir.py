"""IR detection and the mandatory gate.

The gate is the project's most important safety net: naive IR detection was silently
wrong on 2/9 published plastomes, and every position-based analysis
below it inherits the error. Rules:

* identity is always measured on the **full copies**, never on the exact-match core
  (a perfect 22,998 bp core sat inside copies that are only 92.62 % identical);
* the two copies are related by a mirror: position ``i`` in copy A pairs with
  ``C - i`` in copy B (reverse-complemented), with ``C = a1 + b2 - 1`` — this is the
  pairing constant measured;
* **that constant assumes the two copies are the same length.** Annotations do not always
  give them so (`IR_COPY_LENGTH_TOLERANCE`), and when they differ the constant pairs the
  wrong bases and an exact repeat measures as unrelated — H-25. The difference is therefore
  taken off the longer copy at whichever end it sits, and both readings are measured;
  :attr:`IRResult.excess_at` records which one the number came from;
* a failing gate raises :class:`IRGateError`; it is never downgraded to a warning.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace

from Bio.Align import PairwiseAligner

from plastid.core.model import Genome, Interval
from plastid.core.sequence import complement_base, reverse_complement

__all__ = [
    "IR_ALIGN_GAP_BP",
    "IR_ALIGN_KMER",
    "IR_ALIGN_MAX_DIAGONAL",
    "IR_COPY_LENGTH_TOLERANCE",
    "IR_MIN_EVIDENCE_FRACTION",
    "IR_MIN_IDENTITY",
    "IR_MIN_LENGTH",
    "AlignmentStats",
    "CopyMapping",
    "IRDetectionError",
    "IRGateError",
    "IRResult",
    "align_copies",
    "detect_ir",
    "enforce_ir_gate",
    "ir_de_novo",
    "ir_from_annotation",
    "measure_ir",
    "mirror_constant_for",
    "pair_identity",
]

IR_MIN_LENGTH = 5_000
IR_MIN_IDENTITY = 0.999
# Kept for the old annotation-tolerance callers; candidate selection no longer uses it
# (D11-5: unequal copies are legal and are measured, not discarded).
IR_COPY_LENGTH_TOLERANCE = 5

# The alignment (D11-3/D11-4): exact 20-mers seed anchors, a gap larger than 200 bp ends the
# anchor chain, and only the small segments between anchors are aligned globally.
IR_ALIGN_KMER = 20
IR_ALIGN_MAX_DIAGONAL = 200
IR_ALIGN_GAP_BP = 200

# Evidence policy: identity alone can be carried by a tiny known fragment if unknown columns are
# dropped. Unknown columns stay in the denominator, and at least half of the shorter copy must be
# known A/C/G/T on both sides (with a 20-base floor for small constructed cases).
IR_MIN_EVIDENCE_FRACTION = 0.5
_MIN_EVIDENCE_FLOOR = 20

_KMER = IR_ALIGN_KMER  # the validated probe used K=20
_DIAGONALS_CONSIDERED = 12
_DENOVO_LENGTH_WINDOW = (10_000, 35_000)  # validated validity window
_DENOVO_MIN_IDENTITY = 0.95  # validated validity threshold
# A pairwise global alignment larger than this is not attempted: anchors already cover an IR at
# 99 % identity, so a segment this large means the pair diverged beyond the gate. The segment is
# then classified without gaps, which is an upper bound and cannot manufacture a PASS.
_MAX_ALIGN_CELLS = 4_000_000


class IRDetectionError(Exception):
    """No valid IR pair could be determined — position-based analysis must not run."""


class IRGateError(Exception):
    """The two IR copies are less than IR_MIN_IDENTITY identical — stop, do not analyse."""


@dataclass(frozen=True, slots=True)
class IRResult:
    """The IR pair plus the evidence used to trust (or reject) it."""

    ira: Interval
    irb: Interval
    identity: float
    mismatches: int
    compared_length: int
    exact_core: int
    mirror_constant: int
    source: str  # "annotation" | "denovo_diag" | "manual"
    frame_offset: int = 0  # coordinates are in the frame rotated by this much
    denovo_length: int | None = None  # independent de-novo cross-check, when it ran
    # Which end of the longer copy carried the length difference: "none" when the copies are
    # equal, "outer" (the JLB/JLA side) or "inner" (the JSB/JSA side), and "both" when terminal
    # gaps sit at both ends (mapping between the copies is then not supported; M11.11).
    excess_at: str = "none"
    # End-to-end alignment evidence (D11-3). Every column is counted once:
    # matched + mismatched + ambiguous + gap = compared_length.
    matched: int = 0
    known_columns: int = -1  # -1: caller did not record it; derived below from the old fields
    ambiguous_columns: int = 0
    internal_gaps: int = 0
    terminal_gaps: int = 0
    # The annotated copy lengths when the record had exactly two repeat_region copies, even when
    # the sequence's repeat was used instead. Disagreement stays visible (D11-2).
    annotated_lengths: tuple[int, ...] = ()
    annotation_note: str | None = None
    method: str = "mirror"  # "anchored" | "ungapped_fallback" | "mirror"
    # Alignment-derived mapping between the copies, in this result's frame (M11.11). None when the
    # pair was never aligned (a hand-built result) or when the copies were not contiguous in the
    # frame that was measured; a consumer that needs the mapping then refuses by name.
    mapping: CopyMapping | None = None

    def __post_init__(self) -> None:
        if self.known_columns < 0:
            # Compatibility for callers that build an IRResult directly (tests, old manifests):
            # compare_length was the compared copy length, so treat its known columns as such.
            object.__setattr__(
                self, "known_columns", max(0, self.compared_length - self.mismatches)
            )

    @property
    def length(self) -> int:
        return min(self.ira.length, self.irb.length)

    @property
    def minimum_known(self) -> int:
        """The evidence floor for this pair: half the shorter copy, with a floor for tiny cases."""
        return max(_MIN_EVIDENCE_FLOOR, int(self.length * IR_MIN_EVIDENCE_FRACTION))

    @property
    def gate_passed(self) -> bool:
        return self.identity >= IR_MIN_IDENTITY and self.known_columns >= self.minimum_known

    @property
    def denovo_delta(self) -> int | None:
        """Difference between this IR and the de-novo cross-check (None when it did not run)."""
        if self.denovo_length is None:
            return None
        return self.length - self.denovo_length


def pair_identity(sequence: str, interval: Interval, constant: int) -> tuple[float, int, int, int]:
    """Mirror-compare one copy: ``i`` vs ``constant - i`` (reverse-complemented).

    Returns ``(identity, mismatches, compared, exact_core)``. A partner index outside the
    sequence counts as a mismatch — an out-of-range pairing constant must fail, not crash.
    """
    length = len(sequence)
    compared = interval.length
    mismatches = 0
    exact_core = 0
    core_still_open = True
    for offset in range(compared):
        partner = constant - interval.start - offset
        if 0 <= partner < length and sequence[interval.start + offset] == complement_base(
            sequence[partner]
        ):
            if core_still_open:
                exact_core += 1
        else:
            mismatches += 1
            core_still_open = False
    identity = (compared - mismatches) / compared if compared else 0.0
    return identity, mismatches, compared, exact_core


@dataclass(frozen=True, slots=True)
class AlignmentStats:
    """One end-to-end comparison of two copies, every column classified once (D11-3)."""

    identity: float
    matched: int
    mismatches: int
    ambiguous: int
    internal_gaps: int
    terminal_gaps: int
    columns: int
    known: int
    exact_core: int
    excess_at: str
    method: str


@dataclass(frozen=True, slots=True)
class CopyMapping:
    """Position mapping between the two IR copies through the D11-3 alignment (M11.11).

    The copies are legal when they differ in length (D11-5), but every consumer that mirrored a
    position between them assumed equal lengths. This object carries the alignment's own answer:
    for every base of copy A it records where the cursor in copy B's reverse complement stood
    before and after that column, and the same in the other direction. A gap column maps its base
    to the insertion point between two target bases, which is what a boundary on a one-base indel
    is — not a guess and not an off-by-the-difference mirror.

    Coordinates are copy-relative: ``first`` and ``second`` are the two copies in the frame the
    alignment measured (``IRResult.frame_offset`` says how to get back to the record), and
    ``second``-relative positions are on its **reverse complement**, because that is the strand
    the alignment compares. Converting one to GenBank coordinates is the caller's job and happens
    in ``plastid_annotate`` through the two copies' spans.
    """

    first: Interval
    second: Interval
    # b_before[i] / b_after[i]: cursor in B's reverse complement before/after consuming A[i].
    b_before: tuple[int, ...]
    b_after: tuple[int, ...]
    a_before: tuple[int, ...]
    a_after: tuple[int, ...]
    method: str

    def _map(
        self,
        start: int,
        end: int,
        before: tuple[int, ...],
        after: tuple[int, ...],
        length: int,
        source: str,
    ) -> tuple[int, int] | None:
        if not 0 <= start <= end <= length:
            raise ValueError(
                f"{source} interval [{start}, {end}) is outside a copy of {length} bp — a "
                f"mapping is only defined for bases the alignment compared"
            )
        if start == end:
            cursor = before[start]
            return cursor, cursor
        first, last = before[start], after[end - 1]
        if last <= first:
            # The whole interval aligned to gaps in the other copy: there is no counterpart to
            # write, and writing a zero-length feature would be a claim the sequence does not
            # support.
            return None
        return first, last

    def first_to_second(self, start: int, end: int) -> tuple[int, int] | None:
        """Map ``[start, end)`` of copy A to reverse-complement coordinates of copy B."""
        return self._map(start, end, self.b_before, self.b_after, self.first.length, "first")

    def second_to_first(self, start: int, end: int) -> tuple[int, int] | None:
        """Map ``[start, end)`` of copy B's reverse complement to copy A coordinates."""
        return self._map(start, end, self.a_before, self.a_after, self.second.length, "second")


_KNOWN = frozenset("ACGT")
# A 20-mer occurring more than this many times is low-complexity, not an anchor. Indexing it
# would make the hit walk quadratic on a homopolymer or a satellite array; dropping it leaves the
# real anchors intact (two IR copies contribute at most two occurrences) and sends an
# all-repetitive pair to the honest ungapped fallback instead of an unbounded loop.
_MAX_KMER_OCCURRENCES = 32


def _index_kmers(sequence: str, kmer: int) -> dict[str, list[int]]:
    """Known k-mers → positions, excluding words too repetitive to be informative anchors."""
    index: dict[str, list[int]] = {}
    for position in range(len(sequence) - kmer + 1):
        word = sequence[position : position + kmer]
        if _KNOWN.issuperset(word):
            index.setdefault(word, []).append(position)
    return {
        word: positions
        for word, positions in index.items()
        if len(positions) <= _MAX_KMER_OCCURRENCES
    }


def _anchor_chain(a: str, b: str) -> list[tuple[int, int, int]]:
    """Maximal exact runs of ``a`` against ``b``, chained in increasing order.

    ``b`` is already the reverse complement of the second copy, so an anchor at ``(i, j, n)``
    means ``a[i:i+n] == b[j:j+n]``. Hits more than :data:`IR_ALIGN_MAX_DIAGONAL` away from the
    main diagonal are ignored: the copies are the same repeat, not an arbitrary repeat pair.
    """
    kmer = IR_ALIGN_KMER
    index = _index_kmers(b, kmer)

    hits_by_diagonal: dict[int, list[tuple[int, int]]] = {}
    for i in range(len(a) - kmer + 1):
        word = a[i : i + kmer]
        if not _KNOWN.issuperset(word):
            continue
        for j in index.get(word, ()):
            diagonal = j - i
            if abs(diagonal) <= IR_ALIGN_MAX_DIAGONAL:
                hits_by_diagonal.setdefault(diagonal, []).append((i, j))

    anchors: list[tuple[int, int, int]] = []
    for hits in hits_by_diagonal.values():
        hits.sort()
        start_i = start_j = previous_i = previous_j = None
        for i, j in hits:
            if start_i is None:
                start_i, start_j = i, j
            elif (i, j) != (previous_i + 1, previous_j + 1):
                length = previous_i - start_i + kmer
                if length >= kmer:
                    anchors.append((start_i, start_j, length))
                start_i, start_j = i, j
            previous_i, previous_j = i, j
        if start_i is not None:
            length = previous_i - start_i + kmer
            if length >= kmer:
                anchors.append((start_i, start_j, length))

    anchors.sort()
    if not anchors:
        return []
    # Longest increasing chain. Anchors may **overlap** at an indel: the exact run before the
    # indel and the one after it can share a few bases (measured on NC_065154.1, NC_036733.1 and
    # the panel's other unequal pairs). Requiring non-overlap made the chain skip the longer run
    # and then align a large mis-set segment — 86.7 % where a full global alignment gives 99.96 %.
    # Starts must increase; `align_copies` trims the overlap at the next anchor's start.
    n = len(anchors)
    best = [anchor[2] for anchor in anchors]
    previous: list[int] = [-1] * n
    for i, (ai, bi, li) in enumerate(anchors):
        for j in range(i):
            aj, bj, _ = anchors[j]
            if aj < ai and bj < bi and best[j] + li > best[i]:
                best[i] = best[j] + li
                previous[i] = j
    end = max(range(n), key=lambda idx: best[idx])
    chain: list[tuple[int, int, int]] = []
    while end != -1:
        chain.append(anchors[end])
        end = previous[end]
    chain.reverse()
    # Drop a small run that lies wholly inside an earlier run: it adds no coverage and would make
    # the segment between the two real anchors look like a gap.
    kept: list[tuple[int, int, int]] = []
    for anchor in chain:
        if kept:
            last = kept[-1]
            if (
                anchor[0] + anchor[2] <= last[0] + last[2]
                and anchor[1] + anchor[2] <= last[1] + last[2]
            ):
                continue
        kept.append(anchor)
    return kept


def _kind(one: str | None, other: str | None) -> str:
    """One column's class: M (known match), X (known mismatch), G (gap), A (ambiguous)."""
    if one is None or other is None:
        return "G"
    if one in _KNOWN and other in _KNOWN:
        return "M" if one == other else "X"
    return "A"


def _alignment_pairs(sa: str, sb: str) -> tuple[list[tuple[str | None, str | None]], bool]:
    """The alignment's columns as base pairs, and whether every column came from an alignment.

    The pairs carry position, which the counts alone do not: M11.11 maps a gene between unequal
    copies through this path, and a gap column has to know *which* side the gap is on. A segment
    larger than :data:`_MAX_ALIGN_CELLS` is classified without gaps rather than aligned — at this
    point the copies have already diverged far beyond the gate, and an ungapped reading is an
    upper bound that cannot turn a divergent pair into a PASS. That reading is **not** precise,
    so its caller must not build a position mapping from it.
    """
    if not sa:
        return [(None, base) for base in sb], True
    if not sb:
        return [(base, None) for base in sa], True
    if len(sa) * len(sb) > _MAX_ALIGN_CELLS:
        longer, shorter = (sa, sb) if len(sa) >= len(sb) else (sb, sa)
        difference = len(longer) - len(shorter)
        if len(sa) >= len(sb):
            pairs = [(base, None) for base in longer[:difference]] + list(
                zip(sa[difference:], sb, strict=True)
            )
        else:
            pairs = [(None, base) for base in longer[:difference]] + list(
                zip(sa, sb[difference:], strict=True)
            )
        return pairs, False
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    alignment = aligner.align(sa, sb)[0]
    return [
        (None if x == "-" else x, None if y == "-" else y)
        for x, y in zip(alignment[0], alignment[1], strict=True)
    ], True


def _ungapped_pairs(a: str, b: str) -> list[tuple[str | None, str | None]]:
    """The diagonal reading: equal prefix paired, the shorter copy's tail a terminal gap run."""
    return (
        list(zip(a, b, strict=False))
        + [(base, None) for base in a[len(b) :]]
        + [(None, base) for base in b[len(a) :]]
    )


def _stats_from_pairs(
    pairs: list[tuple[str | None, str | None]],
    *,
    method: str,
    exact_core: int = 0,
    known_ends: bool = True,
) -> AlignmentStats:
    """Classify every column once and derive the gate's evidence from the same list.

    ``known_ends=False`` is the diagonal fallback: it placed the length difference at the tail
    because it had to put it somewhere, not because the alignment supports it, so the reported
    ``excess_at`` stays ``"unknown"``.
    """
    kinds = [_kind(one, other) for one, other in pairs]
    matched = kinds.count("M")
    mismatches = kinds.count("X")
    ambiguous = kinds.count("A")
    gaps = kinds.count("G")
    first = next((index for index, kind in enumerate(kinds) if kind != "G"), None)
    last = next((index for index in range(len(kinds) - 1, -1, -1) if kinds[index] != "G"), None)
    if first is None or last is None:
        internal = 0
        terminal = gaps
        excess_at = "unknown"
    else:
        terminal = sum(
            1 for index, kind in enumerate(kinds) if kind == "G" and (index < first or index > last)
        )
        internal = gaps - terminal
        prefix = sum(1 for index in range(first) if kinds[index] == "G")
        suffix = sum(1 for index in range(last + 1, len(kinds)) if kinds[index] == "G")
        if prefix and suffix:
            excess_at = "both"
        elif prefix:
            excess_at = "outer"
        elif suffix:
            excess_at = "inner"
        else:
            excess_at = "none"
        if not known_ends and gaps:
            excess_at = "unknown"

    columns = len(kinds)
    return AlignmentStats(
        identity=matched / columns if columns else 0.0,
        matched=matched,
        mismatches=mismatches,
        ambiguous=ambiguous,
        internal_gaps=internal,
        terminal_gaps=terminal,
        columns=columns,
        known=matched + mismatches,
        exact_core=exact_core,
        excess_at=excess_at,
        method=method,
    )


def _ungapped_stats(a: str, b: str) -> AlignmentStats:
    """No anchor chain: compare on the diagonal and report the length difference as terminal gaps.

    This is not used to accept a repeat — the gate's minimum evidence and identity both reject
    these pairs — it is the honest fallback when no 20-mer chain exists.
    """
    return _stats_from_pairs(_ungapped_pairs(a, b), method="ungapped_fallback", known_ends=False)


def _aligned_pairs(a: str, b: str) -> tuple[list[tuple[str | None, str | None]], str, int]:
    """The anchored path, its method, and the longest exact anchor (0 when there is no chain).

    ``anchored_unverified`` means at least one segment was too large to align and was read
    without gaps. Its identity is an upper bound, and its columns are not evidence for placing a
    gene, so ``measure_ir`` builds no mapping from it.
    """
    chain = _anchor_chain(a, b)
    if not chain:
        return _ungapped_pairs(a, b), "ungapped_fallback", 0

    pairs: list[tuple[str | None, str | None]] = []
    precise = True
    a_cursor = b_cursor = 0
    for index, (a_start, b_start, length) in enumerate(chain):
        # Everything between the previous anchor's end and this anchor's start.
        segment, segment_precise = _alignment_pairs(a[a_cursor:a_start], b[b_cursor:b_start])
        pairs.extend(segment)
        precise = precise and segment_precise
        next_a = chain[index + 1][0] if index + 1 < len(chain) else len(a)
        next_b = chain[index + 1][1] if index + 1 < len(chain) else len(b)
        usable_a = min(a_start + length, next_a) - a_start
        usable_b = min(b_start + length, next_b) - b_start
        # The overlapping runs are exact in both copies; the shorter side's length is real match
        # columns, and a length difference at the overlap is a gap (an indel at the boundary).
        exact = min(usable_a, usable_b)
        pairs.extend((a[a_start + offset], b[b_start + offset]) for offset in range(exact))
        if usable_a > usable_b:
            pairs.extend(
                (a[a_start + exact + offset], None) for offset in range(usable_a - usable_b)
            )
        elif usable_b > usable_a:
            pairs.extend(
                (None, b[b_start + exact + offset]) for offset in range(usable_b - usable_a)
            )
        a_cursor, b_cursor = a_start + usable_a, b_start + usable_b
    segment, segment_precise = _alignment_pairs(a[a_cursor:], b[b_cursor:])
    pairs.extend(segment)
    precise = precise and segment_precise
    return (
        pairs,
        "anchored" if precise else "anchored_unverified",
        max(length for _, _, length in chain),
    )


def align_copies(a: str, b: str) -> AlignmentStats:
    """Compare copy A with the reverse complement of copy B, end to end (D11-3).

    Exact 20-mers seed anchors; the segments between consecutive anchors and the prefix/suffix
    are aligned globally. One SNP leaves two anchors on the same diagonal and one 1-base segment
    between them; one indel shifts the following diagonal, which the chain follows. Columns are
    classified once, and unknown bases are **ambiguous**, never matches.
    """
    pairs, method, exact_core = _aligned_pairs(a, b)
    return _stats_from_pairs(pairs, method=method, exact_core=exact_core)


def _trimmed_pair(
    first: Interval, second: Interval, excess: int, where: str
) -> tuple[Interval, Interval]:
    """The pair with ``abs(excess)`` bases taken off the **longer** copy at the named end.

    ``where`` is ``"outer"`` — the first copy's start or the second copy's end, the JLB/JLA side —
    or ``"inner"``, the first copy's end or the second copy's start, the JSB/JSA side. Trimming the
    longer copy to the shorter's length is what makes the mirror constant valid again: the
    formula assumes the two copies are the same length, and the annotation does not say which end
    the difference sits at.
    """
    if excess == 0:
        return first, second
    if where not in ("outer", "inner"):
        # An IRResult carries "none" only when the copies are equal; saying nothing about where a
        # real difference sits would make the constant depend on a default nobody chose.
        raise ValueError(
            f"the two copies differ by {excess} base(s), so the result must say where the "
            f"difference sits ('outer' or 'inner'); got {where!r}"
        )
    size = abs(excess)
    if excess > 0:  # the first copy is the longer one
        trimmed = (
            Interval(first.start + size, first.end, first.strand)
            if where == "outer"
            else Interval(first.start, first.end - size, first.strand)
        )
        return trimmed, second
    trimmed = (
        Interval(second.start, second.end - size, second.strand)
        if where == "outer"
        else Interval(second.start + size, second.end, second.strand)
    )
    return first, trimmed


def mirror_constant_for(first: Interval, second: Interval, excess_at: str) -> int:
    """The pairing constant for a pair, honouring **where** a length difference sits.

    The one derivation, shared with :func:`plastid.core.normalize.normalize`, which recomputes the
    constant on the rotated intervals. Two derivations would drift, and the drift would be invisible
    until a figure or a table came out wrong on the frame the modules actually use.
    """
    trimmed_first, trimmed_second = _trimmed_pair(
        first, second, first.length - second.length, excess_at
    )
    return trimmed_first.start + trimmed_second.end - 1


def _aligned(
    sequence: str,
    first: Interval,
    second: Interval,
    compared_length: int,
    excess: int,
    where: str,
) -> tuple[float, int, int, int, int, str]:
    """One alignment of the pair: the longer copy trimmed at ``where``, then the mirror constant."""
    trimmed_first, _trimmed_second = _trimmed_pair(first, second, excess, where)
    span = Interval(
        trimmed_first.start, trimmed_first.start + compared_length, trimmed_first.strand
    )
    constant = mirror_constant_for(first, second, where)
    identity, mismatches, compared, exact_core = pair_identity(sequence, span, constant)
    return identity, mismatches, compared, exact_core, constant, where


def _copy_mapping(
    first: Interval,
    second: Interval,
    pairs: list[tuple[str | None, str | None]],
    method: str,
) -> CopyMapping:
    """Turn the alignment's columns into the two-way cursor tables (M11.11).

    ``pairs`` covers copy A (``first``) and copy B's reverse complement (``second``). Each array
    holds, for every base of one copy, where the other copy's cursor stood before and after that
    base's column. A base aligned to a gap is missing from the other copy; its two cursors are
    equal, and an interval made only of such bases has no counterpart.
    """
    b_before = [0] * (first.length + 1)
    b_after = [0] * (first.length + 1)
    a_before = [0] * (second.length + 1)
    a_after = [0] * (second.length + 1)
    one_cursor = other_cursor = 0
    for one, other in pairs:
        before_one = one_cursor
        if one is not None:
            b_before[one_cursor] = other_cursor
            b_after[one_cursor] = other_cursor + (1 if other is not None else 0)
            one_cursor += 1
        if other is not None:
            a_before[other_cursor] = before_one
            a_after[other_cursor] = one_cursor
            other_cursor += 1
    b_before[first.length] = other_cursor
    a_before[second.length] = one_cursor
    if one_cursor != first.length or other_cursor != second.length:
        raise ValueError(
            f"the alignment covers {one_cursor}/{first.length} and {other_cursor}/"
            f"{second.length} bases — a mapping built from a partial alignment would place a "
            f"gene from a coordinate the comparison never reached"
        )
    return CopyMapping(
        first=first,
        second=second,
        b_before=tuple(b_before),
        b_after=tuple(b_after),
        a_before=tuple(a_before),
        a_after=tuple(a_after),
        method=method,
    )


def measure_ir(
    sequence: str,
    ira: Interval,
    irb: Interval,
    *,
    source: str,
    frame_offset: int = 0,
    annotated_lengths: tuple[int, ...] = (),
    annotation_note: str | None = None,
) -> IRResult:
    """Measure the two copies end to end — the one identity used by the gate (D11-3).

    Copy A is compared with the reverse complement of copy B. Exact 20-mers seed anchors; the
    segments between them, and the prefix/suffix, are aligned globally, so a SNP is one mismatch
    rather than a truncated repeat and a one-base indel is one gap rather than a shifted
    comparison. Unknown bases are ambiguous columns, never matches. ``compared_length`` is the
    number of alignment columns (the identity denominator).

    ``mirror_constant`` is retained for the position mapping other modules still use. A pair with
    terminal gaps at both ends has no single mirror constant; the value is 0 and ``excess_at`` says
    ``"both"``. The **mapping** carried in the result is the alignment's own path, so it survives
    any of those shapes (M11.11); it is None when no anchor chain existed, because a diagonal
    reading is not evidence for placing a gene between copies of different length.
    """
    first, second = sorted((ira, irb), key=lambda interval: interval.start)
    if first.end > second.start:
        raise ValueError(f"IR copies must be disjoint and ordered: {first} overlaps {second}")
    if first.length == 0 or second.length == 0:
        raise ValueError("IR copies must not be empty")
    a = sequence[first.start : first.end]
    b = reverse_complement(sequence[second.start : second.end])
    pairs, method, exact_core = _aligned_pairs(a, b)
    stats = _stats_from_pairs(pairs, method=method, exact_core=exact_core)
    if stats.excess_at == "none":
        constant = first.start + second.end - 1
    elif stats.excess_at in ("outer", "inner"):
        constant = mirror_constant_for(first, second, stats.excess_at)
    else:
        constant = 0
    mapping = _copy_mapping(first, second, pairs, method) if method == "anchored" else None
    return IRResult(
        ira=ira,
        irb=irb,
        identity=stats.identity,
        mismatches=stats.mismatches,
        compared_length=stats.columns,
        exact_core=stats.exact_core,
        mirror_constant=constant,
        source=source,
        frame_offset=frame_offset,
        excess_at=stats.excess_at,
        matched=stats.matched,
        known_columns=stats.known,
        ambiguous_columns=stats.ambiguous,
        internal_gaps=stats.internal_gaps,
        terminal_gaps=stats.terminal_gaps,
        annotated_lengths=annotated_lengths,
        annotation_note=annotation_note,
        method=stats.method,
        mapping=mapping,
    )


def _arcs_overlap(first: tuple[int, int], second: tuple[int, int], length: int) -> bool:
    """Whether two circular arcs share a base (``(start, length)`` each)."""
    (first_start, first_length), (second_start, second_length) = first, second
    return (second_start - first_start) % length < first_length or (
        first_start - second_start
    ) % length < second_length


def _annotation_candidates(genome: Genome) -> tuple[tuple[int, int], ...]:
    """The record's repeat copies as circular arcs, one per repeat_region feature.

    Fewer than two long copies means the record makes no IR claim (de-novo is the path); exactly
    two is the IR; anything else is refused by name rather than silently re-read (D11-2).
    """
    candidates = [region for region in genome.repeat_regions if region.length >= IR_MIN_LENGTH]
    if not candidates:
        return ()
    if len(candidates) != 2:
        where = [(region.start + 1, region.end, region.length) for region in candidates]
        raise IRDetectionError(
            f"{genome.accession or genome.name}: expected exactly two repeat_region copies of "
            f"at least {IR_MIN_LENGTH} bp, found {len(candidates)} ({where}) — the record's "
            f"repeat annotation cannot be read as an IR"
        )
    arcs: list[tuple[int, int]] = []
    for region in candidates:
        try:
            arcs.append(region.arc(len(genome)))
        except ValueError as error:
            raise IRDetectionError(
                f"{genome.accession or genome.name}: repeat_region "
                f"{region.start + 1}..{region.end} is not one circular copy — {error}"
            ) from None
    assert len(arcs) == 2
    if _arcs_overlap(arcs[0], arcs[1], len(genome)):
        raise IRDetectionError(
            f"{genome.accession or genome.name}: the two repeat_region copies overlap "
            f"({arcs[0]}, {arcs[1]}) — the annotation cannot be read"
        )
    return tuple(arcs)


def ir_from_annotation(genome: Genome) -> IRResult | None:
    """The IR pair from GenBank ``repeat_region`` annotations, read whole (D11-1/D11-2).

    A joined copy is one copy. When a copy crosses the origin the pair is measured in a frame
    rotated by the length of the piece at the start, so both copies are contiguous;
    ``frame_offset`` records that rotation and ``normalize`` consumes it. Length differences are
    legal and measured end to end — they are not a reason to drop the annotation.
    """
    arcs = _annotation_candidates(genome)
    if not arcs:
        return None
    length = len(genome)
    shift = 0
    for start, arc_length in arcs:
        if start + arc_length > length:
            shift = start + arc_length - length
            break
    measurement_sequence = (
        genome.sequence[shift:] + genome.sequence[:shift] if shift else genome.sequence
    )
    intervals = sorted(
        (
            Interval((start - shift) % length, (start - shift) % length + arc_length)
            for start, arc_length in arcs
        ),
        key=lambda interval: interval.start,
    )
    if intervals[0].end > intervals[1].start:
        raise IRDetectionError(
            f"{genome.accession or genome.name}: the two copies are not disjoint after "
            f"normalizing the frame — the annotation cannot be read"
        )
    return measure_ir(
        measurement_sequence,
        intervals[0],
        intervals[1],
        source="annotation",
        frame_offset=shift,
        annotated_lengths=tuple(interval.length for interval in intervals),
    )


def _mirror_chains(sequence: str, mirror: str) -> list[tuple[int, int, int]]:
    """Exact runs on the strongest mirror diagonals, chained across small gaps.

    Coordinates are ``(position in the sequence, position in the reverse complement, length)``.
    A chain follows a SNP on the same diagonal and a small indel onto a neighbouring diagonal; it
    stops when the next anchor is further than :data:`IR_ALIGN_GAP_BP` away in either copy.
    """
    length = len(sequence)
    index = _index_kmers(sequence, _KMER)

    diagonals: Counter[int] = Counter()
    hits_by_diagonal: dict[int, list[tuple[int, int]]] = {}
    for position in range(length - _KMER + 1):
        word = mirror[position : position + _KMER]
        if not _KNOWN.issuperset(word):
            continue
        for mate in index.get(word, ()):
            diagonal = mate - position
            # (sequence position, reverse-complement position) — the order the chain and geometry
            # use. `mate` comes from the sequence index, `position` walks the reverse complement.
            hits_by_diagonal.setdefault(diagonal, []).append((mate, position))
            diagonals[diagonal] += 1
    if not diagonals:
        return []

    runs: list[tuple[int, int, int]] = []
    for diagonal, _ in diagonals.most_common(_DIAGONALS_CONSIDERED):
        hits = sorted(hits_by_diagonal[diagonal])
        start = previous = None
        for position, mate in hits:
            if start is None:
                start = previous = (position, mate)
            elif (position, mate) == (previous[0] + 1, previous[1] + 1):
                previous = (position, mate)
            else:
                run_length = previous[0] - start[0] + _KMER
                if run_length >= _KMER:
                    runs.append((start[0], start[1], run_length))
                start = previous = (position, mate)
        if start is not None:
            run_length = previous[0] - start[0] + _KMER
            if run_length >= _KMER:
                runs.append((start[0], start[1], run_length))

    runs.sort()
    if not runs:
        return []
    # Chain runs whose positions increase and whose gaps are small enough to be one repeat. The
    # gap limit is what stops the chain at real divergence instead of extending on a global ratio.
    best = [run[2] for run in runs]
    previous_run: list[int] = [-1] * len(runs)
    for i, (ai, bi, li) in enumerate(runs):
        for j in range(i):
            aj, bj, lj = runs[j]
            a_gap = ai - (aj + lj)
            b_gap = bi - (bj + lj)
            if (
                0 <= a_gap <= IR_ALIGN_GAP_BP
                and 0 <= b_gap <= IR_ALIGN_GAP_BP
                and best[j] + li > best[i]
            ):
                best[i] = best[j] + li
                previous_run[i] = j
    end = max(range(len(runs)), key=lambda index: best[index])
    chain: list[tuple[int, int, int]] = []
    while end != -1:
        chain.append(runs[end])
        end = previous_run[end]
    chain.reverse()
    return chain


def _chain_geometry(
    chain: list[tuple[int, int, int]], sequence_length: int
) -> tuple[Interval, Interval] | None:
    """The two copy intervals implied by a chained set of mirror anchors."""
    if not chain:
        return None
    first_a, first_b, _ = chain[0]
    last_a, last_b, last_length = chain[-1]
    a_end = last_a + last_length
    b_end = last_b + last_length
    ira = Interval(first_a, a_end)
    irb = Interval(sequence_length - b_end, sequence_length - first_b)
    if irb.start < 0 or irb.start < ira.end or ira.end > sequence_length:
        return None
    if irb.end > sequence_length or irb.start - ira.end <= 0:
        return None
    return ira, irb


def ir_de_novo(genome: Genome, *, frame_offset: int = 0) -> IRResult | None:
    """Mismatch-tolerant, circular-aware de-novo detection (D11-4).

    Exact 20-mers on the strongest mirror diagonals seed anchors; anchors are chained through a
    SNP (same diagonal) or a small indel (neighbouring diagonal) until the next anchor is more
    than :data:`IR_ALIGN_GAP_BP` away. The full chained extent is then measured by
    :func:`measure_ir`, never the exact core — one SNP cannot truncate the repeat, and extending
    into single-copy sequence lowers the measured identity below the validity threshold.

    The scan runs on the linear record and on a half-rotated record; the longer valid candidate
    wins, so an IR copy crossing the origin is not truncated (A8).
    """
    candidates = [_scan_linear(genome, frame_offset)]
    if genome.topology == "circular":
        shift = len(genome) // 2
        rotated = _rotate_genome(genome, shift)
        candidates.append(_scan_linear(rotated, (frame_offset + shift) % len(genome)))
    valid = [candidate for candidate in candidates if candidate is not None]
    if not valid:
        return None
    return max(valid, key=lambda ir: (ir.length, -ir.frame_offset))


def _scan_linear(genome: Genome, frame_offset: int) -> IRResult | None:
    sequence = genome.sequence
    mirror = reverse_complement(sequence)
    chain = _mirror_chains(sequence, mirror)
    geometry = _chain_geometry(chain, len(sequence))
    if geometry is None:
        return None
    ira, irb = geometry
    ir = measure_ir(sequence, ira, irb, source="denovo_diag", frame_offset=frame_offset)
    if not (_DENOVO_LENGTH_WINDOW[0] <= ir.length <= _DENOVO_LENGTH_WINDOW[1]):
        return None
    if ir.identity < _DENOVO_MIN_IDENTITY:
        return None
    return ir


def _rotate_genome(genome: Genome, shift: int) -> Genome:
    from plastid.core.normalize import rotate_genome

    return rotate_genome(genome, shift)


def detect_ir(genome: Genome, *, source: str = "annotation") -> IRResult:
    """Choose the IR source explicitly; never substitute one path for the other silently (D11-2).

    ``source="annotation"`` uses the record's two copies, or the sequence's repeat when the record
    has no IR annotation. A record that annotates one, three or overlapping copies is refused by
    name, not re-read. ``source="sequence"`` uses the sequence's repeat even when the annotation
    exists, and records the annotation it did not use. A gate failure on the annotation path is
    reported by :func:`enforce_ir_gate` with the way out: ``--set ir_source=sequence``.
    """
    if source not in ("annotation", "sequence"):
        raise ValueError(f"ir_source must be 'annotation' or 'sequence', got {source!r}")
    annotation_error: IRDetectionError | None = None
    annotated: IRResult | None = None
    try:
        annotated = ir_from_annotation(genome)
    except IRDetectionError as error:
        annotation_error = error
    detected = ir_de_novo(genome)

    if source == "sequence":
        if detected is None:
            detail = f" ({annotation_error})" if annotation_error else ""
            raise IRDetectionError(
                f"{genome.accession or genome.name}: no valid sequence repeat found "
                f"(de-novo scan failed){detail}"
            )
        lengths = annotated.annotated_lengths if annotated is not None else ()
        note = str(annotation_error) if annotation_error is not None else None
        return replace(
            detected,
            annotated_lengths=lengths,
            annotation_note=note,
        )

    if annotation_error is not None:
        raise annotation_error
    if annotated is not None:
        if detected is not None:
            annotated = replace(annotated, denovo_length=detected.length)
        return annotated
    if detected is None:
        raise IRDetectionError(
            f"{genome.accession or genome.name}: no valid IR pair found "
            f"(annotation absent; de-novo scan failed)"
        )
    return detected


def enforce_ir_gate(ir: IRResult, *, min_identity: float = IR_MIN_IDENTITY) -> IRResult:
    """Raise :class:`IRGateError` unless the two copies are identical enough, with evidence.

    Identity alone is not enough: at least half of the shorter copy must be known A/C/G/T on both
    sides, so a handful of matching bases surrounded by ``N`` cannot carry a PASS (K-04).
    """
    minimum_known = ir.minimum_known
    if ir.identity < min_identity or ir.known_columns < minimum_known:
        reasons = []
        if ir.identity < min_identity:
            reasons.append(
                f"{ir.identity * 100:.2f} % identity (need >= {min_identity * 100:.1f} %)"
            )
        if ir.known_columns < minimum_known:
            reasons.append(
                f"{ir.known_columns} known columns (need >= {minimum_known}); "
                f"{ir.ambiguous_columns} ambiguous"
            )
        message = (
            f"IR gate failed: the two copies have {' and '.join(reasons)} — "
            f"position-based analysis must not run"
        )
        if ir.source == "annotation":
            message += (
                "; re-run with --set ir_source=sequence to use the sequence's repeat, "
                "or correct the record's repeat_region"
            )
        raise IRGateError(message)
    return ir
