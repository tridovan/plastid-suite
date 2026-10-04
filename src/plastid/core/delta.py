"""MUMmer delta files — the alignment behind the %identity-by-coordinate series (P4-01).

A delta file is what ``nucmer`` writes and what ``delta-filter`` rewrites: the coordinates of
every alignment, plus the *gaps* inside each one. It is small, undocumented in the pipeline's own
files, and decodes only one way — so the walk lives here, as pure functions, with the measurement
that fixed it written down.

Two facts about MUMmer 4's format, both measured 04/10/2026 on the golden pair:

* **one ``>`` header per sequence pair, not per alignment.** Every alignment of a pair is a
  ``7-number`` record terminated by a ``0`` line, all under a single header. A reader that
  expects a header per record silently returns one block for a comparison that has seven.
* **the gap numbers.** After the header, each value ``d`` means: consume ``|d| - 1`` *paired*
  columns (one base from each side); then ``d > 0`` advances the target alone (that base has no
  query partner) and ``d < 0`` advances the query alone. The list ends at a ``0``. Verified on
  all twelve raw alignments of the golden pair: every block lands exactly on its own end
  coordinates and ``mismatches + unpaired bases`` equals the header's own error count.

Coordinates here are **0-based half-open**, like everywhere else in the toolkit; the file's are
1-based inclusive. A query span is stored in the query's forward strand, and the segments carry
*oriented* query coordinates so that ``target_start + k`` pairs with ``query_start + k`` — use
:meth:`DeltaBlock.query_forward` to map a segment position back to the query's own coordinates.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from plastid.core.sequence import complement_base

__all__ = [
    "DeltaBlock",
    "DeltaSegment",
    "WindowIdentity",
    "aligned_totals",
    "iter_aligned_pairs",
    "parse_delta",
    "window_identity",
]

_PROGRAM_LINES = ("NUCMER", "PROMER")
_HEADER_FIELDS = 7
_DNA = frozenset("ACGT")


@dataclass(frozen=True, slots=True)
class DeltaSegment:
    """One gapless stretch: the two sides advance together, ``k`` for ``k``.

    Both spans are 0-based half-open. The query span is *oriented* — measured along the
    alignment — so for a reverse block its coordinates are not the query's own; see
    :meth:`DeltaBlock.query_forward`.
    """

    target_start: int
    target_end: int
    query_start: int
    query_end: int

    @property
    def length(self) -> int:
        """Paired bases: one target base and one query base each."""
        return self.target_end - self.target_start


@dataclass(frozen=True, slots=True)
class DeltaBlock:
    """One alignment: where it sits on both sequences, how wrong it is, and its gapless runs."""

    target_name: str
    query_name: str
    target_length: int
    query_length: int
    target_start: int  # 0-based half-open, target forward strand
    target_end: int
    query_start: int  # 0-based half-open, query forward strand (start < end always)
    query_end: int
    errors: int  # the header's own count: mismatches + unpaired bases
    reverse: bool
    segments: tuple[DeltaSegment, ...]

    @property
    def length(self) -> int:
        """Target bases the alignment consumes."""
        return self.target_end - self.target_start

    def query_forward(self, oriented: int) -> int:
        """Map an oriented query coordinate to the query's own 0-based forward coordinate."""
        if self.reverse:
            return self.query_length - 1 - oriented
        return oriented


def parse_delta(text: str) -> tuple[DeltaBlock, ...]:
    """Read a delta file into blocks; an empty comparison is zero blocks, not an error."""
    lines = text.splitlines()
    start = next(
        (index for index, line in enumerate(lines) if line.strip() in _PROGRAM_LINES), None
    )
    if start is None:
        raise ValueError("not a MUMmer delta file: no NUCMER/PROMER line")

    blocks: list[DeltaBlock] = []
    names: tuple[str, str, int, int] | None = None
    header: tuple[int, ...] | None = None
    gaps: list[int] = []
    for number, line in enumerate(lines[start + 1 :], start=start + 2):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            names = _parse_block_header(stripped, number)
            header, gaps = None, []
            continue
        fields = stripped.split()
        if len(fields) == _HEADER_FIELDS:
            if names is None:
                raise ValueError(f"line {number}: alignment record before any '>' header")
            if header is not None:
                raise ValueError(f"line {number}: alignment record without its terminating 0")
            header = tuple(int(field) for field in fields)
            continue
        if len(fields) == 1:
            if header is None:
                raise ValueError(f"line {number}: gap value outside an alignment record")
            value = int(fields[0])
            if value == 0:
                blocks.append(_build(names, header, gaps))
                header, gaps = None, []
            else:
                gaps.append(value)
            continue
        raise ValueError(f"line {number}: cannot read {stripped!r} as a delta record")

    if header is not None:
        raise ValueError("the last alignment record is not terminated by a 0 line")
    return tuple(blocks)


def _parse_block_header(line: str, number: int) -> tuple[str, str, int, int]:
    fields = line[1:].split()
    if len(fields) != 4:
        raise ValueError(
            f"line {number}: expected '>name name length length' with both lengths "
            f"(this reader needs them to orient reverse alignments), got {line!r}"
        )
    try:
        return fields[0], fields[1], int(fields[2]), int(fields[3])
    except ValueError:
        raise ValueError(f"line {number}: sequence lengths are not integers: {line!r}") from None


def _build(
    names: tuple[str, str, int, int], header: tuple[int, ...], gaps: list[int]
) -> DeltaBlock:
    """One record into a block, walking its gap numbers and checking the walk against the header."""
    target_name, query_name, target_length, query_length = names
    start1, end1, start2, end2, errors, _similar, _non_alphas = header
    reverse = start2 > end2
    if not (1 <= start1 <= end1 <= target_length) or not (
        1 <= min(start2, end2) <= max(start2, end2) <= query_length
    ):
        raise ValueError(
            f"alignment {start1}-{end1} vs {start2}-{end2} lies outside the sequences "
            f"({target_length}, {query_length})"
        )

    target_position = start1 - 1
    if reverse:
        query_position = query_length - start2  # oriented: counts along the alignment
        query_end = query_length - end2 + 1
    else:
        query_position = start2 - 1
        query_end = end2
    target_end = end1

    segments: list[DeltaSegment] = []
    run_target, run_query = target_position, query_position

    def flush(target: int, query: int) -> None:
        if target > run_target:
            segments.append(DeltaSegment(run_target, target, run_query, query))

    for value in gaps:
        step = abs(value) - 1
        target_position += step
        query_position += step
        flush(target_position, query_position)
        if value > 0:
            target_position += 1  # a target base with no query partner
        else:
            query_position += 1  # a query base with no target partner
        run_target, run_query = target_position, query_position
    while target_position < target_end and query_position < query_end:
        target_position += 1
        query_position += 1
    flush(target_position, query_position)

    if target_position != target_end or query_position != query_end:
        raise ValueError(
            f"the gap list does not account for the alignment "
            f"{start1}-{end1} vs {start2}-{end2}: the walk ends at "
            f"({target_position + 1}, {query_position + 1})"
        )
    if reverse:
        # the alignment covers forward positions end2..start2 inclusive
        query_start, query_end_forward = end2 - 1, start2
    else:
        query_start, query_end_forward = start2 - 1, end2
    return DeltaBlock(
        target_name=target_name,
        query_name=query_name,
        target_length=target_length,
        query_length=query_length,
        target_start=start1 - 1,
        target_end=target_end,
        query_start=query_start,
        query_end=query_end_forward,
        errors=errors,
        reverse=reverse,
        segments=tuple(segments),
    )


def iter_aligned_pairs(
    block: DeltaBlock, *, target: str, query: str
) -> Iterator[tuple[int, int, str, str]]:
    """``(target position, query position, target base, query base)`` for every paired column.

    Unpaired bases — the gaps the delta records — are not yielded: they are not comparisons.
    The query base is complemented for a reverse block, so a caller compares the two strings
    directly without knowing the orientation. Positions are 0-based in each sequence's own
    forward strand.
    """
    if len(target) != block.target_length or len(query) != block.query_length:
        raise ValueError(
            f"sequences are {len(target)}/{len(query)} bases but the delta declares "
            f"{block.target_length}/{block.query_length}"
        )
    for segment in block.segments:
        for offset in range(segment.length):
            target_position = segment.target_start + offset
            query_position = block.query_forward(segment.query_start + offset)
            base = query[query_position]
            yield (
                target_position,
                query_position,
                target[target_position],
                complement_base(base) if block.reverse else base,
            )


@dataclass(frozen=True, slots=True)
class WindowIdentity:
    """% identity over one window of the **target**, and the base pairs it was computed from."""

    start_1b: int
    end_1b: int
    compared_bp: int  # aligned base pairs where both sides are unambiguous (ACGT)
    mismatches: int

    @property
    def identity_pct(self) -> float | None:
        """``None`` when the window holds no aligned pair — a value must never be invented."""
        if self.compared_bp == 0:
            return None
        return 100.0 * (1.0 - self.mismatches / self.compared_bp)


def window_identity(
    blocks: Iterable[DeltaBlock],
    *,
    target: str,
    query: str,
    window: int,
    step: int,
) -> tuple[WindowIdentity, ...]:
    """% identity per sliding window of the target, from the blocks' **paired** columns.

    Only base pairs are compared: a base the other sequence does not carry is a gap, and a
    column with an ambiguous base on either side is dropped — the same complete-deletion rule
    :mod:`plastid.core.align` uses. This is why the series and the block table's own ``%IDY``
    differ slightly: the tool counts an indel as an error, the series has no column to count.

    Windows with no aligned pair come back with ``identity_pct`` of ``None``, so an unaligned
    stretch is visible as a hole rather than drawn as a zero.
    """
    if window <= 0 or step <= 0:
        raise ValueError(f"window and step must be positive, got window={window}, step={step}")
    if step > window:
        raise ValueError(f"step {step} exceeds window {window} — columns would be skipped")
    length = len(target)
    compared, mismatches = _per_position(blocks, target=target, query=query)
    # Prefix sums so each window costs the same whatever its length.
    compared_prefix = [0] * (length + 1)
    mismatch_prefix = [0] * (length + 1)
    for index in range(length):
        compared_prefix[index + 1] = compared_prefix[index] + compared[index]
        mismatch_prefix[index + 1] = mismatch_prefix[index] + mismatches[index]

    windows: list[WindowIdentity] = []
    start = 0
    while start < length:
        end = min(start + window, length)
        windows.append(
            WindowIdentity(
                start_1b=start + 1,
                end_1b=end,
                compared_bp=compared_prefix[end] - compared_prefix[start],
                mismatches=mismatch_prefix[end] - mismatch_prefix[start],
            )
        )
        if end == length:
            break
        start += step
    return tuple(windows)


def aligned_totals(blocks: Iterable[DeltaBlock], *, target: str, query: str) -> tuple[int, int]:
    """``(aligned base pairs, mismatches)`` over the whole target.

    Windows overlap, so summing their counts would double-count; this walks the alignment once.
    """
    compared, mismatches = _per_position(blocks, target=target, query=query)
    return sum(compared), sum(mismatches)


def _per_position(
    blocks: Iterable[DeltaBlock], *, target: str, query: str
) -> tuple[bytearray, bytearray]:
    """Per target position: was it an aligned unambiguous pair, and was that pair a mismatch."""
    compared = bytearray(len(target))
    mismatches = bytearray(len(target))
    for block in blocks:
        for target_position, _query_position, target_base, query_base in iter_aligned_pairs(
            block, target=target, query=query
        ):
            if target_base in _DNA and query_base in _DNA:
                compared[target_position] = 1
                if target_base != query_base:
                    mismatches[target_position] = 1
    return compared, mismatches
