"""IR detection and the mandatory gate.

The gate is the project's most important safety net: naive IR detection was silently
wrong on 2/9 published plastomes (PROJECT.md §9A), and every position-based analysis
below it inherits the error. Rules:

* identity is always measured on the **full copies**, never on the exact-match core
  (a perfect 22,998 bp core sat inside copies that are only 92.62 % identical);
* the two copies are related by a mirror: position ``i`` in copy A pairs with
  ``C - i`` in copy B (reverse-complemented), with ``C = a1 + b2 - 1`` — this is the
  pairing constant measured in PROJECT.md §8.3;
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

from plastid.core.model import Genome, Interval
from plastid.core.sequence import complement_base, reverse_complement

__all__ = [
    "IR_COPY_LENGTH_TOLERANCE",
    "IR_MIN_IDENTITY",
    "IR_MIN_LENGTH",
    "IRDetectionError",
    "IRGateError",
    "IRResult",
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
IR_COPY_LENGTH_TOLERANCE = 5

_KMER = 20  # the validated probe used K=20
_DIAGONALS_CONSIDERED = 12
_DENOVO_LENGTH_WINDOW = (10_000, 35_000)  # validated validity window
_DENOVO_MIN_IDENTITY = 0.95  # validated validity threshold


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
    # equal, else "outer" (the JLB/JLA side) or "inner" (the JSB/JSA side). It is recorded because
    # the number alone cannot say which alignment produced it (H-25).
    excess_at: str = "none"

    @property
    def length(self) -> int:
        return min(self.ira.length, self.irb.length)

    @property
    def gate_passed(self) -> bool:
        return self.identity >= IR_MIN_IDENTITY

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


def _trimmed_pair(
    first: Interval, second: Interval, excess: int, where: str
) -> tuple[Interval, Interval]:
    """The pair with ``abs(excess)`` bases taken off the **longer** copy at the named end.

    ``where`` is ``"outer"`` — the first copy's start or the second copy's end, the JLB/JLA side —
    or ``"inner"``, the first copy's end or the second copy's start, the JSB/JSA side. Trimming the
    longer copy to the shorter's length is what makes the mirror constant of §8.3 valid again: the
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
    """One alignment of the pair: the longer copy trimmed at ``where``, then the §8.3 constant."""
    trimmed_first, _trimmed_second = _trimmed_pair(first, second, excess, where)
    span = Interval(
        trimmed_first.start, trimmed_first.start + compared_length, trimmed_first.strand
    )
    constant = mirror_constant_for(first, second, where)
    identity, mismatches, compared, exact_core = pair_identity(sequence, span, constant)
    return identity, mismatches, compared, exact_core, constant, where


def measure_ir(
    sequence: str,
    ira: Interval,
    irb: Interval,
    *,
    source: str,
    frame_offset: int = 0,
) -> IRResult:
    """Measure the identity of an IR pair on the full copies (shorter copy when unequal).

    When the annotated copies differ in length by *d*, the mirror constant that §8.3 states
    (``a1 + b2 - 1``) pairs the wrong bases, and an exact repeat measures as unrelated — the
    identity of two random sequences, ~25 %. Measured 17/10/2026 (H-25): *d* = 1, 2 and 4 on three
    Iridaceae records gave 0.2924 / 0.2636 / 0.2778, and the gate refused all three.

    The annotation does not say **which** end carries the difference, so both alignments it implies
    are measured and the better is kept. That is not a search for a match: a pair that is not a
    repeat scores ~25 % under either alignment, so neither can carry it over the gate.
    """
    first, second = sorted((ira, irb), key=lambda interval: interval.start)
    if first.end > second.start:
        raise ValueError(f"IR copies must be disjoint and ordered: {first} overlaps {second}")
    compared_length = min(first.length, second.length)
    excess = first.length - second.length
    if excess == 0:
        alignments = [_aligned(sequence, first, second, compared_length, 0, "none")]
    else:
        alignments = [
            _aligned(sequence, first, second, compared_length, excess, "inner"),
            _aligned(sequence, first, second, compared_length, excess, "outer"),
        ]
    identity, mismatches, compared, exact_core, constant, excess_at = max(
        alignments, key=lambda alignment: alignment[0]
    )
    return IRResult(
        ira=ira,
        irb=irb,
        identity=identity,
        mismatches=mismatches,
        compared_length=compared,
        exact_core=exact_core,
        mirror_constant=constant,
        source=source,
        frame_offset=frame_offset,
        excess_at=excess_at,
    )


def ir_from_annotation(genome: Genome) -> IRResult | None:
    """The IR pair from GenBank ``repeat_region`` annotations, structurally.

    No substring guessing on qualifier text (the CGAS trap): a candidate is any
    ``repeat_region`` of at least :data:`IR_MIN_LENGTH`; exactly two must remain, they must
    be disjoint, and their lengths may differ by at most :data:`IR_COPY_LENGTH_TOLERANCE`.
    """
    candidates = sorted(
        (region for region in genome.repeat_regions if region.length >= IR_MIN_LENGTH),
        key=lambda interval: interval.start,
    )
    if len(candidates) != 2:
        return None
    first, second = candidates
    if first.end > second.start:
        return None
    if abs(first.length - second.length) > IR_COPY_LENGTH_TOLERANCE:
        return None
    return measure_ir(genome.sequence, first, second, source="annotation")


def ir_de_novo(genome: Genome, *, frame_offset: int = 0) -> IRResult | None:
    """Mirror-diagonal de-novo detection — a port of the validated probe (``strict_best``).

    Build the mirror-diagonal histogram from exact K-mer matches against the reverse
    complement, walk the longest exact run on each of the top diagonals, and accept only
    geometries inside the validated window (10–35 kb, ≥ 95 % identity on the full copies).

    This replaces the mismatch-tolerant extension (``ir_probe_v2.detect_v2``): measured
    2026-10-02, that version over-extends — 27,601 bp against the true 26,217 bp on
    NC_086581 — because it accepts chunks on a global identity ratio. The exact-run version
    is the one behind the §8.3 numbers.

    Genomes whose origin sits inside an IR copy are handled by retrying on a half-rotated
    sequence; the result then carries that rotation in ``frame_offset``.
    """
    result = _scan_linear(genome, frame_offset)
    if result is not None:
        return result
    if genome.topology != "circular":
        return None
    shift = len(genome) // 2
    rotated = _rotate_genome(genome, shift)
    return _scan_linear(rotated, (frame_offset + shift) % len(genome))


def _scan_linear(genome: Genome, frame_offset: int) -> IRResult | None:
    sequence = genome.sequence
    length = len(sequence)
    mirror = reverse_complement(sequence)

    index: dict[str, list[int]] = {}
    for position in range(length - _KMER + 1):
        index.setdefault(sequence[position : position + _KMER], []).append(position)

    diagonals: Counter[int] = Counter()
    for position in range(length - _KMER + 1):
        for mate in index.get(mirror[position : position + _KMER], ()):
            diagonals[mate - position] += 1

    best: IRResult | None = None
    for diagonal, _ in diagonals.most_common(_DIAGONALS_CONSIDERED):
        geometry = _exact_run_geometry(sequence, mirror, diagonal)
        if geometry is None:
            continue
        ira, irb = geometry
        ir = measure_ir(sequence, ira, irb, source="denovo_diag", frame_offset=frame_offset)
        if not (_DENOVO_LENGTH_WINDOW[0] <= ir.length <= _DENOVO_LENGTH_WINDOW[1]):
            continue
        if ir.identity < _DENOVO_MIN_IDENTITY:
            continue
        if best is None or ir.length > best.length:
            best = ir
    return best


def _exact_run_geometry(
    sequence: str, mirror: str, diagonal: int
) -> tuple[Interval, Interval] | None:
    """Longest exact run along one mirror diagonal → the two copy intervals, or None."""
    length = len(sequence)
    if diagonal >= 0:
        n = length - diagonal
        left, right = sequence[diagonal:], mirror[:n]
        start_a, start_r = diagonal, 0
    else:
        n = length + diagonal
        left, right = sequence[:n], mirror[-diagonal:]
        start_a, start_r = 0, -diagonal

    best_run = run = best_start = run_start = 0
    for offset in range(n):
        if left[offset] == right[offset]:
            if run == 0:
                run_start = offset
            run += 1
            if run > best_run:
                best_run, best_start = run, run_start
        else:
            run = 0
    if best_run == 0:
        return None

    ira_start = start_a + best_start
    shifted = start_r + best_start
    ira = Interval(ira_start, ira_start + best_run)
    irb = Interval(length - shifted - best_run, length - shifted)
    if irb.start < 0 or irb.start < ira.end or ira.end > length or irb.end > length:
        return None
    if irb.start - ira.end <= 0:
        return None
    return ira, irb


def _rotate_genome(genome: Genome, shift: int) -> Genome:
    from plastid.core.normalize import rotate_genome

    return rotate_genome(genome, shift)


def detect_ir(genome: Genome) -> IRResult:
    """Annotation first; de-novo detection when there is none. Raises when neither works.

    When both paths produce a result they are cross-checked and the disagreement is stored
    in :attr:`IRResult.denovo_length` — never silently dropped. Measured reason: NC_063920
    annotates a 26,161 bp IR whose copies are 92.6 % identical, while de-novo detection
    finds a perfectly matching 22,998 bp core. Annotation wins, but a 3,163 bp disagreement
    must be visible.
    """
    annotated = ir_from_annotation(genome)
    detected = ir_de_novo(genome)
    if annotated is not None:
        if detected is not None:
            return replace(annotated, denovo_length=detected.length)
        return annotated
    if detected is None:
        raise IRDetectionError(
            f"{genome.accession or genome.name}: no valid IR pair found "
            f"(annotation absent; de-novo scan failed)"
        )
    return detected


def enforce_ir_gate(ir: IRResult, *, min_identity: float = IR_MIN_IDENTITY) -> IRResult:
    """Raise :class:`IRGateError` unless the two copies are identical enough."""
    if ir.identity < min_identity:
        raise IRGateError(
            f"IR gate failed: the two copies are {ir.identity * 100:.2f} % identical "
            f"(need >= {min_identity * 100:.1f} %) — position-based analysis must not run"
        )
    return ir
