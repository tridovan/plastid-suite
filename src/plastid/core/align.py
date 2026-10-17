"""Multiple-alignment maths — pure functions over an alignment, no I/O (P2-07, P2-08).

An alignment is a set of equal-length rows (gaps as ``-``) produced by an external aligner
through the tool door. Everything here is deterministic and offline, so a manuscript number can
be recomputed from the stored alignment without rerunning the aligner — only the alignment
itself needs the tool.

Every measure states how it treats gaps. The default is **complete deletion**: a column holding
any gap or ambiguous base is dropped and counted, never silently scored as a match.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

__all__ = [
    "SUBSTITUTION_TYPES",
    "Alignment",
    "WindowDiversity",
    "codon_frame_breaks",
    "columns_for_span",
    "nucleotide_diversity",
    "reference_positions",
    "substitution_class",
    "substitution_counts",
    "usable_columns",
    "variable_sites",
    "window_diversity",
]

TRANSITIONS = frozenset({("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")})
_DNA = frozenset("ACGT")

# The twelve **directed** changes, named ``<from>_to_<to>``. They are the primitive fact a
# comparison of two rows produces; ``transitions`` and ``transversions`` are sums over them, and
# the six *complementary* pairs a manuscript draws (A↔G, C↔T, A↔C, T↔G, C↔G, A↔T) are sums of two
# directed ones each. Measured 17/10/2026: only the two aggregates existed, so the paper's Figure 7
# — substitution type × species — had no source in this toolkit at all.
SUBSTITUTION_TYPES: tuple[str, ...] = (
    "A_to_G",
    "G_to_A",
    "C_to_T",
    "T_to_C",
    "C_to_G",
    "G_to_C",
    "A_to_T",
    "T_to_A",
    "A_to_C",
    "C_to_A",
    "T_to_G",
    "G_to_T",
)

_TRANSITION_TYPES = frozenset({"A_to_G", "G_to_A", "C_to_T", "T_to_C"})


def substitution_class(type_name: str) -> str:
    """``"Transition"`` or ``"Transversion"`` — the same split :data:`TRANSITIONS` makes."""
    if type_name not in SUBSTITUTION_TYPES:
        raise ValueError(
            f"{type_name!r} is not a substitution type (have: {', '.join(SUBSTITUTION_TYPES)})"
        )
    return "Transition" if type_name in _TRANSITION_TYPES else "Transversion"


@dataclass(frozen=True, slots=True)
class Alignment:
    """Rows of one alignment: equal length, gaps as ``-``, order meaningful (row 0 is reference)."""

    names: tuple[str, ...]
    sequences: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.names) != len(self.sequences):
            raise ValueError(
                f"alignment has {len(self.names)} name(s) but {len(self.sequences)} sequence(s)"
            )
        if len(self.sequences) < 1:
            raise ValueError("alignment has no sequences")
        lengths = {len(sequence) for sequence in self.sequences}
        if len(lengths) != 1:
            raise ValueError(f"alignment rows differ in length: {sorted(lengths)}")
        if not self.sequences[0]:
            raise ValueError("alignment is empty")
        # Fail loud rather than count nothing: a column is classified by ``_DNA``, so a lowercase
        # row would have *every* column skipped and report zero comparable columns instead of an
        # error. Measured 03/10/2026 — MAFFT 7.526 prints lowercase, and the adapters are the
        # boundary that upper-cases it (external/mafft.py).
        if any(character.islower() for sequence in self.sequences for character in sequence):
            raise ValueError(
                "an alignment row carries lowercase bases — the counters read ACGT only, so "
                "every column would be skipped in silence; upper-case it at the boundary"
            )

    @property
    def n_rows(self) -> int:
        return len(self.sequences)

    @property
    def length(self) -> int:
        return len(self.sequences[0])

    def column(self, index: int) -> str:
        return "".join(sequence[index] for sequence in self.sequences)


@dataclass(frozen=True, slots=True)
class WindowDiversity:
    """π over one window of the alignment, with the columns it actually used."""

    start_1b: int
    end_1b: int
    usable_columns: int
    pi: float | None  # None when no column of the window could be used


def _usable(column: str) -> bool:
    """A column counts only when every row carries an unambiguous base — complete deletion."""
    return all(base in _DNA for base in column)


def usable_columns(alignment: Alignment) -> tuple[int, ...]:
    """The columns complete deletion keeps: every row carries an unambiguous base.

    Public because the population-genetic estimators (P4-03) need the same column set as π —
    a statistic computed over a different set of columns than the one beside it is a table
    nobody can reconcile.
    """
    return tuple(index for index in range(alignment.length) if _usable(alignment.column(index)))


def nucleotide_diversity(alignment: Alignment, columns: Iterable[int]) -> tuple[float | None, int]:
    """π over the given columns: ``(value, usable columns)``; ``None`` when nothing is usable.

    ``π = mean over sites of (n/(n−1)) · (1 − Σ p²)``, the Nei & Li (1979) estimator, with
    ``p`` the base frequencies at that site. Gapped columns are dropped, not scored.
    """
    indices = list(columns)
    usable = [index for index in indices if _usable(alignment.column(index))]
    if not usable or alignment.n_rows < 2:
        return None, 0
    n = alignment.n_rows
    total = 0.0
    for index in usable:
        counts = Counter(alignment.column(index))
        homozygosity = sum((count / n) ** 2 for count in counts.values())
        total += n / (n - 1) * (1 - homozygosity)
    return total / len(usable), len(usable)


def window_diversity(
    alignment: Alignment, *, window: int, step: int
) -> tuple[WindowDiversity, ...]:
    """Sliding-window π over alignment columns; the last window may be shorter than ``window``."""
    if window <= 0 or step <= 0:
        raise ValueError(f"window and step must be positive, got window={window}, step={step}")
    if step > window:
        raise ValueError(f"step {step} exceeds window {window} — columns would be skipped")
    windows: list[WindowDiversity] = []
    start = 0
    while start < alignment.length:
        end = min(start + window, alignment.length)
        pi, usable = nucleotide_diversity(alignment, range(start, end))
        windows.append(
            WindowDiversity(start_1b=start + 1, end_1b=end, usable_columns=usable, pi=pi)
        )
        if end == alignment.length:
            break
        start += step
    return tuple(windows)


def columns_for_span(alignment: Alignment, row: int, span: Sequence[int]) -> tuple[int, ...]:
    """Alignment columns holding a reference row's own positions ``[start, end)`` (0-based).

    Gaps in the reference row are skipped: they are not positions of that genome. Used to map an
    annotated feature onto alignment columns, which is the only honest way to report a per-gene
    π from an alignment that contains gaps.
    """
    start, end = int(span[0]), int(span[1])
    if start < 0 or end < start:
        raise ValueError(f"bad span {span!r}")
    columns: list[int] = []
    position = 0
    for index, base in enumerate(alignment.sequences[row]):
        if base == "-":
            continue
        if start <= position < end:
            columns.append(index)
        position += 1
    return tuple(columns)


def reference_positions(alignment: Alignment, row: int = 0) -> tuple[int | None, ...]:
    """For every alignment column, the reference position it holds — ``None`` where it is a gap."""
    positions: list[int | None] = []
    position = 0
    for base in alignment.sequences[row]:
        if base == "-":
            positions.append(None)
        else:
            positions.append(position)
            position += 1
    return tuple(positions)


def variable_sites(alignment: Alignment) -> tuple[int, int]:
    """``(variable, informative)`` column counts — informative means at least two bases twice."""
    variable = 0
    informative = 0
    for index in range(alignment.length):
        column = alignment.column(index)
        if not _usable(column):
            continue
        counts = Counter(column)
        if len(counts) > 1:
            variable += 1
            if sum(1 for count in counts.values() if count >= 2) >= 2:
                informative += 1
    return variable, informative


def codon_frame_breaks(alignment: Alignment) -> tuple[str, ...]:
    """Rows whose gaps break the reading frame — empty when the alignment is codon-aligned.

    A codon model (PAML's ``yn00``) reads every row three bases at a time and has no idea what a
    gap means: a gap that is not a whole codon, or that starts mid-codon, silently shifts the
    frame of that row and every codon after it is read wrong. Measured 04/10/2026 on a synthetic
    pair: a one-base gap made yn00 report a spurious stop codon and exit 255 — a lucky frame, not
    a check. So the check is ours, and it names the row and the position rather than a verdict.
    """
    problems: list[str] = []
    for name, sequence in zip(alignment.names, alignment.sequences, strict=True):
        bases = len(sequence.replace("-", ""))
        if bases % 3:
            problems.append(
                f"row {name!r}: {bases} aligned base(s) is not a whole number of codons"
            )
            continue
        run_start = 0
        position = 0
        while position < len(sequence):
            if sequence[position] == "-":
                run_start = position
                while position < len(sequence) and sequence[position] == "-":
                    position += 1
                run = position - run_start
                before = run_start - sequence[:run_start].count("-")
                if before % 3:
                    problems.append(
                        f"row {name!r}: a gap starts after {before} base(s) — not a codon boundary"
                    )
                elif run % 3:
                    problems.append(f"row {name!r}: a gap of {run} base(s) is not a whole codon")
            else:
                position += 1
    return tuple(problems)


def substitution_counts(reference: str, query: str) -> dict[str, int]:
    """Substitutions between two aligned rows, by type, plus the transition/transversion totals.

    Columns where either row is not an unambiguous base are skipped and counted, so a ratio is
    never computed over an unstated denominator.

    Every one of :data:`SUBSTITUTION_TYPES` is a key of the result, zero included — a reader must
    never have to guess whether an absent key means "none" or "not counted". ``transitions`` and
    ``transversions`` are sums over those twelve, so they can never disagree with them.
    """
    types = dict.fromkeys(SUBSTITUTION_TYPES, 0)
    transitions = 0
    transversions = 0
    skipped = 0
    for ref_base, query_base in zip(reference, query, strict=True):
        if ref_base not in _DNA or query_base not in _DNA:
            skipped += 1
            continue
        if ref_base == query_base:
            continue
        types[f"{ref_base}_to_{query_base}"] += 1
        if (ref_base, query_base) in TRANSITIONS:
            transitions += 1
        else:
            transversions += 1
    return {
        "compared": len(reference) - skipped,
        "skipped": skipped,
        "substitutions": transitions + transversions,
        "transitions": transitions,
        "transversions": transversions,
        **types,
    }
