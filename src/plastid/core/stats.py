"""Statistics that name their test, state their n, and refuse when the data cannot support them.

Decision D5-02: the tests are scipy's, not hand-written. Writing Mann–Whitney's tie correction by
hand is exactly how a plausible wrong p-value gets into a table, and the reason ``popgen`` was
hand-written (D4-03) was that its *alternatives* were a Windows GUI and a VCF-only CLI — a reason
that does not carry here. What is ours is the discipline around it:

* every result carries the **test's name**, the **n** it used and the **assumption** it rests on;
* a statistic the data cannot support is ``NOT_COMPUTABLE`` with a reason — never a number
  (§8.30's rule, applied to every function here);
* every function is pinned by a hand-computed vector in the tests, not by scipy against itself.

scipy is an optional dependency (``pip install 'plastid-suite[stats]'``), imported inside the
function that needs it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from plastid.core.deps import require_module

__all__ = [
    "NOT_COMPUTABLE",
    "OK",
    "Correlation",
    "Descriptives",
    "TestResult",
    "benjamini_hochberg",
    "bonferroni",
    "describe",
    "kruskal_wallis",
    "mann_whitney",
    "pearson",
    "spearman",
]

OK = "ok"
NOT_COMPUTABLE = "not computable"

# Two observations per group is the floor: with one there is nothing to rank against, and the
# test would report p = 1 as if that were a finding.
MIN_PER_GROUP = 2
# Spearman and Pearson need three points before a correlation means anything.
MIN_FOR_CORRELATION = 3


@dataclass(frozen=True, slots=True)
class Descriptives:
    """The five-number summary and the two moments — every field ``None`` when n = 0."""

    n: int
    mean: float | None
    sd: float | None  # sample standard deviation (n − 1), None below two observations
    median: float | None
    q1: float | None
    q3: float | None
    minimum: float | None
    maximum: float | None


@dataclass(frozen=True, slots=True)
class TestResult:
    """A hypothesis test: what it was, on how many, and what it assumes."""

    test: str
    assumption: str
    n: tuple[int, ...]
    statistic: float | None = None
    p: float | None = None
    status: str = OK
    reason: str = ""

    @property
    def computable(self) -> bool:
        return self.status == OK


@dataclass(frozen=True, slots=True)
class Correlation:
    """A correlation coefficient with its n and its p."""

    test: str
    assumption: str
    n: int
    coefficient: float | None = None
    p: float | None = None
    status: str = OK
    reason: str = ""

    @property
    def computable(self) -> bool:
        return self.status == OK


def _not_computable(test: str, assumption: str, n: tuple[int, ...], reason: str) -> TestResult:
    return TestResult(test=test, assumption=assumption, n=n, status=NOT_COMPUTABLE, reason=reason)


def describe(values: Sequence[float]) -> Descriptives:
    """Mean, sample sd and the quartiles; ``None`` wherever the data cannot give one."""
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        return Descriptives(0, None, None, None, None, None, None, None)
    mean = sum(ordered) / n
    sd = math.sqrt(sum((value - mean) ** 2 for value in ordered) / (n - 1)) if n > 1 else None
    return Descriptives(
        n=n,
        mean=mean,
        sd=sd,
        median=_median(ordered),
        q1=_median(ordered[: n // 2]) if n > 1 else None,
        q3=_median(ordered[(n + 1) // 2 :]) if n > 1 else None,
        minimum=ordered[0],
        maximum=ordered[-1],
    )


def _median(ordered: Sequence[float]) -> float | None:
    """The median of an already-sorted sequence — the ordinary midpoint rule."""
    n = len(ordered)
    if n == 0:
        return None
    middle = n // 2
    if n % 2:
        return float(ordered[middle])
    return (ordered[middle - 1] + ordered[middle]) / 2


_MANN_WHITNEY_ASSUMPTION = (
    "two-sided Mann–Whitney U over independent observations; scipy's exact distribution when "
    "there are no ties and both groups are small, otherwise the normal approximation with the "
    "tie correction (scipy's own rule)"
)


def mann_whitney(first: Sequence[float], second: Sequence[float]) -> TestResult:
    """Two groups, ranks — the test to use when normality is not something anyone measured."""
    n = (len(first), len(second))
    if min(n) < MIN_PER_GROUP:
        return _not_computable(
            "Mann–Whitney U",
            _MANN_WHITNEY_ASSUMPTION,
            n,
            f"needs at least {MIN_PER_GROUP} observations per group, has {n[0]} and {n[1]}",
        )
    stats = require_module("scipy.stats", extra="stats", purpose="Statistics")
    outcome = stats.mannwhitneyu(list(first), list(second), alternative="two-sided")
    return TestResult(
        test="Mann–Whitney U",
        assumption=_MANN_WHITNEY_ASSUMPTION,
        n=n,
        statistic=float(outcome.statistic),
        p=float(outcome.pvalue),
    )


_KRUSKAL_ASSUMPTION = (
    "Kruskal–Wallis one-way analysis of variance on ranks, two-sided; the chi-square "
    "approximation, with scipy's tie correction"
)


def kruskal_wallis(*groups: Sequence[float]) -> TestResult:
    """Three or more groups. With two, :func:`mann_whitney` is the sharper test."""
    n = tuple(len(group) for group in groups)
    if len(groups) < 3:
        return _not_computable(
            "Kruskal–Wallis H",
            _KRUSKAL_ASSUMPTION,
            n,
            f"needs at least 3 groups, has {len(groups)} — use Mann–Whitney U for two",
        )
    if min(n) < MIN_PER_GROUP:
        return _not_computable(
            "Kruskal–Wallis H",
            _KRUSKAL_ASSUMPTION,
            n,
            f"needs at least {MIN_PER_GROUP} observations per group, has {list(n)}",
        )
    stats = require_module("scipy.stats", extra="stats", purpose="Statistics")
    outcome = stats.kruskal(*[list(group) for group in groups])
    return TestResult(
        test="Kruskal–Wallis H",
        assumption=_KRUSKAL_ASSUMPTION,
        n=n,
        statistic=float(outcome.statistic),
        p=float(outcome.pvalue),
    )


def _correlation(
    name: str, assumption: str, x: Sequence[float], y: Sequence[float], *, rank: bool
) -> Correlation:
    n = len(x)
    if n != len(y):
        raise ValueError(f"{name}: {n} x-values against {len(y)} y-values")
    if n < MIN_FOR_CORRELATION:
        return Correlation(
            test=name,
            assumption=assumption,
            n=n,
            status=NOT_COMPUTABLE,
            reason=f"needs at least {MIN_FOR_CORRELATION} pairs, has {n}",
        )
    if len(set(x)) < 2 or len(set(y)) < 2:
        return Correlation(
            test=name,
            assumption=assumption,
            n=n,
            status=NOT_COMPUTABLE,
            reason="one of the two series is constant — a correlation with it is undefined",
        )
    stats = require_module("scipy.stats", extra="stats", purpose="Statistics")
    function = stats.spearmanr if rank else stats.pearsonr
    outcome = function(list(x), list(y))
    return Correlation(
        test=name,
        assumption=assumption,
        n=n,
        coefficient=float(outcome.statistic),
        p=float(outcome.pvalue),
    )


def spearman(x: Sequence[float], y: Sequence[float]) -> Correlation:
    """Rank correlation — what to use when the relationship need not be a straight line."""
    return _correlation(
        "Spearman ρ",
        "monotone relationship, no distributional assumption; two-sided",
        x,
        y,
        rank=True,
    )


def pearson(x: Sequence[float], y: Sequence[float]) -> Correlation:
    """Linear correlation — the one that assumes a line, said out loud in ``assumption``."""
    return _correlation(
        "Pearson r",
        "a straight-line relationship between two roughly normal, homoscedastic series; two-sided",
        x,
        y,
        rank=False,
    )


def benjamini_hochberg(pvalues: Sequence[float | None]) -> tuple[float | None, ...]:
    """Benjamini–Hochberg adjusted p-values, monotone by construction.

    ``adjusted(i) = min over j ≥ i of (m · p(j) / j)`` on the sorted p-values — the step-up
    procedure, written out because the monotonicity half is the half that is easy to get wrong.
    """
    positions = [index for index, p in enumerate(pvalues) if p is not None]
    m = len(positions)
    if m == 0:
        return tuple(None for _ in pvalues)
    by_size = sorted(positions, key=lambda index: pvalues[index])
    result: list[float | None] = [None] * len(pvalues)
    running = 1.0
    for rank in range(m, 0, -1):
        value = pvalues[by_size[rank - 1]]
        assert value is not None  # by construction: `positions` holds only the non-blank ones
        running = min(running, value * m / rank)
        result[by_size[rank - 1]] = min(1.0, running)
    return tuple(result)


def bonferroni(pvalues: Sequence[float | None]) -> tuple[float | None, ...]:
    """Bonferroni: ``min(1, m · p)``. The blunter correction, kept for comparison."""
    m = sum(1 for p in pvalues if p is not None)
    if m == 0:
        return tuple(None for _ in pvalues)
    return tuple(None if p is None else min(1.0, p * m) for p in pvalues)
