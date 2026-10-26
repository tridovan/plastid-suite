"""Population-genetic estimators — pure Python, one function per statistic (P4-03; decision D4-03).

DnaSP is a Windows GUI (unscriptable, unusable in Docker) and the scriptable alternatives want
VCF or read data this project does not have, so the estimators are written here from their
published formulas — every one testable against an analytic vector, which a GUI's output is not.
Nothing in this module needs an external tool, so the minimum install runs it (CLAUDE.md).

Each function names its formula and its source, and each returns ``None`` where the statistic is
**not computable** — never a plausible number. The cases are not conventions but arithmetic:

* **Tajima's D is undefined for n ≤ 3.** The variance ``e1·S + e2·S(S−1)`` collapses to exactly
  zero at n = 3 (c1 = c2 = 0), so the statistic is 0/0. Measured here, and the same restriction
  is in EggLib 3.6.1 (``Diversity1.cpp``: ``if (ns > 3)``).
* **Fu & Li's D\\* and F\\* are undefined for n ≤ 3** for the same reason (uD = vD = 0 at n = 3).
* **Fst is undefined when the groups share no differences at all** (H_b = 0 → 1 − H_w/0).

Conventions, stated because a number without them means nothing: columns are kept only when
every sequence carries an unambiguous base (**complete deletion**, as everywhere else in this
toolkit), sites are counted over those columns, and every statistic travels with the n it was
computed from. Fu & Li's statistics need ``η_e`` (external mutations); **without an outgroup
this module uses the number of singletons** — the convention DnaSP uses when no outgroup is
given — and the module says so wherever the number is exported.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

from plastid.core.align import Alignment, nucleotide_diversity, usable_columns

__all__ = [
    "PopgenSummary",
    "fst_gst",
    "fu_li_d_star",
    "fu_li_f_star",
    "haplotype_counts",
    "haplotype_diversity",
    "pairwise_differences",
    "segregating_sites",
    "summarise",
    "tajimas_d",
    "watterson_theta",
]

# n ≥ 4: below it the variance of Tajima's and Fu & Li's statistics is identically zero.
MIN_SEQUENCES_FOR_NEUTRALITY = 4


@dataclass(frozen=True, slots=True)
class PopgenSummary:
    """Every statistic of one sample, each with the n and the site count behind it.

    ``None`` means *not computable* — the reason is in ``notes``. A zero is a measured zero
    (π of a constant alignment is 0, honestly); ``None`` is arithmetic that has no value.
    """

    n: int
    sites: int  # usable columns (complete deletion)
    segregating: int  # S
    singletons: int  # η_e — sites where one base occurs exactly once
    haplotypes: int
    haplotype_diversity: float | None  # Hd
    pi_per_site: float | None  # θ(π)
    theta_w_per_site: float | None  # θ(S)
    tajimas_d: float | None
    fu_li_d_star: float | None
    fu_li_f_star: float | None
    notes: tuple[str, ...] = ()


def haplotype_counts(alignment: Alignment) -> Counter[str]:
    """Distinct sequences (over the usable columns) and how many times each occurs."""
    columns = usable_columns(alignment)
    counts: Counter[str] = Counter()
    for row in alignment.sequences:
        counts["".join(row[index] for index in columns)] += 1
    return counts


def haplotype_diversity(alignment: Alignment) -> float | None:
    """Hd, Nei & Tajima (1981): ``(n/(n−1)) · (1 − Σ p_i²)`` over the distinct haplotypes.

    ``None`` for fewer than two sequences (n/(n−1) has no meaning there); a single haplotype
    among n ≥ 2 sequences is a measured **0**, not ``None`` — every sequence really is the same.
    """
    n = alignment.n_rows
    if n < 2:
        return None
    counts = haplotype_counts(alignment)
    total = sum(counts.values())
    return n / (n - 1) * (1 - sum((count / total) ** 2 for count in counts.values()))


def segregating_sites(alignment: Alignment) -> int:
    """S: usable columns carrying more than one distinct base."""
    return sum(1 for index in usable_columns(alignment) if len(set(alignment.column(index))) > 1)


def _watterson_constant(n: int) -> tuple[float, float]:
    """``a1 = Σ 1/i`` and ``a2 = Σ 1/i²`` for i = 1 … n−1 — the harmonic numbers of the theory."""
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    return a1, a2


def _singletons(alignment: Alignment) -> int:
    """η_e: singletons. Without an outgroup, a singleton is a base occurring exactly once.

    For n = 2 this counts **two** per segregating site — both bases occur once — which is what
    EggLib's ``SiteDiversity::S()`` does and what DnaSP's no-outgroup convention means.
    """
    total = 0
    for index in usable_columns(alignment):
        counts = Counter(alignment.column(index))
        total += sum(1 for count in counts.values() if count == 1)
    return total


def watterson_theta(alignment: Alignment) -> float | None:
    """θ(S), Watterson (1975): ``S / a1`` per site — over the usable columns.

    ``None`` for fewer than two sequences. A sample with no segregating site gives 0.
    """
    n = alignment.n_rows
    if n < 2:
        return None
    sites = len(usable_columns(alignment))
    if sites == 0:
        return None
    a1, _ = _watterson_constant(n)
    return segregating_sites(alignment) / a1 / sites


def tajimas_d(alignment: Alignment) -> float | None:
    """Tajima's D (1989): ``(π_total − θ_total) / sqrt(e1·S + e2·S(S−1))``.

    ``None`` when the variance is zero — which is exactly the case for n ≤ 3 (measured: at
    n = 3, c1 = c2 = 0) and for S = 0. Both are 0/0, not a value.
    """
    n = alignment.n_rows
    if n < MIN_SEQUENCES_FOR_NEUTRALITY:
        return None
    columns = usable_columns(alignment)
    sites = len(columns)
    if sites == 0:
        return None
    pi_per_site, _ = nucleotide_diversity(alignment, columns)
    if pi_per_site is None:
        return None
    s = segregating_sites(alignment)
    a1, a2 = _watterson_constant(n)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    variance = e1 * s + e2 * s * (s - 1)
    if variance <= 0:
        return None
    return (pi_per_site * sites - s / a1) / math.sqrt(variance)


def fu_li_d_star(alignment: Alignment) -> float | None:
    """Fu & Li's D\\* (1993), no-outgroup form: ``((n/(n−1))·η − a1·η_e) / sqrt(uD·η + vD·η²)``.

    ``η`` is the number of mutations (= S under the infinite-sites model this test assumes) and
    ``η_e`` the singletons. The variance constants are Fu & Li (1993) as written by Simonsen,
    Churchill & Aquadro (1995) and implemented in EggLib 3.6.1 — the module's numbers are
    cross-checked against that implementation. ``None`` for n ≤ 3, where
    uD = vD = 0.
    """
    n = alignment.n_rows
    if n < MIN_SEQUENCES_FOR_NEUTRALITY:
        return None
    columns = usable_columns(alignment)
    if not columns:
        return None
    eta = segregating_sites(alignment)
    if eta == 0:
        return None  # the variance is zero; nothing to test
    eta_e = _singletons(alignment)
    a1, a2 = _watterson_constant(n)
    cn = 2 * (n * a1 - 2 * (n - 1)) / ((n - 1) * (n - 2))
    dn = (
        cn
        + (n - 2) / ((n - 1) * (n - 1))
        + (2 / (n - 1)) * (1.5 - (2 * (a1 + 1 / n) - 3) / (n - 2) - 1 / n)
    )
    v_d = (
        (n * n / ((n - 1.0) * (n - 1.0))) * a2
        + a1 * a1 * dn
        - 2 * n * a1 * (a1 + 1) / ((n - 1.0) * (n - 1.0))
    ) / (a1 * a1 + a2)
    u_d = (n / (n - 1.0)) * (a1 - n / (n - 1.0)) - v_d
    variance = u_d * eta + v_d * eta * eta
    if variance <= 0:
        return None
    return ((n / (n - 1.0)) * eta - a1 * eta_e) / math.sqrt(variance)


def fu_li_f_star(alignment: Alignment) -> float | None:
    """Fu & Li's F\\* (1993), no-outgroup form: ``(π_total − (n−1)/n·η_e) / sqrt(uF·η + vF·η²)``.

    Variance constants from Simonsen, Churchill & Aquadro (1995), as in EggLib 3.6.1.
    ``None`` for n ≤ 3 and for a sample with no segregating site.
    """
    n = alignment.n_rows
    if n < MIN_SEQUENCES_FOR_NEUTRALITY:
        return None
    columns = usable_columns(alignment)
    if not columns:
        return None
    eta = segregating_sites(alignment)
    if eta == 0:
        return None
    eta_e = _singletons(alignment)
    a1, a2 = _watterson_constant(n)
    v_f = (
        (2.0 * n * n * n + 110 * n * n - 255 * n + 153) / (9.0 * n * n * (n - 1))
        + (2.0 * (n - 1) * a1) / (n * n)
        - 8.0 * a2 / n
    ) / (a1 * a1 + a2)
    u_f = ((4.0 * n * n + 19 * n + 3 - 12 * (n + 1) * (a1 + 1.0 / n)) / (3 * n * (n - 1))) / (
        a1
    ) - v_f
    variance = u_f * eta + v_f * eta * eta
    if variance <= 0:
        return None
    pi_total, _ = nucleotide_diversity(alignment, columns)
    if pi_total is None:
        return None
    return (pi_total * len(columns) - (n - 1.0) * eta_e / n) / math.sqrt(variance)


def pairwise_differences(left: Alignment, right: Alignment | None = None) -> list[int]:
    """Pairwise differences between sequences, over the usable columns of the pair.

    With ``right=None`` the pairs are within one sample; otherwise every (left, right) pair is
    compared. A gap or ambiguous base in either sequence of a column excludes that column for
    that pair only — the same complete-deletion rule, applied pairwise, which is what the
    between-group statistic needs.
    """
    if right is None:
        pairs = [
            (left.sequences[i], left.sequences[j])
            for i in range(left.n_rows)
            for j in range(i + 1, left.n_rows)
        ]
    else:
        pairs = [(a, b) for a in left.sequences for b in right.sequences]
    differences: list[int] = []
    for first, second in pairs:
        differences.append(
            sum(
                1
                for a, b in zip(first, second, strict=True)
                if a in "ACGT" and b in "ACGT" and a != b
            )
        )
    return differences


def fst_gst(groups: Sequence[Alignment]) -> float | None:
    """Fst between groups as Nei's Gst (1973): ``(H_T − H_S) / H_T``.

    ``H_S`` is the within-group gene diversity (Wright's ``1 − Σ p²`` over the group's haplotype
    frequencies), averaged over groups and weighted by group size; ``H_T`` is the same quantity
    computed on the pooled sample. The two boundary cases are exact, which is why this
    formulation is the one written here: **two groups with the same haplotype composition give
    H_S = H_T and Fst = 0**, and two groups sharing no haplotype give H_S = 0 and Fst = 1.

    The sequence-based estimator (Hudson, Slatkin & Maddison 1992) was implemented first and
    rejected with a measured reason: with two imaginary groups of two sequences each it returns
    **−1** for identical composition, because its within-group mean is a single pair — a
    small-sample artefact that would have gone into a manuscript table as a negative Fst.

    ``H_T = 0`` (every sequence identical) has no value: ``None``, never 0 and never 1. This Gst
    is **uncorrected** — at small n it is biased downward (Nei & Chesser 1983 give the corrected
    estimator, which needs more sequences than this panel has), so the per-group diversities and
    every n are published beside it rather than folded into the ratio.
    """
    samples = [group for group in groups if group.n_rows >= 1]
    if len(samples) < 2:
        return None
    # Every group needs at least two sequences: with one, the within-group diversity is zero by
    # construction, so the estimator returns exactly 1 for any two groups that differ — a number
    # that looks like complete differentiation and is an artefact of the sample size. Measured
    # 04/10/2026 on the mito sketch (one genome per group).
    if any(group.n_rows < 2 for group in samples):
        return None
    pooled = Counter[str]()
    total = 0
    weighted_s = 0.0
    for group in samples:
        counts = haplotype_counts(group)
        size = sum(counts.values())
        diversity = 1 - sum((count / size) ** 2 for count in counts.values())
        weighted_s += size * diversity
        pooled.update(counts)
        total += size
    if total == 0:
        return None
    h_t = 1 - sum((count / total) ** 2 for count in pooled.values())
    if h_t == 0:
        return None
    return (h_t - weighted_s / total) / h_t


def _mean(values: list[int]) -> float | None:
    return sum(values) / len(values) if values else None


def summarise(alignment: Alignment) -> PopgenSummary:
    """Every statistic of one sample, with the reasons for the ones that have no value."""
    n = alignment.n_rows
    columns = usable_columns(alignment)
    notes: list[str] = []
    if n < 2:
        notes.append("fewer than two sequences — nothing is estimable")
    if n < MIN_SEQUENCES_FOR_NEUTRALITY:
        notes.append(
            f"n = {n}: the neutrality tests need at least {MIN_SEQUENCES_FOR_NEUTRALITY} "
            f"sequences — at n ≤ 3 their variance is identically zero"
        )
    if not columns:
        notes.append("no column is usable under complete deletion")
    elif segregating_sites(alignment) == 0:
        notes.append("no segregating site — the neutrality tests have zero variance")
    return PopgenSummary(
        n=n,
        sites=len(columns),
        segregating=segregating_sites(alignment),
        singletons=_singletons(alignment),
        haplotypes=len(haplotype_counts(alignment)),
        haplotype_diversity=haplotype_diversity(alignment),
        pi_per_site=nucleotide_diversity(alignment, columns)[0],
        theta_w_per_site=watterson_theta(alignment),
        tajimas_d=tajimas_d(alignment),
        fu_li_d_star=fu_li_d_star(alignment),
        fu_li_f_star=fu_li_f_star(alignment),
        notes=tuple(notes),
    )
