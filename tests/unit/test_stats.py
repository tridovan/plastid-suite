"""Statistics: every function pinned by a vector computed by hand, not by scipy against itself.

The values below are derived in the comments the way a reader can re-derive them — that is the
point of D5-02's rule. A test that compares scipy with scipy would pass even if both were wrong.
"""

from __future__ import annotations

import importlib.util
import math
import sys

import pytest

from plastid.core.deps import ExportDependencyError
from plastid.core.stats import (
    NOT_COMPUTABLE,
    benjamini_hochberg,
    bonferroni,
    describe,
    kruskal_wallis,
    mann_whitney,
    pearson,
    spearman,
)

requires_stats = pytest.mark.skipif(
    importlib.util.find_spec("scipy") is None,
    reason="scipy is an optional dependency (pip install 'plastid-suite[stats]')",
)


def test_descriptives_are_the_hand_computed_summary() -> None:
    """[1, 2, 3, 4]: mean 2.5, sample sd sqrt(5/3) = 1.2909944…, q1 1.5, q3 3.5."""
    summary = describe([1, 2, 3, 4])
    assert (summary.n, summary.mean, summary.median) == (4, 2.5, 2.5)
    assert summary.sd == pytest.approx(math.sqrt(5 / 3), abs=1e-12)
    assert (summary.q1, summary.q3) == (1.5, 3.5)
    assert (summary.minimum, summary.maximum) == (1, 4)


def test_descriptives_of_nothing_are_blank_not_zero() -> None:
    summary = describe([])
    assert summary.n == 0
    assert summary.mean is None and summary.sd is None and summary.median is None


def test_the_sd_needs_two_observations() -> None:
    assert describe([7.0]).sd is None
    assert describe([7.0]).median == 7.0


@requires_stats
def test_mann_whitney_matches_the_exact_hand_computed_p() -> None:
    """Two groups of three, completely separated: U = 0 and the exact two-sided p = 2/20 = 0.1."""
    outcome = mann_whitney([1, 2, 3], [4, 5, 6])
    assert outcome.n == (3, 3)
    assert outcome.statistic == pytest.approx(0.0)
    assert outcome.p == pytest.approx(0.1, abs=1e-12)
    assert "Mann–Whitney" in outcome.test


@requires_stats
def test_mann_whitney_matches_the_hand_computed_u_statistic() -> None:
    """a = [1, 2, 3, 4] against b = [2, 3, 4, 5] — pooled ranks with ties averaged.

    Pooled: 1 · 2,2 · 3,3 · 4,4 · 5 ⇒ ranks 1 · 2.5,2.5 · 4.5,4.5 · 6.5,6.5 · 8.
    R(a) = 1 + 2.5 + 4.5 + 6.5 = 14.5, so U(a) = 14.5 − 4·5/2 = 4.5.
    """
    outcome = mann_whitney([1, 2, 3, 4], [2, 3, 4, 5])
    assert outcome.statistic == pytest.approx(4.5)
    assert 0.0 < outcome.p < 1.0


def test_a_group_of_one_is_not_computable_and_says_why() -> None:
    """The mito panel's case: one sequence per group has nothing to rank against."""
    outcome = mann_whitney([1.0], [2.0])
    assert not outcome.computable
    assert outcome.status == NOT_COMPUTABLE
    assert outcome.p is None
    assert "at least 2 observations per group" in outcome.reason
    assert outcome.n == (1, 1)


def test_two_groups_are_not_a_kruskal_wallis() -> None:
    outcome = kruskal_wallis([1, 2, 3], [4, 5, 6])
    assert not outcome.computable
    assert "at least 3 groups" in outcome.reason


@requires_stats
def test_kruskal_wallis_of_three_separated_groups_is_the_hand_computed_h() -> None:
    """Groups [1,2,3] [4,5,6] [7,8,9]: R = 6, 15, 24, H = 12/(81·3) · 432 − 27 = 7.2."""
    outcome = kruskal_wallis([1, 2, 3], [4, 5, 6], [7, 8, 9])
    assert outcome.statistic == pytest.approx(7.2, abs=1e-12)
    assert outcome.p == pytest.approx(0.02732372244729253, abs=1e-12)
    assert outcome.n == (3, 3, 3)


@requires_stats
def test_spearman_matches_the_hand_computed_rank_correlation() -> None:
    """x = 1..5, y = [2, 1, 4, 3, 5]: Σd² = 4, ρ = 1 − 6·4/(5·24) = 0.8."""
    outcome = spearman([1, 2, 3, 4, 5], [2, 1, 4, 3, 5])
    assert outcome.coefficient == pytest.approx(0.8, abs=1e-12)
    assert outcome.n == 5


@requires_stats
def test_pearson_of_a_perfect_line_is_one() -> None:
    outcome = pearson([1, 2, 3], [2, 4, 6])
    assert outcome.coefficient == pytest.approx(1.0, abs=1e-12)
    # Three points on a line: the t statistic diverges, so the p-value is at the floor rather
    # than exactly zero — a distinction worth keeping rather than rounding away.
    assert outcome.p < 1e-6


def test_a_constant_series_has_no_correlation_to_report() -> None:
    outcome = spearman([1, 2, 3, 4], [5, 5, 5, 5])
    assert not outcome.computable
    assert "constant" in outcome.reason


def test_a_correlation_needs_three_pairs() -> None:
    outcome = pearson([1, 2], [2, 4])
    assert not outcome.computable
    assert "at least 3 pairs" in outcome.reason


def test_mismatched_correlation_input_is_a_named_error() -> None:
    with pytest.raises(ValueError, match="3 x-values against 2 y-values"):
        spearman([1, 2, 3], [1, 2])


def test_benjamini_hochberg_matches_the_hand_computed_adjustment() -> None:
    """p = [0.01, 0.02, 0.03, 0.04, 0.05], m = 5: min over j ≥ i of m·p(j)/j is 0.05 for all."""
    assert benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05]) == pytest.approx(
        [0.05] * 5, abs=1e-12
    )


def test_benjamini_hochberg_is_monotone_and_keeps_the_input_order() -> None:
    """Unsorted input: 0.04 → 0.04 (rank 2), 0.01 → 0.02 (rank 1), and monotonicity holds."""
    assert benjamini_hochberg([0.04, 0.01]) == pytest.approx([0.04, 0.02], abs=1e-12)
    many = benjamini_hochberg([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074])
    assert many == tuple(sorted(many)), "BH adjusted p-values must not decrease with rank"


def test_benjamini_hochberg_leaves_a_not_computable_blank() -> None:
    adjusted = benjamini_hochberg([0.01, None, 0.02])
    assert adjusted[1] is None
    assert adjusted[0] == pytest.approx(0.02, abs=1e-12)  # m = 2, not 3
    assert adjusted[2] == pytest.approx(0.02, abs=1e-12)


def test_bonferroni_is_the_hand_computed_product_and_caps_at_one() -> None:
    assert bonferroni([0.01, 0.02, 0.5]) == pytest.approx([0.03, 0.06, 1.0], abs=1e-12)
    assert bonferroni([None, None]) == (None, None)


def test_a_missing_stats_extra_is_a_named_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """The minimum install stays pure Python — and says which extra to add, not ImportError."""
    monkeypatch.setitem(sys.modules, "scipy.stats", None)
    with pytest.raises(ExportDependencyError, match="plastid-suite\\[stats\\]"):
        mann_whitney([1, 2, 3], [4, 5, 6])
