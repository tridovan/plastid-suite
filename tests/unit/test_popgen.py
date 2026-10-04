"""The population-genetic estimators against hand-computed vectors (P4-03).

Every expected value here is arithmetic done by hand from the published formula, written out in
the comment beside it — not a number copied from a previous run of this code. The cross-check
against an independent implementation (EggLib 3.6.1) is a separate, recorded measurement
(PROJECT.md §8.28), because a test that asserts what the code already produces proves nothing.

The vector used throughout: four sequences, ten sites, one segregating site where the odd base
appears once.

    s1  AAAAAAAAAA
    s2  AAAAACAAAA
    s3  AAAAACAAAA
    s4  AAAAACAAAA
"""

from __future__ import annotations

import pytest

from plastid.core.align import Alignment
from plastid.core.popgen import (
    fst_gst,
    fu_li_d_star,
    fu_li_f_star,
    haplotype_diversity,
    segregating_sites,
    summarise,
    tajimas_d,
    watterson_theta,
)

SAMPLE = Alignment(
    names=("s1", "s2", "s3", "s4"),
    sequences=(
        "AAAAAAAAAA",
        "AAAAACAAAA",
        "AAAAACAAAA",
        "AAAAACAAAA",
    ),
)


def test_haplotype_diversity_matches_the_hand_computation() -> None:
    # two haplotypes, frequencies 1/4 and 3/4 → Hd = 4/3 · (1 − (1/16 + 9/16)) = 4/3 · 6/16 = 1/2
    assert haplotype_diversity(SAMPLE) == pytest.approx(0.5)


def test_segregating_sites_counts_only_columns_with_two_bases() -> None:
    assert segregating_sites(SAMPLE) == 1
    assert segregating_sites(Alignment(names=("a", "b"), sequences=("AAAA", "AAAA"))) == 0


def test_watterson_theta_matches_the_hand_computation() -> None:
    # n = 4 → a1 = 1 + 1/2 + 1/3 = 11/6; S = 1 → θ per site = 1 / (11/6) / 10 = 6/110
    assert watterson_theta(SAMPLE) == pytest.approx(6 / 110)


def test_tajimas_d_matches_the_hand_computation() -> None:
    # π_total = 0.5 (one site at He = 4/3·2·(3/4)(1/4) = 1/2), θ_total = 6/11,
    # c1 = 5/9 − 6/11 = 1/99, e1 = c1/a1 = 6/1089, V = e1·S = 6/1089
    expected = (0.5 - 6 / 11) / (6 / 1089) ** 0.5
    assert expected == pytest.approx(-0.612372, abs=1e-6)
    assert tajimas_d(SAMPLE) == pytest.approx(expected)


def test_the_neutrality_tests_refuse_n_below_four_and_say_why() -> None:
    """At n = 3 the variance is identically zero — measured, not assumed: with n = 3,
    c1 = b1 − 1/a1 = 2/3 − 2/3 = 0 and c2 = 0 as well, so the statistic is 0/0."""
    three = Alignment(
        names=("s1", "s2", "s3"),
        sequences=("AAAAAAAAAA", "AAAAACAAAA", "AAAAACAAAA"),
    )
    summary = summarise(three)
    assert summary.tajimas_d is None
    assert summary.fu_li_d_star is None
    assert summary.fu_li_f_star is None
    assert summary.haplotype_diversity is not None  # still estimable
    assert any("at least 4 sequences" in note for note in summary.notes)


def test_a_constant_sample_has_zero_diversity_and_no_neutrality_test() -> None:
    constant = Alignment(
        names=("s1", "s2", "s3", "s4"),
        sequences=("ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT"),
    )
    summary = summarise(constant)
    assert summary.segregating == 0
    assert summary.haplotypes == 1
    assert summary.haplotype_diversity == 0  # a measured zero: every sequence is the same
    assert summary.pi_per_site == 0
    assert summary.theta_w_per_site == 0
    assert summary.tajimas_d is None  # zero variance — undefined, and said so
    assert any("no segregating site" in note for note in summary.notes)


def test_fu_and_li_d_star_matches_the_hand_computation() -> None:
    """The same four-sequence vector, computed by hand from Fu & Li (1993) — every step written
    out, with the variance constants of Simonsen, Churchill & Aquadro (1995).

    n = 4 → a1 = 11/6, a2 = 49/36; η = 1 (one mutation), η_e = 1 (it is a singleton).
    c_n = 2(4·11/6 − 6)/(3·2) = 4/9;  d_n = 4/9 + 2/9 + (2/3)(3/2 − (7/6)/2 − 1/4) = 10/9
    v_D = (16/9·49/36 + 121/36·10/9 − 8·11/6·17/6/9) / (121/36 + 49/36) = (83/54)/(85/18) = 83/255
    u_D = 4/3·(11/6 − 4/3) − 83/255 = 2/3 − 83/255 = 29/85
    D* = (4/3·1 − 11/6·1) / sqrt(29/85 + 83/255) = (−1/2)/sqrt(2/3) = −0.612372…
    """
    expected = (-0.5) / (2 / 3) ** 0.5
    assert expected == pytest.approx(-0.6123724356957945)
    assert fu_li_d_star(SAMPLE) == pytest.approx(expected)


def test_fu_and_li_f_star_matches_the_hand_computation() -> None:
    """F* on the same vector: π_total = 1/2, η = η_e = 1, n = 4.

    v_F = (1021/432 + 11/16 − 49/18) / (85/18) = (71/216)·(18/85) = 71/1020
    u_F = (18/36)/a1 − v_F = 3/11 − 71/1020
    F* = (1/2 − 3/4·1) / sqrt(u_F + v_F) = (−1/4)/sqrt(0.272727…) = −0.478713…
    """
    expected = (-0.25) / (3 / 11) ** 0.5
    assert expected == pytest.approx(-0.4787135538781689)
    assert fu_li_f_star(SAMPLE) == pytest.approx(expected)


def test_the_neutrality_statistics_move_with_the_data() -> None:
    """A guard against a statistic that quietly returns a constant: the same mutation carried by
    **two** sequences instead of one is a different frequency spectrum — and, unlike moving the
    mutation to another sequence, it must change the answer."""
    other = Alignment(
        names=("s1", "s2", "s3", "s4"),
        sequences=("ACGTACGTAC", "ACGTACGTTC", "ACGTACGTTC", "ACGTACGTAC"),
    )
    assert tajimas_d(other) != pytest.approx(tajimas_d(SAMPLE))
    assert fu_li_d_star(other) != pytest.approx(fu_li_d_star(SAMPLE))
    assert fu_li_f_star(other) != pytest.approx(fu_li_f_star(SAMPLE))


def test_gap_columns_are_dropped_not_scored() -> None:
    """Complete deletion: a column with a gap in any row is not a site for any statistic."""
    gapped = Alignment(
        names=("s1", "s2", "s3", "s4"),
        sequences=(
            "AAAA-AAAAA",
            "AAAAACAAAA",
            "AAAAACAAAA",
            "AAAAACAAAA",
        ),
    )
    summary = summarise(gapped)
    assert summary.sites == 9  # the gapped column is gone
    assert summary.segregating == 1
    assert summary.haplotype_diversity == pytest.approx(0.5)


def test_fst_is_zero_between_two_groups_with_the_same_composition() -> None:
    """The analytic check the board names: same haplotypes in both groups → no differentiation.

    H_S = 1 − (0.5² + 0.5²) = 0.5 in each group, H_T = the same on the pooled four sequences →
    Fst = (0.5 − 0.5)/0.5 = 0 exactly.
    """
    group_a = Alignment(names=("a1", "a2"), sequences=("AAAAACAAAA", "AAAAAAAAAA"))
    group_b = Alignment(names=("b1", "b2"), sequences=("AAAAACAAAA", "AAAAAAAAAA"))
    assert fst_gst((group_a, group_b)) == pytest.approx(0.0)


def test_fst_is_one_when_the_groups_share_no_haplotype() -> None:
    # H_S = 0 (each group is constant), H_T = 1 − (0.5² + 0.5²) = 0.5 → Fst = 1
    group_a = Alignment(names=("a1", "a2"), sequences=("AAAAAAAAAA", "AAAAAAAAAA"))
    group_b = Alignment(names=("b1", "b2"), sequences=("AAAAACAAAA", "AAAAACAAAA"))
    assert fst_gst((group_a, group_b)) == pytest.approx(1.0)


def test_fst_is_none_when_there_is_no_diversity_to_partition() -> None:
    """Every sequence identical → H_T = 0 → None. It must never come out as 1 or as 0."""
    same = Alignment(names=("a1", "a2"), sequences=("AAAAAAAAAA", "AAAAAAAAAA"))
    assert fst_gst((same, same)) is None


def test_fst_refuses_a_group_with_a_single_sequence() -> None:
    """With one sequence in a group the within-group diversity is zero by construction, so the
    estimator returns exactly 1 for any two groups that differ — a number that looks like
    complete differentiation. It is not computable, and says so."""
    singleton = Alignment(names=("a1",), sequences=("AAAAAAAAAA",))
    other = Alignment(names=("b1",), sequences=("AAAAACAAAA",))
    assert fst_gst((singleton, other)) is None
    mixed = Alignment(names=("b1", "b2"), sequences=("AAAAACAAAA", "AAAAACAAAA"))
    assert fst_gst((singleton, mixed)) is None


def test_fst_needs_at_least_two_groups() -> None:
    one = Alignment(names=("a1", "a2"), sequences=("AAAAACAAAA", "AAAAAAAAAA"))
    assert fst_gst((one,)) is None


def test_summarise_publishes_n_and_the_site_count_beside_every_value() -> None:
    summary = summarise(SAMPLE)
    assert (summary.n, summary.sites) == (4, 10)
    assert summary.singletons == 1  # the odd base occurs once
    assert summary.haplotypes == 2
