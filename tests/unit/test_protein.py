"""Amino-acid composition measures — analytic checks on the Kyte–Doolittle scale (offline)."""

from __future__ import annotations

import pytest

from plastid.core.protein import (
    AROMATIC,
    KYTE_DOOLITTLE,
    RESIDUES,
    aromaticity,
    composition,
    gravy,
    ignored_symbols,
)


def test_the_scale_covers_the_twenty_standard_residues() -> None:
    assert len(RESIDUES) == 20
    assert set(RESIDUES) == set("ACDEFGHIKLMNPQRSTVWY")
    assert KYTE_DOOLITTLE["I"] == 4.5  # the most hydrophobic, per Kyte & Doolittle
    assert KYTE_DOOLITTLE["R"] == -4.5  # the most hydrophilic


def test_composition_counts_only_residues() -> None:
    counts = composition("AA*XG-")
    assert counts["A"] == 2 and counts["G"] == 1
    assert sum(counts.values()) == 3  # stop, X and the gap are not residues


def test_ignored_symbols_are_counted_not_dropped() -> None:
    assert ignored_symbols("AA*XG-") == 3
    assert ignored_symbols("ACDEFG") == 0


def test_gravy_of_a_single_residue_is_its_own_hydropathy() -> None:
    assert gravy({"I": 10}) == pytest.approx(4.5)
    assert gravy({"R": 10}) == pytest.approx(-4.5)


def test_gravy_is_the_mean_over_residues() -> None:
    # (4.5 + (-4.5) + 1.8) / 3
    assert gravy({"I": 1, "R": 1, "A": 1}) == pytest.approx(0.6)


def test_aromaticity_is_the_share_of_f_w_y() -> None:
    assert frozenset({"F", "W", "Y"}) == AROMATIC
    assert aromaticity({"F": 1, "W": 1, "Y": 2, "A": 96}) == pytest.approx(0.04)
    assert aromaticity({"A": 10}) == 0.0


def test_measures_are_none_without_residues() -> None:
    empty = {residue: 0 for residue in RESIDUES}
    assert gravy(empty) is None
    assert aromaticity(empty) is None
