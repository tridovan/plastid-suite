"""Codon mathematics: the code table, RSCU, and the bias indices (offline, pure)."""

from __future__ import annotations

import pytest

from plastid.core.codon import (
    FAMILIES,
    SENSE_CODONS,
    SENSE_FAMILIES,
    STANDARD_CODE,
    TABLE_LABEL,
    cai,
    cbi,
    codon_counts,
    enc,
    fop,
    optimal_codons,
    relative_adaptiveness,
    rscu,
    translate,
    translation_table_label,
)
from plastid.run.defaults import load_defaults


def _uniform(per_codon: int = 1000) -> dict[str, int]:
    return {codon: per_codon for codon in SENSE_CODONS}


def _most_biased(per_codon: int = 1000) -> dict[str, int]:
    return {members[0]: per_codon for members in SENSE_FAMILIES.values()}


def test_the_code_table_is_the_standard_one() -> None:
    assert len(STANDARD_CODE) == 64
    assert sum(1 for amino in STANDARD_CODE.values() if amino == "*") == 3
    assert [STANDARD_CODE[codon] for codon in ("ATG", "TGG", "TAA", "TAG", "TGA")] == [
        "M",
        "W",
        "*",
        "*",
        "*",
    ]
    # degeneracy of the standard code — the family sums check depends on it
    assert len(FAMILIES["L"]) == len(FAMILIES["S"]) == len(FAMILIES["R"]) == 6
    assert len(FAMILIES["M"]) == len(FAMILIES["W"]) == 1
    assert sum(len(members) for members in SENSE_FAMILIES.values()) == 61


def test_translate_marks_stops_and_ignores_a_partial_codon() -> None:
    assert translate("ATGGCTTAA") == "MA*"
    assert translate("ATGG") == "M"  # trailing two bases are not a codon
    assert translate("ATGNNN") == "MX"


def test_codon_counts_drops_stops_and_partials() -> None:
    counts = codon_counts("ATGGCTTAA")
    assert counts["ATG"] == 1 and counts["GCT"] == 1
    assert sum(counts.values()) == 2  # the stop is not a sense codon


def test_rscu_sums_to_the_family_size_for_every_amino_acid() -> None:
    values = rscu(_uniform())
    for members in SENSE_FAMILIES.values():
        assert sum(values[codon] for codon in members) == pytest.approx(len(members))


def test_rscu_is_one_everywhere_under_uniform_usage() -> None:
    values = rscu(_uniform())
    assert all(value == pytest.approx(1.0) for value in values.values())


def test_rscu_of_a_family_with_no_observations_is_zero_not_undefined() -> None:
    values = rscu({codon: 0 for codon in SENSE_CODONS})
    assert all(value == 0.0 for value in values.values())


def test_enc_is_61_without_bias_and_20_at_the_extreme() -> None:
    # Wright 1990: 2 + 9/1/2 + 1/1/3 + 5/1/4 + 3/1/6 == 61 with even usage; 2 + 9 + 1 + 5 + 3 == 20
    assert enc(_uniform()) == pytest.approx(61.0, abs=0.05)
    assert enc(_most_biased()) == pytest.approx(20.0, abs=0.05)


def test_enc_is_none_when_a_degeneracy_class_has_no_family() -> None:
    counts = {codon: 0 for codon in SENSE_CODONS}
    for members in SENSE_FAMILIES.values():
        if len(members) == 2:
            counts[members[0]] = 100
    assert enc(counts) is None  # no four-fold family carries data


def test_adaptiveness_and_optimal_codons_follow_the_reference() -> None:
    reference = {"GCT": 90, "GCC": 10, "GCA": 0, "GCG": 0}
    weights = relative_adaptiveness(reference)
    assert weights["GCT"] == pytest.approx(1.0)
    assert weights["GCC"] == pytest.approx(1 / 9)
    assert weights["GCA"] == 0.0
    assert optimal_codons(reference) == frozenset({"GCT"})

    tied = optimal_codons({"GCT": 50, "GCC": 50, "GCA": 1, "GCG": 0})
    assert tied == frozenset({"GCT", "GCC"})  # a tie is deterministic, not arbitrary


def test_indices_are_one_when_the_sequence_is_its_own_reference() -> None:
    """A sequence adapted to itself: every codon it uses is the best one available."""
    counts = _most_biased()
    assert cai(counts, counts) == pytest.approx(1.0)
    assert fop(counts, counts) == pytest.approx(1.0)
    assert cbi(counts, counts) == pytest.approx(1.0)


def test_cai_skips_codons_the_reference_never_uses() -> None:
    # CAU is absent from the reference, so a sequence using it is not adapted to that reference
    reference = {"CAT": 90, "CAC": 10}
    assert cai({"CAT": 10, "CAC": 0}, reference) == pytest.approx(1.0)
    assert cai({"CAT": 0, "CAC": 10}, reference) == pytest.approx(1 / 9)


def test_indices_are_none_without_data() -> None:
    empty = {codon: 0 for codon in SENSE_CODONS}
    assert cai(empty, empty) is None
    assert fop(empty, empty) is None
    assert cbi(empty, empty) is None


def test_the_recorded_decision_matches_the_compiled_table() -> None:
    """defaults.yaml records D2-03; this keeps the record and core/codon.py from drifting."""
    defaults = load_defaults()["codon"]
    assert defaults["translation_table"] == 1
    assert defaults["reference"] == "self"
    assert TABLE_LABEL.startswith("NCBI table 1")


def test_the_table_label_names_the_tables_the_record_declares() -> None:
    """Measured 05/10/2026 (§9.J): the report said "NCBI table 1" for records whose CDS declare
    11, so the column read as if the toolkit had used the record's own table. The label now
    carries both — and silence (no CDS declares one) changes nothing."""
    from plastid.core.model import Gene, Interval

    def cds(table: int | None) -> Gene:
        return Gene("rbcL", "CDS", (Interval(0, 30),), transl_table=table)

    assert translation_table_label([cds(11), cds(11)]) == (
        "NCBI table 1 (standard); the record's CDS declare 11"
    )
    assert translation_table_label([cds(None)]) == TABLE_LABEL
    assert translation_table_label([]) == TABLE_LABEL
