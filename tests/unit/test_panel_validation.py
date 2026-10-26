"""The IR panel validator's comparison rule — no panel needed, so it runs everywhere."""

from __future__ import annotations

from tools.validate_ir_panel import classify


def _expected(identity: float, length: int = 1_000) -> dict:
    return {"len_a": length, "len_b": length, "identity": identity}


def test_a_matching_record_has_no_finding() -> None:
    assert (
        classify(
            "NC_000001.1",
            pts_length=1_000,
            pts_identity=0.9999,
            pts_passed=True,
            expected=_expected(0.9999),
            scoring_cases={},
        )
        == []
    )


def test_a_length_difference_is_reported() -> None:
    findings = classify(
        "NC_000001.1",
        pts_length=999,
        pts_identity=1.0,
        pts_passed=True,
        expected=_expected(1.0),
        scoring_cases={},
    )
    assert [finding.kind for finding in findings] == ["length"]


def test_a_verdict_difference_is_reported() -> None:
    findings = classify(
        "NC_000001.1",
        pts_length=1_000,
        pts_identity=0.9995,
        pts_passed=True,
        expected=_expected(0.9985),
        scoring_cases={},
    )
    assert [finding.kind for finding in findings] == ["verdict"]


def test_a_documented_scoring_case_is_not_a_finding() -> None:
    findings = classify(
        "NC_035154.1",
        pts_length=1_000,
        pts_identity=0.9995,
        pts_passed=True,
        expected=_expected(0.9985),
        scoring_cases={"NC_035154.1": "14 SNPs vs cheap-gap alignment"},
    )
    assert findings == []


def test_a_broken_expected_record_is_reported() -> None:
    findings = classify(
        "NC_000001.1",
        pts_length=None,
        pts_identity=None,
        pts_passed=None,
        expected={"error": "no two copies"},
        scoring_cases={},
    )
    assert [finding.kind for finding in findings] == ["expected-error"]
