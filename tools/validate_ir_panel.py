"""Validate the pinned IR panel against its independent circular-arc measurements.

The panel lives outside the repository (``.panel/ir-validation/``, gitignored). Its
``measurements/identity.json`` was measured with joined copies concatenated in join-string order;
for a copy that crosses the origin that order is not the circular arc, and it produced two wrong
numbers (NC_054177.1, NC_058533.1). ``measurements/identity_circular.json`` is the corrected
measurement used here.

Run from the product repository root:

    uv run python tools/validate_ir_panel.py
    uv run python tools/validate_ir_panel.py --panel .panel/ir-validation

The script never touches the network and never writes into the panel. Exit codes: 0 when every
record matches its expected length and verdict (scoring-only cases are reported, not failed),
1 on a real disagreement, 2 when the panel or its measurement file is missing.

This is a validation tool, not part of the minimum install. ``parasail`` is not a runtime
dependency; the expected JSON records the independent run that produced it.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, ir_from_annotation
from plastid.core.normalize import normalize

__all__ = ["Finding", "classify", "main", "validate_panel"]

IR_IDENTITY_TOLERANCE = 0.002  # 0.2 percentage points; scoring changes the fourth decimal
PP078674_REGIONS = (("LSC", 83_820), ("IRa", 26_331), ("SSC", 17_195), ("IRb", 26_330))


@dataclass(frozen=True, slots=True)
class Finding:
    accession: str
    kind: str
    detail: str


def classify(
    accession: str,
    *,
    pts_length: int | None,
    pts_identity: float | None,
    pts_passed: bool | None,
    expected: dict,
    scoring_cases: dict,
) -> list[Finding]:
    """Compare one record; scoring-only cases are reported by the caller, not failed here."""
    if "error" in expected:
        return [Finding(accession, "expected-error", str(expected["error"]))]
    findings: list[Finding] = []
    expected_length = min(expected["len_a"], expected["len_b"])
    if pts_length != expected_length:
        findings.append(
            Finding(accession, "length", f"PTS {pts_length} bp, annotation {expected_length} bp")
        )
    delta = float(pts_identity) - float(expected["identity"])
    if abs(delta) > IR_IDENTITY_TOLERANCE and accession not in scoring_cases:
        findings.append(
            Finding(
                accession,
                "identity",
                f"PTS {pts_identity:.6f}, independent {expected['identity']:.6f} "
                f"(delta {delta:+.6f})",
            )
        )
    expected_passed = float(expected["identity"]) >= 0.999
    if pts_passed != expected_passed and accession not in scoring_cases:
        findings.append(
            Finding(
                accession,
                "verdict",
                f"PTS {'pass' if pts_passed else 'fail'}, independent "
                f"{'pass' if expected_passed else 'fail'}",
            )
        )
    return findings


def validate_panel(panel: Path) -> tuple[list[Finding], list[str], dict[str, int]]:
    """Return ``(findings, scoring_notes, counts)``; raise OSError when the panel is unusable."""
    measurement_path = panel / "measurements" / "identity_circular.json"
    payload = json.loads(measurement_path.read_text(encoding="utf-8"))
    records: dict[str, dict] = payload["records"]
    scoring = payload.get("scoring_cases") or {}
    findings: list[Finding] = []
    notes: list[str] = []
    counts = {"records": 0, "passed": 0, "failed": 0, "lengths": 0}

    for accession, expected in sorted(records.items()):
        path = panel / "source" / f"{accession}.gb"
        if not path.is_file():
            findings.append(Finding(accession, "missing", str(path)))
            continue
        counts["records"] += 1
        try:
            genome = read_genbank(path)
            ir = ir_from_annotation(genome)
        except (OSError, ValueError, IRDetectionError) as error:
            findings.append(Finding(accession, "read", f"{type(error).__name__}: {error}"))
            continue
        if ir is None:
            findings.append(Finding(accession, "no-annotation", "no two IR copies found"))
            continue
        if ir.length == min(expected["len_a"], expected["len_b"]):
            counts["lengths"] += 1
        if ir.gate_passed:
            counts["passed"] += 1
        else:
            counts["failed"] += 1
        findings.extend(
            classify(
                accession,
                pts_length=ir.length,
                pts_identity=ir.identity,
                pts_passed=ir.gate_passed,
                expected=expected,
                scoring_cases=scoring,
            )
        )
        if accession in scoring:
            notes.append(f"{accession}: scoring-only — {scoring[accession]}")

    extras = payload.get("extras") or {}
    pp = extras.get("PP078674.1")
    if pp is not None:
        pp_path = panel / "pp078674" / "PP078674.1.gb"
        if not pp_path.is_file():
            findings.append(Finding("PP078674.1", "missing", str(pp_path)))
        else:
            try:
                genome = read_genbank(pp_path)
                ir = ir_from_annotation(genome)
                if ir is None:
                    raise IRDetectionError("PP078674.1 has no two IR copies")
                regions = tuple(
                    (region.name, region.length) for region in normalize(genome, ir).regions
                )
            except (OSError, ValueError, IRDetectionError) as error:
                findings.append(Finding("PP078674.1", "read", f"{type(error).__name__}: {error}"))
            else:
                if regions != PP078674_REGIONS:
                    findings.append(
                        Finding(
                            "PP078674.1",
                            "regions",
                            f"{regions!r}, expected {PP078674_REGIONS!r}",
                        )
                    )
                if ir.length != min(pp["len_a"], pp["len_b"]):
                    findings.append(
                        Finding(
                            "PP078674.1",
                            "length",
                            f"{ir.length} bp, expected {pp['len_a']}/{pp['len_b']}",
                        )
                    )
                if not ir.gate_passed:
                    findings.append(Finding("PP078674.1", "verdict", "IR gate failed"))
    return findings, notes, counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        default=".panel/ir-validation",
        help="pinned panel directory (default: .panel/ir-validation)",
    )
    args = parser.parse_args(argv)
    panel = Path(args.panel)
    if not (panel / "measurements" / "identity_circular.json").is_file():
        reason = (
            f"IR panel: {panel} has no measurements/identity_circular.json — pin the panel first; "
            "this validation is not part of the minimum install"
        )
        if os.environ.get("PTS_REQUIRE_PANEL") == "1":
            print(reason, file=sys.stderr)
            return 2
        print(f"skipped: {reason}")
        return 0
    try:
        findings, notes, counts = validate_panel(panel)
    except (OSError, ValueError, KeyError) as error:
        print(f"IR panel: cannot read {panel}: {error}", file=sys.stderr)
        return 2

    for note in notes:
        print(f"note: {note}")
    for finding in findings:
        print(f"FAIL {finding.accession}: {finding.kind} — {finding.detail}")
    print(
        f"source records: {counts['records']} · pass {counts['passed']} · "
        f"fail {counts['failed']} · length match {counts['lengths']}/{counts['records']}"
    )
    if findings:
        print(f"{len(findings)} disagreement(s)", file=sys.stderr)
        return 1
    print("IR panel: ok (scoring-only cases listed above)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
