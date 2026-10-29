"""``plastid run`` — what a failed step does to the steps behind it (M11.16)."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from plastid.cli import main
from plastid.modules import registry

DATA = Path(__file__).parents[1] / "data"


def test_a_failed_step_marks_the_later_steps_skipped(tmp_path: Path, monkeypatch) -> None:
    """A failure stops the remaining modules; the manifest must say they did not run.

    Measured 27/10/2026: the loop returned on the first failure and the later specs were simply
    absent from ``modules``, so "the run never reached them" and "they were never asked for"
    looked identical in the one file a reader checks.
    """
    first = replace(
        registry.get("genome_stats"),
        run=lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("boom")),
        requires_ir=False,
    )
    second = replace(registry.get("codon_usage"), requires_ir=False)
    monkeypatch.setattr(
        registry,
        "get",
        lambda module_id: second if module_id == "codon_usage" else first,
    )

    code = main(
        [
            "run",
            "-m",
            "genome_stats",
            "-m",
            "codon_usage",
            "-i",
            str(DATA / "NC_050833.1.gb"),
            "-o",
            str(tmp_path / "runs"),
        ]
    )

    assert code == 1
    manifest = json.loads(
        next((tmp_path / "runs").glob("*/*/manifest.json")).read_text(encoding="utf-8")
    )
    assert [module["status"] for module in manifest["modules"]] == ["failed", "skipped"]
    skipped = manifest["modules"][1]
    assert skipped["id"] == "codon_usage"
    assert "genome_stats failed" in skipped["reason"]
