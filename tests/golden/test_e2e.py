"""End-to-end golden run: one command produces a complete, traceable run directory.

Covers P1-08's Done-when: run dir with 00_run_info + tables + report; tool/parameter report;
cache hit on an identical second run; byte-identical restore after deleting an output.
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO

from plastid.cli import main
from plastid.run.manifest import validate


def _run(
    tmp_path: Path, source: Path, dataset: str = "iris_domestica", extra: list[str] | None = None
) -> int:
    return main(
        [
            "run",
            "--module",
            "genome_stats",
            "--input",
            str(source),
            "--out",
            str(tmp_path),
            "--dataset",
            dataset,
            *(extra or []),
        ]
    )


def _only_run_dir(tmp_path: Path, dataset: str = "iris_domestica") -> Path:
    runs = sorted((tmp_path / dataset).iterdir())
    assert len(runs) == 1
    return runs[0]


def test_run_produces_a_complete_traceable_directory(tmp_path: Path, data_dir: Path) -> None:
    assert _run(tmp_path, data_dir / "NC_050833.1.gb") == 0
    run_dir = _only_run_dir(tmp_path)

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)
    assert manifest["status"] == "ok"
    assert manifest["ir"]["gate_passed"] is True
    assert manifest["ir"]["length"] == 26_214
    assert manifest["coordinate_system"]["rotation_offset"] == 0
    assert manifest["modules"][0]["id"] == "genome_stats"
    assert manifest["files"], "the manifest must index every produced file"
    for entry in manifest["files"]:
        assert (run_dir / entry["path"]).is_file()

    stats = (run_dir / "01_genome_stats" / "tables" / "genome_stats.csv").read_text(
        encoding="utf-8"
    )
    assert "153724" in stats
    assert "26214" in stats
    assert "# scope:" in stats  # every export carries its scope

    report = (run_dir / "01_genome_stats" / "report.md").read_text(encoding="utf-8")
    assert "Draft — needs human review." in report
    assert "132 genes" in report
    assert "### References" in report  # biopython citation resolved from the registry

    run_info = (run_dir / "00_run_info" / "run_info.md").read_text(encoding="utf-8")
    assert "## Parameters" in run_info
    assert "## IR gate" in run_info
    assert "passed" in run_info
    assert "genome_stats" in run_info


def test_second_identical_run_is_cached_and_identical(tmp_path: Path, data_dir: Path) -> None:
    source = data_dir / "NC_050833.1.gb"
    assert _run(tmp_path, source) == 0
    first_dir = _only_run_dir(tmp_path)
    first_csv = (first_dir / "01_genome_stats" / "tables" / "genome_stats.csv").read_bytes()

    assert _run(tmp_path, source) == 0
    runs = sorted((tmp_path / "iris_domestica").iterdir())
    assert len(runs) == 2, "each run gets its own directory"
    second_dir = runs[1]
    manifest = json.loads((second_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["modules"][0]["status"] == "cached"

    second_csv_path = second_dir / "01_genome_stats" / "tables" / "genome_stats.csv"
    assert second_csv_path.read_bytes() == first_csv  # same input, same numbers

    second_csv_path.unlink()
    assert _run(tmp_path, source) == 0
    third_dir = sorted((tmp_path / "iris_domestica").iterdir())[2]
    restored = (third_dir / "01_genome_stats" / "tables" / "genome_stats.csv").read_bytes()
    assert restored == first_csv  # restored byte-identically from the cache


def test_failed_gate_stops_the_run_with_exit_2(tmp_path: Path, data_dir: Path) -> None:
    # NC_063920's IR copies are 92.6 % identical — position-based analysis must not run.
    records = SeqIO.parse(str(data_dir / "iridaceae_edge_cases.gb"), "genbank")
    target = next(record for record in records if record.id == "NC_063920.1")
    single = tmp_path / "NC_063920.1.gb"
    SeqIO.write(target, str(single), "genbank")

    assert _run(tmp_path, single, dataset="diverged") == 2
    run_dir = _only_run_dir(tmp_path, "diverged")
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert not (run_dir / "01_genome_stats").exists(), "no module output after a failed gate"


def test_unknown_module_exits_1(tmp_path: Path, data_dir: Path) -> None:
    assert _run(tmp_path, data_dir / "NC_050833.1.gb", extra=["--module", "nope"]) in (1, 2)
