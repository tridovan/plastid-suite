"""A figure that should have been drawn and was not stops a command from reading as a success.

Measured 01/11/2026: `plastid run --figures` lost three R figures, printed nothing about them and
exited 0, so the only trace was a word in the manifest. These tests stand a small script in for
Rscript, so they fail the same way on a machine with R and on one without it.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from plastid.cli import EXIT_ERROR, EXIT_TOOL, main

FAILING_R = """#!/bin/sh
case "$1" in
  --version) echo "Rscript (R) version 4.3.3 (2024-02-29)"; exit 0 ;;
  -e) exit 0 ;;
esac
echo "=== the script's banner, on stdout ==="
echo "Error: the script could not read its input" >&2
echo "Execution halted" >&2
exit 1
"""


@pytest.fixture
def failing_r(tmp_path: Path) -> Path:
    path = tmp_path / "Rscript"
    path.write_text(FAILING_R, encoding="utf-8")
    path.chmod(0o755)
    return path


def _run_dir(out: Path) -> Path:
    return next(out.glob("*/*/manifest.json")).parent


def _ran(out: Path, record: Path) -> Path:
    assert main(["run", "-m", "codon_usage", "-i", str(record), "-o", str(out)]) == 0
    return _run_dir(out)


def test_run_with_figures_names_the_failure_and_says_so_in_its_exit_code(
    tmp_path: Path, data_dir: Path, failing_r: Path, capsys
) -> None:
    pytest.importorskip("matplotlib", reason="the 'figures' extra is not installed")
    out = tmp_path / "runs"
    code = main(
        [
            "run",
            "-m",
            "codon_usage",
            "-i",
            str(data_dir / "NC_050833.1.gb"),
            "-o",
            str(out),
            "--figures",
            "--set",
            f"rscript={failing_r}",
        ]
    )
    assert code == EXIT_TOOL
    said = capsys.readouterr().err
    assert "r_rscu" in said and "could not read its input" in said

    run_dir = _run_dir(out)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "tool_error"
    assert manifest["figures"]["status"] == "failed"
    assert {entry["figure"] for entry in manifest["figures"]["failed"]} == {"r_rscu"}
    # the analysis finished: only the figure is missing, and the tables stay in the run
    assert (run_dir / "01_codon_usage" / "tables" / "rscu.csv").is_file()


def test_figures_exits_with_the_tool_code_when_one_figure_of_several_fails(
    tmp_path: Path, data_dir: Path, failing_r: Path, capsys
) -> None:
    run_dir = _ran(tmp_path / "runs", data_dir / "NC_050833.1.gb")
    capsys.readouterr()
    code = main(["figures", str(run_dir), "--figure", "r_rscu", "--set", f"rscript={failing_r}"])
    assert code == EXIT_TOOL
    assert "could not read its input" in capsys.readouterr().err


def test_a_table_the_run_does_not_hold_is_not_a_tool_failure(
    tmp_path: Path, data_dir: Path, failing_r: Path
) -> None:
    """Nothing was going to be drawn from a table the run does not hold, so R never ran."""
    run_dir = _ran(tmp_path / "runs", data_dir / "NC_050833.1.gb")
    code = main(["figures", str(run_dir), "--figure", "r_kaks", "--set", f"rscript={failing_r}"])
    assert code == EXIT_ERROR


def test_a_batch_that_loses_a_figure_says_so_in_its_exit_code(
    tmp_path: Path, data_dir: Path, failing_r: Path, capsys
) -> None:
    out = tmp_path / "runs"
    code = main(
        [
            "batch",
            "-m",
            "codon_usage",
            "-i",
            str(data_dir / "NC_050833.1.gb"),
            "-i",
            str(data_dir / "NC_024936.1.gb"),
            "-o",
            str(out),
            "--dataset",
            "pair",
            "--figures",
            "--set",
            f"rscript={failing_r}",
        ]
    )
    assert code == EXIT_TOOL
    assert "could not read its input" in capsys.readouterr().err
    payload = json.loads(next((out / "pair").glob("*/batch.json")).read_text(encoding="utf-8"))
    assert payload["figures"]["status"] == "failed"
    assert payload["exit_code"] == EXIT_TOOL


def test_a_missing_plotting_extra_is_not_the_last_analysis_failing(
    tmp_path: Path, data_dir: Path, monkeypatch: pytest.MonkeyPatch, capsys
) -> None:
    """The analyses finished. Their record says so, and the figures say what stopped them."""
    # The extra absent: neither matplotlib nor the module the export draws with can be imported.
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    monkeypatch.setitem(sys.modules, "matplotlib.figure", None)
    out = tmp_path / "runs"
    record = str(data_dir / "NC_050833.1.gb")
    code = main(["run", "-m", "codon_usage", "-i", record, "-o", str(out), "--figures"])
    assert code == EXIT_TOOL
    assert "plastid-suite[figures]" in capsys.readouterr().err

    manifest = json.loads((_run_dir(out) / "manifest.json").read_text(encoding="utf-8"))
    assert [(step["id"], step["status"]) for step in manifest["modules"]] == [("codon_usage", "ok")]
    assert manifest["status"] == "tool_error"
    assert manifest["figures"]["status"] == "failed"
    assert "plastid-suite[figures]" in manifest["figures"]["error"]
