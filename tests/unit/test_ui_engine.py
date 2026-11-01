"""One engine, two ways in: a run started from the interface is a run the command line made."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.ui.app import App, Refused

MODULES = ["genome_stats", "codon_usage"]
LABEL = "NCBI RefSeq, PTS test set"


def _workspace(tmp_path: Path, data_dir: Path) -> App:
    app = App(tmp_path / "work", version="test")
    shutil.copy(data_dir / "NC_050833.1.gb", app.workspace / "records" / "NC_050833.1.gb")
    return app


def _finished(app: App, job_id: str) -> dict:
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        job = app.job(job_id)
        if not job["running"]:
            return job
        time.sleep(0.2)
    raise AssertionError("the job did not finish")


def _form() -> dict:
    return {
        "modules": MODULES,
        "inputs": ["records/NC_050833.1.gb"],
        "dataset": "same",
        "label": LABEL,
    }


def _tables(run_dir: Path) -> dict[str, bytes]:
    return {
        path.relative_to(run_dir).as_posix(): path.read_bytes()
        for path in sorted(run_dir.glob("*/tables/*.csv"))
    }


def test_the_tables_are_the_bytes_a_prompt_writes(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    job = _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    assert job["status"] == "ok"
    from_page = next((app.runs_dir / "same").glob("*/manifest.json")).parent

    out = tmp_path / "prompt"
    argv = ["run", "-m", MODULES[0], "-m", MODULES[1]]
    argv += ["-i", str(data_dir / "NC_050833.1.gb"), "-o", str(out)]
    argv += ["--dataset", "same", "--label", LABEL]
    assert main(argv) == 0
    from_prompt = next((out / "same").glob("*/manifest.json")).parent

    tables = _tables(from_page)
    assert tables and tables == _tables(from_prompt)


def test_the_line_beside_the_form_is_the_line_the_run_records(
    tmp_path: Path, data_dir: Path
) -> None:
    app = _workspace(tmp_path, data_dir)
    shown = app.command({"command": "run", "form": _form()})["args"]
    _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    manifest = next((app.runs_dir / "same").glob("*/manifest.json"))
    recorded = json.loads(manifest.read_text(encoding="utf-8"))["command"]
    assert recorded == ["plastid", *shown]


def test_a_run_from_the_prompt_reads_like_any_other(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    record = str(app.workspace / "records" / "NC_050833.1.gb")
    argv = ["run", "-m", "genome_stats", "-i", record, "-o", str(app.runs_dir), "--dataset", "cli"]
    assert main(argv) == 0
    rows = app.runs()["runs"]
    assert [row["dataset"] for row in rows] == ["cli"] and rows[0]["status"] == "ok"
    detail = app.run("cli", rows[0]["run_id"])
    assert [step["id"] for step in detail["steps"]] == ["genome_stats"]
    table = app.table(detail["steps"][0]["tables"][1]["path"])
    assert table["scope"] and table["columns"]


def test_a_refused_form_starts_nothing(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    with pytest.raises(Refused, match="at least one analysis"):
        app.start_job({"command": "run", "form": {"inputs": ["records/NC_050833.1.gb"]}})
    assert app.jobs() == {"jobs": []}


def test_a_run_that_is_not_running_cannot_be_cancelled(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    run_id = app.runs()["runs"][0]["run_id"]
    with pytest.raises(Refused, match="not running"):
        app.cancel_run({}, "same", run_id)


def test_the_gate_refusal_is_the_one_a_run_gets(tmp_path: Path, data_dir: Path) -> None:
    """`NC_024936.1` carries an annotated pair at 99.80 %, and the page must say what a run says."""
    app = _workspace(tmp_path, data_dir)
    shutil.copy(data_dir / "NC_024936.1.gb", app.workspace / "records" / "NC_024936.1.gb")
    record = app.record("records/NC_024936.1.gb")
    gate = record["gates"]["annotation"]
    assert gate["status"] == "refused" and gate["identity_pct"] == 99.8
    assert gate["annotated_lengths"] == [26220, 26186] and gate["terminal_gaps"] == 52
    assert record["gates"]["sequence"]["status"] == "passed"
    assert record["gates"]["sequence"]["length"] == 26184
