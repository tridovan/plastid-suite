"""``plastid batch`` — one run per sample, one aggregate that copies their cells (M11.6)."""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.modules import registry
from plastid.result.result import Column, Result, Scope, Table


def _stub(monkeypatch: pytest.MonkeyPatch, *, arity: int = 1, broken: str = "") -> None:
    spec = replace(
        registry.get("genome_stats"),
        id="batch_stub",
        arity=arity,
        requires_ir=False,
        reference_param=None,
        run=lambda genome, params, ctx: _result(broken, genome),
    )
    monkeypatch.setattr(registry, "get", lambda module_id: spec)


def _result(broken: str, genome: object) -> Result:
    name = getattr(genome, "name", "")
    if broken and str(name).startswith(broken):
        raise ValueError("this sample is broken on purpose")
    return Result(
        module="batch_stub",
        module_version="0.1.0",
        tables={
            "batch_summary": Table(
                id="batch_summary",
                title="Batch stub",
                columns=(
                    Column("n", "N", dtype="int"),
                    Column("label", "Label"),
                ),
                rows=({"n": 1, "label": str(name)},),
                scope=Scope(label="synthetic"),
            )
        },
    )


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


def _batch_dirs(out: Path) -> list[Path]:
    return sorted((out / "batch").glob("*/batch.json"))


def test_two_samples_get_two_runs_and_one_aggregate(tmp_path: Path, monkeypatch) -> None:
    _stub(monkeypatch)
    first = _write(tmp_path, "a.fa", ">ALPHA.1\nACGTACGTACGTACGT\n")
    second = _write(tmp_path, "b.fa", ">BETA.1\nACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    code = main(["batch", "-m", "batch_stub", "-i", str(first), "-i", str(second), "-o", str(out)])

    assert code == 0
    assert len(list((out / "alpha.1").glob("*/manifest.json"))) == 1
    assert len(list((out / "beta.1").glob("*/manifest.json"))) == 1
    (payload,) = [json.loads(path.read_text(encoding="utf-8")) for path in _batch_dirs(out)]
    assert payload["status"] == "ok"
    assert [sample["sample_id"] for sample in payload["samples"]] == ["ALPHA.1", "BETA.1"]
    assert all(sample["status"] == "ok" for sample in payload["samples"])
    assert all(sample["sha256"] for sample in payload["samples"])
    rows = list(
        csv.DictReader((_batch_dirs(out)[0].parent / "aggregate.csv").open(encoding="utf-8"))
    )
    cells = {(row["sample_id"], row["column"], row["value"]) for row in rows}
    assert ("ALPHA.1", "label", "synthetic alpha") not in cells  # name is the FASTA header
    assert ("ALPHA.1", "label", "ALPHA.1") in cells
    assert ("BETA.1", "n", "1") in cells
    assert all(row["measurement_status"] == "measured" for row in rows)


def test_a_broken_sample_is_named_and_the_other_sample_still_finishes(
    tmp_path: Path, monkeypatch
) -> None:
    _stub(monkeypatch, broken="BROKEN")
    good = _write(tmp_path, "good.fa", ">GOOD.1\nACGTACGTACGTACGT\n")
    broken = _write(tmp_path, "broken.fa", ">BROKEN.1\nACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    code = main(["batch", "-m", "batch_stub", "-i", str(good), "-i", str(broken), "-o", str(out)])

    assert code != 0
    payload = json.loads(_batch_dirs(out)[0].read_text(encoding="utf-8"))
    statuses = {sample["sample_id"]: sample["status"] for sample in payload["samples"]}
    assert statuses == {"GOOD.1": "ok", "BROKEN.1": "failed"}
    assert payload["status"] == "partial"
    # the broken sample contributes no fabricated metric cells
    rows = list(
        csv.DictReader((_batch_dirs(out)[0].parent / "aggregate.csv").open(encoding="utf-8"))
    )
    assert {row["sample_id"] for row in rows} == {"GOOD.1"}


def test_duplicate_sample_names_keep_both_accessions(tmp_path: Path, monkeypatch) -> None:
    _stub(monkeypatch)
    first = _write(tmp_path, "one.fa", ">DUP.1\nACGTACGTACGTACGT\n")
    second = _write(tmp_path, "two.fa", ">DUP.1\nACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    assert (
        main(["batch", "-m", "batch_stub", "-i", str(first), "-i", str(second), "-o", str(out)])
        == 0
    )

    payload = json.loads(_batch_dirs(out)[0].read_text(encoding="utf-8"))
    assert [sample["sample_id"] for sample in payload["samples"]] == ["DUP.1", "DUP.1#2"]
    assert {sample["accession"] for sample in payload["samples"]} == {"DUP.1"}


def test_an_arity_two_module_is_refused_before_anything_runs(tmp_path: Path, monkeypatch) -> None:
    _stub(monkeypatch, arity=2)
    sample = _write(tmp_path, "a.fa", ">ALPHA.1\nACGTACGTACGTACGT\n")
    out = tmp_path / "runs"

    assert main(["batch", "-m", "batch_stub", "-i", str(sample), "-o", str(out)]) == 2
    assert not list(out.glob("**/manifest.json"))


def test_an_aggregate_cell_is_the_single_sample_runs_own_cell(tmp_path: Path) -> None:
    """M11.6: the batch copies the module's cell; it does not recompute it."""
    from plastid.run.figures import read_run_results

    golden = Path(__file__).parents[1] / "data" / "NC_050833.1.gb"
    single = tmp_path / "single"
    batched = tmp_path / "batched"
    assert main(["run", "-m", "genome_stats", "-i", str(golden), "-o", str(single)]) == 0
    assert main(["batch", "-m", "genome_stats", "-i", str(golden), "-o", str(batched)]) == 0

    run_dir = next(single.glob("*/*/manifest.json")).parent
    table = read_run_results(run_dir)["genome_stats"].tables["genome_stats"]
    row = table.rows[0]
    rows = list(
        csv.DictReader((next((batched / "batch").glob("*/aggregate.csv"))).open(encoding="utf-8"))
    )
    cells = {
        (entry["sample_id"], entry["module"], entry["table"], entry["column"]): entry["value"]
        for entry in rows
    }
    for column in table.columns:
        key = ("NC_050833.1", "genome_stats", "genome_stats", column.name)
        assert key in cells
        assert cells[key] == column.format(row[column.name])
