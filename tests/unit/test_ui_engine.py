"""One engine, two ways in: a run started from the interface is a run the command line made."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import pytest
from test_ui_commands import FULL_FORMS, _leaf

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


def test_a_batch_is_found_read_and_cancelled_like_a_run(tmp_path: Path, data_dir: Path) -> None:
    """A batch writes `batch.json` and one run per sample, all under the same process number."""
    app = _workspace(tmp_path, data_dir)
    shutil.copy(data_dir / "NC_024936.1.gb", app.workspace / "records" / "NC_024936.1.gb")
    form = {
        "modules": ["genome_stats"],
        "inputs": ["records/NC_050833.1.gb", "records/NC_024936.1.gb"],
        "dataset": "many",
        "one_by_one": True,
        "sets": {"ir_source": "sequence"},
    }
    job = _finished(app, app.start_job({"command": "run", "form": form})["id"])
    assert job["status"] == "ok" and job["run"]["dataset"] == "many"

    detail = app.run("many", job["run"]["run_id"])
    assert detail["kind"] == "batch"
    assert [sample["accession"] for sample in detail["samples"]] == ["NC_050833.1", "NC_024936.1"]
    assert all(sample["run"]["dataset"] for sample in detail["samples"])
    assert detail["aggregate"].endswith("aggregate.csv")
    assert app.table(detail["aggregate"])["total"] > 0


def test_a_command_on_a_run_leads_back_to_that_run(tmp_path: Path, data_dir: Path) -> None:
    """`figures`, `package` and `qc` make no directory, and the job page still needs a way on."""
    app = _workspace(tmp_path, data_dir)
    _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    run_id = app.runs()["runs"][0]["run_id"]
    form = {"run": f"runs/same/{run_id}"}
    job = _finished(app, app.start_job({"command": "qc", "form": form})["id"])
    assert job["run"] == {"dataset": "same", "run_id": run_id}


def test_the_figures_of_a_run_come_from_the_engine(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    run_id = app.runs()["runs"][0]["run_id"]
    listing = app.figures(f"runs/same/{run_id}")
    can_draw = {figure["id"]: figure["can_draw"] for figure in listing["figures"]}
    assert can_draw["r_rscu"] and not can_draw["pi_window"]
    assert "outgroup" in {option["name"] for option in listing["options"]}


def test_a_folder_that_is_no_run_is_said_in_words(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    (app.runs_dir / "empty" / "x").mkdir(parents=True)
    with pytest.raises(Refused, match="not a run directory"):
        app.figures("runs/empty/x")


def test_a_command_on_a_short_path_still_has_a_page(tmp_path: Path, data_dir: Path) -> None:
    """The engine refuses it by name; the job page must stay up to say so."""
    app = _workspace(tmp_path, data_dir)
    job = _finished(app, app.start_job({"command": "qc", "form": {"run": "x"}})["id"])
    assert job["status"] == "refused" and job["run"] is None
    assert "not a run directory" in job["output"]


def test_a_set_of_references_can_be_listed(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    shutil.copy(data_dir / "NC_024936.1.gb", app.workspace / "references" / "NC_024936.1.gb")
    assert [row["name"] for row in app.records("references")["records"]] == ["NC_024936.1.gb"]
    with pytest.raises(Refused, match="outside the workspace"):
        app.records("../..")


def test_a_run_with_no_dataset_name_can_be_opened(tmp_path: Path, data_dir: Path) -> None:
    """The name is an accession as written, and its folder is that name in lower case."""
    app = _workspace(tmp_path, data_dir)
    form = {"modules": ["genome_stats"], "inputs": ["records/NC_050833.1.gb"]}
    _finished(app, app.start_job({"command": "run", "form": form})["id"])
    row = app.runs()["runs"][0]
    assert row["dataset"] == "NC_050833.1"
    assert app.run(row["folder"], row["run_id"])["status"] == "ok"


@pytest.mark.parametrize("command", sorted(FULL_FORMS))
def test_every_run_lands_where_results_reads(tmp_path: Path, command: str) -> None:
    """Left to the parser's default, a run goes to a hidden folder the page never lists."""
    app = App(tmp_path / "work", version="test")
    args = app.command({"command": command, "form": FULL_FORMS[command]})["args"]
    writes_a_run = {"-o", "--out"} & {
        option for action in _leaf(command)._actions for option in action.option_strings
    }
    if writes_a_run:
        assert args[args.index("-o") + 1] == "runs"


def test_a_record_a_run_wrote_can_join_the_records(tmp_path: Path, data_dir: Path) -> None:
    """Annotate comes before Records on the route, so its output has to be able to get there."""
    app = _workspace(tmp_path, data_dir)
    written = app.runs_dir / "a" / "1" / "01_plastid_annotate" / "plastid_annotation.gb"
    written.parent.mkdir(parents=True)
    shutil.copy(data_dir / "NC_024936.1.gb", written)
    path = app.relative(written)

    assert app.adopt({"path": path, "name": "sample.gb"}) == {"path": "records/sample.gb"}
    assert (app.workspace / "records" / "sample.gb").read_bytes() == written.read_bytes()
    with pytest.raises(Refused, match="not overwritten"):
        app.adopt({"path": path, "name": "sample.gb"})
    with pytest.raises(Refused, match="file name"):
        app.adopt({"path": path, "name": "../sample.gb"})
    with pytest.raises(Refused, match="not a GenBank or FASTA"):
        app.adopt({"path": "runs/a/1/01_plastid_annotate/plastid_annotation.gb", "name": "a.txt"})


def test_a_run_says_what_it_gives_the_next_step(tmp_path: Path, data_dir: Path) -> None:
    app = _workspace(tmp_path, data_dir)
    _finished(app, app.start_job({"command": "run", "form": _form()})["id"])
    run_id = app.runs()["runs"][0]["run_id"]
    assert app.run("same", run_id)["artifacts"] == []
    folder = app.runs_dir / "same" / run_id
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    manifest["artifacts"] = [
        {"path": "01_x/x.gb", "format": "genbank", "sample_id": "S1", "status": "ok"},
        {"path": "01_x/y.fasta", "format": "fasta", "sample_id": "S1", "status": "undetermined"},
    ]
    (folder / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    assert app.run("same", run_id)["artifacts"] == [
        {
            "path": f"runs/same/{run_id}/01_x/x.gb",
            "format": "genbank",
            "sample_id": "S1",
            "status": "ok",
        },
        {
            "path": f"runs/same/{run_id}/01_x/y.fasta",
            "format": "fasta",
            "sample_id": "S1",
            "status": "undetermined",
        },
    ]


def test_the_page_names_the_folder_the_person_knows(tmp_path: Path) -> None:
    """In the image the workspace is /work, which means nothing to someone on Windows."""
    plain = App(tmp_path / "plain", version="test")
    assert plain.about()["workspace"] == str(plain.workspace)
    mounted = App(tmp_path / "work", version="test", shown_as="C:\\Users\\a\\Documents\\PTS")
    assert mounted.about()["workspace"] == "C:\\Users\\a\\Documents\\PTS"


def _kept_in_the_workspace(app: App) -> Path:
    """A copy of CPJSdraw as a person would leave it: the whole folder, under `tools`."""
    folder = app.workspace / "tools" / "cpjsdraw" / "bin"
    folder.mkdir(parents=True)
    script = folder / "CPJSdraw.pl"
    script.write_text('my $version = "0.0.1";\n', encoding="utf-8")
    return script


def test_a_tool_kept_in_the_workspace_is_found(tmp_path: Path) -> None:
    """The image cannot carry CPJSdraw, and someone on Windows has no PATH to put it on."""
    app = App(tmp_path / "work", version="test")
    script = _kept_in_the_workspace(app)
    tools = {tool["name"]: tool for tool in app.tools()["tools"]}
    assert tools["cpjsdraw"]["status"] == "ok"
    assert tools["cpjsdraw"]["found_at"] == str(script)


def test_a_run_finds_the_tool_the_page_found(tmp_path: Path, data_dir: Path) -> None:
    """The stand-in is refused as unpatched, which is how this shows the run reached that copy."""
    app = _workspace(tmp_path, data_dir)
    script = _kept_in_the_workspace(app)
    form = {"inputs": ["records/NC_050833.1.gb"], "dataset": "ir-junction"}
    job = _finished(app, app.start_job({"command": "ir-junction", "form": form})["id"])
    assert job["status"] == "tool_error"
    assert f"{script}: this is not the patched CPJSdraw" in job["output"].replace("\n", " ")


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
