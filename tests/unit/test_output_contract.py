"""Output parity and refresh contracts across the export formats (M11.18).

The same stored result is read through every published format — CSV, JSON, the Excel workbook,
the Word document and Markdown — and every rendered cell must be the *same* string. Scope is
part of the number's meaning, so each format must also state it. The Supplementary package is
then checked as a folder: links resolve, digests match, a rebuild removes what the run no longer
holds, and an absent renderer or workbook is named rather than silently missing.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pytest

from plastid.core.deps import ExportDependencyError
from plastid.report.package import build_package
from plastid.report.qc import render_qc
from plastid.result.docx_export import DocxMeta, export_docx
from plastid.result.export import export_csv, export_json, export_markdown
from plastid.result.figure import FigureRecord
from plastid.result.result import Column, Result, Scope, Table
from plastid.result.xlsx_export import export_xlsx
from plastid.run.batch import SampleRecord, read_batch_results, write_batch


def _result() -> Result:
    summary = Table(
        id="summary",
        title="Genome summary",
        columns=(
            Column("organism", "Organism"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
            Column("ir_source", "IR source"),
        ),
        rows=(
            {
                "organism": "Iris domestica",
                "length_bp": 153_724,
                "gc_pct": 37.85,
                "ir_source": None,  # a value this run did not measure
            },
        ),
        scope=Scope(label="whole plastome · third-party sample", coord_system="genbank"),
    )
    regions = Table(
        id="regions",
        title="Regions",
        columns=(
            Column("region", "Region"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
        ),
        rows=(
            {"region": "LSC", "start_1b": 1, "length_bp": 83_127},
            {"region": "SSC", "start_1b": 109_342, "length_bp": 18_169},
        ),
        scope=Scope(label="normalized frame", coord_system="normalized"),
        key="region",
    )
    return Result(
        module="genome_stats",
        module_version="1.0.0",
        tables={"summary": summary, "regions": regions},
        notes=("a note that belongs to the result",),
        citations=("biopython",),
    )


def _rendered(table: Table) -> tuple[list[str], list[list[str]]]:
    header = [column.header for column in table.columns]
    rows = [[column.format(row.get(column.name)) for column in table.columns] for row in table.rows]
    return header, rows


# ------------------------------------------------------------------- cross-format parity


def test_every_export_renders_the_same_cells_and_the_same_scope() -> None:
    result = _result()

    # CSV: the two comment lines carry the scope, then header and rows.
    for table in result.tables.values():
        lines = export_csv(table).splitlines()
        assert lines[0] == f"# table: {table.title}"
        assert lines[1] == (
            f"# scope: {table.scope.label} | coordinates: {table.scope.coord_system}"
        )
        parsed = list(csv.reader(io.StringIO("\n".join(lines[2:]))))
        header, rows = _rendered(table)
        assert parsed == [header, *rows]

    # JSON: raw values, the same `Column.format`, and the scope as a field.
    payload = json.loads(export_json(result))
    for table_id, table in result.tables.items():
        blob = payload["tables"][table_id]
        header, rows = _rendered(table)
        assert [
            column["label"] + (f" ({column['unit']})" if column["unit"] else "")
            for column in blob["columns"]
        ] == list(header)
        assert [
            [
                table.columns[index].format(row[table.columns[index].name])
                for index in range(len(header))
            ]
            for row in blob["rows"]
        ] == rows
        assert blob["scope"] == {
            "label": table.scope.label,
            "coord_system": table.scope.coord_system,
            "include": [],
            "exclude_regions": [],
        }
    assert payload["citations"] == ["biopython"]

    # Markdown: the same header and cells, with the scope line above each table.
    markdown = export_markdown(result)
    for table in result.tables.values():
        header, rows = _rendered(table)
        assert f"*scope: {table.scope.label} | coordinates: {table.scope.coord_system}*" in markdown
        assert "| " + " | ".join(header) + " |" in markdown
        for row in rows:
            assert "| " + " | ".join(row) + " |" in markdown


def test_the_workbook_holds_the_same_cells_and_names_each_table_s_scope() -> None:
    openpyxl = pytest.importorskip("openpyxl", reason="the 'excel' extra is not installed")
    result = _result()
    book = openpyxl.load_workbook(io.BytesIO(export_xlsx({"genome_stats": result})))

    for table in result.tables.values():
        sheet = book[f"genome_stats.{table.id}"]
        header, rows = _rendered(table)
        assert [cell.value for cell in sheet[1]] == header
        # openpyxl reads an empty cell back as None; the renderer wrote "" (not measured).
        assert [
            [cell.value if cell.value is not None else "" for cell in line]
            for line in sheet.iter_rows(min_row=2)
        ] == rows

    scope = book["plastid.scope"]
    listed = {
        (row[0].value, row[1].value, row[2].value, row[3].value)
        for row in scope.iter_rows(min_row=2)
    }
    assert listed == {
        ("genome_stats", "summary", "whole plastome · third-party sample", "genbank"),
        ("genome_stats", "regions", "normalized frame", "normalized"),
    }


def test_the_word_document_holds_the_same_cells_scope_and_references() -> None:
    docx = pytest.importorskip("docx", reason="the 'docx' extra is not installed")
    result = _result()
    document = docx.Document(
        io.BytesIO(
            export_docx(
                result,
                DocxMeta(
                    title="Genome report",
                    label="NC_050833.1 · third-party sample",
                    references=("Beier S et al. MISA-web. Bioinformatics (2017).",),
                ),
            )
        )
    )

    assert len(document.tables) == len(result.tables)
    for table, word_table in zip(result.tables.values(), document.tables, strict=True):
        header, rows = _rendered(table)
        assert [cell.text for cell in word_table.rows[0].cells] == header
        assert [[cell.text for cell in line.cells] for line in word_table.rows[1:]] == rows

    paragraphs = [paragraph.text for paragraph in document.paragraphs]
    for table in result.tables.values():
        assert any(
            f"scope: {table.scope.label} | coordinates: {table.scope.coord_system}" in text
            for text in paragraphs
        )
    assert any("Beier S et al." in text for text in paragraphs)


# ------------------------------------------------------------------------- the package


def _figure_file(root: Path, name: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / name
    path.write_bytes(b"figure-placeholder")
    return path


def test_the_package_carries_the_figures_it_is_given_and_names_the_ones_it_lacks(
    tmp_path: Path,
) -> None:
    pytest.importorskip("openpyxl", reason="the 'excel' extra is not installed")
    run = tmp_path / "run"
    drawn = _figure_file(tmp_path / "figures", "genome_page.pdf")
    package = build_package(
        run,
        results={"genome_stats": _result()},
        figure_records=(
            FigureRecord("genome_page", "matplotlib", (drawn,), provenance=("summary — scope: x",)),
            FigureRecord("r_rscu", "R", reason="R is not installed"),
        ),
        missing_renderers=(("R", "Rscript not found"),),
        label="panel 80 · third-party",
    )
    index = (package.directory / "index.md").read_text(encoding="utf-8")
    assert "panel 80 · third-party" in index
    assert "[genome_page.pdf](figures/genome_page.pdf)" in index
    assert "summary — scope: x" in index
    assert "`r_rscu` (R) — R is not installed" in index
    assert "Rscript not found" in index
    assert "Biopython" in index  # the result's citation reaches the supplement

    manifest = json.loads((package.directory / "MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["directory"] == "supplementary"
    assert manifest["data_label"] == "panel 80 · third-party"
    assert {entry["figure"] for entry in manifest["figures"]} == {"genome_page", "r_rscu"}
    assert ("R", "Rscript not found") in {
        (entry["output"], entry["reason"]) for entry in manifest["skipped_outputs"]
    }
    for entry in manifest["files"]:
        path = package.directory / entry["path"]
        assert entry["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_a_package_rebuild_removes_what_the_run_no_longer_holds(tmp_path: Path) -> None:
    pytest.importorskip("openpyxl", reason="the 'excel' extra is not installed")
    run = tmp_path / "run"
    first = build_package(run, results={"genome_stats": _result()})
    stale = first.directory / "tables" / "genome_stats__regions.csv"
    assert stale.is_file()

    only_summary = {
        "genome_stats": Result(
            module="genome_stats",
            module_version="1.0.0",
            tables={"summary": _result().tables["summary"]},
        )
    }
    second = build_package(run, results=only_summary)
    assert not stale.exists(), "a table the run no longer holds must not survive a rebuild"
    assert (second.directory / "tables" / "genome_stats__summary.csv").is_file()
    listed = {
        entry["path"]
        for entry in json.loads((second.directory / "MANIFEST.json").read_text(encoding="utf-8"))[
            "files"
        ]
    }
    assert "tables/genome_stats__regions.csv" not in listed


def test_a_failed_rebuild_leaves_the_previous_package_in_place(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A half-written package must never replace a complete one.

    The rebuild stages beside the package and swaps only on success; an export that fails leaves
    the old folder exactly as it was, and the command that called it has failed loudly (M11.18).
    """
    pytest.importorskip("openpyxl", reason="the 'excel' extra is not installed")
    run = tmp_path / "run"
    first = build_package(run, results={"genome_stats": _result()})
    before = (first.directory / "index.md").read_bytes()

    def refuse(table):  # type: ignore[no-untyped-def]
        raise OSError("disk full")

    monkeypatch.setattr("plastid.report.package.export_csv", refuse)
    with pytest.raises(OSError, match="disk full"):
        build_package(run, results={"genome_stats": _result()})

    assert (run / "report" / "supplementary" / "index.md").read_bytes() == before


def test_a_missing_excel_extra_leaves_a_named_gap_rather_than_no_package(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def refuse(results):  # type: ignore[no-untyped-def]
        raise ExportDependencyError(
            "Excel export needs openpyxl — pip install 'plastid-suite[excel]'"
        )

    monkeypatch.setattr("plastid.report.package.export_xlsx", refuse)
    run = tmp_path / "run"
    package = build_package(run, results={"genome_stats": _result()})
    assert (package.directory / "tables" / "genome_stats__summary.csv").is_file()
    assert not (package.directory / "tables.xlsx").exists()
    index = (package.directory / "index.md").read_text(encoding="utf-8")
    assert "`tables.xlsx` — Excel export needs openpyxl" in index
    manifest = json.loads((package.directory / "MANIFEST.json").read_text(encoding="utf-8"))
    assert any(entry["output"] == "tables.xlsx" for entry in manifest["skipped_outputs"])


# ------------------------------------------------------------------------------- QC


def test_the_qc_report_lists_both_catalogues_with_their_status() -> None:
    manifest = {
        "run_id": "r1",
        "command": ["plastid", "run"],
        "package_version": "1.0.0",
        "started_at": "2026-10-27T00:00:00",
        "finished_at": "2026-10-27T00:00:01",
        "status": "ok",
        "params": {},
        "inputs": [],
        "coordinate_system": {"name": "genbank", "rotation_offset": 0},
        "modules": [],
        "files": [],
    }
    text = render_qc(
        manifest,
        figure_records=(
            FigureRecord("pi_window", "matplotlib", reason="needs diversity.diversity_windows"),
            FigureRecord("r_rscu", "R", reason="Rscript not found"),
        ),
        missing_renderers=(("R", "Rscript not found"),),
        data_label="panel 80 · third-party",
    )
    assert "| pi_window | matplotlib | skipped |" in text
    assert "| r_rscu | R | skipped |" in text
    assert "Rscript not found" in text
    assert "data label: panel 80 · third-party" in text


def test_package_and_qc_read_a_run_without_rerunning_analysis(tmp_path: Path) -> None:
    """The two output commands consume the stored run; neither creates a second analysis.

    The stored ``result.json`` is byte-identical after both commands and no run directory is
    added. The package's provenance and the QC's catalogue both name the data label and the IR
    path the run recorded (M11.18).
    """
    from plastid.cli import main

    data = Path(__file__).parents[1] / "data" / "NC_050833.1.gb"
    runs = tmp_path / "runs"
    assert (
        main(
            [
                "run",
                "-m",
                "genome_stats",
                "-i",
                str(data),
                "-o",
                str(runs),
                "--dataset",
                "panel",
                "--label",
                "panel-test",
            ]
        )
        == 0
    )
    run_dir = next((runs / "panel").glob("*/manifest.json")).parent
    result_path = run_dir / "01_genome_stats" / "result.json"
    before = result_path.read_bytes()
    directories = sorted(path.name for path in (runs / "panel").iterdir())

    assert main(["package", str(run_dir)]) == 0
    assert main(["qc", str(run_dir)]) == 0

    assert result_path.read_bytes() == before
    assert sorted(path.name for path in (runs / "panel").iterdir()) == directories
    index = (run_dir / "report" / "supplementary" / "index.md").read_text(encoding="utf-8")
    assert "data label: panel-test" in index
    assert "IR: source annotation" in index
    qc = (run_dir / "report" / "qc.md").read_text(encoding="utf-8")
    assert "data label: panel-test" in qc
    assert "| r_rscu | R |" in qc


def test_a_figure_refresh_drops_a_manifest_entry_whose_file_is_gone(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A manifest that still lists a deleted figure would describe a result that is not there."""
    from types import SimpleNamespace

    from plastid import cli
    from plastid.run.manifest import build_manifest, write_manifest

    run = tmp_path / "runs" / "dataset" / "20261027-000000-abcdef12"
    module_dir = run / "01_genome_stats"
    module_dir.mkdir(parents=True)
    (module_dir / "result.json").write_text(export_json(_result()), encoding="utf-8")
    manifest = build_manifest(
        run_id=run.name,
        command=["plastid", "run"],
        package_version="1.0.0",
        started_at="2026-10-27T00:00:00",
        params={},
        inputs=[],
        coordinate_system={"name": "genbank", "rotation_offset": 0},
    )
    manifest["files"] = [
        {
            "path": "figures/circular_map.svg",
            "module": "figures",
            "sha256": "a" * 64,
            "bytes": 1,
        }
    ]
    write_manifest(run / "manifest.json", manifest)

    monkeypatch.setattr(
        cli,
        "_render_output_figures",
        lambda *args, **kwargs: cli.OutputFigures((), (), (), ()),
    )
    assert cli._cmd_figures(SimpleNamespace(run=str(run), figure=None, sets=[])) == cli.EXIT_ERROR

    written = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    assert written["files"] == []
    assert written["figures"]["status"] == "complete"


def test_a_figure_refresh_updates_the_manifest_digest_of_a_rewritten_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A replaced image must not keep the old digest in the run's file list."""
    from types import SimpleNamespace

    from plastid import cli
    from plastid.run.manifest import build_manifest, write_manifest

    run = tmp_path / "runs" / "dataset" / "20261027-000000-abcdef12"
    module_dir = run / "01_genome_stats"
    module_dir.mkdir(parents=True)
    (module_dir / "result.json").write_text(export_json(_result()), encoding="utf-8")
    figure = run / "figures" / "circular_map.svg"
    figure.parent.mkdir(parents=True)
    figure.write_bytes(b"new drawing")
    manifest = build_manifest(
        run_id=run.name,
        command=["plastid", "run"],
        package_version="1.0.0",
        started_at="2026-10-27T00:00:00",
        params={},
        inputs=[],
        coordinate_system={"name": "genbank", "rotation_offset": 0},
    )
    manifest["files"] = [
        {
            "path": "figures/circular_map.svg",
            "module": "figures",
            "sha256": "a" * 64,
            "bytes": 1,
        }
    ]
    write_manifest(run / "manifest.json", manifest)

    monkeypatch.setattr(
        cli,
        "_render_output_figures",
        lambda *args, **kwargs: cli.OutputFigures((figure,), (), (), ()),
    )
    assert cli._cmd_figures(SimpleNamespace(run=str(run), figure=None, sets=[])) == cli.EXIT_OK

    written = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    (entry,) = written["files"]
    assert entry["sha256"] == hashlib.sha256(b"new drawing").hexdigest()
    assert entry["bytes"] == len(b"new drawing")


# ------------------------------------------------------------------------------ batch


def test_the_batch_aggregate_carries_the_sample_provenance_on_every_row(tmp_path: Path) -> None:
    record = SampleRecord(
        sample_id="NC_050833.1",
        species="Iris domestica",
        accession="NC_050833.1",
        input="panel/NC_050833.1.gb",
        sha256="a" * 64,
        status="ok",
        exit_code=0,
        ir_source="annotation",
        ir_gate="passed",
    )
    write_batch(
        run_dir=tmp_path / "batch",
        dataset="panel",
        records=[record],
        results={"NC_050833.1": {"genome_stats": _result()}},
        params={"data_label": "panel 80 · third-party"},
        command=["plastid", "batch"],
        inputs=["panel/NC_050833.1.gb"],
    )
    path = tmp_path / "batch" / "aggregate.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    assert rows
    assert {row["sample_accession"] for row in rows} == {"NC_050833.1"}
    assert {row["sample_data_label"] for row in rows} == {"panel 80 · third-party"}
    assert {row["sample_ir_source"] for row in rows} == {"annotation"}
    assert {row["sample_ir_gate"] for row in rows} == {"passed"}

    # The provenance columns are not measured cells, so a figure pivot must not add them.
    rebuilt = read_batch_results(tmp_path / "batch")
    names = {
        column.name for table in rebuilt["genome_stats"].tables.values() for column in table.columns
    }
    assert not names & {
        "sample_accession",
        "sample_data_label",
        "sample_ir_source",
        "sample_ir_gate",
    }
