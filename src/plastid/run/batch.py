"""Batch aggregation: one run per sample, one table that copies their cells.

The batch is an **orchestrator**, not a second engine. Every sample is analysed by the same module
the single-sample command runs, in its own run directory with its own manifest; this module reads
the stored results and copies each cell, tagging it with the sample and the scope the module
declared. Nothing is recomputed, averaged or filled in: a failed or not-measured cell is a status,
never a zero (M11.6).
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from plastid.result.result import Column, Result, Scope, Table
from plastid.run.figures import read_run_results

__all__ = [
    "BATCH_SCHEMA_VERSION",
    "SampleRecord",
    "aggregate_rows",
    "read_batch_results",
    "read_sample",
    "write_batch",
]

BATCH_SCHEMA_VERSION = 1

SAMPLES_COLUMNS = (
    "sample_id",
    "species",
    "accession",
    "input",
    "sha256",
    "status",
    "exit_code",
    "run_dir",
    "ir_source",
    "scoring_policy",
    "scoring_policy_version",
    "ir_gate",
    "modules",
)

AGGREGATE_COLUMNS = (
    "sample_id",
    "species",
    "module",
    "table",
    "row_index",
    "row_key",
    "column",
    "value",
    "measurement_status",
    "scope",
    "coord_system",
)


@dataclass
class SampleRecord:
    """One sample's outcome: its run, its status and the provenance the manifest recorded."""

    sample_id: str
    species: str
    accession: str
    input: str
    sha256: str
    status: str  # "ok" | "failed" | "skipped"
    exit_code: int
    run_dir: Path | None = None
    ir_source: str = "not required"
    scoring_policy: str = ""
    scoring_policy_version: int | None = None
    ir_gate: str = ""
    modules: list[dict[str, Any]] = field(default_factory=list)
    error: str = ""


def read_sample(
    *,
    sample_id: str,
    species: str,
    accession: str,
    input_path: str,
    sha256: str,
    exit_code: int,
    run_dir: Path | None,
) -> SampleRecord:
    """Describe one finished sample run from its own manifest — never from the batch's memory."""
    record = SampleRecord(
        sample_id=sample_id,
        species=species,
        accession=accession,
        input=input_path,
        sha256=sha256,
        status="ok" if exit_code == 0 else "failed",
        exit_code=exit_code,
        run_dir=run_dir,
    )
    manifest_path = (run_dir / "manifest.json") if run_dir is not None else None
    if manifest_path is None or not manifest_path.is_file():
        record.status = "failed"
        record.error = "the sample's run wrote no manifest"
        return record
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    record.status = "ok" if manifest.get("status") == "ok" and exit_code == 0 else "failed"
    record.modules = [
        {"id": entry.get("id"), "status": entry.get("status"), "error": entry.get("error", "")}
        for entry in manifest.get("modules", [])
    ]
    ir = manifest.get("ir")
    if isinstance(ir, dict):
        record.ir_source = str(ir.get("source", record.ir_source))
        record.scoring_policy = str(ir.get("scoring_policy", ""))
        version = ir.get("scoring_policy_version")
        record.scoring_policy_version = version if isinstance(version, int) else None
        record.ir_gate = "passed" if ir.get("gate_passed") else "failed"
    else:
        record.ir_gate = "not required"
    failed = [entry for entry in record.modules if entry["status"] == "failed"]
    if failed and not record.error:
        record.error = str(failed[0].get("error", ""))
    return record


def aggregate_rows(
    records: Sequence[SampleRecord], results: Mapping[str, Mapping[str, Result]]
) -> list[dict[str, Any]]:
    """Every stored cell, tagged with the sample, the module, the row key and the scope.

    The value is the table's own rendered cell (the same text the CSV/JSON export carries); a
    ``None`` cell or one the module marked "not measured" is recorded with that status. A cell is
    never defaulted to 0.
    """
    rows: list[dict[str, Any]] = []
    for record in records:
        sample_results = results.get(record.sample_id, {})
        for module_id, result in sample_results.items():
            for table_id in sorted(result.tables):
                table = result.tables[table_id]
                for index, row in enumerate(table.rows):
                    key = row.get(table.key) if table.key is not None else None
                    for column in table.columns:
                        value = row.get(column.name)
                        text = column.format(value)
                        not_measured = value is None or (
                            isinstance(value, str) and value.strip().lower() == "not measured"
                        )
                        status = "not_measured" if not_measured else "measured"
                        rows.append(
                            {
                                "sample_id": record.sample_id,
                                "species": record.species,
                                "module": module_id,
                                "table": table_id,
                                "row_index": index,
                                "row_key": "" if key is None else str(key),
                                "column": column.name,
                                "value": text,
                                "measurement_status": status,
                                "scope": table.scope.label,
                                "coord_system": table.scope.coord_system,
                            }
                        )
    return rows


def _write_csv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_batch(
    *,
    run_dir: Path,
    dataset: str,
    records: Sequence[SampleRecord],
    results: Mapping[str, Mapping[str, Result]],
    params: Mapping[str, Any],
    command: Sequence[str],
    inputs: Sequence[str],
) -> Path:
    """Write ``batch.json``, ``samples.csv`` and ``aggregate.csv`` into one batch run directory."""
    run_dir.mkdir(parents=True, exist_ok=True)
    rows = aggregate_rows(records, results)
    samples = [
        {
            **{
                column: getattr(record, column)
                for column in SAMPLES_COLUMNS
                if hasattr(record, column)
            },
            "run_dir": "" if record.run_dir is None else str(record.run_dir),
            "modules": ";".join(f"{entry['id']}:{entry['status']}" for entry in record.modules),
        }
        for record in records
    ]
    _write_csv(run_dir / "samples.csv", SAMPLES_COLUMNS, samples)
    _write_csv(run_dir / "aggregate.csv", AGGREGATE_COLUMNS, rows)
    ok = sum(1 for record in records if record.status == "ok")
    payload = {
        "schema_version": BATCH_SCHEMA_VERSION,
        "dataset": dataset,
        "status": "ok" if ok == len(records) else ("partial" if ok else "failed"),
        "params": dict(params),
        "command": list(command),
        "inputs": list(inputs),
        "samples": [
            {
                "sample_id": record.sample_id,
                "species": record.species,
                "accession": record.accession,
                "input": record.input,
                "sha256": record.sha256,
                "status": record.status,
                "exit_code": record.exit_code,
                "run_dir": "" if record.run_dir is None else str(record.run_dir),
                "ir_source": record.ir_source,
                "scoring_policy": record.scoring_policy,
                "scoring_policy_version": record.scoring_policy_version,
                "ir_gate": record.ir_gate,
                "modules": record.modules,
                "error": record.error,
            }
            for record in records
        ],
        "cells": len(rows),
        "not_measured_cells": sum(1 for row in rows if row["measurement_status"] != "measured"),
    }
    (run_dir / "batch.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return run_dir


def read_results(record: SampleRecord) -> Mapping[str, Result]:
    """The stored results of one sample run; a failed run has none, never fabricated ones."""
    if record.run_dir is None or not record.run_dir.is_dir():
        return {}
    try:
        return read_run_results(record.run_dir)
    except (OSError, ValueError, KeyError):
        return {}


def read_batch_results(run_dir: str | Path) -> dict[str, Result]:
    """The batch's long table rebuilt as per-module ``Result`` objects with a ``sample`` column.

    A figure pivot needs the tables a single run would have produced, plus the sample dimension.
    Every cell here is the text ``aggregate.csv`` stored — no value is reinterpreted. The tables
    have no key column: one row per (sample, original row), so a key would repeat by design.
    """
    root = Path(run_dir)
    path = root / "aggregate.csv"
    if not path.is_file():
        raise ValueError(f"{root}: no aggregate.csv — is this a batch run directory?")
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["module"], row["table"])
            entry = grouped.setdefault(
                key,
                {"columns": [], "rows": {}, "scope": row["scope"], "coord": row["coord_system"]},
            )
            name = row["column"]
            if name not in entry["columns"]:
                entry["columns"].append(name)
            row_key = (row["sample_id"], row["row_index"])
            cells = entry["rows"].setdefault(
                row_key,
                {
                    "sample": row["sample_id"],
                    "species": row["species"],
                    "row_key": row["row_key"],
                },
            )
            cells[name] = row["value"]

    modules: dict[str, dict[str, Table]] = {}
    for (module_id, table_id), entry in grouped.items():
        columns = (
            Column("sample", "Sample"),
            Column("species", "Species"),
            Column("row_key", "Row key"),
            *(Column(name, name) for name in entry["columns"]),
        )
        rows = tuple(entry["rows"][key] for key in sorted(entry["rows"]))
        modules.setdefault(module_id, {})[table_id] = Table(
            id=table_id,
            title=table_id,
            columns=columns,
            rows=rows,
            scope=Scope(label=entry["scope"], coord_system=entry["coord"]),
        )
    return {
        module_id: Result(module=module_id, module_version="batch", tables=tables)
        for module_id, tables in modules.items()
    }
