"""Read a run directory for the page: what ran, what it wrote, and every figure it did not draw.

The directory is the only source. A run started on the command line and one started here read
the same, because neither leaves anything the other cannot see.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from plastid.run.state import RunState, read_state

__all__ = ["record_of", "run_detail", "state_row", "table"]

_TABLE_ROWS = 500  # rows sent to the page; the file itself holds them all
_IMAGE_SUFFIXES = (".png", ".svg")
_FILE_SUFFIXES = (".png", ".pdf", ".svg", ".eps", ".csv", ".cfg", ".treefile")


def state_row(state: RunState) -> dict[str, Any]:
    """One line of the runs page, with the same words `plastid runs` prints."""
    return {
        "run_id": state.run_id,
        "dataset": state.dataset,
        "folder": state.run_dir.parent.name,  # the address: a run names itself as it was typed
        "status": state.status,
        "kind": state.kind,
        "done": state.done,
        "total": state.total,
        "started_at": state.started_at,
        "finished_at": state.finished_at,
        "reason": state.reason,
        "current": state.current,
        "exit_code": state.exit_code,
    }


def record_of(run_dir: Path) -> dict[str, Any]:
    """What a run says about itself. A batch directory holds `batch.json` and no manifest."""
    for name in ("batch.json", "manifest.json"):
        path = run_dir / name
        if path.is_file():
            return json.loads(path.read_text(encoding="utf-8"))
    return {}


def _step(run_dir: Path, entry: dict[str, Any], root: Path) -> dict[str, Any]:
    folders = sorted(run_dir.glob(f"[0-9]*_{entry.get('id', '')}"))
    folder = folders[0] if folders else None
    tables = sorted((folder / "tables").glob("*.csv")) if folder else []
    report = folder / "report.md" if folder else None
    return {
        "id": entry.get("id"),
        "version": entry.get("version"),
        "status": entry.get("status"),
        "error": entry.get("error"),
        "tables": [
            {"name": path.stem, "path": path.relative_to(root).as_posix()} for path in tables
        ],
        "report": report.relative_to(root).as_posix() if report and report.is_file() else None,
        "tools": [
            {
                "name": tool.get("name"),
                "version": tool.get("version"),
                "rc": tool.get("rc"),
                "argv": tool.get("argv"),
                "duration_s": tool.get("duration_s"),
            }
            for tool in entry.get("tools") or []
        ],
    }


def _figure_id(name: str) -> tuple[str, str]:
    """A figure file's catalogue id and its label: ``r_rscu__Codon.png`` or ``pi_window.svg``."""
    stem = name.split(".", 1)[0]
    return (stem.split("__", 1)[0], stem.split("__", 1)[1]) if "__" in stem else (stem, stem)


def _figures(run_dir: Path, record: dict[str, Any], root: Path) -> dict[str, Any]:
    folder = run_dir / "figures"
    drawn: dict[str, dict[str, Any]] = {}
    if folder.is_dir():
        for path in sorted(folder.iterdir()):
            if not path.is_file() or path.suffix not in _FILE_SUFFIXES:
                continue
            figure_id, label = _figure_id(path.name)
            entry = drawn.setdefault(figure_id, {"id": figure_id, "views": {}, "files": []})
            relative = path.relative_to(root).as_posix()
            entry["files"].append({"name": path.name, "path": relative})
            if path.suffix in _IMAGE_SUFFIXES:
                views = entry["views"]
                # The PNG of a pair is the lighter to show; the SVG of a pair is the same drawing.
                if label not in views or path.suffix == ".png":
                    views[label] = relative
        for figure_id, entry in drawn.items():
            caption = folder / f"{figure_id}.caption.txt"
            if caption.is_file():
                entry["caption"] = caption.relative_to(root).as_posix()
    section = record.get("figures") or {}
    return {
        "status": section.get("status"),
        "error": section.get("error"),
        "drawn": [entry for entry in drawn.values() if entry["views"]],
        "failed": section.get("failed") or [],
        "skipped": [
            entry
            for entry in section.get("skipped") or []
            if isinstance(entry, dict) and entry.get("figure") not in drawn
        ],
        "missing_renderers": section.get("missing_renderers") or [],
    }


def _sample(entry: dict[str, Any], root: Path) -> dict[str, Any]:
    """One row of a batch, with the address of its own run when that lies in the workspace."""
    folder = Path(str(entry.get("run_dir") or ""))
    inside = folder.is_absolute() and folder.is_relative_to(root) and folder.is_dir()
    return {
        "sample_id": entry.get("sample_id"),
        "species": entry.get("species"),
        "accession": entry.get("accession"),
        "status": entry.get("status"),
        "ir_gate": entry.get("ir_gate"),
        "error": entry.get("error"),
        "run": {"dataset": folder.parent.name, "run_id": folder.name} if inside else None,
    }


def _artifact(entry: dict[str, Any], run_dir: Path, root: Path) -> dict[str, Any]:
    """A file a run hands to the next step, with its path from the workspace."""
    return {
        "path": (run_dir / str(entry.get("path"))).relative_to(root).as_posix(),
        "format": entry.get("format"),
        "sample_id": entry.get("sample_id"),
        "status": entry.get("status"),
    }


def run_detail(run_dir: Path, root: Path) -> dict[str, Any]:
    """Everything the run page shows, with paths relative to the workspace."""
    state = read_state(run_dir)
    record = record_of(run_dir)
    aggregate = run_dir / "aggregate.csv"
    return {
        **state_row(state),
        "path": run_dir.relative_to(root).as_posix(),
        "command": record.get("command"),
        "package_version": record.get("package_version"),
        # A batch lists its inputs as the paths typed; their digests are in `samples`.
        "inputs": [
            {"path": entry} if isinstance(entry, str) else entry
            for entry in record.get("inputs") or []
        ],
        "params": record.get("params") or {},
        "data_label": record.get("data_label"),
        "ir": record.get("ir"),
        "reference_sets": record.get("reference_sets") or [],
        "steps": [_step(run_dir, entry, root) for entry in record.get("modules") or []],
        "samples": [_sample(entry, root) for entry in record.get("samples") or []],
        "aggregate": aggregate.relative_to(root).as_posix() if aggregate.is_file() else None,
        "artifacts": [_artifact(entry, run_dir, root) for entry in record.get("artifacts") or []],
        "figures": _figures(run_dir, record, root),
    }


def table(path: Path) -> dict[str, Any]:
    """A table as the run wrote it: its title, its scope label and its cells, unchanged."""
    title = scope = coordinates = ""
    with path.open(encoding="utf-8", newline="") as handle:
        body = []
        for line in handle:
            if line.startswith("# table:"):
                title = line.split(":", 1)[1].strip()
            elif line.startswith("# scope:"):
                scope, _, coordinates = line.split(":", 1)[1].strip().partition(" | coordinates:")
            else:
                body.append(line)
    rows = list(csv.reader(body))
    columns = rows[0] if rows else []
    cells = rows[1:]
    return {
        "title": title,
        "scope": scope,
        "coordinates": coordinates.strip(),
        "columns": columns,
        "rows": cells[:_TABLE_ROWS],
        "total": len(cells),
    }
