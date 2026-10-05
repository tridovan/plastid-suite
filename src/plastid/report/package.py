"""The Supplementary package — regenerated from a run directory, never assembled by hand.

A manuscript's supplement is a folder: the tables, the figures, and a note saying what is in it.
Every file here is written from the run's own stored results, so the package can be rebuilt at any
time and cannot drift from the numbers the run produced. ``index.md`` lists each file with its
row count and its sha256, which is what makes the folder checkable rather than merely tidy.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from plastid.result.export import export_csv
from plastid.result.result import Result
from plastid.result.xlsx_export import export_xlsx

__all__ = ["PACKAGE_DIRNAME", "PackageResult", "build_package"]

PACKAGE_DIRNAME = "supplementary"
_WORKBOOK = "tables.xlsx"


@dataclass(frozen=True, slots=True)
class PackageResult:
    """Where the package landed, and every file it holds with its size and digest.

    The manifest does not list itself: a file cannot carry its own digest, and pretending
    otherwise would mean writing it twice.
    """

    directory: Path
    files: tuple[Path, ...]

    def manifest(self) -> dict:
        return {
            "directory": self.directory.name,
            "files": [
                {
                    "path": path.relative_to(self.directory).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
                for path in self.files
            ],
        }


def _index(
    results: Mapping[str, Result],
    figures: Sequence[Path],
    skipped: Sequence[tuple[str, str]],
    table_files: Sequence[tuple[str, str, Path]],
) -> str:
    """The package's own front page — every file named, with what it holds."""
    lines = [
        "# Supplementary package",
        "",
        "Regenerated from this run's stored results — no analysis was re-run and no file was "
        "assembled by hand. Each table below is also a sheet in `" + _WORKBOOK + "`.",
        "",
        "## Tables",
        "",
        "| Module | Table | Rows | File | sha256 |",
        "|---|---|---|---|---|",
    ]
    for module_id, table_id, path in table_files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows = len(results[module_id].tables[table_id].rows)
        lines.append(f"| {module_id} | {table_id} | {rows} | `{path.name}` | `{digest[:16]}…` |")
    lines += ["", "## Figures", ""]
    if figures:
        for path in figures:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            lines.append(f"- `{path.name}` — {path.stat().st_size} bytes, sha256 `{digest[:16]}…`")
    else:
        lines.append("- none: no figure's source tables are in this run")
    if skipped:
        lines += ["", "Figures not drawn, and why:", ""]
        lines += [f"- `{figure_id}` — {reason}" for figure_id, reason in skipped]
    lines += [
        "",
        "## Provenance",
        "",
        "`MANIFEST.json` carries every file with its size and digest; `../qc.md` carries what ran, "
        "with which versions, on which inputs.",
        "",
    ]
    return "\n".join(lines)


def build_package(
    run_dir: str | Path,
    *,
    results: Mapping[str, Result],
    figures: Sequence[Path] = (),
    skipped: Sequence[tuple[str, str]] = (),
) -> PackageResult:
    """Write ``<run_dir>/report/supplementary/`` from the run's stored results."""
    root = Path(run_dir)
    target = root / "report" / PACKAGE_DIRNAME
    (target / "tables").mkdir(parents=True, exist_ok=True)

    files: list[Path] = []
    table_files: list[tuple[str, str, Path]] = []
    for module_id, result in results.items():
        for table_id, table in result.tables.items():
            path = target / "tables" / f"{module_id}__{table_id}.csv"
            path.write_text(export_csv(table), encoding="utf-8")
            files.append(path)
            table_files.append((module_id, table_id, path))

    workbook = target / _WORKBOOK
    workbook.write_bytes(export_xlsx(results))
    files.append(workbook)

    figure_dir = target / "figures"
    copied: list[Path] = []
    for source in figures:
        figure_dir.mkdir(parents=True, exist_ok=True)
        destination = figure_dir / source.name
        destination.write_bytes(source.read_bytes())
        files.append(destination)
        copied.append(destination)

    index = target / "index.md"
    index.write_text(_index(results, copied, skipped, table_files), encoding="utf-8")
    files.append(index)

    record = PackageResult(directory=target, files=tuple(files))
    manifest = target / "MANIFEST.json"
    manifest.write_text(
        json.dumps(record.manifest(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    files.append(manifest)
    return PackageResult(directory=target, files=tuple(files))
