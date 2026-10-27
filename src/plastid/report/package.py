"""The Supplementary package — regenerated from a run directory, never assembled by hand.

A manuscript's supplement is a folder: the tables, the figures, the draft report and a note
saying what is in it. Every file here is written from the run's own stored results, so the
package can be rebuilt at any time and cannot drift from the numbers the run produced.

**The whole folder is replaced, not edited in place** (M11.18). A rebuild writes into a staging
directory beside the package and only then replaces the old one, so a table or a figure that the
run no longer holds cannot survive as if it were current — and a failed rebuild leaves the old
package in place rather than a half-written one. ``index.md`` lists each file with its row count
and its sha256, links it, and names every output that was not produced (a missing extra or a
missing R) with the reason. A silent absence is what makes a supplement untrustworthy.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from plastid.core.deps import ExportDependencyError
from plastid.report.citations import format_citation, load_citations
from plastid.result.export import export_csv
from plastid.result.figure import FigureRecord
from plastid.result.result import Result
from plastid.result.xlsx_export import export_xlsx

__all__ = ["PACKAGE_DIRNAME", "PackageResult", "build_package"]

PACKAGE_DIRNAME = "supplementary"
_WORKBOOK = "tables.xlsx"


def _figure_entries(records: Sequence[FigureRecord]) -> list[dict]:
    """The manifest's figure inventory — status, copies and the provenance that drew them."""
    return [
        {
            "figure": record.figure_id,
            "renderer": record.renderer,
            "status": record.status,
            "files": [path.name for path in record.files],
            "reason": record.reason,
            "provenance": list(record.provenance),
        }
        for record in records
    ]


@dataclass(frozen=True, slots=True)
class PackageResult:
    """Where the package landed, and every file it holds with its size and digest.

    The manifest does not list itself: a file cannot carry its own digest, and pretending
    otherwise would mean writing it twice.
    """

    directory: Path
    files: tuple[Path, ...]
    label: str = ""
    figures: tuple[FigureRecord, ...] = ()
    skipped_outputs: tuple[tuple[str, str], ...] = ()

    def manifest(self) -> dict:
        return {
            "directory": self.directory.name,
            "data_label": self.label,
            "files": [
                {
                    "path": path.relative_to(self.directory).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
                for path in self.files
            ],
            "figures": _figure_entries(self.figures),
            "skipped_outputs": [
                {"output": output, "reason": reason} for output, reason in self.skipped_outputs
            ],
        }


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _link(path: Path, base: Path) -> str:
    return f"[{path.name}]({path.relative_to(base).as_posix()})"


def _citations(results: Mapping[str, Result]) -> list[str]:
    """Every source a result cites, once, resolved through the shipped registry."""
    keys = list(dict.fromkeys(key for result in results.values() for key in result.citations))
    if not keys:
        return []
    registry = load_citations()
    missing = [key for key in keys if key not in registry]
    if missing:
        raise ValueError(
            f"the run cites {', '.join(missing)}, which the citation registry does not hold — "
            f"a supplement cannot print a reference the toolkit cannot name"
        )
    return [format_citation(registry[key]) for key in keys]


def _index(
    results: Mapping[str, Result],
    figure_records: Sequence[FigureRecord],
    skipped_outputs: Sequence[tuple[str, str]],
    table_files: Sequence[tuple[str, str, Path]],
    reports: Sequence[tuple[str, Path]],
    label: str,
    base: Path,
    provenance: Sequence[str] = (),
) -> str:
    """The package's own front page — every file named, linked and accounted for."""
    lines = [
        "# Supplementary package",
        "",
        "Regenerated from this run's stored results — no analysis was re-run and no file was "
        "assembled by hand. Each table below is also a sheet in `" + _WORKBOOK + "` when the "
        "Excel extra is installed.",
        "",
    ]
    if label:
        lines += [f"> **Data:** {label}", ""]
    lines += [
        "## Tables",
        "",
        "| Module | Table | Rows | File | sha256 |",
        "|---|---|---|---|---|",
    ]
    for module_id, table_id, path in table_files:
        digest = _digest(path)
        rows = len(results[module_id].tables[table_id].rows)
        lines.append(
            f"| {module_id} | {table_id} | {rows} | {_link(path, base)} | `{digest[:16]}…` |"
        )

    lines += ["", "## Figures", ""]
    rendered = [record for record in figure_records if record.files]
    if rendered:
        for record in rendered:
            lines.append(f"### `{record.figure_id}` ({record.renderer})")
            lines.append("")
            for path in record.files:
                lines.append(
                    f"- {_link(path, base)} — {path.stat().st_size} bytes, "
                    f"sha256 `{_digest(path)[:16]}…`"
                )
            lines += [f"  - {line}" for line in record.provenance]
            lines.append("")
    elif not figure_records:
        lines.append("- none: no figure's source tables are in this run")
    else:
        lines.append("- none rendered")

    skipped = [record for record in figure_records if not record.files]
    if skipped:
        lines += ["", "Figures not drawn, and why:", ""]
        lines += [
            f"- `{record.figure_id}` ({record.renderer}) — {record.reason}" for record in skipped
        ]

    if skipped_outputs:
        lines += ["", "Outputs not produced, and why:", ""]
        lines += [f"- `{output}` — {reason}" for output, reason in skipped_outputs]

    if reports:
        lines += ["", "## Reports", ""]
        lines += [
            f"- {_link(path, base)} — {path.stat().st_size} bytes, sha256 `{_digest(path)[:16]}…`"
            for _, path in reports
        ]

    lines += [
        "",
        "## Provenance",
        "",
    ]
    if provenance:
        lines += [f"- {line}" for line in provenance]
    lines += [
        "`MANIFEST.json` carries every file with its size and digest, every figure's status and "
        "provenance, and every output that was not produced. `../qc.md` carries what ran, with "
        "which versions, on which inputs.",
        "",
    ]

    citations = _citations(results)
    if citations:
        lines += ["## References", ""]
        lines += [f"- {citation}" for citation in citations]
        lines.append("")
    return "\n".join(lines)


def _copy(source: Path, destination: Path) -> Path:
    if not source.is_file():
        raise ValueError(
            f"the run records {source} as an output, but it is missing — a package cannot "
            f"carry a file that is not there; re-render or rebuild the run first"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.read_bytes())
    return destination


def build_package(
    run_dir: str | Path,
    *,
    results: Mapping[str, Result],
    figure_records: Sequence[FigureRecord] = (),
    missing_renderers: Sequence[tuple[str, str]] = (),
    reports: Sequence[tuple[str, Path]] = (),
    label: str = "",
    provenance: Sequence[str] = (),
) -> PackageResult:
    """Write ``<run_dir>/report/supplementary/`` from the run's stored results.

    ``figure_records`` is the complete figure inventory — rendered and skipped, from both the
    matplotlib and the R catalogue — so the package never drops a figure merely because another
    command drew it. ``missing_renderers`` names a whole engine that could not run (the plotting
    extra, or R); the tables are still a valid package, and the absence is on the index.
    """
    root = Path(run_dir)
    target = root / "report" / PACKAGE_DIRNAME
    staging = root / "report" / f".{PACKAGE_DIRNAME}.building"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    files: list[Path] = []
    table_files: list[tuple[str, str, Path]] = []
    skipped_outputs: list[tuple[str, str]] = list(missing_renderers)

    (staging / "tables").mkdir()
    for module_id, result in results.items():
        for table_id, table in result.tables.items():
            path = staging / "tables" / f"{module_id}__{table_id}.csv"
            path.write_text(export_csv(table), encoding="utf-8")
            files.append(path)
            table_files.append((module_id, table_id, path))

    try:
        workbook = staging / _WORKBOOK
        workbook.write_bytes(export_xlsx(results))
    except ExportDependencyError as error:
        # The tables are text and remain the source of every number; the workbook is the
        # convenience copy. Its absence is named, never silent (M11.18).
        skipped_outputs.append((_WORKBOOK, str(error)))
    else:
        files.append(workbook)

    figure_dir = staging / "figures"
    copied_records: list[FigureRecord] = []
    for record in figure_records:
        if not record.files:
            copied_records.append(record)
            continue
        copied: list[Path] = []
        for source in record.files:
            destination = figure_dir / source.name
            copied.append(_copy(source, destination))
            files.append(destination)
        copied_records.append(
            FigureRecord(
                record.figure_id,
                record.renderer,
                tuple(copied),
                provenance=record.provenance,
            )
        )

    report_dir = staging / "reports"
    copied_reports: list[tuple[str, Path]] = []
    for name, source in reports:
        destination = _copy(source, report_dir / name)
        files.append(destination)
        copied_reports.append((name, destination))

    index = staging / "index.md"
    index.write_text(
        _index(
            results,
            copied_records,
            skipped_outputs,
            table_files,
            copied_reports,
            label,
            staging,
            provenance,
        ),
        encoding="utf-8",
    )
    files.append(index)

    # The manifest describes the package under its final name, while the files still live in the
    # staging directory — computing it from a `PackageResult` here would name the staging folder.
    manifest = staging / "MANIFEST.json"
    manifest.write_text(
        json.dumps(
            {
                "directory": PACKAGE_DIRNAME,
                "data_label": label,
                "files": [
                    {
                        "path": path.relative_to(staging).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                    for path in files
                ],
                "figures": _figure_entries(copied_records),
                "skipped_outputs": [
                    {"output": output, "reason": reason} for output, reason in skipped_outputs
                ],
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    files.append(manifest)

    # Swap the finished folder into place in one step. The old package is removed only now, so a
    # failure above leaves it untouched (and the command has failed, so nobody reads it as new).
    if target.exists():
        shutil.rmtree(target)
    staging.replace(target)
    placed = tuple(target / path.relative_to(staging) for path in files)
    return PackageResult(
        directory=target,
        files=placed,
        label=label,
        figures=tuple(
            FigureRecord(
                item.figure_id,
                item.renderer,
                tuple(target / path.relative_to(staging) for path in item.files),
                provenance=item.provenance,
                reason=item.reason,
            )
            for item in copied_records
        ),
        skipped_outputs=tuple(skipped_outputs),
    )
