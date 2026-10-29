"""The technical QC report: what ran, with which versions, on which inputs, what was skipped.

Written for a reader who has to decide whether to trust a run's numbers. Everything in it comes
from the run's own ``manifest.json`` — the tool records, the input hashes, the parameters, the IR
gate, the coordinate frame — plus the queryable-number audit over each report, so a run that
carried an orphan number would say so here rather than in a reader's suspicion.
"""

from __future__ import annotations

import json
from pathlib import Path

from plastid.report.audit import audit_run
from plastid.result.figure import FigureRecord

__all__ = ["QC_FILENAME", "render_qc", "render_qc_for_run"]

QC_FILENAME = "qc.md"


def _table(rows: list[tuple[str, ...]], header: tuple[str, ...]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return lines


def _figure_section(
    figure_records: tuple[FigureRecord, ...],
    figures: tuple[str, ...],
    skipped: tuple[tuple[str, str], ...],
    missing_renderers: tuple[tuple[str, str], ...],
) -> list[str]:
    """Every catalogue entry with its status — never just the files that happened to exist."""
    lines = ["", "## Figures", ""]
    if figure_records:
        lines += _table(
            [
                (
                    record.figure_id,
                    record.renderer,
                    record.status,
                    ", ".join(path.name for path in record.files),
                    record.reason,
                )
                for record in figure_records
            ],
            ("Figure", "Renderer", "Status", "Files", "Reason"),
        )
        provenance = [line for record in figure_records for line in record.provenance]
        if provenance:
            lines += ["", "Figure provenance:", ""]
            lines += [f"- {line}" for line in dict.fromkeys(provenance)]
    else:
        lines += [f"- {name}" for name in figures] or ["- none rendered"]
        lines += [f"- skipped `{figure_id}` — {reason}" for figure_id, reason in skipped]
    if missing_renderers:
        lines += ["", "Renderers unavailable:", ""]
        lines += [f"- {renderer} — {reason}" for renderer, reason in missing_renderers]
    return lines


def render_qc(
    manifest: dict,
    *,
    figures: tuple[str, ...] = (),
    skipped: tuple[tuple[str, str], ...] = (),
    audits: tuple[tuple[str, bool, int], ...] = (),
    figure_records: tuple[FigureRecord, ...] = (),
    missing_renderers: tuple[tuple[str, str], ...] = (),
    data_label: str = "",
) -> str:
    """The QC report as markdown — every line traceable to the manifest or to an audit."""
    lines = [
        f"# Technical QC report — run {manifest['run_id']}",
        "",
        f"- command: `{' '.join(manifest['command'])}`",
        f"- package: plastid-suite {manifest['package_version']}",
        f"- started: {manifest['started_at']} · finished: {manifest.get('finished_at')}",
        f"- status: {manifest['status']}",
    ]
    if data_label:
        lines.append(f"- data label: {data_label} (third-party data, not the toolkit's own)")
    lines += [
        "",
        "## Inputs",
        "",
    ]

    lines += _table(
        [
            (entry["path"], entry["sha256"], str(entry["bytes"]))
            for entry in manifest.get("inputs", [])
        ],
        ("Input", "sha256", "bytes"),
    )

    lines += ["", "## Modules", ""]
    lines += _table(
        [
            (
                entry["id"],
                entry["version"],
                entry["status"],
                str(len(entry.get("files", []))),
                entry.get("error", ""),
            )
            for entry in manifest.get("modules", [])
        ],
        ("Module", "Version", "Status", "Files", "Error"),
    )

    tools = [
        (entry["id"], tool["name"], tool["version"], tool["path"], str(tool["exit_code"]))
        for entry in manifest.get("modules", [])
        for tool in entry.get("tools", [])
    ]
    lines += ["", "## Tools", ""]
    lines += (
        _table(tools, ("Module", "Tool", "Version", "Path", "Exit"))
        if tools
        else ["No external tool was called by this run."]
    )

    lines += ["", "## Parameters", ""]
    lines += [f"- `{key}` = `{value}`" for key, value in sorted(manifest.get("params", {}).items())]

    gate = manifest.get("ir") or {}
    lines += ["", "## IR gate and coordinates", ""]
    if gate:
        lines.append(
            f"- IR gate: {'passed' if gate.get('gate_passed') else 'FAILED'} · "
            f"length {gate.get('length')} bp · identity {gate.get('identity')} · "
            f"source {gate.get('source')}"
        )
        if gate.get("denovo_length") is not None:
            lines.append(
                f"- de-novo detector: {gate['denovo_length']} bp (Δ {gate.get('denovo_delta')})"
            )
    else:
        lines.append("- no IR gate: the modules in this run are not position-based")
    coordinates = manifest.get("coordinate_system", {})
    lines.append(
        f"- coordinates: {coordinates.get('name')} · rotation offset "
        f"{coordinates.get('rotation_offset')}"
    )

    lines += _figure_section(figure_records, figures, skipped, missing_renderers)

    lines += ["", "## Queryable-number audit", ""]
    if audits:
        lines += _table(
            [
                (module, "clean" if ok else "FINDINGS", str(numbers))
                for module, ok, numbers in audits
            ],
            ("Module", "Audit", "Numbers checked"),
        )
        lines.append("")
        lines.append(
            "Every number in a report must be a bounded number inside a rendered cell of that "
            "module's own tables; the audit sets aside only the citation registry, the run's data "
            "label and the draft banner."
        )
    else:
        lines.append("No module report was audited (none was found in this run).")

    lines += [
        "",
        "## Files",
        "",
    ]
    lines += _table(
        [
            (entry["path"], entry["module"], entry["sha256"][:16] + "…", str(entry["bytes"]))
            for entry in manifest.get("files", [])
        ],
        ("Path", "Module", "sha256", "bytes"),
    )
    return "\n".join(lines) + "\n"


def render_qc_for_run(
    run_dir: str | Path,
    *,
    figures: tuple[str, ...] = (),
    skipped: tuple[tuple[str, str], ...] = (),
    figure_records: tuple[FigureRecord, ...] = (),
    missing_renderers: tuple[tuple[str, str], ...] = (),
) -> tuple[str, Path]:
    """Read a run's manifest, audit its reports, and write ``report/qc.md``."""
    root = Path(run_dir)
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audits = tuple(
        (report.module, report.audit.ok, report.audit.n_numbers) for report in audit_run(root)
    )
    label = str(manifest.get("params", {}).get("data_label", ""))
    text = render_qc(
        manifest,
        figures=figures,
        skipped=skipped,
        audits=audits,
        figure_records=figure_records,
        missing_renderers=missing_renderers,
        data_label=label,
    )
    target = root / "report" / QC_FILENAME
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return text, target
