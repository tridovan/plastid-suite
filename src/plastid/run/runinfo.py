"""``00_run_info/run_info.md`` — the Materials & Methods page for a run.

Requirement 02/10 (PROJECT.md §2): every run states its tools, versions, parameters, inputs
and the IR verdict, so the methods section can be written from the output folder alone.
"""

from __future__ import annotations

import platform
from collections.abc import Sequence
from typing import Any

from plastid.run.tools import ToolRecord

__all__ = ["render_run_info"]


def render_run_info(manifest: dict[str, Any], tools: Sequence[ToolRecord] = ()) -> str:
    lines: list[str] = [
        f"# Run {manifest['run_id']}",
        "",
        f"- package: plastid-suite {manifest['package_version']}",
        f"- command: {' '.join(manifest['command'])}",
        f"- started: {manifest['started_at']}",
        f"- finished: {manifest.get('finished_at') or '—'}",
        f"- status: {manifest['status']}",
        f"- environment: Python {platform.python_version()} on {platform.platform()}",
        "",
        "## Parameters",
        "",
    ]
    params = manifest.get("params") or {}
    if params:
        lines.extend(f"- {key} = {value}" for key, value in sorted(params.items()))
    else:
        lines.append("- (none)")

    lines.extend(["", "## Inputs", ""])
    inputs = manifest.get("inputs") or []
    if inputs:
        for entry in inputs:
            lines.append(f"- {entry['path']} — sha256 {entry['sha256']} ({entry['bytes']} bytes)")
    else:
        lines.append("- (none)")

    lines.extend(["", "## Tools", ""])
    if tools:
        lines.append("| tool | version | version source | path | exit | seconds |")
        lines.append("|---|---|---|---|---|---|")
        for record in tools:
            lines.append(
                f"| {record.name} | {record.version or '—'} | {record.version_source} | "
                f"{record.path or '—'} | {record.rc} | {record.duration_s} |"
            )
    else:
        lines.append("- (no external tools ran)")

    lines.extend(["", "## IR gate", ""])
    ir = manifest.get("ir")
    if ir:
        verdict = "passed" if ir.get("gate_passed") else "FAILED — position-based analysis blocked"
        lines.append(
            f"- source: {ir.get('source')} · length: {ir.get('length')} bp · "
            f"identity: {float(ir.get('identity', 0)) * 100:.2f} % — {verdict}"
        )
        if ir.get("denovo_length") is not None:
            lines.append(
                f"- de-novo cross-check: {ir['denovo_length']} bp "
                f"(delta {ir.get('denovo_delta')} bp)"
            )
    else:
        lines.append("- (not required by this module)")

    lines.extend(["", "## Coordinates", ""])
    coordinates = manifest.get("coordinate_system") or {}
    lines.append(
        f"- per-genome tables: {coordinates.get('name', 'genbank')} coordinates · "
        f"rotation offset to the normalized frame: {coordinates.get('rotation_offset', 0)}"
    )
    lines.append(
        "- cross-genome comparisons use the normalized frame (LSC → IRa → SSC → IRb) — "
        "every table states its own frame"
    )

    lines.extend(["", "## Modules", ""])
    for entry in manifest.get("modules") or []:
        lines.append(
            f"- {entry['id']} {entry['version']} — {entry['status']} "
            f"({len(entry.get('files', []))} file(s))"
        )
    return "\n".join(lines) + "\n"
