"""``00_run_info/run_info.md`` — the Materials & Methods page for a run.

Requirement 02/10: every run states its tools, versions, parameters, inputs
and the IR verdict, so the methods section can be written from the output folder alone.

**The tool table comes from the manifest, not from one context.** Each module runs with its own
`RunContext`, so a page rendered from a single context named only the tools of whichever module ran
last — a run of `… ssr` said "no external tools ran" while MISA had run (K-08, M11.8). The
manifest already records each module's calls; this page reads those records, so no ordering can
hide a tool. A tool called by several modules appears once, because the question the table answers
is "what ran", not "how many times".
"""

from __future__ import annotations

import platform
from collections.abc import Sequence
from typing import Any

from plastid.run.tools import ToolRecord

__all__ = ["render_run_info"]


def _tools_from_manifest(
    manifest: dict[str, Any], current: Sequence[ToolRecord] = ()
) -> list[ToolRecord]:
    """Every tool call the run recorded, one row per distinct tool, in call order.

    The manifest's per-module ``tools`` are plain dicts (it is written as JSON); the caller may
    also pass the live records. Both are read, so a page rendered before the manifest is final
    and one rendered from a finished manifest agree. Identity is ``(name, version, path)``: the
    same tool at two versions is two rows, because that is two measurements.
    """
    rows: list[ToolRecord] = []
    seen: set[tuple[str, str | None, str | None]] = set()

    def add(record: ToolRecord) -> None:
        key = (record.name, record.version, record.path)
        if key in seen:
            return
        seen.add(key)
        rows.append(record)

    for record in current:
        add(record)
    for module in manifest.get("modules") or []:
        for entry in module.get("tools") or []:
            if isinstance(entry, ToolRecord):
                add(entry)
                continue
            add(
                ToolRecord(
                    name=str(entry.get("name", "?")),
                    kind=str(entry.get("kind", "binary")),
                    version=entry.get("version"),
                    version_source=str(entry.get("version_source", "unknown")),
                    path=entry.get("path"),
                    argv=tuple(entry.get("argv") or ()),
                    cwd=str(entry.get("cwd", "")),
                    rc=int(entry.get("rc", 0)),
                    started_at=str(entry.get("started_at", "")),
                    duration_s=float(entry.get("duration_s", 0.0)),
                    stdout_file=str(entry.get("stdout_file", "")),
                    stderr_file=str(entry.get("stderr_file", "")),
                    stdout_tail=str(entry.get("stdout_tail", "")),
                    interpreter=entry.get("interpreter"),
                )
            )
    return rows


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
    recorded = _tools_from_manifest(manifest, tools)
    if recorded:
        lines.append("| tool | version | version source | path | exit | seconds |")
        lines.append("|---|---|---|---|---|---|")
        for record in recorded:
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
        line = (
            f"- {entry['id']} {entry['version']} — {entry['status']} "
            f"({len(entry.get('files', []))} file(s))"
        )
        # A cache hit is not a computation that ran now. The Tools table above lists only this
        # run's calls, so the original computation is named here — otherwise "no external tools
        # ran" would read as "no tool was involved" (M11.12).
        provenance = entry.get("provenance") or {}
        if entry.get("status") == "cached" and provenance:
            when = provenance.get("computed_at") or "an earlier run"
            tools = provenance.get("tools") or []
            who = (
                ", ".join(
                    f"{tool.get('name')} {tool.get('version') or 'version unknown'}"
                    for tool in tools
                )
                or "no external tools"
            )
            line += f" — payload computed {when} by {who}"
        lines.append(line)
    return "\n".join(lines) + "\n"
