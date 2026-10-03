"""Deterministic exporters: csv · json · markdown.

Byte-identical across runs — fixed column order, fixed float formatting, ``\\n`` endings,
no locale. Every export carries its scope label and coordinate frame.
"""

from __future__ import annotations

import csv
import io
import json

from plastid.result.result import Result, Table

__all__ = ["export_csv", "export_json", "export_markdown"]


def export_csv(table: Table) -> str:
    """One table as CSV; the first two lines carry provenance and are skippable comments."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow([column.header for column in table.columns])
    for row in table.rows:
        writer.writerow([column.format(row.get(column.name)) for column in table.columns])
    comments = (
        f"# table: {table.title}\n"
        f"# scope: {table.scope.label} | coordinates: {table.scope.coord_system}\n"
    )
    return comments + buffer.getvalue()


def export_json(result: Result) -> str:
    """The whole result as JSON, raw values preserved (formatting stays a display concern)."""
    payload = {
        "module": result.module,
        "module_version": result.module_version,
        "tables": {
            table_id: {
                "title": table.title,
                "scope": {
                    "label": table.scope.label,
                    "coord_system": table.scope.coord_system,
                    "include": list(table.scope.include),
                    "exclude_regions": list(table.scope.exclude_regions),
                },
                "columns": [
                    {
                        "name": column.name,
                        "label": column.label,
                        "dtype": column.dtype,
                        "unit": column.unit,
                        "precision": column.precision,
                    }
                    for column in table.columns
                ],
                "key": table.key,
                "rows": [
                    {column.name: row.get(column.name) for column in table.columns}
                    for row in table.rows
                ],
            }
            for table_id, table in result.tables.items()
        },
        "notes": list(result.notes),
        "citations": list(result.citations),
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def export_markdown(result: Result) -> str:
    """The whole result as markdown fragments (for the report stitch-up)."""
    parts: list[str] = []
    for table in result.tables.values():
        parts.append(f"## {table.title}\n")
        parts.append(f"*scope: {table.scope.label} | coordinates: {table.scope.coord_system}*\n")
        parts.append("| " + " | ".join(column.header for column in table.columns) + " |")
        parts.append("|" + "|".join("---" for _ in table.columns) + "|")
        for row in table.rows:
            parts.append(
                "| "
                + " | ".join(column.format(row.get(column.name)) for column in table.columns)
                + " |"
            )
        parts.append("")
    for note in result.notes:
        parts.append(note)
    if result.notes:
        parts.append("")
    for citation in result.citations:
        parts.append(f"- {citation}")
    return "\n".join(parts).rstrip("\n") + "\n"
