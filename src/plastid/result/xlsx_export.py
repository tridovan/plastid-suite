"""Multi-sheet Excel export — one sheet per table, as reproducible as every other export.

The supplement of a manuscript is usually one workbook a reviewer can open; this is that workbook.
It carries the same cells the CSVs carry, through the same :meth:`Column.format`, so a number in
the workbook is the number in the table and not a second rendering of it.

Deterministic for the same reason the Word document is: ``openpyxl`` writes a zip whose entries
carry the wall-clock time and a ``docProps/core.xml`` with created/modified stamps. Both are
pinned here — fixed properties, then the archive rebuilt entry by entry.

``openpyxl`` is an optional dependency (``pip install 'plastid-suite[excel]'``), imported inside
the function.
"""

from __future__ import annotations

import hashlib
import io
import re
from collections.abc import Mapping
from datetime import datetime
from typing import Any

from plastid.core.deps import require_module
from plastid.result.package import deterministic_zip, zip_parts
from plastid.result.result import Result, Table

__all__ = ["SHEET_NAME_LIMIT", "export_xlsx", "sheet_name"]

# Excel's own limits: 31 characters, and none of these characters in a sheet name.
SHEET_NAME_LIMIT = 31
_FORBIDDEN = re.compile(r"[\[\]:*?/\\]")
_EPOCH = datetime(1980, 1, 1, 0, 0, 0)
_CREATOR = "PlastidToolSuite"


def sheet_name(module_id: str, table_id: str, taken: set[str] | None = None) -> str:
    """``<module>.<table>`` as an Excel sheet name — sanitized, truncated, and unique.

    A truncated name gets a short digest of the full one so two long tables cannot collapse onto
    the same sheet and silently overwrite each other.
    """
    raw = _FORBIDDEN.sub("_", f"{module_id}.{table_id}")
    if len(raw) > SHEET_NAME_LIMIT:
        digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:6]  # noqa: S324 — a name, not security
        raw = f"{raw[: SHEET_NAME_LIMIT - 7]}-{digest}"
    name = raw
    suffix = 2
    while taken is not None and name in taken:
        tail = f"-{suffix}"
        name = raw[: SHEET_NAME_LIMIT - len(tail)] + tail
        suffix += 1
    if taken is not None:
        taken.add(name)
    return name


def _write_sheet(sheet: Any, table: Table) -> None:
    sheet.append([column.header for column in table.columns])
    for row in table.rows:
        sheet.append([column.format(row.get(column.name)) for column in table.columns])
    sheet.freeze_panes = "A2"
    for index, column in enumerate(table.columns, start=1):
        width = max(
            len(column.header),
            *(len(column.format(row.get(column.name))) for row in table.rows),
        )
        sheet.column_dimensions[sheet.cell(row=1, column=index).column_letter].width = min(
            48, max(10, width + 2)
        )


def export_xlsx(results: Mapping[str, Result]) -> bytes:
    """Every table of every result in one workbook — a sheet each, in catalogue order."""
    openpyxl = require_module("openpyxl", extra="excel", purpose="Excel export")
    workbook = openpyxl.Workbook()
    workbook.remove(workbook.active)
    taken: set[str] = set()
    for module_id, result in results.items():
        for table in result.tables.values():
            _write_sheet(workbook.create_sheet(title=sheet_name(module_id, table.id, taken)), table)
    if not workbook.sheetnames:  # a workbook with no sheets cannot be saved
        workbook.create_sheet(title="(no tables)")
    workbook.properties.creator = _CREATOR
    workbook.properties.lastModifiedBy = _CREATOR
    workbook.properties.created = _EPOCH
    workbook.properties.modified = _EPOCH
    buffer = io.BytesIO()
    workbook.save(buffer)
    # The properties are pinned, but the archive's own entries still carry "now" — rebuild them.
    return deterministic_zip(zip_parts(buffer.getvalue()))
