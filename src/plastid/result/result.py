"""Normalized result schema (decision §10.7, PROJECT.md §3).

Every number a report shows comes from a table cell; an unresolvable reference raises
:class:`ResultRefError` — it is never printed blank. Cell references:

* ``table.column`` — a single-row table;
* ``table[row].column`` — a multi-row table, addressed through its declared key column.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

__all__ = ["Column", "Result", "ResultRefError", "Scope", "Table"]


class ResultRefError(Exception):
    """A cell reference did not resolve."""


@dataclass(frozen=True, slots=True)
class Scope:
    """What a table covers — attached to every export so no number travels unlabelled."""

    label: str
    coord_system: str = "genbank"  # "genbank" | "normalized"
    include: tuple[str, ...] = ()
    exclude_regions: tuple[str, ...] = ()

    @classmethod
    def cds_no_ir(cls) -> Scope:
        """The fixed Ts/Tv scope (PROJECT.md §6): CDS only, IR excluded."""
        return cls(label="CDS only; IR excluded", include=("CDS",), exclude_regions=("IRa", "IRb"))


@dataclass(frozen=True, slots=True)
class Column:
    """One table column: machine name, human label, and how to render its values."""

    name: str
    label: str
    dtype: str = "str"  # "int" | "float" | "str"
    unit: str | None = None
    precision: int | None = None  # fixed decimals when dtype is float

    def __post_init__(self) -> None:
        if "." in self.name:
            raise ValueError(f"column name must not contain a dot: {self.name!r}")

    @property
    def header(self) -> str:
        return f"{self.label} ({self.unit})" if self.unit else self.label

    def format(self, value: Any) -> str:
        """Deterministic text rendering — no locale, no surprises."""
        if value is None:
            return ""
        if self.dtype == "float" and isinstance(value, float) and self.precision is not None:
            return f"{value:.{self.precision}f}"
        return str(value)


@dataclass(frozen=True, slots=True)
class Table:
    """A table of rows plus the scope it covers. Multi-row tables declare a ``key`` column."""

    id: str
    title: str
    columns: tuple[Column, ...]
    rows: tuple[dict[str, Any], ...]
    scope: Scope
    key: str | None = None

    def __post_init__(self) -> None:
        if "." in self.id:
            raise ValueError(f"table id must not contain a dot: {self.id!r}")
        if not self.rows:
            raise ValueError(f"table {self.id!r} has no rows")
        names = {column.name for column in self.columns}
        if self.key is not None and self.key not in names:
            raise ValueError(f"table {self.id!r}: key column {self.key!r} is not in its columns")
        if self.key is None and len(self.rows) > 1:
            raise ValueError(f"table {self.id!r} has several rows — declare a key column")
        for row in self.rows:
            unknown = set(row) - names
            if unknown:
                raise ValueError(f"table {self.id!r}: row has unknown columns {sorted(unknown)}")

    def value(self, ref: str) -> Any:
        """Resolve ``column`` (single row) or ``rowkey.column`` (multi-row)."""
        if self.key is None:
            if "." in ref:
                raise ResultRefError(
                    f"{self.id}: single-row table takes a bare column, got {ref!r}"
                )
            row_key, column_name = None, ref
        else:
            row_key, separator, column_name = ref.rpartition(".")
            if not separator:
                raise ResultRefError(
                    f"{self.id}: reference {ref!r} needs a row key — use {self.key}.column"
                )
        if column_name not in {column.name for column in self.columns}:
            raise ResultRefError(f"{self.id}: unknown column {column_name!r}")

        if row_key is None:
            row = self.rows[0]
        else:
            matches = [candidate for candidate in self.rows if candidate.get(self.key) == row_key]
            if not matches:
                raise ResultRefError(
                    f"{self.id}: unknown row key {row_key!r} (key column {self.key})"
                )
            row = matches[0]
        if column_name not in row:
            raise ResultRefError(f"{self.id}[{row_key}]: column {column_name!r} has no value")
        return row[column_name]


@dataclass(frozen=True, slots=True)
class Result:
    """Everything one module produced: tables, notes, citations."""

    module: str
    module_version: str
    tables: dict[str, Table] = field(default_factory=dict)
    notes: tuple[str, ...] = ()
    citations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for table_id, table in self.tables.items():
            if table.id != table_id:
                raise ValueError(f"tables[{table_id!r}].id is {table.id!r} — they must match")

    def value(self, ref: str) -> Any:
        """Resolve ``table.column`` or ``table[row].column`` across this result."""
        table_id, separator, rest = ref.partition(".")
        if not separator:
            raise ResultRefError(f"reference {ref!r} needs a table — use table.column")
        table = self.tables.get(table_id)
        if table is None:
            raise ResultRefError(f"unknown table {table_id!r} (have: {sorted(self.tables)})")
        return table.value(rest)
