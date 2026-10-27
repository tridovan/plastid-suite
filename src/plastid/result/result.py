"""Normalized result schema.

Every number a report shows comes from a table cell; an unresolvable reference raises
:class:`ResultRefError` — it is never printed blank. Cell references:

* ``table.column`` — a single-row table;
* ``table.rowkey.column`` — a multi-row table, addressed through its declared key column.
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
        """The fixed Ts/Tv scope: CDS only, IR excluded."""
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
    # Key values to render emphasised — a totals row, in every table that has one. Presentation
    # only: it changes no value and no export other than the Word document's weight.
    emphasis: tuple[str, ...] = ()
    # A module may build a table shape before it knows whether there is anything to put in it,
    # then omit it from the Result when there is not (``repeats`` and ``ssr`` both do). An empty
    # table must never reach a Result — Result.__post_init__ refuses it — because an exported
    # empty table is indistinguishable from a table nobody filled.
    allow_empty: bool = False

    def __post_init__(self) -> None:
        if "." in self.id:
            raise ValueError(f"table id must not contain a dot: {self.id!r}")
        if not self.rows and not self.allow_empty:
            raise ValueError(f"table {self.id!r} has no rows")
        names = {column.name for column in self.columns}
        if self.key is not None and self.key not in names:
            raise ValueError(f"table {self.id!r}: key column {self.key!r} is not in its columns")
        for row in self.rows:
            unknown = set(row) - names
            if unknown:
                raise ValueError(f"table {self.id!r}: row has unknown columns {sorted(unknown)}")
        if self.key is not None:
            values = [row.get(self.key) for row in self.rows]
            repeated = sorted({value for value in values if values.count(value) > 1}, key=str)
            if repeated:
                raise ValueError(
                    f"table {self.id!r}: key column {self.key!r} repeats {repeated} — a "
                    f"reference to a repeated key would silently resolve to one of the rows"
                )
        if self.emphasis and self.key is None:
            raise ValueError(f"table {self.id!r}: emphasis needs a key column to address rows by")
        for value in self.emphasis:
            if self.key is not None and value not in {row.get(self.key) for row in self.rows}:
                raise ValueError(
                    f"table {self.id!r}: emphasis {value!r} is not a value of key column "
                    f"{self.key!r}"
                )

    def _locate(self, ref: str) -> tuple[dict[str, Any], Column]:
        """Resolve ``column`` or ``rowkey.column`` to (row, column)."""
        if self.key is None:
            if len(self.rows) > 1:
                raise ResultRefError(
                    f"{self.id}: {len(self.rows)} rows and no key column — "
                    f"a single value cannot be addressed"
                )
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
        columns = {column.name: column for column in self.columns}
        if column_name not in columns:
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
        return row, columns[column_name]

    def value(self, ref: str) -> Any:
        """The raw value at a reference."""
        row, column = self._locate(ref)
        return row[column.name]

    def format(self, ref: str) -> str:
        """The display string at a reference — the same rendering the exporters use."""
        row, column = self._locate(ref)
        return column.format(row[column.name])


@dataclass(frozen=True, slots=True)
class Result:
    """Everything one module produced: tables, notes, citations, and file artifacts.

    ``artifacts`` are files that are not tables — a submission ``.fsa``/``.tbl``, a tool's
    output — carried as ``{"<name>.<ext>": text}``. They travel through the same export and
    cache path as every other output, so they are written, listed in the manifest, and restored
    on a cache hit exactly like a table. A module never writes a file itself (seam 2).
    """

    module: str
    module_version: str
    tables: dict[str, Table] = field(default_factory=dict)
    notes: tuple[str, ...] = ()
    citations: tuple[str, ...] = ()
    artifacts: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for table_id, table in self.tables.items():
            if table.id != table_id:
                raise ValueError(f"tables[{table_id!r}].id is {table.id!r} — they must match")
            if not table.rows:
                raise ValueError(
                    f"tables[{table_id!r}] is empty — a module that measured nothing omits the "
                    f"table and says so in its notes"
                )

    def _table_for(self, ref: str) -> tuple[Table, str]:
        table_id, separator, rest = ref.partition(".")
        if not separator:
            raise ResultRefError(f"reference {ref!r} needs a table — use table.column")
        table = self.tables.get(table_id)
        if table is None:
            raise ResultRefError(f"unknown table {table_id!r} (have: {sorted(self.tables)})")
        return table, rest

    def value(self, ref: str) -> Any:
        """Resolve ``table.column`` or ``table[row].column`` across this result."""
        table, rest = self._table_for(ref)
        return table.value(rest)

    def format(self, ref: str) -> str:
        """The display string at a result-level reference."""
        table, rest = self._table_for(ref)
        return table.format(rest)
