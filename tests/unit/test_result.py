"""Result schema: tables, cell addressing, deterministic exports.

Decision §10.7 (PROJECT.md §3): every number in a report resolves to exactly one table
cell; an unresolvable reference raises — it never prints blank.
"""

from __future__ import annotations

import csv
import io
import json

import pytest

from plastid.result.export import export_csv, export_json, export_markdown
from plastid.result.result import Column, Result, ResultRefError, Scope, Table


def _single_row_table() -> Table:
    return Table(
        id="genome_stats",
        title="Genome statistics",
        columns=(
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=4),
            Column("accession", "Accession"),
        ),
        rows=({"length_bp": 153_724, "gc_pct": 37.8536, "accession": "NC_050833.1"},),
        scope=Scope(label="whole plastome; both IR copies", coord_system="genbank"),
    )


def _multi_row_table() -> Table:
    return Table(
        id="ssr_by_unit",
        title="SSR by unit",
        columns=(Column("unit", "Unit"), Column("n", "Count", dtype="int")),
        rows=({"unit": "mono", "n": 31}, {"unit": "di", "n": 4}, {"unit": "tri", "n": 2}),
        scope=Scope(label="whole plastome, MISA defaults"),
        key="unit",
    )


def _result() -> Result:
    return Result(
        module="genome_stats",
        module_version="0.1.0",
        tables={"genome_stats": _single_row_table(), "ssr_by_unit": _multi_row_table()},
        notes=["Draft — needs human review."],
    )


def test_a_repeated_key_is_rejected() -> None:
    """A repeated key would make every reference to it silently resolve to one of the rows."""
    with pytest.raises(ValueError, match="repeats"):
        Table(
            id="repeated",
            title="Repeated",
            columns=(Column("k", "K"), Column("v", "V", dtype="int")),
            rows=({"k": "a", "v": 1}, {"k": "a", "v": 2}),
            scope=Scope(label="x"),
            key="k",
        )


def test_emphasis_must_address_a_row_that_exists() -> None:
    """A typo here is invisible in every export except the Word document's weight."""
    table = _multi_row_table()
    emphasised = Table(
        id=table.id,
        title=table.title,
        columns=table.columns,
        rows=table.rows,
        scope=table.scope,
        key="unit",
        emphasis=("tri",),
    )
    assert emphasised.emphasis == ("tri",)

    with pytest.raises(ValueError, match="not a value of key column"):
        Table(
            id=table.id,
            title=table.title,
            columns=table.columns,
            rows=table.rows,
            scope=table.scope,
            key="unit",
            emphasis=("tetra",),
        )


def test_emphasis_needs_a_key_column_to_address_rows_by() -> None:
    with pytest.raises(ValueError, match="needs a key column"):
        Table(
            id="single",
            title="Single",
            columns=(Column("v", "V"),),
            rows=({"v": 1},),
            scope=Scope(label="x"),
            emphasis=("Total",),
        )


def test_single_row_table_addresses_by_column() -> None:
    assert _single_row_table().value("gc_pct") == 37.8536


def test_multi_row_table_addresses_by_row_key() -> None:
    assert _multi_row_table().value("mono.n") == 31
    assert _multi_row_table().value("tri.n") == 2


def test_unknown_references_raise() -> None:
    with pytest.raises(ResultRefError):
        _single_row_table().value("nope")
    with pytest.raises(ResultRefError):  # unknown row key
        _multi_row_table().value("tetra.n")
    with pytest.raises(ResultRefError):  # multi-row table needs a row key
        _multi_row_table().value("n")
    with pytest.raises(ResultRefError):  # unknown table
        _result().value("other.gc_pct")
    with pytest.raises(ResultRefError):  # unknown column inside a known row
        _result().value("ssr_by_unit.mono.nope")


def test_result_addresses_tables() -> None:
    result = _result()
    assert result.value("genome_stats.gc_pct") == 37.8536
    assert result.value("ssr_by_unit.mono.n") == 31


def test_ids_may_not_contain_dots() -> None:
    with pytest.raises(ValueError, match="dot"):
        Table(
            id="a.b",
            title="broken",
            columns=(Column("x", "X"),),
            rows=({"x": 1},),
            scope=Scope(label="s"),
        )


def test_csv_export_is_deterministic_and_carries_scope() -> None:
    first = export_csv(_single_row_table())
    assert first == export_csv(_single_row_table())  # byte-identical re-export

    lines = first.splitlines()
    assert lines[0] == "# table: Genome statistics"
    assert lines[1] == "# scope: whole plastome; both IR copies | coordinates: genbank"
    assert lines[2] == "Length (bp),GC (%),Accession"  # human labels + units

    data = "\n".join(line for line in lines if not line.startswith("#"))
    parsed = list(csv.reader(io.StringIO(data)))
    assert parsed[0] == ["Length (bp)", "GC (%)", "Accession"]
    assert parsed[1] == ["153724", "37.8536", "NC_050833.1"]  # fixed float format, round-trips


def test_json_export_round_trips_with_scope() -> None:
    text = export_json(_result())
    assert text == export_json(_result())
    payload = json.loads(text)
    assert payload["module"] == "genome_stats"
    table = payload["tables"]["genome_stats"]
    assert table["scope"]["label"] == "whole plastome; both IR copies"
    assert table["scope"]["coord_system"] == "genbank"
    assert table["rows"][0]["gc_pct"] == 37.8536
    assert payload["notes"] == ["Draft — needs human review."]


def test_markdown_export_contains_title_scope_and_rows() -> None:
    text = export_markdown(_result())
    assert text == export_markdown(_result())
    assert "## Genome statistics" in text
    assert "scope: whole plastome; both IR copies | coordinates: genbank" in text
    assert "| Length (bp) | GC (%) | Accession |" in text
    assert "| 153724 | 37.8536 | NC_050833.1 |" in text
    assert "Draft — needs human review." in text
