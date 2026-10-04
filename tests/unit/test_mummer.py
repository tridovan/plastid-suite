"""The MUMmer adapter — argv, version refusal, and the two readings of one run (P4-01).

The tools are stubbed (``make_mummer_stub``): what is under test is the adapter's contract, not
MUMmer's algorithm — the real numbers are checked in ``tests/golden/test_tier_b.py`` on a machine
that has the tool.

The one contract worth its own test is :func:`run_nucmer`'s refusal when ``show-coords`` and the
delta file describe different alignments: the module publishes a table and a curve side by side,
and a run where the two disagree would produce a figure nobody can check.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import make_mummer_stub
from plastid.external.mummer import MUMMER_VERSION, parse_coords, run_nucmer
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

# Two rows in the shape show-coords -rcl prints: one forward, one reverse (the IR copy). The
# header row is one line in the tool's output, split here only to keep the source readable.
COORDS = (
    "t.fa q.fa\n"
    "NUCMER\n"
    "\n" + "=" * 63 + "\n"
    "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |"
    "  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n" + "=" * 63 + "\n"
    "       1    20000  |       95    20094  |    20000    20000  |    98.00  |"
    "   100000   100000  |    20.00    20.00  | ref\tqry\n"
    "   80000   100000  |   100000    80000  |    20001    20001  |    99.10  |"
    "   100000   100000  |    20.00    20.00  | ref\tqry\n"
)

DELTA = """t.fa q.fa
NUCMER
>ref qry 100000 100000
1 20000 95 20094 400 400 0
0
80000 100000 100000 80000 180 180 0
0
"""


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "mummer-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="synteny", module_order=1)


def test_parse_coords_reads_the_tool_table_and_skips_its_rulers() -> None:
    rows = parse_coords(COORDS)
    assert len(rows) == 2
    forward, reverse = rows
    assert (forward.target_start, forward.target_end) == (1, 20000)
    assert (forward.query_start, forward.query_end) == (95, 20094)
    assert forward.identity == 98.00
    assert forward.reverse is False
    assert reverse.reverse is True
    assert (reverse.query_start, reverse.query_end) == (100000, 80000)


def test_parse_coords_refuses_a_row_it_cannot_read() -> None:
    broken = COORDS.replace("       1    20000", "       x    20000")
    with pytest.raises(ValueError, match="line 7"):
        parse_coords(broken)


def test_run_nucmer_records_the_pipeline_and_the_version(tmp_path: Path) -> None:
    stub = make_mummer_stub(tmp_path / "tools", delta=DELTA, coords=COORDS)
    ctx = _ctx(tmp_path)
    result = run_nucmer(
        ("ref", "ACGT" * 25000),
        ("qry", "ACGT" * 25000),
        ctx,
        min_match=20,
        min_cluster=65,
        filter_mode="-1",
        directory=stub,
    )
    assert result.version == MUMMER_VERSION
    assert len(result.blocks) == 2
    assert len(result.delta_blocks) == 2
    assert result.query_name == "qry"
    assert [record.name for record in ctx.tools] == ["nucmer", "delta-filter", "show-coords"]
    nucmer = ctx.tools[0]
    argv = list(nucmer.argv)
    assert "--maxmatch" in argv
    assert argv[argv.index("-l") + 1] == "20"
    assert argv[argv.index("-c") + 1] == "65"
    assert argv[argv.index("-t") + 1] == "1"  # threads must not move a manuscript number
    filter_record = ctx.tools[1]
    assert filter_record.argv[1] == "-1"
    assert all(record.version == MUMMER_VERSION for record in ctx.tools)


def test_run_nucmer_refuses_a_different_mummer(tmp_path: Path) -> None:
    stub = make_mummer_stub(tmp_path / "tools", delta=DELTA, coords=COORDS, version="3.23")
    with pytest.raises(ToolError, match="3.23"):
        run_nucmer(
            ("ref", "ACGT" * 25000),
            ("qry", "ACGT" * 25000),
            _ctx(tmp_path),
            min_match=20,
            min_cluster=65,
            filter_mode="-1",
            directory=stub,
        )


def test_run_nucmer_refuses_when_the_two_readings_disagree(tmp_path: Path) -> None:
    """A table and a series drawn from different alignments is the one thing this must not do."""
    short = COORDS.replace(
        "   80000   100000  |   100000    80000  |    20001    20001  |    99.10",
        "   80000    90000  |   100000    90000  |    10001    10001  |    99.10",
    )
    stub = make_mummer_stub(tmp_path / "tools", delta=DELTA, coords=short)
    with pytest.raises(ToolError, match="disagree"):
        run_nucmer(
            ("ref", "ACGT" * 25000),
            ("qry", "ACGT" * 25000),
            _ctx(tmp_path),
            min_match=20,
            min_cluster=65,
            filter_mode="-1",
            directory=stub,
        )
