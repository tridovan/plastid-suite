"""Genome structure: skew conventions, region invariants, and parameter guardrails (offline)."""

from __future__ import annotations

from itertools import pairwise

import pytest

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome
from plastid.core.sequence import gc_skew
from plastid.modules import structure
from plastid.run.context import RunContext


def _ctx(genome: Genome) -> RunContext:
    ir = detect_ir(genome)
    enforce_ir_gate(ir)
    # the module reads only params and the IR — no layout is needed for these tests
    return RunContext(layout=None, params={}, module_id="structure", module_order=1, ir=ir)  # type: ignore[arg-type]


def _rows(genome: Genome, params: dict | None = None) -> dict[str, dict]:
    result = structure.run(genome, params or {}, _ctx(genome))
    return {row["region"]: row for row in result.tables["regions"].rows}


def test_gc_skew_sign_convention() -> None:
    # (G − C)/(G + C), stated in the column header — a flip here is the phase's classic error
    assert gc_skew("GGGC") == pytest.approx(0.5)
    assert gc_skew("GCCC") == pytest.approx(-0.5)
    assert gc_skew("GCGC") == pytest.approx(0.0)


def test_gc_skew_is_undefined_without_g_and_c() -> None:
    with pytest.raises(ValueError, match="no G and no C"):
        gc_skew("ATAT")
    with pytest.raises(ValueError, match="empty"):
        gc_skew("")


def test_the_two_ir_copies_have_opposite_skew(golden_genome: Genome) -> None:
    """The IR copies are reverse complements of each other — their skew must mirror exactly."""
    rows = _rows(golden_genome)
    assert rows["IRa"]["gc_skew"] == pytest.approx(-rows["IRb"]["gc_skew"], abs=1e-9)
    assert rows["IRa"]["gc_pct"] == pytest.approx(rows["IRb"]["gc_pct"], abs=0.01)


def test_regions_are_contiguous_and_ordered(golden_genome: Genome) -> None:
    rows = _rows(golden_genome)
    order = ["LSC", "IRa", "SSC", "IRb"]
    for earlier, later in pairwise(order):
        assert rows[later]["start_1b"] == rows[earlier]["end_1b"] + 1
    assert rows["LSC"]["start_1b"] == 1
    assert rows["IRb"]["end_1b"] == len(golden_genome)


def test_cumulative_skew_starts_anew_in_each_region(golden_genome: Genome) -> None:
    table = structure.run(golden_genome, {}, _ctx(golden_genome)).tables["gc_skew_windows"]
    seen: set[str] = set()
    for row in table.rows:
        if row["region"] not in seen:
            seen.add(row["region"])
            assert row["cumulative_gc_skew"] == pytest.approx(row["gc_skew"], abs=1e-9)


def test_window_parameters_reach_the_table_and_the_scope(golden_genome: Genome) -> None:
    result = structure.run(golden_genome, {"window": 2000, "step": 500}, _ctx(golden_genome))
    summary = result.tables["structure_summary"].rows[0]
    assert (summary["window_bp"], summary["step_bp"]) == (2000, 500)
    assert "2000 bp window / 500 bp step" in result.tables["gc_skew_windows"].scope.label


def test_step_larger_than_window_is_rejected(golden_genome: Genome) -> None:
    """A step wider than the window would silently skip sequence between windows."""
    with pytest.raises(ValueError, match="exceeds window"):
        structure.run(golden_genome, {"window": 100, "step": 200}, _ctx(golden_genome))


def test_non_positive_window_is_rejected(golden_genome: Genome) -> None:
    with pytest.raises(ValueError, match="must be positive"):
        structure.run(golden_genome, {"window": 0}, _ctx(golden_genome))
