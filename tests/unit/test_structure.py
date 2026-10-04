"""Genome structure: skew conventions, region invariants, and parameter guardrails (offline)."""

from __future__ import annotations

from itertools import pairwise

import pytest

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, Interval
from plastid.core.normalize import Region
from plastid.core.sequence import gc_skew, reverse_complement
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


def _genome_with_an_at_tail(ir_length: int = 6_000) -> Genome:
    """A minimal plastome whose LSC ends in a window with no G and no C.

    The windowing walks the region in ``step`` increments, so the last window is
    ``region.length mod step`` bases long — here three, and they are ``TTA``.
    """
    lsc = "GC" * 500 + "TTA"  # 1,003 bp: the trailing 3 bp hold no G and no C
    ir = "ACGT" * (ir_length // 4)
    ssc = "AACCGGTT" * 62 + "AC"  # 498 bp
    sequence = lsc + ir + ssc + reverse_complement(ir)
    return Genome(
        name="synthetic",
        accession="SYNTH_TAIL",
        sequence=sequence,
        topology="circular",
        repeat_regions=(
            Interval(start=len(lsc), end=len(lsc) + len(ir)),
            Interval(
                start=len(lsc) + len(ir) + len(ssc),
                end=len(lsc) + len(ir) + len(ssc) + len(ir),
            ),
        ),
    )


def test_a_window_with_no_g_and_no_c_is_blank_rather_than_fatal() -> None:
    """Measured 04/10/2026 on the corpus (P4a-03): *Chasmanthe floribunda* NC_072562.1 has an
    81,003 bp LSC, so its 1,000 bp windowing ends in a **3 bp** window — ``TTA``. That window has
    no G and no C, ``gc_skew`` is undefined for it, and the module raised: the whole run for that
    genome died, 1 of the corpus's 34.

    An undefined window is a value to leave blank, not a reason to abandon the other 33,000. The
    trailing window is *always* this shape when a region's length is not a multiple of the step,
    so which genomes survive was luck — the golden genome's LSC ends in a 127 bp window that
    happens to hold G and C.
    """
    region = Region(name="LSC", start=0, end=1_003)
    rows = structure._window_rows(region, "GC" * 500 + "TTA", window=1_000, step=1_000)

    assert [row["window_length_bp"] for row in rows] == [1_000, 3]
    assert rows[0]["gc_skew"] == 0.0
    assert rows[1]["gc_skew"] is None, "undefined is not 0.0 and not an error"
    assert rows[1]["cumulative_gc_skew"] is None, (
        "a running sum with an undefined term is undefined"
    )


def test_a_region_with_no_g_and_no_c_still_stops_the_run() -> None:
    """The line between the two: a *region* without G and C is a broken record, not a short
    window — that one keeps raising, because a whole region of pure A/T is not a rounding
    artifact of the windowing but a claim about the data."""
    with pytest.raises(ValueError, match="no G and no C"):
        gc_skew("AT" * 5_000)


def test_the_whole_genome_runs_when_its_tail_window_has_no_gc() -> None:
    """The module end-to-end on such a genome: a result, not an exception — and the summary
    counts the undefined window instead of leaving a blank nobody can explain."""
    genome = _genome_with_an_at_tail()
    result = structure.run(genome, {}, _ctx(genome))

    windows = result.tables["gc_skew_windows"].rows
    lsc = [row for row in windows if row["region"] == "LSC"]
    assert len(lsc) == 2
    assert lsc[-1]["gc_skew"] is None

    summary = result.tables["structure_summary"].rows[0]
    assert summary["n_undefined_windows"] == 1
    assert any("no G and no C" in note for note in result.notes)
