"""The ``repeats`` module — the class names, the region of each copy, and what is refused.

The tool is stubbed: what is under test is the module's contract (the two classes as the *tool*
marks them, both copies' regions, longest-first ranking, the floor in the label, an empty result
that is still a result), not MUMmer's algorithm — the real numbers are pinned in
``tests/golden/test_tier_b.py``.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRGateError, IRResult, detect_ir, enforce_ir_gate
from plastid.core.model import Genome, Interval
from plastid.external.repeatmatch import REPEAT_MATCH_VERSION
from plastid.modules import repeats
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

DATA = Path(__file__).parents[1] / "data"
# Two matches on the golden, in the shape the real tool printed for this record: the IR
# pair itself, and the 42 bp match inside IRa.
STUB_OUTPUT = "  83128  153724r  26214\n  88397   88418   42\n"


def _stub(directory: Path, output: str = STUB_OUTPUT) -> Path:
    """A ``repeat-match`` printing canned output, beside the ``mummer`` that names its release."""
    directory.mkdir(parents=True, exist_ok=True)
    for name, body in (
        ("repeat-match", f"printf '%s' '{output}'\n"),
        ("mummer", f'echo "{REPEAT_MATCH_VERSION}"\n'),
    ):
        path = directory / name
        path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return directory


def _run(tmp_path: Path, ir: IRResult | None, tools: Path, genome: Genome | None = None):
    layout = RunLayout.create(tmp_path / "run", "unit", fingerprint="test")
    layout.make_dirs()
    params = {"repeats": str(tools)}
    context = RunContext(layout=layout, params=params, module_id="repeats", module_order=1, ir=ir)
    return repeats.run(genome or read_genbank(DATA / "NC_050833.1.gb"), params, context)


@pytest.fixture(scope="module")
def golden_ir() -> IRResult:
    ir = detect_ir(read_genbank(DATA / "NC_050833.1.gb"))
    enforce_ir_gate(ir)
    return ir


def test_each_copy_carries_its_own_region(tmp_path: Path, golden_ir: IRResult) -> None:
    result = _run(tmp_path, golden_ir, _stub(tmp_path / "bin"))
    rows = result.tables["repeats"].rows
    assert [(row["repeat_id"], row["class"], row["length_bp"]) for row in rows] == [
        ("R01", "palindromic", 26_214),
        ("R02", "dispersed", 42),
    ]
    # the IR pair: copy 1 is IRa, copy 2 is IRb — the regions a manuscript table prints
    assert (rows[0]["region1"], rows[0]["region2"]) == ("IRa", "IRb")
    assert (rows[0]["start1_1b"], rows[0]["end1_1b"]) == (83_128, 109_341)
    assert (rows[0]["start2_1b"], rows[0]["end2_1b"]) == (127_511, 153_724)
    # the 42 bp match sits inside IRa, twice over
    assert (rows[1]["region1"], rows[1]["region2"]) == ("IRa", "IRa")


def test_the_summary_counts_the_two_classes_and_names_the_floor(
    tmp_path: Path, golden_ir: IRResult
) -> None:
    result = _run(tmp_path, golden_ir, _stub(tmp_path / "bin"))
    summary = result.tables["repeat_summary"].rows[0]
    assert (summary["n_matches"], summary["n_dispersed"], summary["n_palindromic"]) == (2, 1, 1)
    assert summary["min_length_bp"] == 30  # the floor travels with the count
    assert summary["tool"] == f"repeat-match (MUMmer {REPEAT_MATCH_VERSION})"
    label = result.tables["repeat_summary"].scope.label
    assert "exact" in label and "30 bp" in label and "tandem" in label
    assert any("inverted-repeat pair itself" in note for note in result.notes)


def test_the_floor_reaches_the_tool_and_the_label(tmp_path: Path, golden_ir: IRResult) -> None:
    """A run asked for 50 bp: the argv carries it, the label prints it, the count obeys it."""
    layout = RunLayout.create(tmp_path / "run", "unit", fingerprint="test")
    layout.make_dirs()
    stub = _stub(tmp_path / "bin", output="  83128  153724r  26214\n")
    params = {"repeats": str(stub), "min_length_bp": 50}
    context = RunContext(
        layout=layout, params=params, module_id="repeats", module_order=1, ir=golden_ir
    )
    result = repeats.run(read_genbank(DATA / "NC_050833.1.gb"), params, context)
    summary = result.tables["repeat_summary"].rows[0]
    assert (summary["min_length_bp"], summary["n_matches"]) == (50, 1)
    assert "50 bp" in result.tables["repeat_summary"].scope.label
    assert "-n" in context.tools[0].argv and "50" in context.tools[0].argv


def test_a_record_with_no_ir_cannot_have_its_regions_labelled(
    tmp_path: Path, golden_ir: IRResult
) -> None:
    """``requires_ir`` is not decoration: without the gate there is no LSC/IRa/SSC/IRb to print."""
    with pytest.raises(IRGateError, match="repeats"):
        _run(tmp_path, None, _stub(tmp_path / "bin"))


def test_finding_nothing_is_a_result_with_a_zero_and_a_note(tmp_path: Path) -> None:
    """A run that found nothing must not look like a run that did not happen."""
    stub = _stub(tmp_path / "empty", output="")
    # the canonical frame of tests/unit/test_spans.py: LSC [0,30) · IRa [30,50) · SSC [50,70) ·
    # IRb [70,90) — the LSC longer than the SSC, which is what makes it canonical
    genome = Genome(name="synthetic", sequence="ACGT" * 22 + "AC")
    ir = IRResult(
        ira=Interval(30, 50),
        irb=Interval(70, 90),
        identity=1.0,
        mismatches=0,
        compared_length=20,
        exact_core=20,
        mirror_constant=119,
        source="manual",
    )
    result = _run(tmp_path, ir, stub, genome=genome)
    assert result.tables["repeat_summary"].rows[0]["n_matches"] == 0
    assert "repeats" not in result.tables  # absent rather than empty (Table refuses zero rows)
    assert any("no exact match" in note for note in result.notes)
