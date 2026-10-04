"""The synteny module with stubbed MUMmer binaries — wiring, labels, and the loud failures.

The alignments are canned; the *sequences* are real to the module, so the % identity series is
computed from them exactly as in a real run (the stub's delta says which positions are paired,
the module compares the bases). That is the split worth testing here: the tool decides the
alignment, the toolkit decides what the numbers mean.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import make_mummer_stub
from plastid.core.ir import IRResult
from plastid.core.model import Genome, Interval
from plastid.core.sequence import reverse_complement
from plastid.modules.synteny import run
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

LENGTH = 10_500  # LSC 6000 · IRa 2000 · SSC 500 · IRb 2000, the IR copies mirrored


def _reference() -> Genome:
    lsc = "ACGT" * 1_500  # 6,000
    ira = "ACGTTGCA" * 250  # 2,000
    ssc = "TTGACCAG" * 62 + "TTGA"  # 500
    irb = reverse_complement(ira)
    sequence = lsc + ira + ssc + irb
    assert len(sequence) == LENGTH
    return Genome(name="synthetic reference", accession="REF.1", sequence=sequence)


def _ir() -> IRResult:
    return IRResult(
        ira=Interval(6_000, 8_000),
        irb=Interval(8_500, 10_500),
        identity=1.0,
        mismatches=0,
        compared_length=2_000,
        exact_core=2_000,
        mirror_constant=6_000 + 10_500 - 1,
        source="annotation",
    )


def _query(reference: Genome) -> Genome:
    """The reference with three substitutions — every one of them a real mismatch to count."""
    sequence = list(reference.sequence)
    for position in (100, 5_000, 9_000):
        sequence[position] = "A" if sequence[position] != "A" else "C"
    return Genome(name="synthetic query", accession="QRY.1", sequence="".join(sequence))


DELTA = """t.fa q.fa
NUCMER
>REF.1 QRY.1 10500 10500
1 10500 1 10500 3 3 0
0
"""

COORDS = (
    "t.fa q.fa\n"
    "NUCMER\n"
    "\n" + "=" * 63 + "\n"
    "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |"
    "  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n" + "=" * 63 + "\n"
    "       1    10500  |        1    10500  |    10500    10500  |    99.97  |"
    "    10500    10500  |   100.00   100.00  | REF.1\tQRY.1\n"
)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "synteny-test", fingerprint="test")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="synteny", module_order=1)
    ctx.ir = _ir()
    return ctx


def _run_with(tmp_path: Path, *, delta: str = DELTA, coords: str = COORDS):
    reference = _reference()
    query = _query(reference)
    stub = make_mummer_stub(tmp_path / "tools", delta=delta, coords=coords)
    result = run((reference, query), {"mummer": str(stub)}, _ctx(tmp_path))
    return result


def test_the_series_is_computed_from_the_sequences_not_taken_from_the_tool(tmp_path: Path) -> None:
    result = _run_with(tmp_path)
    genomes = result.tables["synteny_genomes"]
    row = genomes.rows[0]
    assert row["genome"] == "QRY.1"
    assert row["n_blocks"] == 1
    assert row["aligned_bp"] == LENGTH  # every base paired
    assert row["coverage_pct"] == pytest.approx(100.0)
    assert row["min_identity"] == row["max_identity"] == row["mean_identity"] == 99.97
    # … while the identity over base pairs comes from the three substitutions planted above
    assert row["aligned_identity"] == pytest.approx(100 * (1 - 3 / LENGTH))


def test_the_series_has_one_value_per_window_and_no_zero_for_a_hole(tmp_path: Path) -> None:
    result = _run_with(tmp_path)
    series = result.tables["synteny_identity"]
    first = series.rows[0]
    assert (first["window_start"], first["window_end"]) == (1, 600)
    assert first["compared_bp"] == 600
    assert first["identity"] == pytest.approx(100 * (1 - 1 / 600))  # the substitution at 100
    assert series.rows[-1]["window_end"] == LENGTH


def test_the_summary_states_the_window_and_the_reference(tmp_path: Path) -> None:
    result = _run_with(tmp_path)
    summary = result.tables["synteny_summary"]
    assert summary.rows[0]["accession"] == "REF.1"
    assert summary.rows[0]["reference_length"] == LENGTH
    assert summary.rows[0]["n_genomes"] == 1
    assert summary.rows[0]["window_bp"] == 600
    assert summary.rows[0]["step_bp"] == 200


def test_the_raw_tool_outputs_travel_with_the_result(tmp_path: Path) -> None:
    result = _run_with(tmp_path)
    assert set(result.artifacts) == {"mummer_QRY.1.coords", "mummer_QRY.1.delta"}
    assert result.artifacts["mummer_QRY.1.delta"].startswith("t.fa q.fa\nNUCMER")


def test_every_table_carries_a_scope_label(tmp_path: Path) -> None:
    result = _run_with(tmp_path)
    for table in result.tables.values():
        assert table.scope.label
        assert table.scope.coord_system == "normalized"


def test_an_empty_comparison_fails_loudly_rather_than_reporting_nothing(tmp_path: Path) -> None:
    empty_coords = COORDS.split("=" * 63 + "\n")[0]
    with pytest.raises(ValueError, match="no alignment block"):
        _run_with(tmp_path, delta="t.fa q.fa\nNUCMER\n", coords=empty_coords)
