"""M11.15 on the validation data: "no intron" is not "could not assess".

Two independent sets, both gitignored (``.panel/``), so the check skips with a named reason
without them:

* **80 RefSeq plastomes outside Iridaceae** — 5 fail the IR gate, 31 carry ambiguous bases, 3 are
  linear and 2 are reduced. A gate failure must stop as ``IRGateError``; the module never reports
  a failed record as having zero introns.
* **30 PGA annotation variants** — every one passes the IR gate; two carry a malformed tRNA join
  whose parts overlap by 5 bp. That gap is named "not assessed" in the notes and is not counted
  as an intron, while the rest of the record is still analysed.

Counts are measured 26/10/2026 on the pinned panel and asserted so a rebuilt panel that no longer
exercises these cases fails here rather than passing on a smaller sample. The panel is validation
data, not a replacement for the golden.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, detect_ir, enforce_ir_gate
from plastid.core.model import Genome
from plastid.modules import introns
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"
VARIANTS = ROOT / ".panel" / "ir-validation" / "annotation-variants"

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not PANEL.is_dir() or not VARIANTS.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]

# The panel's documented shape and the measured intron facts (README, 26/10/2026).
PANEL_RECORDS = 80
PANEL_GATE_FAILURES = 5
PANEL_TRNA_INTRONS = 576
PANEL_WITHOUT_TRNA_INTRON = 1
VARIANTS_RECORDS = 30
VARIANTS_TRNA_INTRONS = 213
VARIANTS_UNASSESSED = 2


def _context(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "intron-panel", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="introns", module_order=1)


def _records(directory: Path) -> list[Genome]:
    return [read_genbank(path) for path in sorted(directory.glob("*.gb"))]


@pytest.fixture(scope="module")
def panel() -> dict[str, object]:
    records = _records(PANEL)
    gate_failures: list[Genome] = []
    analysed: list[tuple[Genome, object]] = []
    for genome in records:
        try:
            ir = enforce_ir_gate(detect_ir(genome))
        except (IRGateError, IRDetectionError):
            gate_failures.append(genome)
            continue
        analysed.append((genome, ir))
    return {"records": records, "gate_failures": gate_failures, "analysed": analysed}


def test_the_panel_still_has_the_shapes_this_milestone_is_about(panel) -> None:
    assert len(panel["records"]) == PANEL_RECORDS
    assert len(panel["gate_failures"]) == PANEL_GATE_FAILURES


def test_a_gate_failure_is_named_and_never_reported_as_zero_introns(panel, tmp_path: Path) -> None:
    ctx = _context(tmp_path)

    for genome in panel["gate_failures"]:
        with pytest.raises(IRGateError):
            introns.run(genome, {}, ctx)


def test_every_analysed_panel_record_measures_the_class_it_has(panel, tmp_path: Path) -> None:
    """All 75 gate-passing records annotate CDS and tRNA, so both counts are measurements."""
    ctx = _context(tmp_path)
    trna_introns = 0
    without = 0

    for genome, ir in panel["analysed"]:
        ctx.ir = ir
        result = introns.run(genome, {}, ctx)
        summary = result.tables["intron_summary"].rows[0]
        assert isinstance(summary["n_cds_introns"], int), genome.accession
        assert isinstance(summary["n_trna_introns"], int), genome.accession
        assert summary["assessment"].startswith("CDS and tRNA introns measured")
        trna_introns += summary["n_trna_introns"]
        if summary["n_trna_introns"] == 0:
            without += 1

    assert trna_introns == PANEL_TRNA_INTRONS
    assert without == PANEL_WITHOUT_TRNA_INTRON


def test_the_pga_variants_keep_the_rest_of_the_record_when_one_join_is_malformed(
    tmp_path: Path,
) -> None:
    """Two variants write a tRNA join whose parts overlap; that gap is named, not published."""
    ctx = _context(tmp_path)
    trna_introns = 0
    unassessed = 0

    for genome in _records(VARIANTS):
        ctx.ir = enforce_ir_gate(detect_ir(genome))
        result = introns.run(genome, {}, ctx)
        summary = result.tables["intron_summary"].rows[0]
        if any("overlapping join" in note for note in result.notes):
            unassessed += 1
            assert summary["assessment"].endswith("1 tRNA gap(s) not assessed")
            assert "not assessed" in " ".join(result.notes)
        trna_introns += summary["n_trna_introns"]

    assert len(_records(VARIANTS)) == VARIANTS_RECORDS
    assert trna_introns == VARIANTS_TRNA_INTRONS
    assert unassessed == VARIANTS_UNASSESSED
