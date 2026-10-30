"""M11.5 on the validation panel: the data a module actually needs decides whether it runs.

The panel is 80 RefSeq plastomes outside Iridaceae. It contains the four shapes that separate a
sequence-based analysis from a position-based one: **5** records whose IR fails the 0.999 gate,
**31** with ambiguous bases, **3** deposited as linear and **2** reduced plastomes. Codon usage
and amino-acid composition read CDS sequences and must run on all of them; a position-based
module must still stop at the gate.

The panel is gitignored; without it the check skips with a named reason. It is validation data,
not a replacement for the golden and not a reason to change any gate threshold.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, detect_ir, enforce_ir_gate
from plastid.core.model import Genome, cds_features
from plastid.modules import amino_acid, codon_usage, genome_stats
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.panel,
    pytest.mark.skipif(
        not PANEL.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]

# The panel's documented shape (README, 26/10/2026). Asserted, so a rebuilt panel that no longer
# exercises these shapes fails here rather than passing on a smaller sample.
GATE_FAILURES = 5
AMBIGUOUS = 31
LINEAR = 3
REDUCED = 2


def _context(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "module-needs", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="test", module_order=1, ir=None)


@pytest.fixture(scope="module")
def panel() -> dict[str, list[Genome]]:
    """The panel split by the shapes this milestone is about."""
    gate_failures: list[Genome] = []
    ambiguous: list[Genome] = []
    linear: list[Genome] = []
    reduced: list[Genome] = []
    for path in sorted(PANEL.glob("*.gb")):
        genome = read_genbank(path)
        if set(genome.sequence.upper()) - set("ACGT"):
            ambiguous.append(genome)
        if genome.topology == "linear":
            linear.append(genome)
        if len(genome) < 100_000:
            reduced.append(genome)
        try:
            enforce_ir_gate(detect_ir(genome))
        except (IRGateError, IRDetectionError):
            gate_failures.append(genome)
    return {
        "gate_failures": gate_failures,
        "ambiguous": ambiguous,
        "linear": linear,
        "reduced": reduced,
    }


def test_the_panel_still_has_the_shapes_this_milestone_is_about(panel) -> None:
    assert len(panel["gate_failures"]) == GATE_FAILURES
    assert len(panel["ambiguous"]) == AMBIGUOUS
    assert len(panel["linear"]) == LINEAR
    assert len(panel["reduced"]) == REDUCED


@pytest.mark.parametrize("shape", ["gate_failures", "ambiguous", "linear", "reduced"])
def test_sequence_based_modules_run_on_every_shape(
    panel: dict[str, list[Genome]], shape: str, tmp_path: Path
) -> None:
    """No IR result is passed: that is the run these modules have when nothing else needs one."""
    ctx = _context(tmp_path)
    for genome in panel[shape]:
        cdss = cds_features(genome)
        assert cdss, f"{genome.accession}: the fixture has no CDS to count"
        codon = codon_usage.run(genome, {}, ctx)
        amino = amino_acid.run(genome, {}, ctx)
        assert codon.tables["codon_usage"].rows[0]["n_cds"] == len(cdss)
        assert amino.tables["amino_acid"].rows[0]["n_cds"] == len(cdss)


@pytest.mark.parametrize("shape", ["gate_failures", "ambiguous", "linear", "reduced"])
def test_a_position_based_module_is_still_blocked_on_the_same_records(
    panel: dict[str, list[Genome]], shape: str, tmp_path: Path
) -> None:
    ctx = _context(tmp_path)
    for genome in panel[shape]:
        with pytest.raises(IRGateError):
            genome_stats.run(genome, {}, ctx)
