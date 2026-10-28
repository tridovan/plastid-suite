"""One inserted or deleted base inside an IR copy can halve the de-novo IR and still pass the gate.

**What this file holds, and why it is here rather than in `test_ir.py`.** The defect is not in the
measured identity — that number is right for the pair the detector chose. It is in *which pair the
detector chooses*: a one-base indel moves the second half of the repeat onto a neighbouring
mirror diagonal, the chain builder starts a fresh run there, and the half that still matches on
the first diagonal becomes a valid-looking 12–13 kb "IR" at 100 %. `_DENOVO_LENGTH_WINDOW` begins
at 10 kb, so the half passes every validity test the detector applies, and `enforce_ir_gate` then
accepts it.

**Why it is not fixed here.** `M11.8` is acceptance, not repair, and the obvious repair was
measured and rejected: allowing a small negative gap in the chain (so the two halves join across
the indel) fixes this record and creates a *new* half-length pass on `NC_057683.1` and changes
`NC_068894.1`'s detected length from 22,652 to 25,121 bp. A correct fix has to decide which of
two overlapping candidate pairs is the repeat, and that is a detector design change with its own
milestone. Until then the test below is `xfail(strict=True)`: it records the defect, it fails if
the defect ever disappears, and it is never deleted.

**The construction is the simplest possible: it uses the record's own annotated copies.** Take a
published plastome, delete one base from the middle of copy B, and ask the de-novo detector — which
never reads the annotation — what it finds. The expected IR is the copies' own length minus one;
the observed value is about half of it. A fully random synthetic genome does **not** reproduce
this (measured: 12/20/26 kb random IRs all lose the pair or return it whole), because the collapse
depends on the repeat's real k-mer spectrum; that is why this test reads the pinned panel and
skips, with a named reason, where `.panel/` is absent.
"""

from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import detect_ir, ir_de_novo, ir_from_annotation

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / ".panel" / "ir-validation" / "source" / "NC_028027.1.gb"


def _require_source() -> Path:
    if SOURCE.is_file():
        return SOURCE
    reason = f"the IR panel record {SOURCE.name} is not in this tree (.panel/ is gitignored)"
    if os.environ.get("PTS_REQUIRE_PANEL") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _deleted_mid_copy_b():
    """Copy B's middle base removed; the IR should then be ``len(copy) - 1`` at ~100 %."""
    genome = read_genbank(_require_source())
    sequence = genome.sequence
    copies = ir_from_annotation(genome)
    if copies is None:
        pytest.skip(f"{genome.accession}: the record annotates no IR pair")
    middle = (copies.irb.start + copies.irb.end) // 2
    deleted = sequence[:middle] + sequence[middle + 1 :]
    bare = replace(genome, genes=(), repeat_regions=(), sequence=deleted)
    return genome, copies, deleted, bare


@pytest.mark.xfail(
    strict=True,
    reason=(
        "K-07: a 1 bp indel splits the repeat across two mirror diagonals; the chain builder "
        "keeps one half and the 10 kb validity floor accepts it as an IR at 100 %"
    ),
)
def test_one_deleted_base_does_not_halve_the_de_novo_inverted_repeat() -> None:
    genome, copies, deleted, bare = _deleted_mid_copy_b()
    found = ir_de_novo(bare)
    assert found is not None, "one deleted base lost the IR pair entirely"
    # Independent expectation: the annotated copies are exact (identity 1.0), so removing one base
    # from one copy makes the pair shorter by one and identical on every remaining column.
    assert found.length == copies.length - 1, (
        f"{genome.accession}: deleted 1 bp from a {copies.length} bp copy and the detector found "
        f"a {found.length} bp IR at {found.identity:.4f} — the half-repeat, not the repeat"
    )
    assert found.identity >= 0.999


def test_the_reproduction_is_the_construction_it_claims_to_be() -> None:
    """The negative control for the xfail above: the *unedited* record must detect its full IR.

    Without this, the xfail could pass because de-novo detection is broken on the record for an
    unrelated reason, and the test would report the wrong cause for the rest of its life.
    """
    genome, copies, deleted, bare = _deleted_mid_copy_b()
    intact = replace(genome, genes=(), repeat_regions=())
    found = ir_de_novo(intact)
    assert found is not None and found.length == copies.length
    assert len(deleted) == len(genome) - 1  # the edit is one base, as the docstring says


def test_the_halved_pair_would_pass_the_gate() -> None:
    """The reason this matters: the half-repeat is not refused, it is *accepted*.

    An identity of ~100 % over a pair inside the detector's own length window satisfies every
    check, so a user running `plastid run` on such a sequence gets tables built on an IR that is
    half its real size — and nothing in the output says so. `detect_ir` is what the CLI calls, so
    the check goes through the same entry point a run uses.
    """
    genome, copies, deleted, bare = _deleted_mid_copy_b()
    found = detect_ir(bare)
    if found.length == copies.length - 1:
        pytest.skip("the detector now finds the full repeat — K-07 is fixed; delete this skip")
    assert found.gate_passed, (
        "the halved pair is refused after all — re-measure before reading further"
    )
    assert found.length < copies.length * 0.75
