"""One inserted or deleted base inside an IR copy must not halve the de-novo inverted repeat.

**What this file holds, and why it is here rather than in `test_ir.py`.** The defect was never in
the measured identity — that number was right for the pair the detector chose. It was in *which
pair the detector chooses*: a one-base indel moved the second half of the repeat onto a
neighbouring mirror diagonal, the chain builder started a fresh run there, and the half that still
matched on the first diagonal became a valid-looking 12–13 kb "IR" at 100 %. `_DENOVO_LENGTH_WINDOW`
begins at 10 kb, so the half passed every validity test the detector applied, and `enforce_ir_gate`
then accepted it.

**Fixed in the chain rule (K-07, 28/10/2026).** Adjacent exact runs overlap at an indel, because
both are seeded by the same 20-mers around it: a one-base indel leaves a net overlap of 19 bases
between the run before it and the run after it, on both mirror diagonals. The chain now allows an
overlap as long as neither run lies wholly inside the other, and scores a chain by the bases it
covers rather than by the sum of its run lengths — the same rule `_anchor_chain` already applied
on the alignment path. Measured: 132 of 132 constructions (22 positions × {delete, insert} × three
records) now return the full repeat, where the non-overlap rule returned 66.

**The construction is the simplest possible: it uses the record's own annotated copies.** Take a
published plastome, delete one base from the middle of copy B, and ask the de-novo detector — which
never reads the annotation — what it finds. The expected IR is the copies' own length minus one.
A fully random synthetic genome does **not** reproduce this (measured: 12/20/26 kb random IRs all
lose the pair or return it whole), because the collapse depends on the repeat's real k-mer
spectrum; that is why this test reads the pinned panel and skips, with a named reason, where
`.panel/` is absent.
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

pytestmark = pytest.mark.panel


def _require_source() -> Path:
    if SOURCE.is_file():
        return SOURCE
    reason = f"the IR panel record {SOURCE.name} is not in this tree (.panel/ is gitignored)"
    if os.environ.get("PTS_REQUIRE_PANEL") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _edited_mid_copy_b(kind: str):
    """Copy B's middle base removed (``delete``) or one base inserted (``insert``).

    The expectation is arithmetic: an exact pair of length ``n`` stays ``n`` after an insertion
    (one gap column) and becomes ``n - 1`` after a deletion. ``sequence`` is the edited molecule.
    """
    if kind not in ("delete", "insert"):
        raise ValueError(f"kind must be 'delete' or 'insert', got {kind!r}")
    genome = read_genbank(_require_source())
    sequence = genome.sequence
    copies = ir_from_annotation(genome)
    if copies is None:
        pytest.skip(f"{genome.accession}: the record annotates no IR pair")
    middle = (copies.irb.start + copies.irb.end) // 2
    if kind == "delete":
        edited = sequence[:middle] + sequence[middle + 1 :]
    else:
        edited = sequence[:middle] + "A" + sequence[middle:]
    bare = replace(genome, genes=(), repeat_regions=(), sequence=edited)
    return genome, copies, edited, bare


@pytest.mark.parametrize("kind", ["delete", "insert"])
def test_a_single_base_indel_does_not_halve_the_de_novo_inverted_repeat(kind: str) -> None:
    genome, copies, edited, bare = _edited_mid_copy_b(kind)
    found = ir_de_novo(bare)
    assert found is not None, f"one {kind}d base lost the IR pair entirely"
    # Independent expectation: the annotated copies are exact (identity 1.0). A deleted base makes
    # the pair one base shorter and identical on every remaining column; an inserted base leaves
    # the pair 25,195 columns with one gap column, so its length is unchanged.
    expected = copies.length - (1 if kind == "delete" else 0)
    assert found.length == expected, (
        f"{genome.accession}: {kind}d 1 bp from a {copies.length} bp copy and the detector found "
        f"a {found.length} bp IR — the half-repeat, not the repeat"
    )
    assert found.identity >= 0.999


def test_the_reproduction_is_the_construction_it_claims_to_be() -> None:
    """The negative control for the test above: the *unedited* record must detect its full IR.

    Without this, the regression could pass because de-novo detection is broken on the record for
    an unrelated reason, and the test would report the wrong cause for the rest of its life.
    """
    genome = read_genbank(_require_source())
    copies = ir_from_annotation(genome)
    if copies is None:
        pytest.skip(f"{genome.accession}: the record annotates no IR pair")
    intact = replace(genome, genes=(), repeat_regions=())
    found = ir_de_novo(intact)
    assert found is not None and found.length == copies.length


def test_the_halved_pair_would_still_pass_the_gate() -> None:
    """Why the collapse mattered: the half-repeat is not refused, it is *accepted*.

    An identity of ~100 % over a pair inside the detector's own length window satisfies every
    check, so nothing downstream could have noticed. The check therefore goes through `detect_ir`,
    the entry point a run uses, and asserts the full pair — not merely that the detector changed.
    """
    genome, copies, edited, bare = _edited_mid_copy_b("delete")
    found = detect_ir(bare)
    assert found.length == copies.length - 1
    assert found.gate_passed
    assert found.length > copies.length * 0.9
