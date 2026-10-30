"""The unequal-copy mapping, measured on the panel that was collected for it (M11.11).

`.panel/ir-validation/source` is 80 RefSeq plastomes outside Iridaceae; 32 of them annotate IR
copies that differ in length by 1–12 bp. Those are the cases the old mirror arithmetic could not
express: it paired every base after the indel with the wrong partner, so a gene's second copy was
placed off by the length difference. The mapping now follows the gate's own alignment.

The check is independent of the mapping code in the part that matters: every mapped gene is
compared with the **record's own annotated copy** of the same name in the other repeat, part for
part. A mapping that was off by an indel cannot land there. Five comparisons disagree for reasons
that are properties of the deposited record — a copy annotated in the wrong orientation, or a
second copy annotated truncated — and they are named below rather than hidden by a looser
tolerance.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRGateError, detect_ir, enforce_ir_gate
from plastid.core.model import Gene, Interval
from plastid.core.normalize import rotate_gene
from plastid.modules.plastid_annotate import _ir_copies

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

# The panel's documented shape: 32 unequal pairs, of which three fail the 0.999 gate on their
# own alignment and therefore never reach the annotator.
UNEQUAL_RECORDS = 32
REFUSED = {"NC_047256.1", "NC_057683.1", "NC_084414.1"}

# Mapped gene vs the record's own second copy. These five disagree because the deposited
# annotation does, not because the mapping does; each was read by hand on 26/10/2026:
#   * the two NC_044733 rRNAs and both NC_049129 rRNAs: the record annotates both IR copies on
#     the **same** strand, which an inverted repeat cannot have; the mapped coordinates are exact
#     and the mapping's strand is the inversion's.
#   * NC_044733 23S: the two copies are annotated 2,961 and 2,810 bp — a 151 bp disagreement
#     between the record's own copies.
#   * NC_062492 ycf2: the second copy is annotated 201 bp against 6,843 bp, the truncated
#     duplicate the reference path already refuses to search (`TRUNCATED_DUPLICATE_FRACTION`).
KNOWN_RECORD_DISAGREEMENTS: dict[tuple[str, str], str] = {
    ("NC_044733.1", "16S ribosomal RNA"): "second copy annotated on the same strand",
    ("NC_044733.1", "23S ribosomal RNA"): "the record's copies are 2,961 and 2,810 bp",
    ("NC_049129.1", "rrn4.5"): "second copy annotated on the same strand",
    ("NC_049129.1", "rrn5"): "second copy annotated on the same strand",
    ("NC_062492.1", "ycf2"): "second copy annotated 201 bp (truncated duplicate)",
}

# The largest disagreement left is 10 bp, at an annotation boundary moved by the copy's indel;
# 15 leaves room for one more base on either side without admitting an indel-sized error (the
# panel's largest copy difference is 12 bp).
TOLERANCE = 15

# The measured histogram, re-measured 26/10/2026. Every mapped copy lands inside the other copy;
# 265 of the 272 have a deposited copy of the same name, kind and part count to compare with
# (seven do not). 254 are exact, two are one base out, seven are 4–10 bp out, and two are the
# named record disagreements above (152 and 6,682 bp). The three same-strand disagreements are
# coordinate-exact and are already inside the 254.
EXPECTED_DELTAS = {0: 254, 1: 2, 4: 2, 6: 1, 9: 2, 10: 2, 152: 1, 6682: 1}


def _in(gene: Gene, region: Interval) -> bool:
    return all(region.start <= part.start and part.end <= region.end for part in gene.parts)


def _envelope(gene: Gene) -> tuple[int, int]:
    return (min(part.start for part in gene.parts), max(part.end for part in gene.parts))


@pytest.fixture(scope="module")
def measured() -> dict[str, object]:
    """One pass over the panel: unequal pairs, gate verdicts, and every mapped gene's agreement."""
    unequal = 0
    refused: set[str] = set()
    checked = 0
    outside_target: list[str] = []
    disagreements: dict[tuple[str, str], str] = {}
    deltas: Counter[int] = Counter()
    for path in sorted(PANEL.glob("*.gb")):
        genome = read_genbank(path)
        accession = genome.accession or genome.name
        ir = detect_ir(genome)
        if ir.ira.length == ir.irb.length:
            continue
        unequal += 1
        try:
            enforce_ir_gate(ir)
        except IRGateError:
            refused.add(accession)
            continue
        frame = ir.frame_offset
        first, second = sorted((ir.ira, ir.irb), key=lambda interval: interval.start)
        in_frame = [rotate_gene(gene, frame, len(genome)) for gene in genome.genes]
        sources = [gene for gene in in_frame if gene.kind in ("CDS", "rRNA") and _in(gene, first)]
        published = [
            gene for gene in in_frame if gene.kind in ("CDS", "rRNA") and _in(gene, second)
        ]
        for gene in sources:
            mapped_copies = [copy for copy in _ir_copies([gene], ir) if _in(copy, second)]
            if not mapped_copies:
                continue
            copy = mapped_copies[0]
            checked += 1
            if not all(
                part.strand == -source.strand
                for source, part in zip(gene.parts, copy.parts, strict=True)
            ):
                outside_target.append(f"{accession}:{gene.name}: strand not inverted")
            candidates = [
                candidate
                for candidate in published
                if candidate.name == copy.name
                and candidate.kind == copy.kind
                and len(candidate.parts) == len(copy.parts)
            ]
            if not candidates:
                continue
            mapped_low, mapped_high = _envelope(copy)
            best = min(
                candidates,
                key=lambda candidate: (
                    abs(_envelope(candidate)[0] - mapped_low)
                    + abs(_envelope(candidate)[1] - mapped_high)
                ),
            )
            published_low, published_high = _envelope(best)
            delta = max(abs(published_low - mapped_low), abs(published_high - mapped_high))
            deltas[delta] += 1
            if delta > TOLERANCE or copy.strand != best.strand:
                disagreements[(accession, copy.name)] = (
                    f"{delta} bp / strands {copy.strand}/{best.strand}"
                )
    return {
        "unequal": unequal,
        "refused": refused,
        "checked": checked,
        "outside_target": outside_target,
        "disagreements": disagreements,
        "deltas": deltas,
    }


def test_the_panel_still_holds_the_shape_this_mapping_was_built_for(measured) -> None:
    """The data claim, not a test count: 32 unequal pairs, three refused by their own alignment."""
    assert measured["unequal"] == UNEQUAL_RECORDS
    assert measured["refused"] == REFUSED


def test_every_mapped_gene_lands_inside_the_other_copy_on_the_opposite_strand(measured) -> None:
    assert measured["checked"] >= 250, "the panel stopped exercising the mapping"
    assert measured["outside_target"] == []


def test_the_mapping_agrees_with_the_records_own_second_copy(measured) -> None:
    """The independent answer: each record annotates the other copy, part for part.

    265 of the 272 mapped genes have a deposited copy of the same name, kind and part count:
    254 land exactly, two are one base out, seven are 4–10 bp out, and the two large differences
    are the named record disagreements. Any change to this distribution fails here rather than
    being absorbed by the tolerance.

    Both sides are compared in the IR's frame. Comparing the mapped interval against an
    un-rotated deposited gene on a record whose ``frame_offset`` is 1 reports a false one-base
    error — the panel has two such records with nine genes each, which is why the frame matters.
    """
    assert measured["disagreements"].keys() == KNOWN_RECORD_DISAGREEMENTS.keys(), (
        "the set of records whose own annotation disagrees with the mapping changed: "
        f"{measured['disagreements']}"
    )
    assert dict(measured["deltas"]) == EXPECTED_DELTAS, (
        f"the mapped-gene disagreement distribution changed: {dict(measured['deltas'])}"
    )
    assert max(delta for delta in measured["deltas"] if delta <= TOLERANCE) <= 10
