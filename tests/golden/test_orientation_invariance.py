"""A genome written on the other strand describes the same plastome.

The rotation cases move the origin. This file flips the strand: a record is reverse-complemented
together with its `repeat_region` copies, so it is the same circular molecule read the other way.
A whole-plastome alignment does not depend on which strand the file spells, and the normalized
frame is re-derived from the reference's own IR, so the window and region tables — and π — cannot
move. The pair is synthetic so the case costs seconds; the claim is about the frame, not about a
particular record.
"""

from __future__ import annotations

import random
from dataclasses import replace
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Gene, Genome, Interval, RepeatRegion
from plastid.core.normalize import normalize, rotate_genome
from plastid.core.sequence import reverse_complement
from plastid.modules import barcodes, diversity
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

LSC_BP, IR_BP, SSC_BP = 8_000, 6_000, 4_000

pytestmark = pytest.mark.tier_b


def _random_sequence(length: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


def _pair() -> tuple[Genome, Genome]:
    """Reference and query over the same molecule; the query differs at every 101st base."""
    lsc = _random_sequence(LSC_BP, 1)
    ir = _random_sequence(IR_BP, 2)
    ssc = _random_sequence(SSC_BP, 3)
    sequence = lsc + ir + ssc + reverse_complement(ir)
    query = list(sequence)
    for index in range(0, len(query), 101):
        query[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[query[index]]
    copies = (
        RepeatRegion((Interval(LSC_BP, LSC_BP + IR_BP),)),
        RepeatRegion((Interval(LSC_BP + IR_BP + SSC_BP, len(sequence)),)),
    )
    # Two CDS give the module its gene and spacer tables; their sequence is irrelevant here —
    # the claim is about the frame, not about a protein.
    genes = (
        Gene(name="rbcL", kind="CDS", parts=(Interval(100, 400),), transl_table=11),
        Gene(
            name="psbA",
            kind="CDS",
            parts=(Interval(LSC_BP - 400, LSC_BP - 100),),
            transl_table=11,
        ),
    )
    reference = Genome(
        name="Testus reference",
        accession="T.REF",
        sequence=sequence,
        genes=genes,
        repeat_regions=copies,
    )
    return reference, Genome(
        name="Testus query",
        accession="T.QRY",
        sequence="".join(query),
        genes=genes,
        repeat_regions=copies,
    )


def _flip(interval: Interval, length: int) -> Interval:
    return Interval(
        length - interval.end,
        length - interval.start,
        -interval.strand,
        partial_start=interval.partial_end,
        partial_end=interval.partial_start,
    )


def _reverse(genome: Genome) -> Genome:
    """The same molecule written 5'→3' on the other strand, features flipped with it."""
    length = len(genome)
    return replace(
        genome,
        sequence=reverse_complement(genome.sequence),
        genes=tuple(
            replace(gene, parts=tuple(_flip(part, length) for part in reversed(gene.parts)))
            for gene in genome.genes
        ),
        repeat_regions=tuple(
            RepeatRegion(tuple(_flip(part, length) for part in region.parts))
            for region in genome.repeat_regions
        ),
    )


def _rowset(rows) -> list[str]:
    """Rows as an order-free set (a mirrored annotation lists its loci in the other order)."""
    return sorted(repr(dict(row)) for row in rows)


def _ctx(tmp_path: Path, module: str, genome: Genome) -> RunContext:
    layout = RunLayout.create(tmp_path, module, fingerprint="orientation")
    layout.make_dirs()
    context = RunContext(layout=layout, params={}, module_id=module, module_order=1)
    context.ir = enforce_ir_gate(detect_ir(genome))
    return context


def test_the_reverse_complement_is_the_same_molecule() -> None:
    reference, query = _pair()
    flipped = _reverse(reference)
    assert len(flipped) == len(reference)
    assert flipped.sequence == reverse_complement(reference.sequence)
    assert [region.length for region in flipped.repeat_regions] == [
        region.length for region in reference.repeat_regions
    ]
    for genome in (reference, flipped):
        ir = enforce_ir_gate(detect_ir(genome))
        assert ir.identity >= 0.999
    assert flipped.sequence != reference.sequence


def _to_reference_frame(reference: Genome, query: Genome) -> Genome:
    """The query written from the base the reference's normalized frame starts at.

    Only the reference is normalized by the modules; a query is aligned as deposited. A record
    written on the other strand starts at the mirror of its old last base, so without this the
    linear aligner sees a circular permutation and inserts a whole-genome end gap (measured:
    32,000 columns, 16,000 usable). Bringing the query to the same start isolates the strand from
    the origin — the two claims the rotation cases separate.
    """
    offset = normalize(reference, enforce_ir_gate(detect_ir(reference))).rotation_offset
    return rotate_genome(query, offset)


def test_diversity_whole_genome_tables_do_not_depend_on_the_strand(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    require_tools(tools, "mafft")
    reference, query = _pair()
    flipped_reference, flipped_query = _reverse(reference), _reverse(query)
    forward = diversity.run(
        [reference, query], {}, _ctx(tmp_path / "forward", "diversity", reference)
    )
    reversed_run = diversity.run(
        [flipped_reference, _to_reference_frame(flipped_reference, flipped_query)],
        {},
        _ctx(tmp_path / "reversed", "diversity", flipped_reference),
    )

    left = forward.tables["diversity_summary"].rows[0]
    right = reversed_run.tables["diversity_summary"].rows[0]
    for key in ("alignment_length", "usable_columns", "pi", "window_bp", "step_bp"):
        assert left[key] == right[key], key
    # π per region is an aggregate over every usable column, so mirroring the molecule cannot
    # move it.
    assert forward.tables["diversity_regions"].rows == reversed_run.tables["diversity_regions"].rows

    forward_windows = forward.tables["diversity_windows"].rows
    reversed_windows = reversed_run.tables["diversity_windows"].rows
    assert len(forward_windows) == len(reversed_windows)
    assert [row["usable_columns"] for row in forward_windows] == [
        row["usable_columns"] for row in reversed_windows
    ]
    # Per-window π is *not* asserted equal: a mirrored frame puts a different set of bases under
    # the same fixed grid. The values permute (the multiset is equal), which is the grid's
    # property, not the molecule's.

    # Per-gene π is an aggregate over a gene's own columns and cannot move. The intergenic-spacer
    # table is not compared: this synthetic genome has two genes on a circle, so *both* arcs
    # between them are adjacent in annotation order, and which one is the spacer flips with the
    # strand. That is a property of a two-gene circle, not of the frame.
    assert _rowset(forward.tables["diversity_genes"].rows) == _rowset(
        reversed_run.tables["diversity_genes"].rows
    )


def test_barcodes_sequence_tables_do_not_depend_on_the_strand(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    require_tools(tools, "mafft")
    reference, query = _pair()
    flipped_reference, flipped_query = _reverse(reference), _reverse(query)
    forward = barcodes.run(
        [reference, query], {}, _ctx(tmp_path / "forward", "barcodes", reference)
    )
    reversed_run = barcodes.run(
        [flipped_reference, _to_reference_frame(flipped_reference, flipped_query)],
        {},
        _ctx(tmp_path / "reversed", "barcodes", flipped_reference),
    )

    left = forward.tables["barcodes_summary"].rows[0]
    right = reversed_run.tables["barcodes_summary"].rows[0]
    for key in ("n_windows", "threshold_rule", "n_regions", "n_candidates", "n_found"):
        assert left[key] == right[key], key
    for key in ("mean_pi", "threshold_pi"):
        # means over the same multiset of window values; the last digit follows summation order
        assert left[key] == pytest.approx(right[key], rel=1e-12), key

    forward_windows = forward.tables["barcodes_windows"].rows
    reversed_windows = reversed_run.tables["barcodes_windows"].rows
    assert len(forward_windows) == len(reversed_windows)
    assert sorted(row["pi"] for row in forward_windows if row["pi"] is not None) == sorted(
        row["pi"] for row in reversed_windows if row["pi"] is not None
    )
    if "barcodes_regions" in forward.tables:
        assert len(forward.tables["barcodes_regions"].rows) == len(
            reversed_run.tables["barcodes_regions"].rows
        )

    forward_candidates = {
        row["locus"]: (row["found"], row["span_bp"], row["pi"])
        for row in forward.tables["barcodes_candidates"].rows
    }
    reversed_candidates = {
        row["locus"]: (row["found"], row["span_bp"], row["pi"])
        for row in reversed_run.tables["barcodes_candidates"].rows
    }
    assert forward_candidates == reversed_candidates
