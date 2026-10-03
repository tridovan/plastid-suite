"""IR detection and the mandatory gate — the highest-risk component (PROJECT.md §9A).

The measured failure modes this file guards (each one stays forever):
  * naive detection was silently wrong on 2/9 published plastomes;
  * identity measured on the exact-match core hides divergence (22,998 bp core at 100 %
    while the true 26,162 bp copies are only 92.62 % identical);
  * an off-by-one pairing constant collapses identity to ~29 % — it must fail, loudly.
"""

from __future__ import annotations

import pytest

from plastid.core.ir import (
    IRGateError,
    detect_ir,
    enforce_ir_gate,
    ir_de_novo,
    ir_from_annotation,
    measure_ir,
    pair_identity,
)
from plastid.core.model import Genome


def test_annotation_ir_on_the_golden_genome(golden_genome: Genome) -> None:
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    assert ir.source == "annotation"
    assert ir.length == 26_214
    assert ir.identity == 1.0
    assert ir.mismatches == 0
    assert ir.mirror_constant == 236_850  # = a1 + b2 - 1, matches PROJECT.md §8.3
    enforce_ir_gate(ir)  # must not raise


def test_off_by_one_pairing_constant_collapses_identity(golden_genome: Genome) -> None:
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    shifted_identity, _, _, _ = pair_identity(
        golden_genome.sequence, ir.ira, ir.mirror_constant + 1
    )
    assert shifted_identity < 0.5  # measured ~29 %; never silently accepted


def test_identity_is_measured_on_full_copies_not_the_exact_core(
    edge_genomes: dict[str, Genome],
) -> None:
    # NC_063920: a perfect 22,998 bp exact core, but the full 26,162 bp copies are 92.62 %
    # identical (PROJECT.md §8.3). The gate must use the full copies.
    ir = ir_from_annotation(edge_genomes["NC_063920.1"])
    assert ir is not None
    assert ir.exact_core == 22_998
    # 1,929 mismatches over the 26,161 compared bases (the shorter copy; the two differ by
    # 1 bp). PROJECT.md §8.3 rounds this to 92.63 %.
    assert ir.compared_length == 26_161
    assert ir.mismatches == 1_929
    assert ir.identity == pytest.approx(1 - 1_929 / 26_161)
    with pytest.raises(IRGateError):
        enforce_ir_gate(ir)


def test_gate_flags_the_second_diverged_genome(edge_genomes: dict[str, Genome]) -> None:
    # Re-measured 2026-10-02 with the structural pairing (a1 + b2 - 1): 98.87 %, not the
    # 53.65 % recorded in PROJECT.md §8.3. Verdict unchanged: below 99.9 % → gate fails.
    ir = ir_from_annotation(edge_genomes["NC_063921.1"])
    assert ir is not None
    assert ir.identity == pytest.approx(0.988695, abs=1e-5)
    with pytest.raises(IRGateError):
        enforce_ir_gate(ir)


def test_short_repeat_regions_are_ignored(edge_genomes: dict[str, Genome]) -> None:
    # NC_063921 annotates five repeat_regions, three of them spurious (1 bp, 393 bp, 2 bp).
    # Structural filtering (length >= 5 kb, exactly two candidates) must pick the real pair.
    genome = edge_genomes["NC_063921.1"]
    assert len(genome.repeat_regions) == 5
    ir = ir_from_annotation(genome)
    assert ir is not None
    assert ir.length == 26_183


def test_copies_of_unequal_length_are_measured_over_the_shorter(
    edge_genomes: dict[str, Genome],
) -> None:
    # NC_063920 annotates 26,162 bp and 26,161 bp copies — a 1 bp difference.
    genome = edge_genomes["NC_063920.1"]
    ir = ir_from_annotation(genome)
    assert ir is not None
    assert ir.compared_length == 26_161
    assert ir.length == 26_161


def test_de_novo_reproduces_the_two_known_failures(edge_genomes: dict[str, Genome]) -> None:
    # PROJECT.md §8.3: the validated detector matches GenBank on 7/9 annotated genomes and
    # fails on exactly these two — one loudly, one silently. Both must stay visible.
    loud = ir_de_novo(edge_genomes["NC_063921.1"])
    assert loud is None  # no candidate passes the validity window — loud failure

    silent = ir_de_novo(edge_genomes["NC_063920.1"])
    assert silent is not None
    assert silent.length == 22_998  # the exact-match core, shorter than the true 26,161


def test_detect_ir_records_the_de_novo_disagreement(edge_genomes: dict[str, Genome]) -> None:
    ir = detect_ir(edge_genomes["NC_063920.1"])
    assert ir.source == "annotation"  # annotation wins ...
    assert ir.length == 26_161
    assert ir.denovo_length == 22_998
    assert ir.denovo_delta == 3_163  # ... but the disagreement is stored, never dropped


def test_de_novo_agrees_with_annotation_on_the_golden(golden_genome: Genome) -> None:
    ir = detect_ir(golden_genome)
    assert ir.denovo_length == 26_214
    assert ir.denovo_delta == 0


def test_de_novo_on_a_genome_without_any_annotation(edge_genomes: dict[str, Genome]) -> None:
    genome = edge_genomes["NC_085783.1"]
    assert genome.repeat_regions == ()
    ir = detect_ir(genome)
    assert ir.source == "denovo_diag"
    assert ir.length == 21_148
    assert ir.identity == 1.0
    enforce_ir_gate(ir)


def test_measure_ir_rejects_intervals_that_are_not_a_pair() -> None:
    from plastid.core.model import Interval

    with pytest.raises(ValueError, match="disjoint"):
        measure_ir("ACGT" * 10, Interval(0, 10), Interval(5, 15), source="manual")
