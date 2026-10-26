"""Canonical origin and the LSC → IRa → SSC → IRb split (decision 2026-10-02)."""

from __future__ import annotations

from plastid.core.ir import IRResult, detect_ir, enforce_ir_gate, ir_from_annotation
from plastid.core.model import Genome, Interval
from plastid.core.normalize import _shift_interval, normalize, rotate_genome


def _synthetic_ir(a: tuple[int, int], b: tuple[int, int]) -> IRResult:
    ira, irb = Interval(*a), Interval(*b)
    return IRResult(
        ira=ira,
        irb=irb,
        identity=1.0,
        mismatches=0,
        compared_length=ira.length,
        exact_core=ira.length,
        mirror_constant=ira.start + irb.end - 1,
        source="manual",
    )


def test_golden_genome_is_already_canonical(golden_genome: Genome) -> None:
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    frame = normalize(golden_genome, ir)
    assert frame.rotation_offset == 0
    assert [(region.name, region.length) for region in frame.regions] == [
        ("LSC", 83_127),
        ("IRa", 26_214),
        ("SSC", 18_169),
        ("IRb", 26_214),
    ]
    assert sum(region.length for region in frame.regions) == len(frame.genome) == 153_724
    assert frame.regions[0].start == 0
    assert frame.regions[-1].end == len(frame.genome)


def test_origin_picks_the_long_single_copy_gap() -> None:
    # 100 bp circle, copies [10,20) and [50,60): gaps are 30 (between) and 50 (wrapping),
    # so the wrapping gap is the LSC and the origin moves to position 60.
    genome = Genome(name="synthetic", sequence="ACGT" * 25)
    frame = normalize(genome, _synthetic_ir((10, 20), (50, 60)))
    assert frame.rotation_offset == 60
    assert [(region.name, region.length) for region in frame.regions] == [
        ("LSC", 50),
        ("IRa", 10),
        ("SSC", 30),
        ("IRb", 10),
    ]
    assert frame.genome.sequence == genome.sequence[60:] + genome.sequence[:60]


def test_origin_keeps_the_between_gap_as_lsc_when_it_is_longer() -> None:
    # copies [10,20) and [70,80): gaps are 50 (between) and 30 (wrapping) → origin at 20,
    # and the copy that follows the LSC in the new frame is the one at [70,80).
    genome = Genome(name="synthetic", sequence="ACGT" * 25)
    frame = normalize(genome, _synthetic_ir((10, 20), (70, 80)))
    assert frame.rotation_offset == 20
    assert [(region.name, region.length) for region in frame.regions] == [
        ("LSC", 50),
        ("IRa", 10),
        ("SSC", 30),
        ("IRb", 10),
    ]
    assert frame.ir.ira.start == 50  # the copy adjacent to the LSC's 3' end
    assert frame.ir.irb.end == len(genome)  # the other copy ends at the origin


def test_shift_splits_a_feature_that_crosses_the_origin() -> None:
    parts = _shift_interval(Interval(90, 120), -10, 100)
    assert [(part.start, part.end) for part in parts] == [(80, 100), (0, 10)]


def test_normalize_a_genome_that_starts_inside_an_ir(edge_genomes: dict[str, Genome]) -> None:
    # NC_085783.1 has no repeat_region annotation and its origin sits inside an IR copy —
    # the case that forces the detector's half-rotation retry.
    genome = edge_genomes["NC_085783.1"]
    ir = detect_ir(genome)
    enforce_ir_gate(ir)
    frame = normalize(genome, ir)
    assert frame.regions[0].name == "LSC"
    assert frame.regions[0].start == 0
    assert frame.regions[-1].end == len(frame.genome)
    assert sum(region.length for region in frame.regions) == len(genome) == 153_162
    assert frame.ir.length == 21_148


def test_rotate_genome_moves_genes_and_ir_consistently(golden_genome: Genome) -> None:
    rotated = rotate_genome(golden_genome, 12_345)
    assert rotated.sequence == golden_genome.sequence[12_345:] + golden_genome.sequence[:12_345]
    assert len(rotated.genes) == len(golden_genome.genes)
    # total gene length is invariant under rotation (wrap splits add parts, not bases)
    assert sum(gene.length for gene in rotated.genes) == sum(
        gene.length for gene in golden_genome.genes
    )
