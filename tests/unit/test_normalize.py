"""Canonical origin and the LSC → IRa → SSC → IRb split (decision 2026-10-02)."""

from __future__ import annotations

import pytest

from plastid.core.genbank import read_genbank, to_genbank
from plastid.core.ir import IRResult, detect_ir, enforce_ir_gate, ir_from_annotation
from plastid.core.model import Gene, Genome, Interval, coding_sequence, gene_qualifiers
from plastid.core.normalize import _shift_interval, normalize, rotate_genome
from plastid.modules.structure import run as structure_run
from plastid.run.context import RunContext


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
    # The old detector returned a 21,148 bp exact core here; the chained de-novo detector
    # (M11.4) finds the full 26,184 bp repeat. Independent check on the rotated sequence: the
    # two copies are exact reverse complements over their whole length, and extending either by
    # one base breaks the equality — so 26,184 is the extent, not an over-extension.
    assert frame.ir.length == 26_184


def test_rotate_genome_moves_genes_and_ir_consistently(golden_genome: Genome) -> None:
    rotated = rotate_genome(golden_genome, 12_345)
    assert rotated.sequence == golden_genome.sequence[12_345:] + golden_genome.sequence[:12_345]
    assert len(rotated.genes) == len(golden_genome.genes)
    # total gene length is invariant under rotation (wrap splits add parts, not bases)
    assert sum(gene.length for gene in rotated.genes) == sum(
        gene.length for gene in golden_genome.genes
    )


@pytest.mark.parametrize(
    ("strand", "expected"),
    [
        # (start, end, partial_start, partial_end), in transcript order
        (1, [(94, 100, True, False), (0, 4, False, True)]),
        (-1, [(0, 4, False, True), (94, 100, True, False)]),
    ],
)
def test_rotation_keeps_partial_flags_on_the_endpoints_they_describe(
    strand: int, expected: list[tuple[int, int, bool, bool]]
) -> None:
    """K-05: a partial flag is about a coordinate, not about the new record's low/high end.

    The interval [95, 105) on a 100 bp circle is already a wrap in genomic terms; rotated by one
    base it splits at the new origin, and both flags stay on the same physical endpoints. A
    minus-strand feature is returned in transcript order (low fragment first), because
    ``coding_sequence`` reads the parts in the order it is given.
    """
    genome = Genome(
        "partial",
        "ACGT" * 25,
        genes=(
            Gene(
                "g",
                "CDS",
                (Interval(95, 105, strand=strand, partial_start=True, partial_end=True),),
            ),
        ),
    )
    rotated = rotate_genome(genome, 1)
    parts = rotated.genes[0].parts
    assert [
        (part.start, part.end, part.partial_start, part.partial_end) for part in parts
    ] == expected


def test_a_split_introduced_by_rotation_is_not_partial() -> None:
    """A coordinate this toolkit invented must not travel as if the record had marked it."""
    genome = Genome(
        "plain",
        "ACGT" * 25,
        genes=(Gene("g", "CDS", (Interval(90, 120),)),),
    )
    rotated = rotate_genome(genome, 13)
    parts = rotated.genes[0].parts
    assert [(part.start, part.end) for part in parts] == [(77, 100), (0, 7)]
    assert not any(part.partial_start or part.partial_end for part in parts)


@pytest.mark.parametrize("shift", [1, 7, 12_345, 83_127, 140_000])
def test_rotation_preserves_every_coding_sequence(golden_genome: Genome, shift: int) -> None:
    """Only the origin moves: a CDS's bases and their transcript order do not.

    This includes the golden's trans-spliced, mixed-strand rps12 copies, so a fixed
    high-then-low fragment order or a re-sort of parts would change a protein and fail here.
    """
    rotated = rotate_genome(golden_genome, shift)
    for original, moved in zip(golden_genome.genes, rotated.genes, strict=True):
        assert coding_sequence(golden_genome, original) == coding_sequence(rotated, moved), (
            f"{original.name} changed coding sequence after a {shift} bp rotation"
        )
        assert original.length == moved.length


def test_rotation_preserves_a_multipart_trans_spliced_gene_across_the_origin() -> None:
    """Mixed strands and transcript order survive a rotation that splits one of the parts."""
    genome = Genome(
        "trans-spliced",
        "ACGT" * 25,
        genes=(
            Gene(
                "rps12",
                "CDS",
                (Interval(90, 100, strand=-1), Interval(0, 10, strand=1)),
                trans_splicing=True,
            ),
        ),
    )
    before = coding_sequence(genome, genome.genes[0])
    for shift in (1, 5, 13, 95):
        rotated = rotate_genome(genome, shift)
        after = coding_sequence(rotated, rotated.genes[0])
        assert after == before, f"transcript order changed at shift {shift}"
        assert rotated.genes[0].trans_splicing is True
        assert sum(part.length for part in rotated.genes[0].parts) == 20


def test_read_rotate_write_read_keeps_the_annotation_model(tmp_path) -> None:
    """K-05's acceptance path: file → rotate → file → read, not just intervals in memory.

    The feature carries both partial markers, a qualifier of every kind the model writes, and a
    mixed-strand join; the repeat region carries a partial marker too. Total bases per feature
    must not move, and the coding sequences must be identical.
    """
    genome = Genome(
        name="synthetic_circle",
        accession="SYN0001",
        sequence="ACGT" * 25,
        genes=(
            Gene(
                "partial_plus",
                "CDS",
                (
                    Interval(90, 100, partial_start=True),
                    Interval(0, 10, partial_end=True),
                ),
                product="partial protein",
                transl_table=11,
                codon_start=1,
                exception="RNA editing",
            ),
            Gene(
                "partial_minus",
                "tRNA",
                (
                    Interval(0, 5, strand=-1, partial_end=True),
                    Interval(80, 100, strand=-1, partial_start=True),
                ),
                product="tRNA-Gly",
            ),
            Gene(
                "trans",
                "CDS",
                (Interval(20, 30, strand=-1), Interval(60, 70, strand=1)),
                product="trans-spliced protein",
                trans_splicing=True,
                transl_table=1,
                translation="MKV" * 3,
            ),
            Gene("pseudo", "CDS", (Interval(40, 50, strand=-1),), pseudo=True),
        ),
        repeat_regions=(Interval(0, 10, partial_end=True),),
        source={"organism": "Synthetic testus", "mol_type": "genomic DNA"},
    )
    first = read_genbank(_write(tmp_path / "first.gb", genome))
    rotated = rotate_genome(first, 7)
    second = read_genbank(_write(tmp_path / "second.gb", rotated))

    assert second.accession == first.accession
    assert second.topology == first.topology
    assert second.source == first.source
    assert second.sequence == rotated.sequence
    assert len(second.genes) == len(first.genes)
    for before, after in zip(rotated.genes, second.genes, strict=True):
        assert after.name == before.name
        assert after.kind == before.kind
        assert after.product == before.product
        assert after.pseudo is before.pseudo
        assert after.trans_splicing is before.trans_splicing
        assert after.transl_table == before.transl_table
        assert after.codon_start == before.codon_start
        assert after.exception == before.exception
        assert after.translation == before.translation
        assert gene_qualifiers(after) == gene_qualifiers(before)
        assert after.length == before.length
        assert coding_sequence(second, after) == coding_sequence(rotated, before)
    assert second.repeat_regions == rotated.repeat_regions


def _write(path, genome: Genome):
    path.write_text(to_genbank(genome, translations=True), encoding="utf-8")
    return path


def test_a_module_using_the_normalized_frame_is_rotation_invariant(
    golden_genome: Genome,
) -> None:
    """The module path, not just the helper: ``structure`` reads the normalized frame.

    The same record at two origins must produce the same four regions with the same lengths, GC
    and skew. Only ``rotation_offset`` — the map back to GenBank — may differ.
    """

    def regions(genome: Genome) -> dict[str, dict]:
        ir = ir_from_annotation(genome)
        assert ir is not None
        ctx = RunContext(
            layout=None,  # type: ignore[arg-type]
            params={},
            module_id="structure",
            module_order=1,
            ir=ir,
        )
        params = {"window": 100_000, "step": 100_000}
        result = structure_run(genome, params, ctx)
        return {row["region"]: row for row in result.tables["regions"].rows}

    original = regions(golden_genome)
    for shift in (1, 7, 12_345, 83_127):
        rotated = rotate_genome(golden_genome, shift)
        assert regions(rotated) == original, f"normalized regions moved at shift {shift}"
