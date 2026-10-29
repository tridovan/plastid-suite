"""The annotated record written from five different origins.

`test_rotation_structure.py` covers the 30 unannotated rotated FASTA with `structure`. The
gene-based modules cannot run on those files — they carry no annotation — so this file rotates the
**annotated** panel records by the same offsets `fasta/EXPECTED.tsv` uses (`rotate_genome`, which
moves the sequence and every feature together) and asks the same question: a rotation renames
positions and changes nothing else, so no number may move.

Three modules read genes but no alignment — `codon_usage`, `amino_acid` and `stats` — and their
invariance is exact: rotate a record and every codon, residue and per-locus measurement is the
same number. `stats` rotates all four genomes at once, each by its own offset, because its output
is a comparison across genomes, not a property of one.

`barcodes`, `diversity` and `substitutions` align whole plastomes. Their reference defines the
coordinates every table reports, so the reference is rotated and must still give the same tables.
These three run on the **golden pair** rather than on two panel plastomes: the panel's six rotated
records are from six different families, and a two-way whole-plastome MAFFT of them takes ~3
minutes per alignment (measured 28/10/2026), while the golden pair aligns in ~7 s. The claim being
checked — that the frame the reference is analyzed in follows its own annotation, not the byte it
starts at — does not depend on how far apart the two genomes are.

The `diversity_igs` table is compared as a keyed set, not as an ordered list: its rows follow the
genomic order, and rotation changes which spacer comes first. That is presentation, not a number.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.normalize import rotate_genome
from plastid.modules import amino_acid, barcodes, codon_usage, diversity, stats, substitutions
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation"
SOURCES = PANEL / "source"
EXPECTED = PANEL / "fasta" / "EXPECTED.tsv"
GOLDEN = ROOT / "tests" / "data"

# The four panel plastomes the FASTA set also rotates; `stats` needs at least two per group, so
# its case uses all four and the single-genome modules use the first three.
PANEL_RECORDS = ("NC_008641.1", "NC_016708.2", "NC_028027.1", "NC_016736.1")

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not EXPECTED.is_file() or not SOURCES.is_dir(),
        reason="the rotated FASTA panel is not in this tree (.panel/ is gitignored)",
    ),
]


def _offsets() -> dict[str, list[int]]:
    """The five rotation offsets of each panel plastome, in the order the panel wrote them."""
    lines = EXPECTED.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    offsets: dict[str, list[int]] = {}
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t"), strict=True))
        offsets.setdefault(row["source"], []).append(int(row["rotation_offset"]))
    return offsets


def _ctx(root: Path, module: str) -> RunContext:
    layout = RunLayout.create(root, module, fingerprint="m11-17")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module, module_order=1)


def _rows(result) -> dict[str, list[tuple[tuple[str, object], ...]]]:
    return {
        table_id: sorted(tuple(sorted(row.items())) for row in table.rows)
        for table_id, table in result.tables.items()
    }


def _rotations(accession: str):
    genome = read_genbank(SOURCES / f"{accession}.gb")
    return [(offset, rotate_genome(genome, offset)) for offset in _offsets()[accession]]


def test_codon_usage_does_not_depend_on_the_origin(tmp_path: Path) -> None:
    """RSCU, ENC, CAI, CBI and Fop are properties of the CDS set, and the CDS set is unchanged."""
    for accession in PANEL_RECORDS[:3]:
        tables = [
            _rows(codon_usage.run(genome, {}, _ctx(tmp_path / f"{accession}-{offset}", "codon")))
            for offset, genome in _rotations(accession)
        ]
        assert all(table == tables[0] for table in tables[1:]), accession
        assert "rscu" in tables[0] and "rscu_families" in tables[0]


def test_amino_acid_composition_does_not_depend_on_the_origin(tmp_path: Path) -> None:
    """GRAVY and aromaticity are sums over residues; a rotation adds none."""
    for accession in PANEL_RECORDS[:3]:
        tables = [
            _rows(amino_acid.run(genome, {}, _ctx(tmp_path / f"{accession}-{offset}", "amino")))
            for offset, genome in _rotations(accession)
        ]
        assert all(table == tables[0] for table in tables[1:]), accession
        assert "amino_acid_composition" in tables[0]


def test_group_statistics_do_not_depend_on_the_origin(tmp_path: Path) -> None:
    """Each genome is measured from its own CDS; moving all four origins moves no measurement."""
    groups = {accession: "A" if index < 2 else "B" for index, accession in enumerate(PANEL_RECORDS)}
    per_record = {accession: _rotations(accession) for accession in PANEL_RECORDS}
    tables = []
    for index in range(len(_offsets()[PANEL_RECORDS[0]])):
        genomes = [per_record[accession][index][1] for accession in PANEL_RECORDS]
        result = stats.run(genomes, {"groups": groups}, _ctx(tmp_path / f"stats-{index}", "stats"))
        tables.append(_rows(result))
    assert all(table == tables[0] for table in tables[1:])
    assert "stats_tests" in tables[0] and "stats_descriptives" in tables[0]


def test_the_alignment_modules_follow_the_reference_annotation_not_its_origin(
    tmp_path: Path,
) -> None:
    """Rotating the reference moves no number in `barcodes`, `diversity` or `substitutions`.

    The reference is the genome every coordinate is reported against, and all three modules now
    align it in its **normalized** frame (the frame its regions and genes are named in). Before
    28/10/2026 `diversity` and `substitutions` aligned the record as deposited: the same two
    genomes then gave 230,663 alignment columns instead of 154,827, 386 substitutions instead of
    1,094, and π `None` on nearly every gene row. The panel case is the measurement that found it;
    this test is the regression.
    """
    reference = read_genbank(GOLDEN / "NC_050833.1.gb")
    query = read_genbank(GOLDEN / "NC_024936.1.gb")
    shifted = rotate_genome(reference, len(reference) // 2)

    for module in (barcodes, diversity, substitutions):
        results = []
        for index, genome in enumerate((reference, shifted)):
            ctx = _ctx(tmp_path / f"{module.MODULE_ID}-{index}", module.MODULE_ID)
            ctx.ir = enforce_ir_gate(detect_ir(genome, source="sequence"))
            results.append(_rows(module.run([genome, query], {}, ctx)))
        assert results[0] == results[1], (
            f"{module.MODULE_ID}: the tables moved when the reference was rotated — "
            f"{sorted(set(results[0]) ^ set(results[1]))}"
        )
