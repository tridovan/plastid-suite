"""Tier A golden numbers (PROJECT.md §6) — pure Python, no external tools, no network.

Golden sample: Iris domestica NC_050833.1.
"""

from __future__ import annotations

import pytest

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, coding_sequence, count_genes
from plastid.core.normalize import normalize
from plastid.core.sequence import gc_percent
from plastid.modules import amino_acid, codon_usage, introns, structure
from plastid.run.context import RunContext


def _ctx(genome: Genome, module_id: str) -> RunContext:
    ir = detect_ir(genome)
    enforce_ir_gate(ir)
    # these modules read only params and the IR — no layout is involved
    return RunContext(layout=None, params={}, module_id=module_id, module_order=1, ir=ir)  # type: ignore[arg-type]


def test_length_and_gc(golden_genome: Genome) -> None:
    assert len(golden_genome) == 153_724
    assert gc_percent(golden_genome.sequence) == pytest.approx(37.8536, abs=0.001)


def test_gene_counts_follow_the_fixed_rule(golden_genome: Genome) -> None:
    counts = count_genes(golden_genome)
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)
    assert counts.total == 132  # never report a bare "132" without the breakdown


def test_feature_gene_count_is_a_diagnostic_only(golden_genome: Genome) -> None:
    # PROJECT.md §6: `gene` features = 133 (one truncated ycf1 copy has no CDS);
    # unique names = 112; names repeated in the IR = 21.
    counts = count_genes(golden_genome)
    assert counts.gene_features == 133
    assert counts.unique_names == 112
    assert counts.repeated_names == 21


def test_truncated_ycf1_copy(golden_genome: Genome) -> None:
    # PROJECT.md §9.9: ycf1 has two `gene` features — a truncated pseudo copy at
    # 108,448–109,347 (gene + misc_feature, no CDS) and the real copy at
    # 122,984–128,404, which is the only one carrying a CDS.
    copies = [gene for gene in golden_genome.genes if gene.name == "ycf1"]
    gene_features = [gene for gene in copies if gene.kind == "gene"]
    cds = [gene for gene in copies if gene.kind == "CDS"]

    truncated = [gene for gene in gene_features if gene.pseudo]
    assert len(gene_features) == 2
    assert len(truncated) == 1
    assert truncated[0].span == (108_448, 109_347)  # 1-based inclusive

    assert len(cds) == 1
    assert cds[0].span == (122_984, 128_404)
    assert not any(gene.span[1] < 110_000 for gene in cds)  # nothing CDS-like in the truncated copy


def test_trans_spliced_rps12(golden_genome: Genome) -> None:
    cdss = [gene for gene in golden_genome.genes if gene.name == "rps12" and gene.kind == "CDS"]
    assert len(cdss) == 2  # one per IR copy
    assert all(len(gene.parts) == 3 for gene in cdss)  # trans-spliced exons
    assert any(gene.strand == 0 for gene in cdss)  # the IRb copy mixes strands


def test_golden_genome_identity(golden_genome: Genome) -> None:
    assert golden_genome.accession == "NC_050833.1"
    assert golden_genome.name == "Iris domestica"
    assert golden_genome.topology == "circular"


def test_ir_and_regions(golden_genome: Genome) -> None:
    ir = detect_ir(golden_genome)
    enforce_ir_gate(ir)
    assert ir.length == 26_214
    assert ir.identity == 1.0
    assert ir.source == "annotation"
    assert ir.denovo_delta == 0  # the independent de-novo scan agrees exactly

    frame = normalize(golden_genome, ir)
    assert frame.rotation_offset == 0  # the golden genome already starts at the LSC
    assert [(region.name, region.length) for region in frame.regions] == [
        ("LSC", 83_127),
        ("IRa", 26_214),
        ("SSC", 18_169),
        ("IRb", 26_214),
    ]
    assert sum(region.length for region in frame.regions) == 153_724


def test_coding_sequence_reads_each_part_on_its_own_strand(golden_genome: Genome) -> None:
    """Measured 2026-10-03: only this reading of trans-spliced rps12 yields a canonical ORF.

    Re-sorting the parts gives AAT…TAT; reading them all on the gene's strand gives ATA…ATT.
    """
    cdss = [gene for gene in golden_genome.genes if gene.kind == "CDS" and gene.name == "rps12"]
    assert len(cdss) == 2  # one per IR copy
    for gene in cdss:
        sequence = coding_sequence(golden_genome, gene)
        assert len(sequence) % 3 == 0
        assert sequence.startswith("ATG")
        assert sequence.endswith("TAA")


def test_codon_usage_covers_the_counted_cds_set(golden_genome: Genome) -> None:
    """P2-03's Done-when: the CDS set is the counted one and every family sum checks out."""
    result = codon_usage.run(golden_genome, {}, _ctx(golden_genome, "codon_usage"))

    summary = result.tables["codon_usage"].rows[0]
    assert summary["n_cds"] == count_genes(golden_genome).cds == 86
    assert summary["internal_stops"] == 0
    assert summary["translation_table"].startswith("NCBI table 1")
    assert summary["enc"] == pytest.approx(49.69, abs=0.01)
    assert summary["cai"] == pytest.approx(0.5375, abs=0.0001)
    assert summary["fop"] == pytest.approx(0.5118, abs=0.0001)

    families = result.tables["rscu_families"].rows
    assert len(families) == 20  # every amino acid, no family silently skipped
    assert [row["amino_acid"] for row in families if row["status"] != "ok"] == []


def test_amino_acid_totals_reconcile_with_codon_usage(golden_genome: Genome) -> None:
    """P2-04's Done-when: the same CDS set and the same totals as the P2-03 codon table."""
    codons = codon_usage.run(golden_genome, {}, _ctx(golden_genome, "codon_usage"))
    amino = amino_acid.run(golden_genome, {}, _ctx(golden_genome, "amino_acid"))

    codon_summary = codons.tables["codon_usage"].rows[0]
    summary = amino.tables["amino_acid"].rows[0]
    assert summary["n_cds"] == codon_summary["n_cds"] == 86
    assert summary["n_residues"] == codon_summary["n_codons"]  # every sense codon is one residue
    assert summary["ignored_symbols"] == 86  # the terminal stop codon of each CDS
    assert summary["gravy"] == pytest.approx(-0.1034, abs=0.0001)
    assert summary["aromaticity"] == pytest.approx(0.1118, abs=0.0001)

    composition = amino.tables["amino_acid_composition"].rows
    assert len(composition) == 20
    assert sum(row["count"] for row in composition) == summary["n_residues"]


def test_introns_cross_check_the_annotation_part_boundaries(golden_genome: Genome) -> None:
    """P2-05's Done-when: the intron list reconciles with the CDS part boundaries."""
    result = introns.run(golden_genome, {}, _ctx(golden_genome, "introns"))
    summary = result.tables["intron_summary"].rows[0]

    assert summary["n_introns"] == 17
    assert summary["n_genes_with_introns"] == 12
    assert summary["total_length_bp"] == 12_547
    assert summary["longest_bp"] == 1_055
    # the two rps12 junctions cross a region boundary: trans-splicing, not introns
    assert summary["n_trans_splicing"] == 2
    assert summary["annotation_check"].startswith("agrees")

    rows = result.tables["introns"].rows
    assert len(rows) == 17
    assert len({row["intron_id"] for row in rows}) == 17  # both IR copies named apart
    assert all(0 <= row["phase"] <= 2 for row in rows)
    ndh_a = next(row for row in rows if row["gene"] == "ndhA")
    assert (ndh_a["start_1b"], ndh_a["end_1b"]) == (119_328, 120_382)  # from the CDS parts
    assert ndh_a["phase"] == 1


def test_structure_module_tiles_the_golden_genome(golden_genome: Genome) -> None:
    """P2-02's Done-when: LSC 83,127 · SSC 18,169 · IRa = IRb = 26,214 · exact tiling."""
    result = structure.run(golden_genome, {}, _ctx(golden_genome, "structure"))

    rows = {row["region"]: row for row in result.tables["regions"].rows}
    assert rows["LSC"]["length_bp"] == 83_127
    assert rows["SSC"]["length_bp"] == 18_169
    assert rows["IRa"]["length_bp"] == rows["IRb"]["length_bp"] == 26_214
    assert sum(row["length_bp"] for row in result.tables["regions"].rows) == 153_724
    assert result.tables["structure_summary"].rows[0]["n_windows"] == 157
