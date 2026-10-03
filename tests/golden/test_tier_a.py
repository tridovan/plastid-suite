"""Tier A golden numbers (PROJECT.md §6) — pure Python, no external tools, no network.

Golden sample: Iris domestica NC_050833.1.
"""

from __future__ import annotations

import pytest

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome, count_genes
from plastid.core.normalize import normalize
from plastid.core.sequence import gc_percent


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
