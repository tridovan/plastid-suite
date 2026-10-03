"""GenBank reader — real fixtures plus error handling. Never touches the network."""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Genome, count_genes


def test_reads_a_second_genome(data_dir: Path) -> None:
    genome = read_genbank(data_dir / "NC_024936.1.gb")
    assert genome.accession == "NC_024936.1"
    assert genome.name == "Iris gatesii"
    assert genome.topology == "circular"
    assert len(genome) == 153_441
    counts = count_genes(genome)
    assert counts.total == 132
    # No truncated ycf1 copy here — the diagnostic differs from the golden genome (133).
    assert counts.gene_features == 132


def test_rejects_multi_record_files(tmp_path: Path, data_dir: Path) -> None:
    text = (data_dir / "NC_050833.1.gb").read_text(encoding="utf-8")
    two_records = tmp_path / "two_records.gb"
    two_records.write_text(text + text, encoding="utf-8")
    with pytest.raises(ValueError, match="exactly one record"):
        read_genbank(two_records)


def test_only_gene_features_become_genes(golden_genome: Genome) -> None:
    kinds = {gene.kind for gene in golden_genome.genes}
    # No source / repeat_region / misc_feature entries.
    assert kinds == {"gene", "CDS", "tRNA", "rRNA"}
