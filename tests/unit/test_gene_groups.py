"""Gene grouping: the loader's guards, prefix claiming, and the nothing-is-dropped rule."""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.gene_groups import (
    GeneGroupError,
    assign,
    load_gene_groups,
)
from plastid.core.model import Gene as Feature
from plastid.core.model import Genome, Interval, count_genes


def _write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "groups.yaml"
    path.write_text(text, encoding="utf-8")
    return path


def _genome(*features: Feature) -> Genome:
    return Genome(name="synthetic", sequence="ATGC" * 300, genes=features)


def test_shipped_file_loads_and_is_not_empty() -> None:
    groups = load_gene_groups()
    assert [group.category for group in groups.groups][:2] == ["Photosynthesis", "Photosynthesis"]
    assert {group.name for group in groups.groups} >= {
        "Subunits of photosystem I",
        "Transfer RNA genes",
        "Ribosomal RNA genes",
    }


def test_a_prefix_group_claims_a_name_no_list_mentions(tmp_path: Path) -> None:
    """The whole point of matching tRNA by prefix: another species brings another anticodon."""
    path = _write(
        tmp_path,
        "categories:\n  - name: Self-replication\n    groups:\n"
        "      - name: Transfer RNA genes\n        prefix: trn\n",
    )
    groups = load_gene_groups(path)
    assert groups.claim("trnL-TAA") is not None
    assert groups.claim("trnA-UGC") is not None
    assert groups.claim("rbcL") is None


def test_a_gene_no_group_claims_comes_back_unclassified(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "categories:\n  - name: Photosynthesis\n    groups:\n"
        "      - name: Large subunit of rubisco\n        genes: [rbcL]\n",
    )
    groups = load_gene_groups(path)
    genome = _genome(
        Feature("rbcL", "CDS", (Interval(0, 100),)),
        Feature("mystery1", "CDS", (Interval(200, 300),)),
    )
    rows = assign(genome, groups)
    assert [(row.gene.name, row.group and row.group.name) for row in rows] == [
        ("rbcL", "Large subunit of rubisco"),
        ("mystery1", None),  # reported, not dropped and not folded into a neighbour
    ]


def test_a_name_shadowed_by_an_earlier_prefix_rule_is_rejected(tmp_path: Path) -> None:
    """Silent by construction: the gene still lands somewhere, just not where the file says."""
    path = _write(
        tmp_path,
        "categories:\n  - name: One\n    groups:\n"
        "      - name: Transfer RNA genes\n        prefix: trn\n"
        "      - name: Maturase\n        genes: [trnK-UUU]\n",
    )
    with pytest.raises(GeneGroupError, match="already claims it"):
        load_gene_groups(path)


def test_a_group_must_claim_by_exactly_one_rule(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "categories:\n  - name: One\n    groups:\n"
        "      - name: Both\n        genes: [rbcL]\n        prefix: rbc\n",
    )
    with pytest.raises(GeneGroupError, match="exactly one"):
        load_gene_groups(path)


def test_a_typo_in_a_key_is_reported_not_ignored(tmp_path: Path) -> None:
    path = _write(
        tmp_path,
        "categories:\n  - name: One\n    groups:\n      - name: Rubisco\n        gene: [rbcL]\n",
    )
    with pytest.raises(GeneGroupError, match="unknown key"):
        load_gene_groups(path)


def test_only_counted_kinds_are_assigned(tmp_path: Path) -> None:
    """`gene` features are a diagnostic (CLAUDE.md); the table counts CDS/tRNA/rRNA only."""
    groups = load_gene_groups()
    genome = _genome(
        Feature("rbcL", "gene", (Interval(0, 100),)),
        Feature("rbcL", "CDS", (Interval(0, 100),)),
    )
    assert [row.gene.kind for row in assign(genome, groups)] == ["CDS"]


def test_every_gene_of_the_golden_genome_is_classified(data_dir: Path) -> None:
    """The claim the publication table rests on: the shipped file covers a real plastome."""
    from plastid.core.genbank import read_genbank

    genome = read_genbank(data_dir / "NC_050833.1.gb")
    groups = load_gene_groups()
    rows = assign(genome, groups)
    assert len(rows) == count_genes(genome).total == 132
    assert [row.gene.name for row in rows if row.group is None] == []
