"""GFF3 reading: how an annotator's own feature file becomes the same ``Gene`` objects a
downloaded GenBank record produces (P3-03).

The rules under test are the ones that decide a gene count, so they are pinned here rather than
discovered in a golden number later: which rows make one feature, in what order their parts come
back, and which malformed rows are refused by name instead of quietly dropped.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from plastid.core.gff import GffError, genes_from_gff, read_gff_rows

_HEADER = "##gff-version 3\n##source-version Chloe 1.0\n"


def _row(kind: str, start: int, end: int, strand: str, attributes: str) -> str:
    return f"NC_050833.1\t.\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"


# A minus-strand CDS with one intron, written the way Chloe writes it: the exons are separate rows
# sharing a parent, in transcript order (5'→3', so descending coordinates on the minus strand).
_INTRON_CDS = _HEADER + (
    _row("gene", 4697, 5808, "-", "ID=g1;gene=rps16")
    + _row("CDS", 5769, 5808, "-", "ID=c1;parent=g1;gene=rps16")
    + _row("CDS", 4697, 4908, "-", "ID=c1;parent=g1;gene=rps16")
    + _row("intron", 4909, 5768, "-", "ID=i1;parent=g1;gene=rps16")
)


def test_rows_sharing_a_parent_become_one_feature_in_file_order() -> None:
    genes = genes_from_gff(_INTRON_CDS)
    assert [gene.kind for gene in genes] == ["gene", "CDS"]
    cds = genes[1]
    assert cds.name == "rps16"
    # file order, not sorted: it is the transcript order the reading frame depends on
    assert [(part.start, part.end) for part in cds.parts] == [(5768, 5808), (4696, 4908)]
    assert cds.strand == -1
    assert cds.length == 40 + 212


def test_a_gene_row_and_its_child_rows_stay_two_features() -> None:
    """``gene`` and ``CDS`` share an id space in some files and not others; the kind is part of
    the identity, so a parent row and its child never collapse into one feature."""
    genes = genes_from_gff(_INTRON_CDS)
    assert [(gene.kind, gene.name, gene.length) for gene in genes] == [
        ("gene", "rps16", 1112),  # the row's own span, intron included
        ("CDS", "rps16", 252),  # exons only
    ]


def test_other_kinds_are_left_out_not_mistaken_for_genes() -> None:
    kinds = {gene.kind for gene in genes_from_gff(_INTRON_CDS)}
    assert "intron" not in kinds  # the toolkit derives introns from the exon gaps itself


def test_the_name_comes_from_the_parent_when_the_child_has_none() -> None:
    text = _HEADER + (
        _row("gene", 100, 200, "+", "ID=g1;gene=psbA")
        + _row("CDS", 100, 200, "+", "ID=c1;parent=g1")
    )
    assert [gene.name for gene in genes_from_gff(text)] == ["psbA", "psbA"]


def test_pseudo_and_translation_table_are_read_when_the_file_declares_them() -> None:
    text = _HEADER + _row("CDS", 100, 200, "+", "ID=c1;gene=ycf1;pseudo=true;transl_table=11")
    gene = genes_from_gff(text)[0]
    assert gene.pseudo is True
    assert gene.transl_table == 11


def test_a_file_without_a_translation_table_says_so_rather_than_assuming_one() -> None:
    text = _HEADER + _row("CDS", 100, 200, "+", "ID=c1;gene=psbA")
    assert genes_from_gff(text)[0].transl_table is None


def test_a_row_with_neither_id_nor_parent_is_a_named_error() -> None:
    text = _HEADER + _row("CDS", 100, 200, "+", "gene=psbA")
    with pytest.raises(GffError, match="neither ID nor parent"):
        genes_from_gff(text)


def test_a_feature_with_no_name_is_a_named_error() -> None:
    """An unnamed feature would be counted as a gene and dropped by every name-based analysis."""
    text = _HEADER + _row("CDS", 100, 200, "+", "ID=c1")
    with pytest.raises(GffError, match="no gene or name"):
        genes_from_gff(text)


def test_an_unknown_strand_is_refused_because_the_frame_depends_on_it() -> None:
    text = _HEADER + _row("CDS", 100, 200, ".", "ID=c1;gene=psbA")
    with pytest.raises(GffError, match="strand"):
        genes_from_gff(text)


@pytest.mark.parametrize(
    ("line", "match"),
    [
        ("NC_050833.1\t.\tCDS\t100\t.\tgene=psbA", "expected 9 tab-separated fields"),
        ("NC_050833.1\t.\tCDS\tx\t200\t.\t+\t.\tID=c1;gene=psbA", "not an integer"),
        (_row("CDS", 100, 200, "+", "ID"), "attribute"),
    ],
)
def test_a_malformed_row_is_refused_with_its_line_number(line: str, match: str) -> None:
    with pytest.raises(GffError, match=match) as caught:
        read_gff_rows(_HEADER + line)
    assert "line 3" in str(caught.value)  # the two header lines come first


def test_a_translation_table_that_is_not_a_number_is_refused() -> None:
    text = _HEADER + _row("CDS", 100, 200, "+", "ID=c1;gene=psbA;transl_table=x")
    with pytest.raises(GffError, match="transl_table='x' is not an integer"):
        genes_from_gff(text)


def test_two_rows_of_one_feature_may_not_declare_two_translation_tables() -> None:
    text = _HEADER + (
        _row("CDS", 100, 150, "+", "ID=c1;gene=rps16;transl_table=11")
        + _row("CDS", 160, 200, "+", "ID=c1;gene=rps16;transl_table=1")
    )
    with pytest.raises(GffError, match="different translation tables"):
        genes_from_gff(text)


def test_coordinates_are_read_1_based_inclusive_like_every_other_input() -> None:
    row = read_gff_rows(_HEADER + _row("CDS", 100, 200, "+", "ID=c1;gene=psbA"))[0]
    assert (row.start, row.end) == (100, 200)
    gene = genes_from_gff(_HEADER + _row("CDS", 100, 200, "+", "ID=c1;gene=psbA"))[0]
    assert gene.span == (100, 200)
    assert gene.parts[0].to_1based() == (100, 200)


def test_reading_stops_at_the_fasta_directive() -> None:
    """A GFF may carry the sequence after ``##FASTA``; those lines are not features."""
    text = _INTRON_CDS + "##FASTA\n>NC_050833.1\nACGT\tACGT\tACGT\n"
    assert len(read_gff_rows(text)) == 4


def test_an_empty_file_is_an_empty_feature_list_not_an_error() -> None:
    assert genes_from_gff(_HEADER) == ()


def test_an_annotation_with_no_genes_is_visible_as_such() -> None:
    """The reader does not invent a feature to avoid returning nothing: a caller that needs genes
    checks for them (PGA's silent-empty failure mode is why)."""
    text = _HEADER + _row("intron", 100, 200, "+", "ID=i1;gene=psbA")
    assert genes_from_gff(text) == ()


# The real thing: what Chloe 0.1.18 wrote for the golden genome (tests/data/PROVENANCE.md).
def test_a_real_annotator_run_reads_back_as_the_golden_gene_counts(data_dir: Path) -> None:
    """The golden contract, reached from an annotator's own file: **132 = 86 CDS + 38 tRNA +
    8 rRNA**. This is the check a gene count needs and a gene count alone cannot replace — the
    same file read without --no-transform also counts 132 while placing a third of the genes in
    another frame (§8.17)."""
    from plastid.core.genbank import read_genbank
    from plastid.core.model import Genome, count_genes

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    genes = genes_from_gff((data_dir / "chloe_NC_050833.1.gff").read_text(encoding="utf-8"))
    counts = count_genes(Genome(name=golden.name, sequence=golden.sequence, genes=genes))
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)
    assert counts.total == 132
    # every feature came back named and every part with a direction: nothing was skipped, and
    # nothing was given a made-up strand, to reach that count. (A gene's own strand is 0 when its
    # parts disagree — the trans-spliced rps12 copies cross from the LSC into an inverted repeat.)
    assert all(gene.name for gene in genes)
    assert all(part.strand in (1, -1) for gene in genes for part in gene.parts)


def test_a_real_annotator_runs_cds_all_translate(data_dir: Path) -> None:
    """The coordinates in that file are this sequence's own: every CDS translates cleanly, with
    the same 80 + 6 split the golden record's own annotation gives (81 + 5 — Chloe calls one
    start differently, and neither call fails the check)."""
    from plastid.core.genbank import read_genbank
    from plastid.core.model import Genome
    from plastid.core.orf import check_genome

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    genes = genes_from_gff((data_dir / "chloe_NC_050833.1.gff").read_text(encoding="utf-8"))
    checks = check_genome(
        Genome(name=golden.name, sequence=golden.sequence, genes=genes), default_table=11
    )
    statuses = Counter(check.status for check in checks)
    assert statuses["ok"] + statuses["alternative_start"] == 86
    assert not [check for check in checks if check.is_failure]


def test_a_real_annotator_run_keeps_the_trans_spliced_gene_in_transcript_order(
    data_dir: Path,
) -> None:
    """rps12 is trans-spliced across the LSC and both inverted repeats. Its exon rows must come
    back in the file's order — the order that spells the protein — and not sorted by coordinate:
    sorted, the second exon's 3 bp shift would slide the frame of the whole gene."""
    from Bio.Seq import Seq

    from plastid.core.genbank import read_genbank
    from plastid.core.model import Genome, coding_sequence

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    genes = genes_from_gff((data_dir / "chloe_NC_050833.1.gff").read_text(encoding="utf-8"))
    genome = Genome(name=golden.name, sequence=golden.sequence, genes=genes)
    copies = [gene for gene in genes if gene.kind == "CDS" and gene.name == "rps12"]
    assert len(copies) == 2  # one per inverted-repeat copy, as the record has it
    proteins = {str(Seq(coding_sequence(genome, copy)).translate(table=11)) for copy in copies}
    assert len(proteins) == 1  # the same protein from both copies
    assert proteins.pop().startswith("MPTIKQLIRNTRQPIRNVTKSPALRGCPQRRGTCTRVYTI")
