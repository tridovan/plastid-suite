"""GenBank reader and writer — real fixtures plus error handling. Never touches the network."""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank, to_genbank
from plastid.core.model import Gene, Genome, Interval, count_genes


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


# --- writing back out (P3A-01) ----------------------------------------------------------------
#
# The writer exists so an annotator's output can become a file our own reader reads back
# identically, and a submission can be built from. Every test here is therefore a round trip —
# read · write · read again · compare — and the negative case is a feature the sequence cannot
# hold. What the writer deliberately does NOT reproduce is stated in its docstring, not hidden.


def _assert_same_feature(before: Gene, after: Gene) -> None:
    assert after.kind == before.kind
    assert after.name == before.name
    assert after.product == before.product
    assert after.pseudo == before.pseudo
    assert after.trans_splicing == before.trans_splicing
    assert after.transl_table == before.transl_table
    assert after.codon_start == before.codon_start
    assert after.exception == before.exception
    # Interval is frozen, so this compares coordinates *and* the partial flags
    assert after.parts == before.parts


def test_the_golden_record_round_trips(golden_genome: Genome, tmp_path: Path) -> None:
    """265 features come back with the same kind, name, span, product, table and partial flags."""
    written = tmp_path / "round_trip.gb"
    written.write_text(to_genbank(golden_genome), encoding="utf-8")
    back = read_genbank(written)

    assert back.sequence == golden_genome.sequence
    assert back.accession == golden_genome.accession
    assert back.name == golden_genome.name
    assert back.topology == golden_genome.topology
    assert back.source == golden_genome.source
    assert back.repeat_regions == golden_genome.repeat_regions
    assert len(back.genes) == len(golden_genome.genes) == 265
    for before, after in zip(golden_genome.genes, back.genes, strict=True):
        _assert_same_feature(before, after)


def _round_trip(tmp_path: Path, genome: Genome) -> Genome:
    """Write a genome out and read it straight back — the shape every test below needs."""
    written = tmp_path / "round_trip.gb"
    written.write_text(to_genbank(genome), encoding="utf-8")
    return read_genbank(written)


def test_a_joined_minus_strand_cds_keeps_its_parts_in_transcript_order(tmp_path: Path) -> None:
    """rps12's shape: two parts on the minus strand, and the order of the parts is the message.

    Re-sorting them would silently produce a different protein (measured on both golden copies,
    see ``model.coding_sequence``), so the round trip has to preserve the order, not just the set.
    """
    gene = Gene(
        name="rps12",
        kind="CDS",
        product="ribosomal protein S12",
        parts=(Interval(20, 30, strand=-1), Interval(2, 10, strand=-1)),
        trans_splicing=True,
        transl_table=11,
    )
    text = to_genbank(Genome(name="Testus testus", sequence="ACGT" * 20, genes=(gene,)))
    assert "complement(join(" in text

    back = _round_trip(tmp_path, Genome(name="Testus testus", sequence="ACGT" * 20, genes=(gene,)))
    assert back.genes[0].parts == gene.parts
    assert back.genes[0].strand == -1


def test_a_partial_coordinate_keeps_its_marker(tmp_path: Path) -> None:
    """``<`` / ``>`` say the coordinate is incomplete — dropping one turns "runs off the contig"
    into "starts here", which is a different claim in a submitted record."""
    gene = Gene(
        name="ycf1",
        kind="CDS",
        parts=(Interval(0, 12, strand=1, partial_start=True, partial_end=True),),
        transl_table=11,
    )
    genome = Genome(name="Testus testus", sequence="ATGAAATAGCCCTAA", genes=(gene,))
    assert "<1..>12" in to_genbank(genome)

    back = _round_trip(tmp_path, genome)
    assert back.genes[0].parts[0].partial_start is True
    assert back.genes[0].parts[0].partial_end is True


def test_a_name_that_is_its_own_product_is_written_once(tmp_path: Path) -> None:
    """A record with no ``/gene`` names its features by their product sentence — measured on the
    stand-in mitogenome ``OX459128.1``. Writing both would claim a gene *symbol* that is a
    sentence, so ``/gene`` is written only when it says something ``/product`` does not."""
    gene = Gene(
        name="NADH dehydrogenase subunit 1",
        kind="CDS",
        product="NADH dehydrogenase subunit 1",
        parts=(Interval(0, 9),),
    )
    genome = Genome(name="Testus testus", sequence="ATGAAATAG", genes=(gene,))
    text = to_genbank(genome)
    assert "/gene=" not in text
    assert '/product="NADH dehydrogenase subunit 1"' in text

    _assert_same_feature(gene, _round_trip(tmp_path, genome).genes[0])


def test_bare_qualifiers_are_written_without_a_value(tmp_path: Path) -> None:
    """``/pseudo`` and ``/trans_splicing`` carry no value in every real record, including the
    golden one — BioPython writes ``=""`` for an empty string, so the value has to be ``None``."""
    gene = Gene(
        name="rps12",
        kind="CDS",
        parts=(Interval(0, 9),),
        pseudo=True,
        trans_splicing=True,
    )
    genome = Genome(name="Testus testus", sequence="ATGAAATAG", genes=(gene,))
    text = to_genbank(genome)
    assert "/pseudo\n" in text
    assert "/trans_splicing\n" in text
    assert '=""' not in text

    back = _round_trip(tmp_path, genome)
    assert back.genes[0].pseudo is True
    assert back.genes[0].trans_splicing is True


def test_the_source_feature_carries_the_records_own_qualifiers(
    golden_genome: Genome, tmp_path: Path
) -> None:
    """The source feature is where a submission states what the sample *is* — organism,
    organelle, mol_type, cultivar. It is written from ``Genome.source`` and read back whole."""
    text = to_genbank(golden_genome)
    assert '/organism="Iris domestica"' in text
    assert '/organelle="plastid:chloroplast"' in text
    assert '/mol_type="genomic DNA"' in text
    assert _round_trip(tmp_path, golden_genome).source == golden_genome.source


def test_a_feature_past_the_sequence_end_is_refused() -> None:
    """Fail loud: writing ``1..30`` for a 9 bp sequence produces a file that looks fine and is
    a lie. An origin-spanning feature on a circular genome is a join, not an out-of-range part."""
    genome = Genome(
        name="Testus testus",
        sequence="ATGAAATAG",
        genes=(Gene(name="atp1", kind="CDS", parts=(Interval(0, 30),)),),
    )
    with pytest.raises(ValueError, match="past the end"):
        to_genbank(genome)
