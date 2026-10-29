"""A CDS is measured on its transcript — strand, part order and frame — not on a forward slice.

The codon and amino-acid modules read every CDS through `core.model.coding_sequence`, which reads
each part on its own strand in annotation order. A minus-strand CDS read as a forward slice gives
a different codon composition and a different protein, so the test builds a genome whose only CDS
is on the minus strand and compares the module's table with the transcript's own counts.
"""

from __future__ import annotations

from pathlib import Path

from plastid.core.codon import codon_counts, translate
from plastid.core.model import Gene, Genome, Interval, coding_sequence
from plastid.core.protein import composition
from plastid.core.sequence import reverse_complement
from plastid.modules import amino_acid, codon_usage
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

TRANSCRIPT = "ATGGCTTAA" * 3  # nine codons, the last one a stop
FLANK = "ACGTACGTAC" * 5


def _minus_strand_genome() -> Genome:
    sequence = FLANK + reverse_complement(TRANSCRIPT) + FLANK
    start = len(FLANK)
    gene = Gene(
        name="psbA",
        kind="CDS",
        parts=(Interval(start, start + len(TRANSCRIPT), strand=-1),),
        transl_table=11,
    )
    return Genome(name="Testus", accession="T.2", sequence=sequence, genes=(gene,))


def _ctx(tmp_path: Path, module_id: str) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", module_id, fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module_id, module_order=1)


def test_the_fixture_reads_the_minus_strand_cds_as_its_transcript() -> None:
    genome = _minus_strand_genome()
    assert coding_sequence(genome, genome.genes[0]) == TRANSCRIPT
    assert genome.subseq(genome.genes[0].parts[0]) != TRANSCRIPT


def test_codon_usage_measures_a_minus_strand_cds_on_its_transcript(tmp_path: Path) -> None:
    genome = _minus_strand_genome()
    result = codon_usage.run(genome, {}, _ctx(tmp_path, "codon_usage"))
    summary = result.tables["codon_usage"].rows[0]
    expected = {codon: count for codon, count in codon_counts(TRANSCRIPT).items() if count}
    observed = {row["codon"]: row["count"] for row in result.tables["rscu"].rows if row["count"]}
    assert observed == expected
    assert summary["n_codons"] == sum(expected.values())
    assert summary["internal_stops"] == translate(TRANSCRIPT)[:-1].count("*")


def test_amino_acid_measures_a_minus_strand_cds_on_its_transcript(tmp_path: Path) -> None:
    genome = _minus_strand_genome()
    result = amino_acid.run(genome, {}, _ctx(tmp_path, "amino_acid"))
    summary = result.tables["amino_acid"].rows[0]
    expected = composition(translate(TRANSCRIPT))
    observed = {
        row["residue"]: row["count"]
        for row in result.tables["amino_acid_composition"].rows
        if row["count"]
    }
    assert observed == {residue: count for residue, count in expected.items() if count}
    assert summary["n_residues"] == sum(expected.values())
