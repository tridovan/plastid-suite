"""ORF checks (P3-01): the classes, the negative cases, and what is refused rather than guessed.

The rule these tests encode is measured, not assumed: 5 of the golden genome's 86 CDS do not
start with ATG (PROJECT.md §6, TASKS.md pitfalls), so "must be ATG" would fail 3.7 % of real
plastid CDS. Every negative case here is a silent-wrong mode the check exists to make loud.
"""

from __future__ import annotations

import pytest

from plastid.core.model import Gene, Genome, Interval
from plastid.core.orf import (
    STATUS_ALTERNATIVE_START,
    STATUS_INTERNAL_STOP,
    STATUS_NO_START,
    STATUS_NO_STOP,
    STATUS_NOT_MULTIPLE_OF_THREE,
    STATUS_OK,
    STATUS_PSEUDOGENE,
    OrfError,
    annotated_not_coding,
    check_cds,
    check_genome,
)

_PAD = "A" * 12


def _genome(codons: str, name: str = "gene1", **gene_kwargs: object) -> Genome:
    """A linear genome holding one CDS at a known offset; padding keeps coordinates honest."""
    sequence = _PAD + codons + _PAD
    gene = Gene(
        name,
        "CDS",
        (Interval(len(_PAD), len(_PAD) + len(codons)),),
        **gene_kwargs,  # type: ignore[arg-type]
    )
    return Genome(name="synthetic", sequence=sequence, genes=(gene,), topology="linear")


def _check(codons: str, *, default_table: int = 11, **gene_kwargs: object):
    genome = _genome(codons, **gene_kwargs)
    return check_cds(genome, genome.genes[0], default_table=default_table)


def test_a_clean_cds_is_ok() -> None:
    check = _check("ATG" + "TTC" * 5 + "TAA")  # 21 bp = 7 codons
    assert check.status == STATUS_OK
    assert (check.start_codon, check.stop_codon) == ("ATG", "TAA")
    assert check.internal_stop_positions == ()
    assert check.length_bp == 21
    assert check.length_is_multiple_of_three


def test_a_table11_alternative_start_is_flagged_not_failed() -> None:
    """rps19 GTG / ndhD ATC / rpl33 ATA are real: measured across the panel, not invented."""
    for codon in ("GTG", "ATC", "ATA", "TTG", "CTG", "ATT"):
        check = _check(codon + "TTC" * 5 + "TAA")
        assert check.status == STATUS_ALTERNATIVE_START, codon
        assert check.start_class == "alternative", codon
        assert check.start_codon == codon


def test_acg_is_an_alternative_start_but_marked_non_canonical() -> None:
    """ACG is not in NCBI table 11's initiator list; it is used for rpl2/ndhD and is normally
    edited to AUG in the transcript. Reporting it as ``no_start`` would be a false alarm on two
    of the golden genome's 86 CDS, so it is flagged separately instead of hidden or excused."""
    check = _check("ACG" + "TTC" * 5 + "TAA")
    assert check.status == STATUS_ALTERNATIVE_START
    assert check.start_class == "non_canonical"


def test_a_codon_that_cannot_initiate_is_no_start() -> None:
    check = _check("AAA" + "TTC" * 5 + "TAA")
    assert check.status == STATUS_NO_START
    assert check.start_class == "invalid"


def test_a_missing_stop_is_reported() -> None:
    check = _check("ATG" + "TTC" * 5 + "GGG")
    assert check.status == STATUS_NO_STOP
    assert check.stop_codon == "GGG"


def test_an_internal_stop_is_reported_with_its_codon_position() -> None:
    check = _check("ATG" + "TTC" * 2 + "TGA" + "TTC" * 2 + "TAA")
    assert check.status == STATUS_INTERNAL_STOP
    # 1-based codon index inside the CDS: ATG(1) TTC(2) TTC(3) TGA(4)
    assert check.internal_stop_positions == (4,)


def test_a_length_not_divisible_by_three_is_reported() -> None:
    check = _check("ATG" + "TTC" * 5 + "TA")
    assert check.status == STATUS_NOT_MULTIPLE_OF_THREE
    assert not check.length_is_multiple_of_three


def test_the_frame_problem_outranks_the_others() -> None:
    """A shifted frame makes the later questions meaningless; report the frame, keep the facts."""
    check = _check("AAA" + "TTC" * 5 + "GG")
    assert check.status == STATUS_NOT_MULTIPLE_OF_THREE
    assert check.start_class == "invalid"  # still visible, just not the headline


def test_a_broken_pseudogene_is_marked_not_counted_as_a_failure() -> None:
    """The internal stop is still visible in the facts; the headline says what the copy is."""
    check = _check("ATG" + "TTC" * 2 + "TGA" + "TTC" * 2 + "TAA", pseudo=True)
    assert check.gene.pseudo
    assert check.status == STATUS_PSEUDOGENE
    assert check.is_failure is False
    assert check.internal_stop_positions == (4,)  # not hidden, just not the headline


def test_a_clean_pseudogene_is_still_reported_as_a_pseudogene() -> None:
    """It translates cleanly, but it is *annotated* as a pseudogene — that is what is reported."""
    check = _check("ATG" + "TTC" * 5 + "TAA", pseudo=True)
    assert check.status == STATUS_PSEUDOGENE
    assert check.is_failure is False


def test_codon_start_other_than_one_is_refused() -> None:
    """A ``/codon_start=2`` shifts the frame; translating it as 1 would answer a question the
    record did not ask. Refuse rather than translate."""
    with pytest.raises(OrfError, match="codon_start"):
        _check("ATG" + "TTC" * 5 + "TAA", codon_start=2)


def test_an_unknown_translation_table_is_refused() -> None:
    with pytest.raises(OrfError, match="translation table"):
        _check("ATG" + "TTC" * 5 + "TAA", transl_table=4)


def test_the_records_own_table_is_used_and_the_default_applies_when_absent() -> None:
    alt = _check("GTG" + "TTC" * 5 + "TAA", transl_table=11)
    assert alt.transl_table == 11 and alt.status == STATUS_ALTERNATIVE_START

    # table 1 does not list GTG as an initiator: the same sequence, read under the fallback,
    # is a start problem — the record's declaration is what decides, not our preference
    fallback = _check("GTG" + "TTC" * 5 + "TAA", default_table=1)
    assert fallback.transl_table == 1 and fallback.status == STATUS_NO_START


def test_a_minus_strand_cds_is_read_on_its_own_strand() -> None:
    codons = "ATG" + "TTC" * 5 + "TAA"
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse = "".join(complement[base] for base in reversed(codons))
    sequence = _PAD + reverse + _PAD
    gene = Gene(
        "gene1",
        "CDS",
        (Interval(len(_PAD), len(_PAD) + len(reverse), -1),),
    )
    genome = Genome(name="synthetic", sequence=sequence, genes=(gene,), topology="linear")
    check = check_cds(genome, gene, default_table=11)
    assert check.status == STATUS_OK
    assert check.start_codon == "ATG"


def test_a_gene_feature_with_no_cds_over_it_is_annotated_not_coding() -> None:
    """The truncated-ycf1 shape: a `gene` feature with no CDS anywhere under it."""
    sequence = _PAD + "ATG" + "TTC" * 5 + "TAA" + _PAD
    genome = Genome(
        name="synthetic",
        sequence=sequence,
        genes=(
            Gene("ycf1", "gene", (Interval(0, 5),), pseudo=True),
            Gene("ycf1", "CDS", (Interval(len(_PAD), len(_PAD) + 18),)),
        ),
        topology="linear",
    )
    orphans = annotated_not_coding(genome)
    assert [gene.name for gene in orphans] == ["ycf1"]
    assert orphans[0].kind == "gene" and orphans[0].pseudo


def test_a_gene_feature_covered_by_a_cds_is_not_reported() -> None:
    sequence = _PAD + "ATG" + "TTC" * 5 + "TAA" + _PAD
    genome = Genome(
        name="synthetic",
        sequence=sequence,
        genes=(
            Gene("rbcL", "gene", (Interval(len(_PAD), len(_PAD) + 21),)),
            Gene("rbcL", "CDS", (Interval(len(_PAD), len(_PAD) + 21),)),
        ),
        topology="linear",
    )
    assert annotated_not_coding(genome) == ()


def test_a_neighbouring_cds_over_the_tail_does_not_excuse_the_claim() -> None:
    """Measured on the golden genome: the truncated ycf1 overlaps ndhF's first 21 bp. A rule
    that accepts *any* CDS overlap would call this copy coding and hide it."""
    sequence = _PAD + "ATG" + "TTC" * 5 + "TAA" + _PAD
    genome = Genome(
        name="synthetic",
        sequence=sequence,
        genes=(
            Gene("ycf1", "gene", (Interval(len(_PAD), len(_PAD) + 21),), pseudo=True),
            Gene("ndhF", "CDS", (Interval(len(_PAD) + 18, len(_PAD) + 21),)),
        ),
        topology="linear",
    )
    assert [gene.name for gene in annotated_not_coding(genome)] == ["ycf1"]


def test_a_trna_gene_feature_is_implemented_by_its_own_feature() -> None:
    """A tRNA gene has no CDS and is not a defect: the tRNA feature is what implements it."""
    sequence = _PAD + "TTT" * 7 + _PAD
    genome = Genome(
        name="synthetic",
        sequence=sequence,
        genes=(
            Gene("trnK-UUU", "gene", (Interval(len(_PAD), len(_PAD) + 21),)),
            Gene("trnK-UUU", "tRNA", (Interval(len(_PAD), len(_PAD) + 21),)),
            Gene("trnI", "gene", (Interval(0, 3),), pseudo=True),
        ),
        topology="linear",
    )
    assert [gene.name for gene in annotated_not_coding(genome)] == ["trnI"]


def test_check_genome_covers_every_cds_feature_including_pseudogenes() -> None:
    first, second = "ATG" + "TTC" * 5 + "TAA", "ATG" + "TTC" * 5 + "TAG"  # 21 bp each
    genome = Genome(
        name="synthetic",
        sequence=_PAD + first + second + _PAD,
        genes=(
            Gene("a", "CDS", (Interval(len(_PAD), len(_PAD) + len(first)),)),
            Gene(
                "b",
                "CDS",
                (Interval(len(_PAD) + len(first), len(_PAD) + len(first) + len(second)),),
                pseudo=True,
            ),
            Gene("a", "gene", (Interval(len(_PAD), len(_PAD) + len(first)),)),
        ),
        topology="linear",
    )
    checks = check_genome(genome, default_table=11)
    assert [check.gene.name for check in checks] == ["a", "b"]
    assert [check.status for check in checks] == [STATUS_OK, STATUS_PSEUDOGENE]
