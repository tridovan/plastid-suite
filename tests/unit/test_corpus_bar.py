"""The classification rule P3B-10's bar is read from — tested like any other claim in this repo.

The phase's Done-when compares two readings of 34 published plastomes and requires ours to differ
from the published annotation in fewer conclusion-level places than PGA's. A counting rule that
quietly favours one side would make that comparison meaningless, so the rule is pinned here: one
locus in one place matches, a different spelling is vocabulary, and anything else — a locus missing,
a locus invented, a gene not covered — is a conclusion.

The cases below are the ones the corpus actually produced (§8.50): the golden's `psbZ`/`lhbA`
pair, the IR's second copy, the trans-spliced `rps12`, and a tRNA whose ends moved by a base.
"""

from __future__ import annotations

import pytest
from tools.corpus_bar import Difference, by_gene, compare, features, summarise

from plastid.core.model import Gene, Interval


def gene(name: str, start: int, end: int, *, kind: str = "CDS", strand: int = 1) -> Gene:
    """A one-exon feature, written in 1-based inclusive coordinates like the spans it is read as."""
    return Gene(name=name, kind=kind, parts=(Interval(start - 1, end, strand),))


def spliced(name: str, exons: list[tuple[int, int]], *, kind: str = "CDS") -> Gene:
    return Gene(
        name=name, kind=kind, parts=tuple(Interval(start - 1, end, 1) for start, end in exons)
    )


def test_the_same_gene_in_the_same_place_is_no_difference() -> None:
    assert compare([gene("psbA", 100, 300)], [gene("psbA", 100, 300)]) == []


def test_another_spelling_of_one_locus_is_a_vocabulary_difference() -> None:
    """D5-09's whole subject: the same gene under the group's dialect, not a wrong call."""
    differences = compare([gene("psbZ", 100, 300)], [gene("lhbA", 100, 300)])

    assert len(differences) == 1
    assert (differences[0].verdict, differences[0].detail) == ("vocabulary", "name")


def test_a_gene_the_call_does_not_cover_is_a_conclusion_difference() -> None:
    """The golden's `petB`, measured: published 74233-75632, our call 74985-75632. The intron's
    first exon is 6 bp and `blastx` cannot align two codons across it, so the call starts at the
    second exon — 752 bp of the gene uncovered, and counted, not smoothed over (§8.50)."""
    differences = compare([gene("petB", 74_233, 75_632)], [gene("petB", 74_985, 75_632)])

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]
    assert differences[0].note == "-752 bp"


def test_one_envelope_two_different_exon_splits_is_a_conclusion() -> None:
    """The golden's `rpoC1`, measured: both readings span 20179-22962 and their inner exon
    boundaries differ by 27 bp. Comparing only the outer ends calls that pair identical — and the
    first version of this rule did exactly that, which is why the case is pinned here."""
    published = [spliced("rpoC1", [(20_179, 21_801), (22_531, 22_962)])]
    called = [spliced("rpoC1", [(20_179, 21_813), (22_516, 22_962)])]

    assert (published[0].parts[0].start, published[0].parts[-1].end) == (
        called[0].parts[0].start,
        called[0].parts[-1].end,
    )
    differences = compare(published, called)

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]


def test_an_end_that_moved_one_base_is_reported_as_a_boundary_not_a_conclusion() -> None:
    """The corpus's dominant class, measured on 719 of 1,264 matched pairs: ARAGORN ends a tRNA
    one base outside where the record does. Same gene, same anticodon, same place — reported, and
    not counted as a conclusion (§8.51). Kept to one base on purpose: the tighter setting, since a
    tolerance chosen to make a total come out right would measure nothing."""
    differences = compare(
        [gene("trnD-GUC", 30_541, 30_614, kind="tRNA", strand=-1)],
        [gene("trnD-GUC", 30_540, 30_615, kind="tRNA", strand=-1)],
    )

    assert [(d.verdict, d.detail) for d in differences] == [("boundary", "extent")]
    assert summarise(differences)["conclusion"] == 0


def test_an_end_that_moved_two_bases_is_a_conclusion() -> None:
    """One base is the tolerance, and this is the other side of it."""
    differences = compare(
        [gene("trnD-GUC", 30_541, 30_614, kind="tRNA", strand=-1)],
        [gene("trnD-GUC", 30_539, 30_615, kind="tRNA", strand=-1)],
    )

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]


def test_a_call_on_the_other_strand_is_never_within_tolerance() -> None:
    """The corpus carries one: our `trnS-UGA` on the minus strand where the record has plus.
    Identical coordinates on the other strand is a different call, not a boundary."""
    differences = compare(
        [gene("trnS-UGA", 35_472, 35_564, kind="tRNA", strand=1)],
        [gene("trnS-UGA", 35_472, 35_564, kind="tRNA", strand=-1)],
    )

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]


def test_a_different_number_of_exons_is_never_within_tolerance() -> None:
    """The golden's `petB`: the record has the 6 bp first exon, the call has only the second."""
    published = [spliced("petB", [(74_233, 74_238), (74_985, 75_632)])]
    called = [spliced("petB", [(74_985, 75_632)])]

    differences = compare(published, called)

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]


def test_an_overlap_under_half_is_two_differences_not_one() -> None:
    """Under half of the shorter feature is not the same gene in the same place — so the record's
    gene is missing *and* the call is an extra, which is what the reader has to be told.

    Half and above still matches (the `psbJ` case below): the rule is a threshold, and both sides
    of it are pinned here so a later change of the threshold cannot pass unnoticed.
    """
    under = compare([gene("rps19", 1000, 1200)], [gene("rps19", 1150, 1350)])
    assert sorted(d.detail for d in under) == ["extra", "missing"]

    over = compare([gene("rps19", 1000, 1200)], [gene("rps19", 1100, 1350)])
    assert sorted(d.detail for d in over) == ["extent"]  # 101 shared of 201 = just over half


def test_a_second_copy_in_the_inverted_repeat_is_matched_as_the_second_copy() -> None:
    """The IR puts every repeat gene twice; the second copy must pair with the second, never be
    reported as an extra call beside a missing one (which is what P3B-10's first corpus run
    measured, before the mirror was fixed — 12 records × 6 loci)."""
    published = [gene("rpl2", 1000, 2000), gene("rpl2", 30_000, 31_000)]
    called = [gene("rpl2", 1000, 2000), gene("rpl2", 30_000, 31_000)]

    assert compare(published, called) == []


def test_the_trans_spliced_rps12_is_compared_by_exon_not_by_envelope() -> None:
    """Measured: the published `rps12` envelope spans 68,976–139,634. An envelope comparison
    "overlaps" half the plastome and matches it against `ycf2`; by exon it is one gene."""
    published = [spliced("rps12", [(68_976, 68_990), (98_000, 98_400), (139_600, 139_634)])]
    called = [spliced("rps12", [(68_976, 68_990), (98_000, 98_400), (139_600, 139_634)])]

    assert compare(published, called) == []
    assert len(features(published, "CDS")[0].spans) == 3


def test_an_overlap_under_half_of_the_longer_still_matches_when_it_is_half_the_shorter() -> None:
    """The rule is stated on the shorter feature: a 60 bp call inside a 200 bp gene is half its
    own length and is the same gene, placed short — an extent difference, not a missing gene."""
    differences = compare([gene("psbJ", 1000, 1199)], [gene("psbJ", 1000, 1059)])

    assert [(d.verdict, d.detail) for d in differences] == [("conclusion", "extent")]


def test_the_counts_separate_what_a_reader_needs_to_tell_apart() -> None:
    published = [gene("psbZ", 100, 300), gene("rps7", 500, 700), gene("ndhB", 900, 1200)]
    called = [gene("lhbA", 100, 300), gene("rps7", 500, 660)]

    counts = summarise(compare(published, called))

    assert counts["vocabulary"] == 1
    assert counts["conclusion"] == 2  # rps7 short, ndhB not called
    assert counts["conclusion:missing"] == 1
    assert counts["conclusion:extent"] == 1
    assert counts["CDS:conclusion"] == 2
    assert counts["tRNA:conclusion"] == 0


def test_a_locus_called_twice_under_two_names_is_one_extra() -> None:
    """The measured shape of H-18: our caller emitted both `psbZ` and `lhbA` at 35895-36083.
    The first match takes the published feature; the second is an extra call, and the count
    says so rather than hiding the duplicate."""
    published = [gene("psbZ", 35_895, 36_083)]
    called = [gene("psbZ", 35_895, 36_083), gene("lhbA", 35_895, 36_083)]

    differences = compare(published, called)

    assert [(d.detail, d.called.split()[0]) for d in differences] == [("extra", "lhbA")]


def test_the_composition_helper_names_the_loci_behind_a_count() -> None:
    published = [gene("rps7", 500, 700), gene("rps7", 30_000, 30_200)]
    called = [gene("rps7", 500, 700)]

    differences = compare(published, called)

    assert by_gene(differences, "missing") == {"rps7": 1}
    assert by_gene(differences, "extra") == {}


@pytest.mark.parametrize("kind", ["CDS", "tRNA", "rRNA"])
def test_only_the_features_of_the_kind_asked_for_are_compared(kind: str) -> None:
    everything = [
        gene("psbA", 100, 300, kind="CDS"),
        gene("trnA-UGC", 400, 480, kind="tRNA"),
        gene("rrn16", 500, 1900, kind="rRNA"),
    ]

    compared = features(everything, kind)

    assert [feature.name for feature in compared] == [
        {"CDS": "psbA", "tRNA": "trnA-UGC", "rRNA": "rrn16"}[kind]
    ]


def test_a_difference_line_carries_both_spans_and_both_names() -> None:
    """The Done-when asks for every difference "kind + span, nothing summarized away"."""
    difference = Difference(
        "CDS",
        "conclusion",
        "extent",
        "petB 74233-75632(+1)",
        "petB 74985-75632(+1)",
        "-752 bp",
    )

    assert difference.line() == "CDS extent petB 74233-75632(+1) -> petB 74985-75632(+1) (-752 bp)"
