"""The plastid calling layer (P3B-07): its reference index, its floors, its repeat reading.

The rules themselves live in ``core/caller.py`` and are guarded by ``test_mito_caller.py``, which
P3B-07's extraction left untouched — that suite is the evidence the machinery did not move. What
is checked here is what the plastid path states *for itself*: which sequence a reference entry
carries, which names a call may print, and the pair of floors no other path may hand it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.caller import Hsp, ReferenceEntry
from plastid.core.genbank import read_genbank
from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, Interval
from plastid.core.plastid_caller import (
    DEFAULT_MIN_COVERAGE,
    DEFAULT_MIN_IDENTITY,
    PLASTID_TRANSL_TABLE,
    call_cds,
    call_rrna,
    read_inverted_repeat,
    reference_entries,
    reference_index,
    repeat_agreement,
    vocabulary_of,
)

DATA = Path(__file__).parents[1] / "data"


def _genome() -> Genome:
    """Three genes: a CDS with a protein, one without, and an rRNA."""
    return Genome(
        name="synthetic plastome",
        accession="SYN_PLASTID.1",
        sequence="ATGGCC" * 200 + "ACGT" * 50,
        genes=(
            Gene(name="psbA", kind="CDS", parts=(Interval(0, 30),), translation="MADV"),
            Gene(name="rps3", kind="CDS", parts=(Interval(40, 70),)),  # no /translation
            Gene(name="rrn16", kind="rRNA", parts=(Interval(90, 120),)),
        ),
        source={"organelle": "plastid"},
    )


def test_a_reference_entry_carries_the_records_own_translation() -> None:
    """Never a self-translation (D5-08 §8.21) — and no protein at all means no entry at all."""
    entries = reference_entries([_genome()], kind="CDS", vocabulary=("psbA",))
    assert [entry.sequence for entry in entries] == ["MADV"]
    assert entries[0].source_gene == "psbA"
    assert entries[0].key == "SYN_PLASTID.1|psbA|psbA"


def test_an_rrna_entry_carries_nucleotides_not_a_translation() -> None:
    entries = reference_entries([_genome()], kind="rRNA", vocabulary=("rrn16",))
    assert len(entries) == 1
    assert set(entries[0].sequence) <= set("ACGT")


def test_a_name_the_vocabulary_does_not_hold_is_kept_but_unnameable() -> None:
    """H-06's rule: the gene is reported, never renamed into a vocabulary that lacks it."""
    index = reference_index([_genome()], kind="CDS", vocabulary=("psbA",))
    unnameable = index.unnameable()
    assert [entry.source_gene for entry in unnameable] == []
    other = Genome(
        name="other",
        accession="SYN2",
        sequence="ACGT" * 10,
        genes=(Gene(name="ycf99", kind="CDS", parts=(Interval(0, 9),), translation="MADV"),),
        source={"organelle": "plastid"},
    )
    index = reference_index([_genome(), other], kind="CDS", vocabulary=("psbA",))
    assert [entry.source_gene for entry in index.unnameable()] == ["ycf99"]
    assert index.entry_for("SYN2|ycf99|-") is not None


def test_the_vocabulary_is_the_reference_sets_own_spelling() -> None:
    """D5-09 puts the naming step in front of the annotator, so A is what the set already says."""
    assert vocabulary_of([_genome()], kind="CDS") == ("psbA", "rps3")
    assert vocabulary_of([_genome()], kind="rRNA") == ("rrn16",)


def test_the_two_kinds_are_two_indexes_never_one_pool() -> None:
    """A blastn percentage and a blastx one are not one scale (D5-08 §3)."""
    genome = _genome()
    cds = reference_index([genome], kind="CDS", vocabulary=("psbA",))
    rrna = reference_index([genome], kind="rRNA", vocabulary=("rrn16",))
    assert {entry.source_gene for entry in cds.entries} == {"psbA"}
    assert {entry.source_gene for entry in rrna.entries} == {"rrn16"}
    with pytest.raises(ValueError, match="kind must be"):
        reference_index([genome], kind="tRNA", vocabulary=())


def test_the_plastid_callers_state_their_own_floors() -> None:
    """Measured on the answer key (§8.46); the mito path's coverage of 0.5 is not enough here."""
    assert (DEFAULT_MIN_IDENTITY, DEFAULT_MIN_COVERAGE) == (50.0, 0.7)
    assert PLASTID_TRANSL_TABLE == 11
    # the two wrappers must send their own pair on, not the shared core's required arguments
    empty = reference_index([], kind="CDS", vocabulary=())
    result = call_cds(_genome(), empty, ())
    assert (result.min_identity, result.min_coverage) == (50.0, 0.7)
    rrna = call_rrna(_genome(), reference_index([], kind="rRNA", vocabulary=()), ())
    assert (rrna.min_identity, rrna.min_coverage) == (50.0, 0.7)


def test_a_name_with_a_space_gets_a_header_a_blast_hit_can_come_back_as() -> None:
    """Measured 08/10/2026: NC_042827.1 names its rRNA ``23S ribosomal RNA``.

    BLAST cuts ``sseqid`` at the first whitespace, so the raw header came back as
    ``NC_042827.1|23S`` and the parser refused it — rightly: a key that can never match the entry
    it came from is worse than an error. The fix is in ``key`` itself, so the FASTA writer and
    the lookup cannot disagree.
    """
    entry = ReferenceEntry(
        accession="NC_042827.1",
        organism="Iris missouriensis",
        source_gene="23S ribosomal RNA",
        vocabulary="rrn23",
        sequence="ACGT",
    )
    assert entry.key == "NC_042827.1|23S_ribosomal_RNA|rrn23"
    assert " " not in entry.key


# ------------------------------------------------------------------ the inverted repeat reading


def _self_hit(qlo: int, qhi: int, slo: int, shi: int, *, pident: float, length: int) -> Hsp:
    """A blastn hit on the minus strand — what an inverted copy looks like to a self-search."""
    return Hsp(
        subject_key="NC_000000.1|repeat|-",
        pident=pident,
        length=length,
        qstart=qlo,
        qend=qhi,
        sstart=shi,
        send=slo,
        evalue=0.0,
        bitscore=1000.0,
        frame=None,
    )


def test_the_diagonal_is_not_a_repeat() -> None:
    """A self-search's longest, most identical hit is the sequence matching itself."""
    diagonal = Hsp(
        subject_key="self",
        pident=100.0,
        length=26000,
        qstart=1,
        qend=26000,
        sstart=1,
        send=26000,
        evalue=0.0,
        bitscore=50000.0,
        frame=None,
    )
    assert read_inverted_repeat([diagonal]) is None


def test_the_anti_diagonal_is_read_as_the_two_copies() -> None:
    hit = _self_hit(83128, 109341, 127511, 153724, pident=99.9, length=26214)
    repeat = read_inverted_repeat([hit])
    assert repeat is not None
    assert repeat.first == (83128, 109341)
    assert repeat.second == (127511, 153724)
    assert repeat.identity == 99.9


def test_a_short_or_sloppy_duplication_is_not_called_the_repeat() -> None:
    short = _self_hit(100, 900, 5000, 5800, pident=100.0, length=800)
    sloppy = _self_hit(100, 3000, 5000, 7900, pident=90.0, length=2900)
    assert read_inverted_repeat([short, sloppy]) is None


def test_the_two_readings_of_the_repeat_are_compared_not_merged() -> None:
    """The structural frame governs; the search is a second opinion, a difference is a finding."""
    structural = IRResult(
        ira=Interval(83127, 109341),
        irb=Interval(127510, 153724),
        identity=1.0,
        mismatches=0,
        compared_length=26214,
        exact_core=26214,
        mirror_constant=0,
        source="annotation",
    )
    agree = read_inverted_repeat(
        [_self_hit(83128, 109341, 127511, 153724, pident=99.9, length=26214)]
    )
    assert "agrees with the structural reading" in repeat_agreement(structural, agree)
    differ = read_inverted_repeat(
        [_self_hit(83128, 101341, 137511, 155724, pident=99.9, length=18214)]
    )
    message = repeat_agreement(structural, differ)
    assert (
        "differ" in message
        and "26214" in message
        and "the structural frame is the one kept" in message
    )
    assert "no anti-diagonal" in repeat_agreement(structural, None)
    assert "was not asked" in repeat_agreement(None, agree)


@pytest.mark.skipif(
    not (DATA / "NC_050833.1.gb").is_file(), reason="the golden record is not vendored"
)
def test_the_goldens_reference_index_names_every_cds_locus() -> None:
    """The pinned panel is the real shape: 86 CDS → one entry each, all of them nameable."""
    golden = read_genbank(DATA / "NC_050833.1.gb")
    vocabulary = vocabulary_of([golden], kind="CDS")
    index = reference_index([golden], kind="CDS", vocabulary=vocabulary)
    assert len(index.entries) == 86
    assert index.unnameable() == ()
    assert len(vocabulary) == 79  # 86 features, 79 loci: the inverted repeat's copies
    assert set(vocabulary_of([golden], kind="rRNA")) == {"rrn16", "rrn23", "rrn4.5", "rrn5"}
