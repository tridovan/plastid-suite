"""The plastid calling layer (P3B-07): its reference index, its floors, its repeat reading.

The rules themselves live in ``core/caller.py`` and are guarded by ``test_mito_caller.py``, which
P3B-07's extraction left untouched — that suite is the evidence the machinery did not move. What
is checked here is what the plastid path states *for itself*: which sequence a reference entry
carries, which names a call may print, and the pair of floors no other path may hand it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core import caller as core
from plastid.core.caller import Hsp, ReferenceEntry, ReferenceIndex
from plastid.core.genbank import read_genbank
from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, Interval
from plastid.core.plastid_caller import (
    DEFAULT_MIN_COVERAGE,
    DEFAULT_MIN_IDENTITY,
    PLASTID_TRANSL_TABLE,
    TRUNCATED_DUPLICATE_FRACTION,
    call_cds,
    call_rrna,
    read_inverted_repeat,
    reference_entries,
    reference_index,
    repeat_agreement,
    truncated_annotations,
    vocabulary_of,
)
from plastid.core.sequence import reverse_complement

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
    """Never a self-translation (D5-08) — and no protein at all means no entry at all."""
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
    """A blastn percentage and a blastx one are not one scale (D5-08)."""
    genome = _genome()
    cds = reference_index([genome], kind="CDS", vocabulary=("psbA",))
    rrna = reference_index([genome], kind="rRNA", vocabulary=("rrn16",))
    assert {entry.source_gene for entry in cds.entries} == {"psbA"}
    assert {entry.source_gene for entry in rrna.entries} == {"rrn16"}
    with pytest.raises(ValueError, match="kind must be"):
        reference_index([genome], kind="tRNA", vocabulary=())


def test_the_plastid_callers_state_their_own_floors() -> None:
    """Measured on the answer key; the mito path's coverage of 0.5 is not enough here."""
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


# ---------------------------------------------------------------- one codon before the alignment


def _one_codon_genome(upstream: str, *, reverse: bool = False) -> Genome:
    """38 codons and a stop: ``upstream`` + 37 × GCT + TAA, plus filler."""
    coding = upstream + "GCT" * 37 + "TAA"
    sequence = reverse_complement(coding) if reverse else coding
    return Genome(
        name="one-codon fixture",
        accession="SYN_ONE_CODON.1",
        sequence=sequence + "ATGC" * 20,
        genes=(),
        source={"organelle": "plastid"},
    )


def _psbl_reference() -> ReferenceEntry:
    """A 38-residue reference; only its length and key matter to the caller."""
    return ReferenceEntry(
        accession="NC_000000.1",
        organism="Testus testus",
        source_gene="psbL",
        vocabulary="psbL",
        sequence="M" * 38,
    )


def _missing_first_codon(entry: ReferenceEntry, *, frame: int = 1) -> Hsp:
    """The alignment starts at codon 2 and ends at codon 38, so the 5' end is not reached."""
    return Hsp(
        subject_key=entry.key,
        pident=100.0,
        length=37,
        qstart=4,
        qend=114,
        sstart=2,
        send=38,
        evalue=0.0,
        bitscore=200.0,
        frame=frame,
    )


def test_an_initiator_one_codon_before_the_alignment_is_read() -> None:
    """Measured 27/10/2026 on 39 Crataegus records: psbL's alignment starts at codon 2.

    The codon before it is ACG, an accepted plastid initiator. Reading that one codon makes the
    call 117 bp — the reference's own length — with a complete ORF; leaving it partial made 27
    calls that no ORF check accepts. One codon is read, never a walk to a distant ATG.
    """
    entry = _psbl_reference()
    genome = _one_codon_genome("ACG")

    index = ReferenceIndex.from_entries([entry], vocabulary=("psbL",))
    call = call_cds(genome, index, (_missing_first_codon(entry),)).one("psbL")

    assert call.parts[0].to_1based() == (1, 117)
    assert not call.partial_start
    assert call.parts[0].start == 0  # the upstream ACG is now the first codon
    assert any("one codon is read" in note for note in call.notes)


def test_a_non_initiator_one_codon_before_the_alignment_stays_partial() -> None:
    """ACA is not an initiator in table 11, so the partial start stands; nothing is invented."""
    entry = _psbl_reference()
    genome = _one_codon_genome("ACA")

    index = ReferenceIndex.from_entries([entry], vocabulary=("psbL",))
    call = call_cds(genome, index, (_missing_first_codon(entry),)).one("psbL")

    assert call.parts[0].to_1based() == (4, 117)
    assert call.partial_start
    assert any("not moved to the next ATG" in note for note in call.notes)


def test_the_same_read_works_on_the_minus_strand() -> None:
    """psbL is on the minus strand in the measured records; the 5' end is the high coordinate."""
    entry = _psbl_reference()
    genome = _one_codon_genome("ACG", reverse=True)

    index = ReferenceIndex.from_entries([entry], vocabulary=("psbL",))
    call = call_cds(genome, index, (_missing_first_codon(entry, frame=-1),)).one("psbL")

    assert call.parts[0].to_1based() == (1, 117)
    assert call.parts[0].strand == -1
    assert not call.partial_start


def test_the_mito_default_does_not_read_the_upstream_codon() -> None:
    """The one-codon read is opt-in: the frozen mito path keeps its previous partial start."""
    from plastid.core import caller as core

    entry = _psbl_reference()
    genome = _one_codon_genome("ACG")

    call = core.call_cds(
        genome,
        ReferenceIndex.from_entries([entry], vocabulary=("psbL",)),
        (_missing_first_codon(entry),),
        min_identity=50.0,
        min_coverage=0.7,
        transl_table=11,
    ).one("psbL")

    assert call.parts[0].to_1based() == (4, 117)
    assert call.partial_start


# ---------------------------------------------------------------- a locus a record annotated twice


def _twice(first: str, second: str) -> Genome:
    """One record carrying two CDS features under one name — the measured ycf1 shape."""
    return Genome(
        name="doubly annotated",
        accession="SYN_DUP.1",
        sequence="ATGGCC" * 900,
        genes=(
            Gene(name="ycf1", kind="CDS", parts=(Interval(0, 30),), translation=first),
            Gene(name="ycf1", kind="CDS", parts=(Interval(40, 70),), translation=second),
        ),
        source={"organelle": "plastid"},
    )


def test_a_locus_a_record_annotated_twice_is_searched_once() -> None:
    """Measured 08/10/2026 (H-19): 24 of the 34 pinned records carry a **truncated** ``ycf1``
    beside the complete one, and the truncated copy scores *higher* on identity — so the golden's
    5,421 bp gene was called as 900 bp. Structure decides first (D5-08), and the feature the
    run does not search is **named**, not dropped in silence.
    """
    long_protein = "M" + "A" * 1799
    genome = _twice(long_protein, long_protein[:300])
    entries = reference_entries([genome], kind="CDS")
    assert [len(entry.sequence) for entry in entries] == [1800]
    reported = truncated_annotations([genome], kind="CDS")
    assert len(reported) == 1
    assert (reported[0].locus, reported[0].kept, reported[0].dropped) == ("ycf1", 1800, 300)
    assert reported[0].ratio == pytest.approx(300 / 1800)


def test_a_two_residue_difference_is_not_a_truncated_duplicate() -> None:
    """The measured boundary: the pinned set's 243 within-record duplicate pairs sit **either**
    below 0.2992 **or** at 0.9927 and above — nothing in between. A rule that called every unequal
    pair a truncation would drop ``rpl2``'s own 271-against-273-residue pair.
    """
    assert TRUNCATED_DUPLICATE_FRACTION == 0.5
    genome = _twice("M" + "A" * 272, "M" + "A" * 270)
    assert len(reference_entries([genome], kind="CDS")) == 2
    assert truncated_annotations([genome], kind="CDS") == ()


def test_the_inverted_repeats_two_copies_are_two_entries() -> None:
    """Equal length is the repeat's shape (216 of the 243 measured pairs). Those two features are
    one gene in two copies, and both are searched — ``core/caller`` assigns an HSP to the smallest
    reference that contains it."""
    protein = "M" + "A" * 299
    genome = _twice(protein, protein)
    assert len(reference_entries([genome], kind="CDS")) == 2
    assert truncated_annotations([genome], kind="CDS") == ()


def test_a_lone_short_feature_is_still_searchable() -> None:
    """The rule fires on a *comparison*, and one feature offers none — so it does not fire here.

    A record that annotates a locus once, briefly, may be carrying a short ortholog; a partial CDS
    is the audit's own case and it blocks on the coordinates themselves (the audit blocked
    ``NC_024936.1`` that way). Dropping this one would take a locus out of the search on no
    evidence at all — which is why the audit, not this rule, is where such a record is answered.
    """
    genome = Genome(
        name="annotated once, briefly",
        accession="SYN_ONCE.1",
        sequence="ATGGCC" * 900,
        genes=(
            Gene(name="ycf1", kind="CDS", parts=(Interval(0, 30),), translation="M" + "A" * 299),
        ),
        source={"organelle": "plastid"},
    )
    assert [len(entry.sequence) for entry in reference_entries([genome], kind="CDS")] == [300]
    assert truncated_annotations([genome], kind="CDS") == ()


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


# ------------------------------------------------- structure decides before the score (H-20)


def _entry(name: str, protein: str, accession: str) -> ReferenceEntry:
    return ReferenceEntry(
        accession=accession,
        organism="Testus testus",
        source_gene=name,
        vocabulary=name,
        sequence=protein,
    )


def _one_hsp(entry: ReferenceEntry, *, qstart: int, qend: int, pident: float, bitscore: float):
    from plastid.core.caller import Hsp

    return Hsp(
        subject_key=entry.key,
        pident=pident,
        length=len(entry.sequence),
        qstart=qstart,
        qend=qend,
        sstart=1,
        send=len(entry.sequence),
        evalue=1e-30,
        bitscore=bitscore,
        frame=1,
    )


def test_a_record_that_starts_a_gene_late_does_not_win_on_identity() -> None:
    """Measured on the golden's ``rps19``: three references spell it 225 bp / 74 aa, and
    that protein is the **suffix** of the 92 aa protein the target's own record carries — the same
    gene annotated from an inner ATG. They score 100.00 % against 98.91 % for the complete ones,
    so identity alone transfers the truncation and the call lands 54 bp short of the published CDS.

    Structure decides before the score (D5-08): H-19 keeps a truncated duplicate out of the
    *index*; this is the same rule between candidates.
    """
    protein = "MK" + "A" * 90
    complete = _entry("rps19", protein, "NC_000001.1")
    late = _entry("rps19", protein[18:], "NC_000002.1")  # the complete protein's last 74
    index = ReferenceIndex.from_entries((complete, late), vocabulary=("rps19",))
    hsps = (
        _one_hsp(complete, qstart=1, qend=279, pident=98.9, bitscore=189.0),
        _one_hsp(late, qstart=19, qend=240, pident=100.0, bitscore=154.0),
    )

    plain = core.call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
    structured = call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)

    assert plain.one("rps19").reference.accession == "NC_000002.1"  # identity picks the tail
    assert structured.one("rps19").reference.accession == "NC_000001.1"  # the plastid path does not
    # and it is a *choice*, not a silence: the loser is still reported with its own numbers
    assert {row.reference_accession: row.chosen for row in structured.competing} == {
        "NC_000001.1": True,
        "NC_000002.1": False,
    }


def test_a_gene_cut_short_is_not_read_as_a_late_start() -> None:
    """A prefix is a different claim from a suffix, and the rule deliberately leaves it alone:
    no measurement yet says which of the two a record meant (see ``caller._not_truncated``)."""
    protein = "MK" + "A" * 90
    complete = _entry("rps19", protein, "NC_000001.1")
    short = _entry("rps19", protein[:20], "NC_000002.1")  # a prefix, not a suffix
    index = ReferenceIndex.from_entries((complete, short), vocabulary=("rps19",))
    hsps = (
        _one_hsp(complete, qstart=1, qend=279, pident=98.9, bitscore=189.0),
        _one_hsp(short, qstart=1, qend=60, pident=100.0, bitscore=60.0),
    )

    assert (
        call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
        .one("rps19")
        .reference.accession
        == "NC_000002.1"
    )


# ----------------------------------------------- the reference model decides (Crataegus 27/10/2026)


def _model_hsp(entry: ReferenceEntry, *, pident: float, residues: int, bitscore: float):
    """An HSP covering ``residues`` codons of ``entry``, query 1..3*residues."""
    from plastid.core.caller import Hsp

    return Hsp(
        subject_key=entry.key,
        pident=pident,
        length=residues,
        qstart=1,
        qend=3 * residues,
        sstart=1,
        send=residues,
        evalue=1e-30,
        bitscore=bitscore,
        frame=1,
    )


def test_a_minority_reference_model_does_not_win_on_identity() -> None:
    """Measured 27/10/2026 on 39 Crataegus records.

    Three of 39 references carry extended small-CDS models: psbM 270 bp where 36 records say
    105, psbK 195 where 36 say 186, psaJ 318 where 38 say 135, psbZ 279 where 38 say 189,
    psbJ 144 where 36 say 123, accD 1,527 where 36 say 1,476, rpoB 3,258 where 38 say 3,207.
    The caller keeps the highest identity, so those three won at every locus and every genome
    inherited the extension. Removing them from the reference set raised CDS exact agreement
    from 76.7 % to 85.7 % and CDS within 1 bp from 80.8 % to 89.9 % — the model, not the gene
    finding, was the limit.

    Every one of those seven outlier proteins is the consensus protein **plus residues before
    its start**, so the plastid path displaces such a rare model when the consensus model is a
    suffix of it. The mito prototype keeps identity-only.

    This also pins the order: the consensus protein is a *suffix* of the winner, which is what
    ``_not_truncated`` drops, so the guard has to run before the structure rule or it never sees
    the candidate it needs.
    """
    modal = "M" + "A" * 61  # 62 aa = 186 bp, the model 36 of 39 references carry
    extended = "GGG" + modal  # 65 aa = 195 bp, the model 3 of 39 carry: 3 residues before the start
    entries = tuple(_entry("psbK", modal, f"NC_0000{i:02d}.1") for i in range(36)) + tuple(
        _entry("psbK", extended, f"NC_0001{i:02d}.1") for i in range(3)
    )
    index = ReferenceIndex.from_entries(entries, vocabulary=("psbK",))
    modal_entry = next(entry for entry in index.entries if entry.sequence == modal)
    extended_entry = next(entry for entry in index.entries if entry.sequence == extended)
    hsps = (
        _model_hsp(extended_entry, pident=100.0, residues=65, bitscore=300.0),
        _model_hsp(modal_entry, pident=98.4, residues=62, bitscore=290.0),
    )

    plain = core.call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
    guarded = call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)

    assert plain.one("psbK").reference.sequence == extended  # identity picks the 3-of-39 model
    assert plain.one("psbK").span == (1, 195)
    assert guarded.one("psbK").reference.sequence == modal  # the plastid path follows the majority
    assert guarded.one("psbK").span == (1, 186)  # 36 of 39 records say 186 bp, not 195
    # a choice, not a silence: the identity winner is still reported with its own numbers
    assert {row.reference_accession: row.chosen for row in guarded.competing} == {
        modal_entry.accession: True,
        extended_entry.accession: False,
    }


def test_a_three_prime_extension_is_not_an_annotation_outlier() -> None:
    """Measured on the same 39 Crataegus records: five ``rps18`` and two ``atpA`` models are
    three or four residues longer than the consensus at the **3'** end, and they are the real
    models. A plain extension rule shortened those calls by 9 and 12 bp and left them without a
    stop codon (``no_stop``); exact agreement fell by seven calls. Only a 5'-over-long model is
    displaced — the end the start annotation can be stretched past.
    """
    modal = "M" + "A" * 61  # 62 aa
    extended = modal + "GGG"  # 65 aa: residues after the consensus end, not before its start
    entries = tuple(_entry("rps18", modal, f"NC_0000{i:02d}.1") for i in range(36)) + tuple(
        _entry("rps18", extended, f"NC_0001{i:02d}.1") for i in range(3)
    )
    index = ReferenceIndex.from_entries(entries, vocabulary=("rps18",))
    extended_entry = next(entry for entry in index.entries if entry.sequence == extended)
    modal_entry = next(entry for entry in index.entries if entry.sequence == modal)
    hsps = (
        _model_hsp(extended_entry, pident=100.0, residues=65, bitscore=300.0),
        _model_hsp(modal_entry, pident=98.4, residues=62, bitscore=290.0),
    )

    assert (
        call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
        .one("rps18")
        .reference.sequence
        == extended
    )


def test_the_model_guard_does_nothing_without_a_majority_model() -> None:
    """4 + 4 + 1 references over three lengths: the identity winner is a 5' extension of the
    65 aa model, but 4 of 9 references is not a majority, so nothing is transferred."""
    short = "M" + "A" * 61  # 62 aa
    middle = "GGG" + short  # 65 aa
    long = "CCC" + middle  # 68 aa: a 5' extension of `middle`
    entries = (
        tuple(_entry("psbK", short, f"NC_0000{i:02d}.1") for i in range(4))
        + tuple(_entry("psbK", middle, f"NC_0001{i:02d}.1") for i in range(4))
        + (_entry("psbK", long, "NC_000200.1"),)
    )
    index = ReferenceIndex.from_entries(entries, vocabulary=("psbK",))
    chosen = next(entry for entry in index.entries if entry.sequence == long)
    middle_entry = next(entry for entry in index.entries if entry.sequence == middle)
    hsps = (
        _model_hsp(chosen, pident=100.0, residues=68, bitscore=320.0),
        _model_hsp(middle_entry, pident=99.0, residues=65, bitscore=300.0),
    )

    assert (
        call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
        .one("psbK")
        .reference.sequence
        == long
    )


def test_a_majority_model_that_did_not_align_cannot_displace_the_one_that_did() -> None:
    """The guard asks for a majority model **among the passing candidates**. Eight references
    say 186 bp but none of them aligned above the floors; only the 195 bp reference did, and
    the genome is then best explained by what actually aligned."""
    modal = "M" + "A" * 61
    extended = "GGG" + modal
    entries = tuple(_entry("psbK", modal, f"NC_0000{i:02d}.1") for i in range(8)) + (
        _entry("psbK", extended, "NC_000100.1"),
    )
    index = ReferenceIndex.from_entries(entries, vocabulary=("psbK",))
    extended_entry = next(entry for entry in index.entries if entry.sequence == extended)
    hsps = (_model_hsp(extended_entry, pident=100.0, residues=65, bitscore=300.0),)

    assert (
        call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
        .one("psbK")
        .reference.sequence
        == extended
    )


def test_a_rare_model_that_is_not_an_extension_keeps_identitys_choice() -> None:
    """The golden's own ``rpl2`` (measured 27/10/2026) is the negative case for the guard.

    Across the 33-record pinned set, 19 references spell 271 aa from an inner start codon, 11
    spell 273, 2 spell 274 and 4 spell 244; the target's published protein is 273. The 274 model
    is rare (2 of 33), but the modal 271 protein is not the 274 protein with residues removed
    from one end — it is a different start. A plain majority rule displaced the identity winner
    there and pulled the call 6 bp off the published 5' end. A rare model is displaced only by a
    majority model it **extends at the 5' end**.
    """
    modal = "M" + "A" * 61  # 62 aa
    different_start = "M" + "D" + "C" * 60  # 62 aa, and not a suffix of `rare`
    rare = "M" + "C" * 67  # 68 aa — rare, and no model is a 5' extension of it
    entries = (
        tuple(_entry("rpl2", modal, f"NC_0000{i:02d}.1") for i in range(30))
        + (_entry("rpl2", different_start, "NC_000100.1"),)
        + (_entry("rpl2", rare, "NC_000200.1"),)
    )
    index = ReferenceIndex.from_entries(entries, vocabulary=("rpl2",))
    rare_entry = next(entry for entry in index.entries if entry.sequence == rare)
    modal_entry = next(entry for entry in index.entries if entry.sequence == modal)
    other_entry = next(entry for entry in index.entries if entry.sequence == different_start)
    hsps = (
        _model_hsp(rare_entry, pident=100.0, residues=68, bitscore=340.0),
        _model_hsp(modal_entry, pident=98.5, residues=62, bitscore=310.0),
        _model_hsp(other_entry, pident=97.0, residues=62, bitscore=300.0),
    )

    assert (
        call_cds(_genome(), index, hsps, min_identity=50.0, min_coverage=0.5)
        .one("rpl2")
        .reference.sequence
        == rare
    )
