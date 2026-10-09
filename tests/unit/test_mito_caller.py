"""The mito caller (P3A-03) — the rules that decide what a BLAST hit is allowed to become.

Every rule here comes from a measurement on the answer-key record; the tests are where those
measurements stop being folklore:

* BLAST truncates ``sseqid`` at the first whitespace, so a subject is identified by the first
  three ``|`` fields and ``-`` stands for "no vocabulary entry";
* one reference header can stand for two different sequences (measured: Arabidopsis ``atp6``
  385/349 aa, Psychotria ``cox2`` 257/390 aa), so an HSP belongs to the **smallest** reference
  length that contains its subject range — and two sequences under one key are never chained;
* real HSP sets contain redundant alignments (a 27 aa 100 %-identity atp8 fragment at 8,557), a
  second genomic copy of a gene's part, and two frames over one query region (``sdh4``, ``mttB``)
  — each is dropped **with a reason**, never silently;
* the published record's own three ``no_stop`` CDS end flush against the protein, so the stop
  codon is taken only when the sequence confirms it — one codon is read, never a search
  (D3A-02 forbids ORF-walking).
"""

from __future__ import annotations

import pytest

from plastid.core.mito_caller import (
    CallResult,
    GeneCall,
    Hsp,
    ReferenceIndex,
    call_cds,
    call_rrna,
)
from plastid.core.mito_reference import ReferenceEntry
from plastid.core.model import Genome
from plastid.core.sequence import reverse_complement

# A synthetic mitogenome whose every codon is readable at a glance. Positions are 1-based
# inclusive, the way BLAST and GenBank report them.
_PLUS_CDS = "ATG" + "GCT" * 8 + "TAA"  # 1..30: nine codons and a stop, plus strand
_NO_STOP = "GCT" * 9 + "TTA"  # 31..60: nothing a plus-strand call can stop on…
_MINUS_CDS = reverse_complement(_PLUS_CDS)  # 61..90: the same gene on the minus strand
_FILLER = "ATGC" * 40  # 91..250
GENOME = Genome(
    name="synthetic mito",
    accession="SYN_MITO.1",
    sequence=_PLUS_CDS + _NO_STOP + _MINUS_CDS + _FILLER,
    topology="circular",
)
# …except its last codon (58..60), which reads TAA on the minus strand: that is what the
# minus-strand stop extension below must find.
assert reverse_complement(GENOME.sequence[57:60]) == "TAA"


_SAME_AS_GENE = object()  # "no vocabulary override" — None means None, and must be spellable


def _ref(
    gene: str,
    length: int,
    *,
    accession: str = "NC_000001.1",
    vocabulary: str | None | object = _SAME_AS_GENE,
) -> ReferenceEntry:
    """A reference entry of ``length`` residues, all methionine (only its length matters here)."""
    return ReferenceEntry(
        accession=accession,
        organism="Testus testus",
        source_gene=gene,
        vocabulary=gene if vocabulary is _SAME_AS_GENE else vocabulary,  # type: ignore[arg-type]
        sequence="M" * length,
    )


def _index(*entries: ReferenceEntry, vocabulary: tuple[str, ...] | None = None) -> ReferenceIndex:
    return ReferenceIndex.from_entries(
        entries,
        vocabulary=vocabulary or tuple(e.vocabulary for e in entries if e.vocabulary),
    )


def _hsp(
    entry: ReferenceEntry,
    *,
    qstart: int,
    qend: int,
    sstart: int = 1,
    send: int | None = None,
    frame: int = 1,
    pident: float = 95.0,
    length: int | None = None,
    bitscore: float = 100.0,
) -> Hsp:
    """One protein HSP against ``entry``, with BLAST's 1-based inclusive coordinates.

    A minus-frame hit is reported by BLAST with ``qstart > qend``; the helper does that swap so
    every test can state the query interval the way a reader thinks of it.
    """
    send = send if send is not None else len(entry.sequence)
    if frame < 0:
        qstart, qend = max(qstart, qend), min(qstart, qend)
    return Hsp(
        subject_key=entry.key,
        pident=pident,
        length=length if length is not None else abs(qend - qstart) // 3 + 1,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        evalue=1e-30,
        bitscore=bitscore,
        frame=frame,
    )


def _rrna_hsp(entry: ReferenceEntry, *, qstart: int, qend: int, sstart: int, send: int) -> Hsp:
    """A nucleotide HSP: no frame, so its strand has to come from the subject's direction."""
    return Hsp(
        subject_key=entry.key,
        pident=97.0,
        length=abs(qend - qstart) + 1,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        evalue=0.0,
        bitscore=500.0,
        frame=None,
    )


# --------------------------------------------------------------------------- one HSP, one gene


def test_a_plus_frame_hsp_becomes_a_cds_on_the_plus_strand() -> None:
    """The simplest transfer: the alignment spans the reference, the query interval is the call."""
    entry = _ref("ccmB", 10)
    result = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=1, qend=30),))

    call = result.one("ccmB")
    assert call.kind == "CDS"
    assert [(p.to_1based(), p.strand) for p in call.parts] == [((1, 30), 1)]


def test_a_minus_frame_hsp_becomes_a_complement_interval() -> None:
    """BLAST reports a minus-frame hit with qstart > qend. Getting this backwards puts the whole
    gene on the wrong strand while every count stays right — this project's definition of a
    silent wrong.

    The 3' end is 58, one codon below the alignment's own end: this fixture's minus strand reads
    TAA at 58–60, and the stop codon is read whether or not the alignment reached the reference's
    end (``test_a_stop_after_an_alignment_that_stopped_short_is_still_the_end``). The strand and
    the orientation are what this test is about.
    """
    entry = _ref("cox3", 10)
    hsps = (_hsp(entry, qstart=90, qend=61, frame=-2, send=9),)
    part = call_cds(GENOME, _index(entry), hsps).one("cox3").parts[0]

    assert (part.to_1based(), part.strand) == ((58, 90), -1)


def test_a_first_codon_that_is_not_an_initiator_is_a_partial_start_never_a_shift() -> None:
    """Measured on the answer key: 3 of its 31 CDS start with a codon no table calls an initiator
    (rps4 ACG, mttB CGC, rps10 ATA). The coordinate is reported as incomplete; it is never moved
    to the next ATG, which would invent a record NCBI did not publish."""
    entry = _ref("rps4", 10)
    call = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=31, qend=60, send=10),)).one("rps4")

    assert call.partial_start is True
    assert call.parts[0].to_1based() == (31, 60)  # unmoved
    assert any("initiator" in note for note in call.notes)


# --------------------------------------------------------------------------- several exons


def test_exons_are_ordered_by_the_reference_not_by_position() -> None:
    """A trans-spliced gene's exons sit in an order the protein knows and the genome does not.
    Ordering by the alignment is evidence; ordering by coordinate is a guess that looks tidy."""
    entry = _ref("nad5", 40)
    hsps = (
        _hsp(entry, qstart=200, qend=229, sstart=30, send=40, bitscore=50),
        _hsp(entry, qstart=1, qend=30, sstart=1, send=10, bitscore=90),
        _hsp(entry, qstart=100, qend=139, sstart=11, send=29, bitscore=70),
    )
    call = call_cds(GENOME, _index(entry), hsps).one("nad5")

    assert [p.to_1based() for p in call.parts] == [(1, 30), (100, 139), (200, 229)]


def test_mixed_strand_exons_keep_their_own_strand() -> None:
    """Measured: nad1 and nad5 align on both strands of one mitogenome. Each exon carries the
    strand it aligned on, and the gene reports mixed rather than picking one."""
    entry = _ref("nad1", 30)
    hsps = (
        _hsp(entry, qstart=61, qend=90, sstart=1, send=10, frame=-1, bitscore=90),
        _hsp(entry, qstart=1, qend=30, sstart=11, send=30, frame=1, bitscore=80),
    )
    call = call_cds(GENOME, _index(entry), hsps).one("nad1")

    assert [(p.to_1based(), p.strand) for p in call.parts] == [((61, 90), -1), ((1, 30), 1)]
    assert call.to_gene().strand == 0


# --------------------------------------------------------------------------- what gets dropped


def test_a_redundant_hsp_adds_no_exon_and_is_reported() -> None:
    """Measured: a 27 aa 100 %-identity atp8 fragment at 8,557 against a full-length hit. It adds
    no residue the kept alignment does not cover — dropping it in silence is how a gene ends up
    with a second, invented exon."""
    entry = _ref("atp8", 20)
    kept = _hsp(entry, qstart=1, qend=60, sstart=1, send=20, bitscore=200)
    redundant = _hsp(entry, qstart=200, qend=220, sstart=1, send=7, bitscore=60)
    result = call_cds(GENOME, _index(entry), (kept, redundant))

    assert [p.to_1based() for p in result.one("atp8").parts] == [(1, 60)]
    assert [row.hsp for row in result.dropped] == [redundant]
    assert "second copy" in result.dropped[0].reason


def test_a_redundant_hsp_far_from_the_call_is_a_second_copy_not_a_silent_drop() -> None:
    """The same rule, stated as the finding it is: the fragment at 200-220 is a copy of the
    gene's 5' part sitting somewhere else in the genome."""
    entry = _ref("nad3", 20)
    kept = _hsp(entry, qstart=1, qend=60, sstart=1, send=20, bitscore=200)
    copy = _hsp(entry, qstart=150, qend=170, sstart=5, send=11, bitscore=70)
    result = call_cds(GENOME, _index(entry), (kept, copy))

    assert copy.qlo >= 150
    assert any("copy" in row.reason for row in result.dropped)


def test_an_hsp_in_another_frame_over_the_same_region_is_a_conflict_not_an_exon() -> None:
    """Measured on sdh4: one query region aligns at 50 % in one frame and 75 % in another. The
    better alignment is the call; the other is reported, because "two frames over one region" is
    something a reader has to see."""
    entry = _ref("sdh4", 20)
    kept = _hsp(entry, qstart=20, qend=79, frame=-1, sstart=1, send=20, bitscore=200)
    conflict = _hsp(entry, qstart=30, qend=59, frame=-3, sstart=5, send=15, bitscore=120)
    result = call_cds(GENOME, _index(entry), (kept, conflict))

    assert [p.to_1based() for p in result.one("sdh4").parts] == [(20, 79)]
    assert [row.hsp for row in result.dropped] == [conflict]
    assert "frame" in result.dropped[0].reason


def test_every_vocabulary_locus_appears_in_called_or_not_called() -> None:
    """Nothing may vanish: the tables account for the whole vocabulary."""
    found, absent = _ref("ccmB", 10), _ref("cox1", 10, accession="NC_000007.1")
    result = call_cds(GENOME, _index(found, absent), (_hsp(found, qstart=1, qend=30),))

    assert {call.symbol for call in result.calls} == {"ccmB"}
    assert [row.symbol for row in result.not_called] == ["cox1"]
    assert "no alignment" in result.not_called[0].reason


# --------------------------------------------------------------------------- floors


def test_a_hit_below_the_identity_floor_is_not_called_and_says_why() -> None:
    entry = _ref("nad7", 20)
    hsps = (_hsp(entry, qstart=1, qend=60, send=20, pident=41.3),)
    result = call_cds(GENOME, _index(entry), hsps, min_identity=50.0)

    assert result.calls == ()
    row = result.not_called[0]
    assert (row.symbol, row.best_identity) == ("nad7", pytest.approx(41.3))
    assert [dropped.hsp for dropped in result.dropped] == list(hsps)


def test_the_floors_are_applied_per_reference_and_a_passing_runner_up_wins() -> None:
    """Measured on atp6: the highest-identity reference can fail the coverage floor while another
    passes. Applying the floors to the winner would drop a locus that has a good call."""
    better = _ref("atp6", 10, accession="NC_000008.1")  # 93 %, full coverage
    closest = _ref("atp6", 100, accession="NC_000009.1")  # 99 %, 10 % coverage
    hsps = (
        _hsp(better, qstart=1, qend=30, send=10, pident=93.0, bitscore=200),
        _hsp(closest, qstart=1, qend=30, send=10, pident=99.0, bitscore=210),
    )
    result = call_cds(GENOME, _index(better, closest), hsps)

    assert result.one("atp6").reference.accession == "NC_000008.1"
    assert {row.reference_accession: row.chosen for row in result.competing} == {
        "NC_000008.1": True,
        "NC_000009.1": False,
    }
    assert [row.passes_floors for row in result.competing] == [True, False]


def test_a_locus_whose_best_reference_also_fails_reports_its_numbers() -> None:
    entry = _ref("rps1", 100)
    result = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=1, qend=30, send=10),))

    assert result.calls == ()
    row = result.not_called[0]
    assert row.best_coverage == pytest.approx(0.1)
    assert "coverage" in row.reason


# --------------------------------------------------------------------------- reference indexing


def test_an_hsp_belongs_to_the_smallest_reference_length_that_contains_it() -> None:
    """Measured: two sequences under one header key (Arabidopsis atp6 385/349 aa). An alignment
    ending at subject position 380 can only have come from the 385 aa copy — the 349 aa one must
    not claim it, which is what keeps a truncated reference out of a call's denominator."""
    short = _ref("atp6", 349, accession="NC_000002.1")
    long = _ref("atp6", 385, accession="NC_000002.1")
    assert short.key == long.key, "the two copies share one BLAST subject key"

    hsps = (_hsp(long, qstart=1, qend=300, sstart=1, send=380, length=100),)
    result = call_cds(GENOME, _index(short, long), hsps)

    assert [row.reference_length for row in result.competing] == [385]


def test_two_sequences_under_one_key_are_never_chained_into_a_chimera() -> None:
    short = _ref("cox2", 10, accession="NC_000003.1")
    long = _ref("cox2", 20, accession="NC_000003.1")
    hsps = (
        _hsp(short, qstart=1, qend=30, sstart=1, send=10, pident=99.0, bitscore=200),
        _hsp(long, qstart=100, qend=159, sstart=11, send=20, pident=90.0, bitscore=210),
    )
    call = call_cds(GENOME, _index(short, long), hsps).one("cox2")

    assert call.reference.sequence == short.sequence
    assert [p.to_1based() for p in call.parts] == [(1, 30)], "the other sequence's exon joined it"


def test_an_hit_to_a_gene_the_vocabulary_cannot_name_is_reported_not_called() -> None:
    """H-06: the reference set carries ten standard genes the vocabulary lacks (rps12, rpl2 …). A
    run must be able to *say* it found one — renaming it here would change naming for every mito
    record as a side effect of an annotator task."""
    named = _ref("ccmB", 10, vocabulary="ccmB")
    unnamed = _ref("rps12", 10, accession="NC_000004.1", vocabulary=None)
    hsps = (
        _hsp(named, qstart=1, qend=30, send=10),
        _hsp(unnamed, qstart=31, qend=60, send=10),
    )
    result = call_cds(GENOME, _index(named, unnamed, vocabulary=("ccmB",)), hsps)

    assert {call.symbol for call in result.calls} == {"ccmB"}
    assert [row.symbol for row in result.unnameable] == ["rps12"]
    assert result.unnameable[0].parts[0].to_1based() == (31, 60)


def test_the_subject_key_of_an_unnameable_reference_is_the_one_blast_prints() -> None:
    """The header's vocabulary field is ``-`` when the toolkit cannot name the gene, and BLAST
    prints that field like any other: a key that used ``None`` would never match a hit."""
    entry = _ref("rpl16", 10, vocabulary=None)
    assert entry.key == "NC_000001.1|rpl16|-"


# --------------------------------------------------------------------------- the stop codon


def test_the_stop_codon_is_taken_when_the_sequence_confirms_it() -> None:
    """The reference protein has no stop, so an alignment reaching its C-terminus is extended by
    exactly one codon — and only when that codon is a stop. Measured: ccmFc and nad6 are
    annotated in the record *without* their stop, which sits immediately after."""
    entry = _ref("nad6", 9)
    call = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=1, qend=27, send=9),)).one("nad6")

    assert call.parts[0].to_1based() == (1, 30)
    assert call.partial_end is False


def test_a_next_codon_that_is_not_a_stop_is_not_searched_for() -> None:
    """D3A-02: boundaries come from the alignment. One codon is read; if it is not a stop the end
    is left where the alignment put it and reported incomplete — walking on to the next stop is
    ORF-walking, and it would invent an end the evidence does not show."""
    entry = _ref("cox2", 9)
    call = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=31, qend=57, send=9),)).one("cox2")

    assert call.parts[0].to_1based() == (31, 57), "the call walked past the alignment's end"
    assert call.partial_end is True
    assert any("not a stop" in note for note in call.notes)


def test_the_stop_extension_moves_the_transcript_end_on_the_minus_strand() -> None:
    """On the minus strand the 3' end is the *lower* coordinate; extending the wrong end is a
    silent three-base error in every downstream position-based analysis."""
    entry = _ref("cox3", 10)
    hsps = (_hsp(entry, qstart=90, qend=61, frame=-2, send=10),)
    call = call_cds(GENOME, _index(entry), hsps).one("cox3")

    assert call.parts[0].to_1based() == (58, 90)
    assert call.partial_end is False


def test_a_stop_after_an_alignment_that_stopped_short_is_still_the_end() -> None:
    """Measured (P3B-10, §8.55): the stop codon is read whether or not the alignment reached the
    reference's 3' end. A target whose gene is shorter than the reference's presents an alignment
    that simply ends — blastx stops translating at a stop and cannot align past one — and
    requiring the alignment to reach the reference's end left every such call one codon short of
    the published CDS. 65 of the corpus's remaining differences were exactly this, ``psaJ`` (30
    records) at the top of them."""
    entry = _ref("cox3", 10)
    hsps = (_hsp(entry, qstart=90, qend=61, frame=-2, send=9),)
    call = call_cds(GENOME, _index(entry), hsps).one("cox3")

    assert call.parts[0].to_1based() == (58, 90)
    assert call.partial_end is False
    assert any("one codon shorter than the reference" in note for note in call.notes)


def test_an_alignment_that_stops_short_of_the_reference_end_is_partial() -> None:
    entry = _ref("rps3", 10)
    call = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=1, qend=30, send=9),)).one("rps3")

    assert (call.partial_start, call.partial_end) == (False, True)


def test_an_extension_that_would_cross_the_record_end_is_refused_and_reported() -> None:
    """A circular record wraps, so extending past its end would quietly read the *other* end of
    the genome — and on a linear record it would read past the data. Neither is an end the
    evidence shows, so the call keeps the alignment's end and the run says why."""
    entry = _ref("nad9", 10)
    hsps = (_hsp(entry, qstart=len(GENOME) - 29, qend=len(GENOME), send=10),)
    result = call_cds(GENOME, _index(entry), hsps)
    call = result.one("nad9")

    assert call.parts[0].to_1based() == (len(GENOME) - 29, len(GENOME))
    assert call.partial_end is True
    assert any("end of the record" in note for note in call.notes)


# --------------------------------------------------------------------------- the call as data


def test_a_call_becomes_a_gene_with_the_translation_table_it_was_given() -> None:
    entry = _ref("nad9", 10)
    call = call_cds(GENOME, _index(entry), (_hsp(entry, qstart=1, qend=30),), transl_table=1)
    gene = call.one("nad9").to_gene()

    assert (gene.name, gene.kind, gene.transl_table) == ("nad9", "CDS", 1)
    assert gene.product is None, "product text is the submitter's; the vocabulary is not a source"


def test_coverage_is_the_merged_exon_bases_over_three_times_the_reference_length() -> None:
    entry = _ref("ccmB", 20)
    call = call_cds(
        GENOME, _index(entry), (_hsp(entry, qstart=1, qend=30, send=10),), min_coverage=0.4
    )
    assert call.one("ccmB").coverage == pytest.approx(30 / (3 * 20))


def test_identity_is_weighted_by_the_alignment_length() -> None:
    entry = _ref("nad2", 30)
    hsps = (
        _hsp(entry, qstart=1, qend=30, sstart=1, send=10, length=10, pident=100.0, bitscore=90),
        _hsp(entry, qstart=100, qend=159, sstart=11, send=30, length=20, pident=70.0, bitscore=80),
    )
    call = call_cds(GENOME, _index(entry), hsps).one("nad2")

    assert call.identity == pytest.approx((100.0 * 10 + 70.0 * 20) / 30)


def test_the_result_records_the_floors_it_ran_with() -> None:
    entry = _ref("ccmB", 10)
    result = call_cds(
        GENOME,
        _index(entry),
        (_hsp(entry, qstart=1, qend=30),),
        min_identity=45.0,
        min_coverage=0.4,
    )

    assert isinstance(result, CallResult)
    assert isinstance(result.calls[0], GeneCall)
    assert (result.min_identity, result.min_coverage) == (45.0, 0.4)


# --------------------------------------------------------------------------- rRNA


def test_the_rrna_strand_comes_from_the_subject_not_the_query() -> None:
    """The reference rRNA is stored 5'→3', so ``sstart > send`` means the gene lies on the query's
    minus strand — the inverse of the blastx convention. Measured against the answer key: rrn18
    and rrn5 are minus, rrn26 is plus, and all three spans agree to within 4 bp."""
    plus = _ref("rrn26", 30)
    minus = _ref("rrn18", 30, accession="NC_000005.1")
    hsps = (
        _rrna_hsp(plus, qstart=1, qend=30, sstart=1, send=30),
        _rrna_hsp(minus, qstart=61, qend=90, sstart=30, send=1),
    )
    result = call_rrna(GENOME, _index(plus, minus, vocabulary=("rrn26", "rrn18")), hsps)

    assert (result.one("rrn26").parts[0].to_1based(), result.one("rrn26").parts[0].strand) == (
        (1, 30),
        1,
    )
    assert (result.one("rrn18").parts[0].to_1based(), result.one("rrn18").parts[0].strand) == (
        (61, 90),
        -1,
    )


def test_an_rrna_call_is_not_refined_by_codons() -> None:
    """rRNA has no reading frame: no stop extension, no initiator test, and no codon multiple."""
    entry = _ref("rrn5", 100)
    call = call_rrna(
        GENOME,
        _index(entry, vocabulary=("rrn5",)),
        (_rrna_hsp(entry, qstart=31, qend=120, sstart=1, send=90),),
    ).one("rrn5")

    assert call.parts[0].to_1based() == (31, 120)
    assert call.kind == "rRNA"
    assert (call.partial_start, call.partial_end) == (False, True)


def test_an_rrna_alignment_that_misses_the_reference_start_is_a_partial_start() -> None:
    entry = _ref("rrn26", 100)
    call = call_rrna(
        GENOME,
        _index(entry, vocabulary=("rrn26",)),
        (_rrna_hsp(entry, qstart=1, qend=99, sstart=11, send=99),),
    ).one("rrn26")

    assert (call.partial_start, call.partial_end) == (True, True)


# --------------------------------------------------------------------------- findings


def test_two_calls_that_overlap_each_other_are_a_finding() -> None:
    """A repeat-rich mitogenome can hand two loci the same region. Both calls stand — each is the
    best evidence for its own gene — and the overlap is reported, because a silent overlap is how
    a chimeric annotation gets published."""
    left = _ref("ccmC", 10)
    right = _ref("ccmB", 10, accession="NC_000006.1")
    hsps = (_hsp(left, qstart=1, qend=30), _hsp(right, qstart=20, qend=50))
    result = call_cds(GENOME, _index(left, right), hsps)

    assert len(result.calls) == 2
    assert any("overlap" in finding for finding in result.findings)


def test_an_hsp_set_that_saturates_the_cap_is_reported() -> None:
    """-max_hsps truncates silently: a gene whose exon set was cut off looks like a shorter gene,
    and its coverage quietly drops. A pair that returns exactly the cap is reported."""
    entry = _ref("nad4", 200)
    hsps = tuple(
        _hsp(
            entry,
            qstart=1 + 10 * i,
            qend=10 + 10 * i,
            sstart=1 + 10 * i,
            send=10 + 10 * i,
            bitscore=100 - i,
        )
        for i in range(5)
    )
    result = call_cds(GENOME, _index(entry), hsps, max_hsps_per_subject=5)

    assert any("cap" in finding for finding in result.findings)


def test_a_call_whose_exons_meet_at_the_record_ends_is_flagged_as_possibly_circular() -> None:
    """A gene crossing the origin of a circular record arrives as two HSPs, one at each end.
    Reading them as two exons is a guess; the run reports the possibility instead."""
    entry = _ref("nad5", 20)
    hsps = (
        _hsp(entry, qstart=1, qend=30, sstart=1, send=10, bitscore=90),
        _hsp(entry, qstart=len(GENOME) - 29, qend=len(GENOME), sstart=11, send=20, bitscore=80),
    )
    result = call_cds(GENOME, _index(entry), hsps)

    assert any("origin" in finding for finding in result.findings)
