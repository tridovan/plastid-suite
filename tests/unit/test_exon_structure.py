"""The reference's own exon structure, placed on the target — the rule, on synthetic DNA.

Every case here is one the corpus measured (P3B-10): an exon the protein search cannot
seed on (`petB`'s 6 bp), a boundary an alignment ran straight through (`rpoC1`), a call exon the
transfer must **not** lose (`rps12`'s, across the trans-splice), a record whose exons are not
colinear at all (`NC_072562.1`), and a hit that lands somewhere else (the inverted repeat's other
copy). The alignments are fakes — what is under test is the arithmetic and the refusals, not BLAST.
"""

from __future__ import annotations

import pytest

from plastid.core.caller import (
    NOT_REACHED_FIVE_PRIME,
    CallResult,
    GeneCall,
    ReferenceEntry,
    transfer_structure,
)
from plastid.core.exon_structure import (
    MAX_INTRON,
    Junction,
    Transfer,
    junctions,
    place,
    render_junctions,
)
from plastid.core.model import Gene, Genome, Interval


def gene(name: str, *spans: tuple[int, int], strand: int = 1) -> Gene:
    """A feature from 1-based inclusive exons, in transcript order, the way a record states one."""
    return Gene(name=name, kind="CDS", parts=tuple(Interval(s - 1, e, strand) for s, e in spans))


def genome(sequence: str) -> Genome:
    return Genome(name="synthetic", sequence=sequence, topology="circular")


class Hit:
    """An alignment the way ``TraceHit`` answers: junction offset → target position.

    ``pairs`` is the walk itself — ``(junction position, target position)`` for every aligned
    column, in order — so a test can state a gap exactly where it wants one.
    """

    def __init__(self, pairs: list[tuple[int, int]], bitscore: float = 100.0) -> None:
        self.bitscore = bitscore
        self._table = dict(pairs)

    def map_subject(self, position: int) -> int | None:
        return self._table.get(position)


# --------------------------------------------------------------------------------- the structure


def test_a_boundary_is_read_on_the_strand_its_own_two_exons_share() -> None:
    """Per pair, never per gene: a trans-spliced gene's exons may lie on both strands, and the
    two exons of one junction are still colinear (measured: NC_072562.1's ``rps12``)."""
    feature = Gene(
        name="rps12",
        kind="CDS",
        parts=(Interval(99, 113, -1), Interval(1099, 1128, 1), Interval(1300, 1325, 1)),
    )

    (only,) = junctions(genome("ACGT" * 400), feature, "x|rps12|j")

    assert only.index == 1  # the pair across the trans-splice has no junction; this one does
    assert (only.left, only.right, only.intron) == ((1, 29), (202, 226), 172)
    assert len(only.sequence) == 226


def test_one_junction_per_internal_boundary_carries_both_exons_and_the_intron() -> None:
    """The junction is exon k, the whole intron, exon k+1 — the three pieces the search needs."""
    feature = gene("clpP", (1, 71), (100, 391), (500, 751))
    (first, second) = junctions(genome("ACGT" * 200), feature, "NC_1|clpP|j")

    assert first.key == "NC_1|clpP|j0"
    assert (first.left, first.boundary, first.right) == ((1, 71), 72, (100, 391))
    assert first.intron == 28  # 72..99
    assert (first.index, first.exons, first.first, first.last) == (0, 3, True, False)
    assert (second.index, second.first, second.last) == (1, False, True)


def test_a_gap_that_is_not_an_intron_gets_no_junction() -> None:
    """``rps12``'s 69 kb trans-splice: the DNA between those exons is not an intron to search."""
    feature = gene("rps12", (1, 114), (200_000, 200_228), strand=-1)

    assert junctions(genome("ACGT" * 60_000), feature, "x|rps12|j") == ()
    assert MAX_INTRON < 69_056  # the measured gap this refuses, with room to spare


def test_a_single_exon_gene_has_no_boundary_to_place() -> None:
    assert junctions(genome("ACGT" * 10), gene("psbA", (1, 30)), "x|psbA|j") == ()


def test_junctions_render_as_a_fasta_a_hit_can_be_read_back_from() -> None:
    """The id is the key the search returns, and the subject key rule keeps three fields."""
    (only,) = junctions(genome("ACGT" * 200), gene("petD", (1, 8), (100, 111)), "NC_9|petD|j")

    text = render_junctions([only])

    lines = text.splitlines()
    assert lines[0] == ">NC_9|petD|j0"
    assert "".join(lines[1:]) == only.sequence  # wrapped, as every FASTA writer here wraps


# ------------------------------------------------------------------------------------ placing


def mapped(only: Junction, left_at: int, right_at: int) -> Hit:
    """A hit that walks a junction onto the target: exon k lands at ``left_at``, exon k+1 after it.

    Stated from the junction's own offsets, so a test says where each exon goes and nothing else.
    """
    pairs = [
        (offset, left_at + offset - only.left[0])
        for offset in range(only.left[0], only.left[1] + 1)
    ]
    intron_at = left_at + only.left[1] - only.left[0] + 1
    pairs += [
        (offset, intron_at + offset - only.boundary)
        for offset in range(only.boundary, only.right[0])
    ]
    pairs += [
        (offset, right_at + offset - only.right[0])
        for offset in range(only.right[0], only.right[1] + 1)
    ]
    return Hit(pairs)


def test_a_first_exon_the_protein_search_cannot_see_is_placed_in_front_of_the_call() -> None:
    """`petB` on the golden: a 6 bp exon, and an alignment that simply started after it.

    The junction walks exon (6 bp) → intron → exon, and the call's single exon is the second one;
    the placed 6 bp exon goes *in front* of it, which is where the published record has it.
    """
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="petB", kind="CDS", parts=(Interval(0, 6, 1), Interval(694, 1000, 1))),
        "x|petB|j",
    )
    call = Transfer(
        symbol="petB", strand=1, spans=((74985, 75632),), partial_start=True, partial_end=False
    )
    hit = mapped(only, left_at=74233, right_at=74985)

    placed = place(call, [only], {only.key: [hit]})

    assert placed is not None
    assert placed.spans == ((74233, 74238), (74985, 75632))
    assert (placed.partial_start, placed.partial_end) == (False, False)


def test_a_boundary_inside_one_alignment_is_moved_to_where_the_reference_puts_it() -> None:
    """`rpoC1` measured 15 and 12 bp away: the protein alignment ran into the intron and stopped."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rpoC1", kind="CDS", parts=(Interval(0, 100, 1), Interval(300, 400, 1))),
        "x|rpoC1|j",
    )
    call = Transfer(
        symbol="rpoC1",
        strand=1,
        spans=((5000, 5105), (5300, 5400)),
        partial_start=False,
        partial_end=False,
    )
    hit = mapped(only, left_at=4999, right_at=5300)

    placed = place(call, [only], {only.key: [hit]})

    assert placed is not None
    assert placed.spans == ((5000, 5098), (5300, 5400))


def test_the_calls_own_end_is_kept_where_the_call_established_it() -> None:
    """Its initiator and stop tests are evidence about *this* genome; the placement only moves
    what the alignment could not see. An established start stays where the caller put it."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rpoC1", kind="CDS", parts=(Interval(0, 100, 1), Interval(300, 400, 1))),
        "x|rpoC1|j",
    )
    call = Transfer(
        symbol="rpoC1",
        strand=1,
        spans=((5000, 5105), (5300, 5400)),
        partial_start=False,
        partial_end=False,
    )
    hit = mapped(only, left_at=4999, right_at=5300)

    placed = place(call, [only], {only.key: [hit]})

    assert placed is not None
    assert placed.spans[0][0] == 5000  # not 4999: the caller's start codon stands


def test_a_call_exon_across_the_trans_splice_survives_the_transfer() -> None:
    """`rps12`'s first exon has no junction (69 kb away) — dropping it would delete a real exon."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rps12", kind="CDS", parts=(Interval(0, 229, 1), Interval(300, 326, 1))),
        "x|rps12|j",
    )
    call = Transfer(
        symbol="rps12",
        strand=1,
        spans=((100, 213), (5000, 5229), (6000, 6025)),
        partial_start=False,
        partial_end=True,
    )
    hit = mapped(only, left_at=5000, right_at=6000)

    placed = place(call, [only], {only.key: [hit]})

    assert placed is not None
    assert placed.spans == ((100, 213), (5000, 5228), (6000, 6025))


def test_a_hit_that_lands_somewhere_else_is_rejected() -> None:
    """The inverted repeat's other copy, or a paralog: a junction is only used when its exons
    land on the call. Measured — without this, `NC_072562.1`'s misoriented junction anchored on
    26 shared bases and would have moved the golden's `rps12` 700 bp."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rpl2", kind="CDS", parts=(Interval(0, 100, 1), Interval(300, 400, 1))),
        "x|rpl2|j",
    )
    call = Transfer(
        symbol="rpl2",
        strand=1,
        spans=((5000, 5100), (5300, 5400)),
        partial_start=False,
        partial_end=False,
    )
    elsewhere = mapped(only, left_at=20_000, right_at=21_000)

    assert place(call, [only], {only.key: [elsewhere]}) is None


def test_a_junction_that_did_not_align_keeps_the_alignment_boundary_and_says_so() -> None:
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="petD", kind="CDS", parts=(Interval(0, 8, 1), Interval(300, 400, 1))),
        "x|petD|j",
    )
    call = Transfer(
        symbol="petD",
        strand=1,
        spans=((5000, 5008), (5300, 5400)),
        partial_start=False,
        partial_end=False,
    )

    placed = place(call, [only], {})

    assert placed is None  # nothing placed at all: the call stands exactly as the caller made it


def test_an_exon_placed_on_the_other_side_of_the_intron_is_refused() -> None:
    """The two exons must sit in gene order, an intron apart: an alignment that says otherwise
    is not describing these two exons."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="petB", kind="CDS", parts=(Interval(0, 6, 1), Interval(300, 360, 1))),
        "x|petB|j",
    )
    call = Transfer(
        symbol="petB", strand=1, spans=((5000, 5360),), partial_start=True, partial_end=False
    )
    overlapping = Hit([(offset, 5000 + offset) for offset in range(1, 61)])  # both exons overlap

    assert place(call, [only], {only.key: [overlapping]}) is None


def test_a_minus_strand_gene_is_placed_in_its_own_direction() -> None:
    """The junction is written 5'→3' of the gene, so on the genome its exons run backwards:
    exon k lands on the *higher* coordinates and exon k+1 on the lower ones."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rpl16", kind="CDS", parts=(Interval(391, 400, -1), Interval(0, 100, -1))),
        "x|rpl16|j",
    )
    call = Transfer(
        symbol="rpl16", strand=-1, spans=((80539, 80940),), partial_start=True, partial_end=False
    )
    hit = mapped(only, left_at=81865, right_at=80539)

    placed = place(call, [only], {only.key: [hit]})

    assert placed is not None
    assert placed.spans == ((81865, 81873), (80539, 80638))


def _split_hits(only: Junction) -> tuple[Hit, Hit]:
    """The two alignments the search returned for MW201730's rps16 junction, in miniature.

    One starts at the first exon and runs into the intron, the other starts inside the intron and
    runs through the second exon; neither carries both exons.
    """
    left = Hit([(offset, 7416 + offset - only.left[0]) for offset in range(only.left[0], 101)])
    right = Hit(
        [(offset, 5250 + offset - only.right[0]) for offset in range(900, only.right[1] + 1)]
    )
    return left, right


def test_a_junction_the_search_splits_into_two_alignments_still_places_both_exons() -> None:
    """Measured 27/10/2026 on ``MW201730``'s ``rps16``: the 1,134 bp junction came back as two
    alignments (766 bp + 421 bp) because the intron's own length diverged, and each one carried a
    single exon. Requiring one alignment for both deleted the gene's 40 bp first exon: the call
    was written as the 3' exon alone, partial 5'.

    Both exons are still mapped end to end, each by an alignment of the junction's own DNA, so
    the reference's structure is evidence about both. The pair is used only where it orders and
    anchors exactly like a single hit.
    """
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rps16", kind="CDS", parts=(Interval(6399, 6439, -1), Interval(5249, 5477, -1))),
        "x|rps16|j",
    )
    call = Transfer(
        symbol="rps16", strand=-1, spans=((5250, 5477),), partial_start=True, partial_end=False
    )
    left, right = _split_hits(only)

    placed = place(call, [only], {only.key: [left, right]})

    assert placed is not None
    assert placed.spans == ((7416, 7455), (5250, 5477))
    assert (placed.partial_start, placed.partial_end) == (False, False)


def test_two_alignments_that_do_not_sit_in_gene_order_place_nothing() -> None:
    """The pair rule is the single-hit rule: the second exon still has to sit an intron away, in
    the gene's direction, or the two alignments are not describing this gene."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="rps16", kind="CDS", parts=(Interval(6399, 6439, -1), Interval(5249, 5477, -1))),
        "x|rps16|j",
    )
    call = Transfer(
        symbol="rps16", strand=-1, spans=((5250, 5477),), partial_start=True, partial_end=False
    )
    left, right = _split_hits(only)
    above = Hit(
        [(offset, 7500 + offset - only.right[0]) for offset in range(900, only.right[1] + 1)]
    )

    assert place(call, [only], {only.key: [left, above]}) is None
    assert right  # the ordered pair is the one the first test used


@pytest.mark.parametrize("offset", [0, 2_000])
def test_a_position_outside_the_alignment_has_no_place_on_the_target(offset: int) -> None:
    from plastid.external.blast import TraceHit

    hit = TraceHit(
        subject_key="x|g|j0",
        pident=100.0,
        length=3,
        qstart=100,
        qend=102,
        sstart=1,
        send=3,
        evalue=0.0,
        bitscore=50.0,
        qseq="ACG",
        sseq="ACG",
    )

    assert hit.map_subject(offset) is None


def test_a_gap_in_the_alignment_is_not_a_position() -> None:
    """An insertion in the subject has no base on the target, and must not be given one."""
    from plastid.external.blast import TraceHit

    hit = TraceHit(
        subject_key="x|g|j0",
        pident=100.0,
        length=4,
        qstart=10,
        qend=12,
        sstart=1,
        send=4,
        evalue=0.0,
        bitscore=50.0,
        qseq="AC-G",
        sseq="ACGG",
    )

    assert [hit.map_subject(position) for position in (1, 2, 3, 4)] == [10, 11, None, 12]


def test_a_traced_search_needs_both_halves_of_the_alignment() -> None:
    from plastid.external.blast import parse_trace_hits

    line = "seq\tacc|g|j0\t99.0\t3\t0\t0\t1\t3\t1\t3\t1e-9\t40\tACG\tACG"

    assert len(parse_trace_hits(line, seq_id="seq")) == 1
    with pytest.raises(ValueError, match="alignment whose halves do not line up"):
        parse_trace_hits(line.replace("\tACG\tACG", "\tACG\tAC"), seq_id="seq")


# ------------------------------------------------------------------------------ rebuilding calls


def reference(name: str, *spans: tuple[int, int], junctions_of_it=()) -> ReferenceEntry:
    """A reference entry carrying a record's own exon structure, the way the index builds one."""
    return ReferenceEntry(
        accession="NC_1",
        organism="Iris synthetica",
        source_gene=name,
        vocabulary=name,
        sequence="M" * 100,
        junctions=junctions_of_it,
    )


def called(entry: ReferenceEntry, *spans: tuple[int, int], notes: tuple[str, ...] = ()) -> GeneCall:
    return GeneCall(
        symbol=entry.vocabulary or entry.source_gene,
        kind="CDS",
        parts=tuple(Interval(start - 1, end, 1) for start, end in spans),
        coverage=1.0,
        identity=100.0,
        bitscore=100.0,
        reference=entry,
        transl_table=11,
        notes=notes,
    )


def test_a_call_with_nothing_to_place_is_returned_untouched() -> None:
    """A single-exon gene, an rRNA, a trans-spliced one: this is a transfer of what the reference
    states, not a second opinion about every call."""
    entry = reference("psbA")
    result = CallResult(calls=(called(entry, (100, 400)),))

    assert transfer_structure(result, {}) is result


def test_a_call_whose_junction_did_not_align_keeps_the_callers_own_boundaries() -> None:
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="petD", kind="CDS", parts=(Interval(0, 8, 1), Interval(300, 400, 1))),
        "NC_1|petD|j",
    )
    entry = reference("petD", junctions_of_it=(only,))
    result = CallResult(calls=(called(entry, (5000, 5400)),))

    assert transfer_structure(result, {}) is result


def test_a_placed_start_takes_the_note_that_said_the_alignment_missed_it_with_it() -> None:
    """A table that kept saying "partial" beside an established start would contradict itself."""
    (only,) = junctions(
        genome("A" * 50),
        Gene(name="petB", kind="CDS", parts=(Interval(0, 6, 1), Interval(300, 360, 1))),
        "NC_1|petB|j",
    )
    entry = reference("petB", junctions_of_it=(only,))
    call = called(entry, (5000, 5060), notes=(NOT_REACHED_FIVE_PRIME + " — the start is partial",))
    result = CallResult(calls=(call,))

    rebuilt = transfer_structure(result, {only.key: [mapped(only, left_at=4994, right_at=5000)]})

    (only_call,) = rebuilt.calls
    assert [part.start + 1 for part in only_call.parts] == [4994, 5000]
    assert only_call.notes == ()
    assert not only_call.partial_start
    assert any("exon structure was placed" in note for note in rebuilt.findings)
