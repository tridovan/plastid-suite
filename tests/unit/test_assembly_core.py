"""Reading an assembly: circular identity and read support, on answers known by construction.

The read-support tests use reads placed by hand on a small circle, with the expected counts worked
out from the definition. `test_profile_matches_a_direct_count` adds a second oracle: it counts every
boundary by walking every read, a different computation from the running sums under test.
"""

from __future__ import annotations

import random

import pytest

from plastid.core.assembly import (
    SUPPORT_SUPPORTED,
    SUPPORT_UNDETERMINED,
    AlignedRead,
    AssemblyError,
    SupportRule,
    build_profile,
    call_boundary,
    compare_to_reference,
    ir_junctions,
    orient_to,
    parse_sam,
    same_circular,
    weak_regions,
)
from plastid.core.model import Interval
from plastid.core.sequence import reverse_complement

RULE = SupportRule(min_spanning=2, clip_window=3)


def _random_sequence(length: int, seed: int) -> str:
    source = random.Random(seed)
    return "".join(source.choice("ACGT") for _ in range(length))


def _read(start: int, end: int, *, left: int = 0, right: int = 0) -> AlignedRead:
    return AlignedRead(
        start=start,
        end=end,
        left_clip=left,
        right_clip=right,
        matched=end - start,
        mismatched=0,
        inserted=0,
        deleted=0,
        blocks=((start, end),),
        segments=((start, end),),
    )


# --- circular identity ---------------------------------------------------------------------------


def test_same_circular_reads_any_origin_and_either_strand() -> None:
    sequence = _random_sequence(500, 1)
    rotated = sequence[137:] + sequence[:137]
    assert same_circular(sequence, rotated)
    assert same_circular(sequence, reverse_complement(rotated))


def test_same_circular_refuses_a_substitution_and_a_length_difference() -> None:
    sequence = _random_sequence(500, 2)
    changed = sequence[:250] + ("A" if sequence[250] != "A" else "C") + sequence[251:]
    assert not same_circular(sequence, changed)
    assert not same_circular(sequence, sequence + "A")


def test_orient_to_recovers_the_reference_frame_from_a_rotated_reverse_complement() -> None:
    reference = _random_sequence(8_000, 3)
    query = reverse_complement(reference[3_100:] + reference[:3_100])
    assert orient_to(query, reference) == reference


def test_orient_to_is_none_when_no_anchor_exists() -> None:
    assert orient_to(_random_sequence(8_000, 4), _random_sequence(8_000, 5)) is None


def test_comparison_of_one_circle_says_identical_whatever_the_frame() -> None:
    reference = _random_sequence(6_000, 6)
    query = reverse_complement(reference[1_000:] + reference[:1_000])
    comparison = compare_to_reference(query, reference)
    assert comparison.identical
    assert comparison.identity == 1.0
    assert comparison.length_delta == 0


def test_comparison_counts_a_known_substitution_and_a_known_deletion() -> None:
    reference = _random_sequence(6_000, 7)
    substituted = "A" if reference[3_000] != "A" else "C"
    # One substitution at 3000, and 40 bases deleted at 4500: 5960 bases of the reference remain,
    # 5959 of them equal and one different.
    query = reference[:3_000] + substituted + reference[3_001:4_500] + reference[4_540:]
    comparison = compare_to_reference(query[700:] + query[:700], reference)
    assert not comparison.identical
    assert comparison.length_delta == -40
    assert comparison.mismatches == 1
    assert comparison.gap_columns == 40
    assert comparison.origin_gap_columns == 0
    assert comparison.matched == 5_959
    assert comparison.identity == pytest.approx(5_959 / 6_000)
    [first, second] = comparison.differences
    assert (first.reference_start, first.reference_end, first.mismatches) == (3_000, 3_001, 1)
    assert (second.reference_start, second.reference_end) == (4_500, 4_540)
    assert second.gap_columns == 40


def test_comparison_places_a_long_insertion_on_either_side() -> None:
    reference = _random_sequence(9_000, 12)
    insert = _random_sequence(429, 13)
    query = reference[:5_000] + insert + reference[5_000:]
    gained = compare_to_reference(query, reference)
    assert gained.gap_columns == 429
    assert gained.mismatches == 0
    assert [item.reference_start for item in gained.differences] == [5_000]
    lost = compare_to_reference(reference, query)
    assert lost.gap_columns == 429
    assert [item.reference_start for item in lost.differences] == [5_000]
    assert lost.differences[0].reference_end == 5_429


def test_comparison_counts_bases_missing_at_the_join_of_the_circle() -> None:
    # A circle lacking 139 bases at the reference's own origin: a real difference on a circle, and
    # no rotation can make it an edge effect.
    reference = _random_sequence(8_000, 14)
    query = reference[:-139]
    comparison = compare_to_reference(query, reference)
    assert not comparison.identical
    assert comparison.gap_columns == 139
    assert comparison.origin_gap_columns == 139
    assert comparison.identity == pytest.approx((8_000 - 139) / 8_000)


def test_comparison_gives_no_identity_for_a_region_it_could_only_read_ungapped() -> None:
    # The other isomer of a plastome: one long region reversed, flanked by an exact repeat. Reading
    # it column by column would report most of the region as mismatches.
    source = random.Random(11)
    lsc = "".join(source.choice("ACGT") for _ in range(12_000))
    repeat = "".join(source.choice("ACGT") for _ in range(4_000))
    ssc = "".join(source.choice("ACGT") for _ in range(9_000))
    reference = lsc + repeat + ssc + reverse_complement(repeat)
    isomer = lsc + repeat + reverse_complement(ssc) + reverse_complement(repeat)
    comparison = compare_to_reference(isomer, reference)
    assert not comparison.identical
    assert comparison.identity is None
    assert comparison.mismatches is None


def test_comparison_reports_no_number_when_the_sequences_cannot_be_placed() -> None:
    comparison = compare_to_reference(_random_sequence(6_000, 8), _random_sequence(6_100, 9))
    assert not comparison.identical
    assert comparison.identity is None
    assert comparison.method == "no anchor"


# --- SAM -----------------------------------------------------------------------------------------


def _sam(flag: int, position: int, cigar: str) -> str:
    return f"r\t{flag}\tseq\t{position}\t60\t{cigar}\t*\t0\t0\tACGT\tIIII\n"


def test_parse_sam_reads_a_span_clips_and_operations_from_the_cigar() -> None:
    [read] = parse_sam(["@HD\tVN:1.6\n", _sam(0, 11, "5S10=1X4=2I3=1D6=7S")])
    assert (read.start, read.left_clip, read.right_clip) == (10, 5, 7)
    # Reference consumed: 10 + 1 + 4 + 3 + 1 + 6 = 25, so the span ends at 35.
    assert read.end == 35
    assert (read.matched, read.mismatched, read.inserted, read.deleted) == (23, 1, 2, 1)
    assert read.blocks == ((10, 20), (20, 21), (21, 25), (25, 28), (29, 35))
    assert read.segments == ((10, 35),)


def test_a_long_indel_cuts_the_read_into_two_segments() -> None:
    [deletion] = parse_sam([_sam(0, 1, "40=139D40=")])
    assert (deletion.start, deletion.end) == (0, 219)
    assert deletion.segments == ((0, 40), (179, 219))
    [insertion] = parse_sam([_sam(0, 1, "40=12I40=")])
    assert insertion.segments == ((0, 40), (40, 80))
    [short] = parse_sam([_sam(0, 1, "40=9D40=")])
    assert short.segments == ((0, 89),)


@pytest.mark.parametrize("flag", [0x4, 0x100, 0x800, 0x100 | 0x10])
def test_parse_sam_drops_everything_but_primary_alignments(flag: int) -> None:
    assert list(parse_sam([_sam(flag, 1, "4=")])) == []


def test_parse_sam_refuses_a_cigar_that_cannot_count_mismatches() -> None:
    with pytest.raises(AssemblyError, match="--eqx"):
        list(parse_sam([_sam(0, 1, "50M")]))


# --- read support --------------------------------------------------------------------------------


def test_a_read_across_the_origin_folds_back_onto_the_circle() -> None:
    # On a circle of 100 with 30 bases of extension: a read over 85..115 covers bases 85-99 and
    # 0-14. With flank 10 it spans boundaries 95..105, which fold to 95-99 and 0-5.
    profile = build_profile([_read(85, 115)], length=100, extension=30, flank=10, min_clip=5)
    spanned = [b for b in range(100) if profile.spanning[b]]
    assert spanned == [0, 1, 2, 3, 4, 5, 95, 96, 97, 98, 99]
    covered = [i for i in range(100) if profile.depth[i]]
    assert covered == list(range(15)) + list(range(85, 100))


def test_a_read_in_the_extension_counts_at_the_start_of_the_circle() -> None:
    profile = build_profile([_read(102, 122)], length=100, extension=30, flank=5, min_clip=5)
    assert [b for b in range(100) if profile.spanning[b]] == list(range(7, 18))
    assert profile.reads == 1


def test_profile_matches_a_direct_count() -> None:
    source = random.Random(10)
    length, extension, flank, min_clip = 200, 40, 7, 4
    reads = []
    for _ in range(400):
        start = source.randrange(0, length + extension - 30)
        end = start + source.randrange(20, 31)
        reads.append(
            _read(start, end, left=source.choice((0, 3, 9)), right=source.choice((0, 2, 12)))
        )
    profile = build_profile(
        reads, length=length, extension=extension, flank=flank, min_clip=min_clip
    )
    expected_span = [0] * length
    expected_clip = [0] * length
    expected_depth = [0] * length
    for read in reads:
        for boundary in range(length + extension + 1):
            if read.start + flank <= boundary <= read.end - flank:
                expected_span[boundary % length] += 1
        for base in range(read.start, read.end):
            expected_depth[base % length] += 1
        if read.left_clip >= min_clip:
            expected_clip[read.start % length] += 1
        if read.right_clip >= min_clip:
            expected_clip[read.end % length] += 1
    assert list(profile.spanning) == expected_span
    assert list(profile.clipped) == expected_clip
    assert list(profile.depth) == expected_depth


def test_an_alignment_past_the_extended_sequence_is_refused() -> None:
    with pytest.raises(AssemblyError, match="past the extended sequence"):
        build_profile([_read(120, 140)], length=100, extension=30, flank=5, min_clip=5)


def test_error_rate_is_undefined_without_a_matched_base() -> None:
    profile = build_profile([], length=100, extension=10, flank=5, min_clip=5)
    assert profile.error_rate is None
    assert profile.mean_depth == 0


def test_a_join_crossed_by_reads_is_supported() -> None:
    # A read [30 + o, 80 + o) keeps 10 bases left of boundary 55 only while o <= 15: 16 reads.
    reads = [_read(30 + offset, 80 + offset) for offset in range(20)]
    profile = build_profile(reads, length=200, extension=50, flank=10, min_clip=5)
    call = call_boundary(profile, RULE, name="JLB", boundary=55)
    assert call.call == SUPPORT_SUPPORTED
    assert call.spanning == 16
    assert call.clipped == 0
    assert call.fraction == 1.0


def test_a_read_does_not_cross_what_it_jumps_over() -> None:
    # 139 bases of the target that the molecule lacks: the aligner bridges them with a deletion, so
    # the span runs from 0 to 219 and not a single read has a base inside the bridged stretch.
    reads = list(parse_sam([_sam(0, 1, "40=139D40=") for _ in range(50)]))
    profile = build_profile(reads, length=500, extension=50, flank=10, min_clip=5)
    assert profile.spanning[100] == 0
    assert profile.spanning[20] == 50
    call = call_boundary(profile, RULE, name="origin", boundary=100)
    assert call.call == SUPPORT_UNDETERMINED


def test_a_join_no_read_crosses_is_undetermined_however_deep_the_flanks() -> None:
    # Reads on both sides stop at 100 with sequence left over: the signature of a wrong join.
    reads = [_read(40 + offset, 100, right=30) for offset in range(30)]
    reads += [_read(100, 160 - offset, left=30) for offset in range(30)]
    profile = build_profile(reads, length=300, extension=50, flank=10, min_clip=5)
    call = call_boundary(profile, RULE, name="JLB", boundary=100)
    assert call.spanning == 0
    assert call.clipped == 60
    assert call.call == SUPPORT_UNDETERMINED


def test_reads_stopping_are_reported_and_do_not_decide() -> None:
    # Half the reads stop at a motif on correct sequence; those that cross still prove the join.
    reads = [_read(40 + offset, 100, right=30) for offset in range(30)]
    reads += [_read(60 + offset, 120 + offset) for offset in range(12)]
    profile = build_profile(reads, length=300, extension=50, flank=10, min_clip=5)
    call = call_boundary(profile, RULE, name="JLB", boundary=100)
    assert call.spanning == 12
    assert call.clipped == 30
    assert call.fraction == pytest.approx(12 / 42)
    assert call.call == SUPPORT_SUPPORTED


def test_a_join_nobody_read_is_undetermined_and_has_no_fraction() -> None:
    profile = build_profile([], length=100, extension=10, flank=5, min_clip=5)
    call = call_boundary(profile, RULE, name="origin", boundary=0)
    assert call.fraction is None
    assert call.call == SUPPORT_UNDETERMINED


def test_clips_a_few_bases_from_the_boundary_still_count_against_it() -> None:
    reads = [_read(20, 58, right=20) for _ in range(30)]
    profile = build_profile(reads, length=200, extension=50, flank=10, min_clip=5)
    call = call_boundary(profile, RULE, name="x", boundary=60)
    assert call.clipped == 30


def test_weak_regions_finds_a_gap_and_keeps_a_run_across_the_origin_whole() -> None:
    # Reads cover boundaries 20..80 well and leave the rest of the circle of 100 unread.
    reads = [_read(10, 90) for _ in range(5)]
    profile = build_profile(reads, length=100, extension=20, flank=10, min_clip=5)
    [region] = weak_regions(profile, RULE)
    assert (region.start, region.end, region.length, region.lowest) == (81, 19, 39, 0)


def test_weak_regions_is_empty_on_a_well_covered_circle() -> None:
    reads = [_read(start, start + 60) for start in range(0, 100, 5) for _ in range(4)]
    profile = build_profile(reads, length=100, extension=60, flank=10, min_clip=5)
    assert weak_regions(profile, RULE) == ()


def test_weak_regions_names_an_unread_circle_once() -> None:
    profile = build_profile([], length=100, extension=10, flank=5, min_clip=5)
    [region] = weak_regions(profile, RULE)
    assert (region.start, region.end, region.length) == (0, 99, 100)


def test_a_correct_low_depth_sample_is_not_a_finding() -> None:
    # Three reads cross everywhere, as at 22×: an absolute count of ten reads calls this all weak.
    reads = [_read(start, start + 60) for start in range(0, 200, 30) for _ in range(3)]
    profile = build_profile(reads, length=200, extension=60, flank=10, min_clip=5)
    assert min(profile.spanning) >= 2
    assert weak_regions(profile, RULE) == ()


# --- junction names ------------------------------------------------------------------------------


def test_junctions_on_the_usual_frame_lsc_first() -> None:
    # LSC 0-80, IRb 80-100, SSC 100-120, IRa 120-140, on a circle of 140.
    first, second = Interval(80, 100), Interval(120, 140)
    assert ir_junctions(first, second, length=140) == (
        ("JLB", 80),
        ("JSB", 100),
        ("JSA", 120),
        ("JLA", 0),
    )


def test_junctions_when_the_small_region_is_at_the_origin() -> None:
    # SSC 0-20, IRa 20-40, LSC 40-120, IRb 120-140: the longer arc lies between the copies.
    first, second = Interval(20, 40), Interval(120, 140)
    assert ir_junctions(first, second, length=140) == (
        ("JSA", 20),
        ("JLA", 40),
        ("JLB", 120),
        ("JSB", 0),
    )


def test_junctions_come_back_on_the_sequence_the_frame_was_rotated_from() -> None:
    first, second = Interval(80, 100), Interval(120, 140)
    names = dict(ir_junctions(first, second, length=140, frame_offset=30))
    assert names == {"JLB": 110, "JSB": 130, "JSA": 10, "JLA": 30}
