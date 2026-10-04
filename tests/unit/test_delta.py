"""The MUMmer delta reader — the alignment behind the %identity-by-coordinate series (P4-01).

The format is small and its failure mode is silent, so the tests pin the two things that were
measured, not assumed:

* **one ``>`` header per sequence pair, not per alignment** — MUMmer 4 writes every alignment of
  a pair under a single header, each record terminated by a ``0`` line (measured 04/10/2026 on the
  Iris pair: 7 alignments under one header).
* **the gap numbers only decode one way.** The walk was fixed against the real file: each value
  ``d`` consumes ``|d| - 1`` paired columns; ``d > 0`` then advances the *target* index (the query
  has a gap), ``d < 0`` advances the *query* index. Verified on the golden pair: all seven blocks
  land exactly on their end coordinates and ``mismatches + gap bases`` equals the header's own
  error count (block 1: 152 + 152? — see the vendored fixture test, which checks the real file).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.delta import DeltaBlock, iter_aligned_pairs, parse_delta
from plastid.core.sequence import reverse_complement

# One perfect forward block.
SIMPLE = """/tmp/t.fa /tmp/q.fa
NUCMER
>t q 10 10
1 10 1 10 0 0 0
0
"""

# One forward block with a single-base gap in the query: target 12 bases, query 11.
# d = 4 ⇒ 3 paired columns, then one target base with no query partner, then the tail.
GAPPED = """/tmp/t.fa /tmp/q.fa
NUCMER
>t q 12 11
1 12 1 11 2 2 0
4
0
"""

# A reverse block: the query runs backwards (s2 > e2) — what the IR copies produce.
REVERSE = """/tmp/t.fa /tmp/q.fa
NUCMER
>t q 10 10
1 10 10 1 0 0 0
0
"""

TWO_BLOCKS_ONE_HEADER = """/tmp/t.fa /tmp/q.fa
NUCMER
>t q 40 40
1 20 1 20 0 0 0
0
21 40 21 40 0 0 0
0
"""


def test_one_header_holds_every_alignment_of_a_pair() -> None:
    blocks = parse_delta(TWO_BLOCKS_ONE_HEADER)
    assert len(blocks) == 2
    assert [(block.target_start, block.target_end) for block in blocks] == [(0, 20), (20, 40)]


def test_coordinates_are_converted_to_zero_based_half_open() -> None:
    (block,) = parse_delta(SIMPLE)
    assert (block.target_start, block.target_end) == (0, 10)
    assert (block.query_start, block.query_end) == (0, 10)
    assert block.target_length == 10
    assert block.query_length == 10
    assert block.reverse is False


def test_a_gap_value_consumes_paired_columns_and_then_advances_one_side() -> None:
    (block,) = parse_delta(GAPPED)
    # 3 paired columns, then a target base the query does not have, then the remaining 8 pairs:
    # the target consumes 12 bases in two gapless segments, the query 11.
    assert [(segment.target_start, segment.target_end) for segment in block.segments] == [
        (0, 3),
        (4, 12),
    ]
    assert [(segment.query_start, segment.query_end) for segment in block.segments] == [
        (0, 3),
        (3, 11),
    ]
    assert sum(segment.length for segment in block.segments) == 11  # paired bases only


def test_a_reverse_block_is_marked_and_its_query_coordinates_are_oriented() -> None:
    (block,) = parse_delta(REVERSE)
    assert block.reverse is True
    # forward coordinates keep start < end …
    assert (block.query_start, block.query_end) == (0, 10)
    # … and the segments carry oriented query coordinates, so the pairing is 1:1 in k
    (segment,) = block.segments
    assert (segment.target_start, segment.target_end) == (0, 10)
    assert (segment.query_start, segment.query_end) == (0, 10)


def test_aligned_pairs_walk_both_strands_and_count_matches() -> None:
    target = "ACGTACGTAA"
    query = "ACGTACCTAA"  # one mismatch at position 6
    (block,) = parse_delta(SIMPLE)
    pairs = list(iter_aligned_pairs(block, target=target, query=query))
    assert len(pairs) == 10
    assert sum(1 for _, _, t, q in pairs if t != q) == 1
    assert pairs[6][2:] == ("G", "C")

    # the reverse block pairs the target's first base with the query's last, complemented
    (reverse_block,) = parse_delta(REVERSE)
    reverse_query = reverse_complement(target)  # a perfect reverse match
    pairs = list(iter_aligned_pairs(reverse_block, target=target, query=reverse_query))
    assert pairs[0] == (0, 9, "A", "A")  # query base T, complemented
    assert pairs[1] == (1, 8, "C", "C")
    assert all(t == q for _, _, t, q in pairs)


def test_a_gap_base_is_skipped_rather_than_scored() -> None:
    (block,) = parse_delta(GAPPED)
    target = "ACGTACGTACGT"  # index 3 is the base the query does not carry
    query = "ACGACGTACGT"
    pairs = list(iter_aligned_pairs(block, target=target, query=query))
    assert len(pairs) == 11  # paired bases only — the gap is not a row
    assert all(t == q for _, _, t, q in pairs)
    assert 3 not in [target_position for target_position, _, _, _ in pairs]


def test_window_identity_counts_only_paired_unambiguous_bases() -> None:
    from plastid.core.delta import window_identity

    target = "ACGTACGTAA"
    query = "ACGTACCTAA"  # mismatch at position 6
    (block,) = parse_delta(SIMPLE)
    windows = window_identity((block,), target=target, query=query, window=5, step=5)
    assert [(w.start_1b, w.end_1b) for w in windows] == [(1, 5), (6, 10)]
    assert windows[0].compared_bp == 5
    assert windows[0].identity_pct == 100.0
    assert windows[1].mismatches == 1
    assert windows[1].identity_pct == 80.0

    # a window with no alignment has no value — never a zero
    empty = window_identity((), target=target, query=query, window=5, step=5)
    assert [w.identity_pct for w in empty] == [None, None]
    assert [w.compared_bp for w in empty] == [0, 0]


def test_window_identity_does_not_score_gap_bases() -> None:
    from plastid.core.delta import window_identity

    (block,) = parse_delta(GAPPED)
    target = "ACGTACGTACGT"
    query = "ACGACGTACGT"  # 11 pairs; target position 3 has no partner
    (window,) = window_identity((block,), target=target, query=query, window=12, step=12)
    assert window.compared_bp == 11
    assert window.identity_pct == 100.0


def test_aligned_totals_do_not_double_count_overlapping_windows() -> None:
    from plastid.core.delta import aligned_totals

    target = "ACGTACGTAA"
    query = "ACGTACCTAA"
    (block,) = parse_delta(SIMPLE)
    assert aligned_totals((block,), target=target, query=query) == (10, 1)
    # the same numbers the windows must sum to when they do not overlap
    from plastid.core.delta import window_identity

    windows = window_identity((block,), target=target, query=query, window=5, step=5)
    assert sum(w.compared_bp for w in windows) == 10
    assert sum(w.mismatches for w in windows) == 1


def test_window_identity_rejects_impossible_windows() -> None:
    from plastid.core.delta import window_identity

    (block,) = parse_delta(SIMPLE)
    with pytest.raises(ValueError):
        window_identity((block,), target="ACGTACGTAA", query="ACGTACGTAA", window=0, step=5)
    with pytest.raises(ValueError):
        window_identity((block,), target="ACGTACGTAA", query="ACGTACGTAA", window=4, step=5)


def test_a_malformed_record_raises_rather_than_guessing() -> None:
    with pytest.raises(ValueError):
        parse_delta("/tmp/t.fa /tmp/q.fa\nNUCMER\n>t q 10 10\n1 10 1 10 0 0\n0\n")
    with pytest.raises(ValueError):
        parse_delta("/tmp/t.fa /tmp/q.fa\nNUCMER\n>t q 10 10\n1 10 1 10 0 0 0\n5\n")  # no 0
    with pytest.raises(ValueError):
        parse_delta("not a delta file\n")
    with pytest.raises(ValueError):
        parse_delta("/tmp/t.fa /tmp/q.fa\nNUCMER\n>t q\n1 10 1 10 0 0 0\n0\n")  # no lengths


def test_an_empty_delta_is_no_blocks_not_an_error() -> None:
    """A pair with nothing in common is a legitimate result: zero blocks, said plainly."""
    assert parse_delta("/tmp/t.fa /tmp/q.fa\nNUCMER\n") == ()


def test_the_walk_lands_on_the_headers_own_end_coordinates() -> None:
    """The parser's internal consistency: every byte the header claims is consumed once."""
    (block,) = parse_delta(GAPPED)
    assert block.segments[-1].target_end == block.target_end
    assert block.segments[-1].query_end == block.query_end
    assert isinstance(block, DeltaBlock)


def test_a_real_delta_fragment_parses_to_its_measured_spans() -> None:
    """The vendored fragment is two alignment records of the measured Iris-pair run.

    It is here because a synthetic delta cannot catch a wrong walk: the synthetic cases above
    were written by the same hand as the parser. This one was produced by MUMmer 4.0.1 on the
    golden pair (2026-10-04) and its spans, orientation and totals are what the tool printed.
    """
    path = Path(__file__).parents[1] / "data" / "iris_pair_sample.delta"
    blocks = parse_delta(path.read_text(encoding="utf-8"))
    assert len(blocks) == 2
    forward, reverse = blocks
    assert (forward.target_name, forward.query_name) == ("NC_050833.1", "NC_024936.1")
    assert (forward.target_length, forward.query_length) == (153724, 153441)

    # block 1 of the real run: 1..13576 vs 95..13614, 272 errors
    assert (forward.target_start, forward.target_end) == (0, 13576)
    assert (forward.query_start, forward.query_end) == (94, 13614)
    assert forward.errors == 272
    assert forward.reverse is False

    # the reverse block that carries an IR copy: 83128..109341 vs 153441..127260
    assert reverse.reverse is True
    assert (reverse.target_start, reverse.target_end) == (83127, 109341)
    assert (reverse.query_start, reverse.query_end) == (127259, 153441)

    # every base the header claims is consumed exactly once, on both sides
    for block in blocks:
        last = block.segments[-1]
        assert last.target_end == block.target_end
        # the walk's last paired base is the block's last base on the query too — the forward
        # end when the block runs with the query, the forward start when it runs against it
        expected = block.query_start if block.reverse else block.query_end - 1
        assert block.query_forward(last.query_end - 1) == expected
