"""Alignment maths — analytic checks on hand-built alignments (offline, no aligner)."""

from __future__ import annotations

import pytest

from plastid.core.align import (
    SUBSTITUTION_TYPES,
    Alignment,
    columns_for_span,
    nucleotide_diversity,
    reference_positions,
    substitution_class,
    substitution_counts,
    variable_sites,
    window_diversity,
)


def _two(second: str) -> Alignment:
    return Alignment(names=("ref", "other"), sequences=("ACGTACGT", second))


def test_rows_must_be_equal_length() -> None:
    with pytest.raises(ValueError, match="differ in length"):
        Alignment(names=("a", "b"), sequences=("ACGT", "ACG"))


def test_diversity_is_zero_for_identical_rows() -> None:
    pi, usable = nucleotide_diversity(_two("ACGTACGT"), range(8))
    assert pi == 0.0
    assert usable == 8


def test_diversity_of_one_difference_in_four_sites() -> None:
    # one differing site in four: its site value is 2/1 * (1 - 2*(1/2)^2) = 1, so π = 1/4
    alignment = Alignment(names=("a", "b"), sequences=("AAAA", "AATA"))
    pi, usable = nucleotide_diversity(alignment, range(4))
    assert usable == 4
    assert pi == pytest.approx(1 / 4)


def test_gapped_columns_are_dropped_and_counted() -> None:
    alignment = Alignment(names=("a", "b"), sequences=("ACGT", "A-GT"))
    pi, usable = nucleotide_diversity(alignment, range(4))
    assert usable == 3  # the gap column is not scored as a match
    assert pi == 0.0


def test_diversity_is_none_when_nothing_is_usable() -> None:
    alignment = Alignment(names=("a", "b"), sequences=("-N", "-N"))
    assert nucleotide_diversity(alignment, range(2)) == (None, 0)


def test_three_rows_weight_the_majority_allele() -> None:
    # one odd base out of three: 3/2 * (1 - (2/3)^2 - (1/3)^2) = 3/2 * 4/9 = 2/3
    alignment = Alignment(names=("a", "b", "c"), sequences=("A", "A", "T"))
    pi, usable = nucleotide_diversity(alignment, range(1))
    assert (pi, usable) == (pytest.approx(2 / 3), 1)


def test_window_step_is_validated() -> None:
    alignment = _two("ACGTACGT")
    with pytest.raises(ValueError, match="exceeds window"):
        window_diversity(alignment, window=4, step=5)
    with pytest.raises(ValueError, match="must be positive"):
        window_diversity(alignment, window=0, step=1)


def test_windows_cover_the_alignment_and_report_usable_columns() -> None:
    alignment = _two("ACGTACGT")
    windows = window_diversity(alignment, window=4, step=4)
    assert [(w.start_1b, w.end_1b) for w in windows] == [(1, 4), (5, 8)]
    assert all(w.usable_columns == 4 for w in windows)


def test_reference_positions_skip_gaps() -> None:
    alignment = Alignment(names=("a", "b"), sequences=("A-CGT", "ATCGT"))
    assert reference_positions(alignment) == (0, None, 1, 2, 3)


def test_span_to_columns_follows_the_reference_row() -> None:
    alignment = Alignment(names=("a", "b"), sequences=("A-CGT", "ATCGT"))
    # reference positions 1..2 (0-based half-open) sit at columns 2 and 3
    assert columns_for_span(alignment, 0, (1, 3)) == (2, 3)


def test_substitution_counts_split_transitions_from_transversions() -> None:
    counts = substitution_counts("AAAA", "GGGG")  # A>G four times: all transitions
    assert counts["transitions"] == 4 and counts["transversions"] == 0
    counts = substitution_counts("AAAA", "CCCC")  # A>C: a transversion each
    assert counts["transitions"] == 0 and counts["transversions"] == 4
    assert counts["compared"] == 4 and counts["skipped"] == 0


def test_substitution_counts_skip_gaps_and_ambiguous_bases() -> None:
    counts = substitution_counts("ACGT", "A-GN")
    assert counts["compared"] == 2  # only A and T were comparable
    assert counts["skipped"] == 2
    assert counts["substitutions"] == 0


def test_every_substitution_is_counted_under_its_own_directed_type() -> None:
    """The twelve directed pairs are the primitive fact; transitions/transversions aggregate them.

    The manuscript's Figure 7 draws the six *complementary* pairs, so the toolkit has to hold the
    directed counts they are summed from — measured 17/10/2026, it held only the two aggregates.
    """
    counts = substitution_counts("AGCT", "GATC")  # one of each transition, no transversion
    assert counts["A_to_G"] == 1
    assert counts["G_to_A"] == 1
    assert counts["C_to_T"] == 1
    assert counts["T_to_C"] == 1
    assert counts["transversions"] == 0
    assert sum(counts[name] for name in SUBSTITUTION_TYPES) == counts["substitutions"] == 4


def test_the_type_counters_partition_into_transitions_and_transversions() -> None:
    counts = substitution_counts("ACGTACGT", "GCTAGCAT")
    transitions = sum(
        counts[name] for name in SUBSTITUTION_TYPES if substitution_class(name) == "Transition"
    )
    transversions = sum(
        counts[name] for name in SUBSTITUTION_TYPES if substitution_class(name) == "Transversion"
    )
    assert transitions == counts["transitions"]
    assert transversions == counts["transversions"]
    assert transitions + transversions == counts["substitutions"]


def test_a_skipped_column_contributes_to_no_substitution_type() -> None:
    """An ambiguous base is not a type: it is counted as skipped, never as a substitution."""
    counts = substitution_counts("ACGT", "A-GN")
    assert sum(counts[name] for name in SUBSTITUTION_TYPES) == 0
    # every type is present as a key, so a reader never has to guess whether a zero means "none"
    assert set(SUBSTITUTION_TYPES) <= set(counts)


def test_the_type_names_are_the_twelve_directed_pairs() -> None:
    assert len(SUBSTITUTION_TYPES) == 12
    assert set(SUBSTITUTION_TYPES) == {f"{a}_to_{b}" for a in "ACGT" for b in "ACGT" if a != b}
    assert {name for name in SUBSTITUTION_TYPES if substitution_class(name) == "Transition"} == {
        "A_to_G",
        "G_to_A",
        "C_to_T",
        "T_to_C",
    }


def test_the_recorded_window_decision_matches_the_module() -> None:
    """defaults.yaml records D2-02; this keeps the record and the module from drifting."""
    from plastid.modules import diversity
    from plastid.run.defaults import load_defaults

    recorded = load_defaults()["diversity"]
    assert recorded["window_bp"] == diversity.DEFAULT_WINDOW_BP
    assert recorded["step_bp"] == diversity.DEFAULT_STEP_BP


def test_variable_and_informative_sites() -> None:
    # a site is informative only when two alleles each appear at least twice — four rows, not three
    three = Alignment(names=("a", "b", "c"), sequences=("AAT", "AAG", "AAT"))
    assert variable_sites(three) == (1, 0)
    four = Alignment(names=("a", "b", "c", "d"), sequences=("AAT", "AAG", "AAT", "AAG"))
    assert variable_sites(four) == (1, 1)


def test_an_alignment_with_lowercase_bases_is_refused() -> None:
    """Measured 03/10/2026: MAFFT 7.526 prints lowercase, and a lowercase row made every column
    count as "skipped" — the first real-data Ts/Tv run reported zero comparable columns instead
    of an error. The boundary upper-cases (external/mafft.py); the model refuses the rest."""
    with pytest.raises(ValueError, match="lowercase"):
        Alignment(names=("ref", "other"), sequences=("ACGT", "acgt"))


def test_the_mafft_parser_upper_cases_what_the_aligner_prints() -> None:
    """The sample is verbatim from the measured run: MAFFT 7.526's own output, lowercase."""
    from plastid.external.mafft import _parse_fasta

    aligned = _parse_fasta(">NC_050833.1\n------cccgatatcttact\n>NC_024936.1\nACGTACGTACGTACGTAC\n")

    assert aligned["NC_050833.1"] == "------CCCGATATCTTACT"
    assert aligned["NC_024936.1"] == "ACGTACGTACGTACGTAC"


def test_codon_frame_breaks_flags_a_gap_that_shifts_the_reading_frame() -> None:
    """yn00 reads an alignment three bases at a time, so a gap that is not a whole codon does not
    fail there — it silently reads the wrong frame. Measured 04/10/2026: a one-base gap in a
    synthetic pair made yn00 report a spurious stop codon (exit 255) only because that frame
    happened to hit one; a frame without a stop would have been scored as if it were real. The
    check is therefore ours, and it names the row and the position."""
    from plastid.core.align import codon_frame_breaks

    # 3-base gap after 6 bases: every codon stays whole — no breaks
    clean = Alignment(
        names=("a", "b"),
        sequences=("ATGGCT---TCTAAACCCGGGTTT", "ATGGCTTCTAAACCCGGGTTT---"),
    )
    assert codon_frame_breaks(clean) == ()

    # a 1-base gap: a frameshift, and the row that carries it is named
    broken = Alignment(
        names=("a", "b"), sequences=("ATGGCT-TCTAAACCCGGGTTT-", "ATGGCTTCTAAACCCGGGTTT--")
    )
    breaks = codon_frame_breaks(broken)
    assert breaks
    assert any("row 'a'" in message for message in breaks)

    # a 3-base gap that starts mid-codon: also a frameshift
    shifted = Alignment(
        names=("a", "b"), sequences=("A---TGGCTTCTAAACCCGGGTTT", "ATGGCTTCTAAACCCGGGTTT---")
    )
    breaks = codon_frame_breaks(shifted)
    assert len(breaks) == 1
    assert "row 'a'" in breaks[0]
    assert "codon boundary" in breaks[0]


def test_codon_frame_breaks_flags_a_row_whose_length_is_not_whole_codons() -> None:
    from plastid.core.align import codon_frame_breaks

    alignment = Alignment(names=("a", "b"), sequences=("ATGGCTT", "ATGGCTT"))
    breaks = codon_frame_breaks(alignment)
    assert len(breaks) == 2  # both rows, each named
    assert all("codon" in message for message in breaks)


# ---------------------------------------------------------------- MAFFT strategy (M11.20)


def test_every_strategy_name_maps_to_a_mafft_flag() -> None:
    """The mapping is data, so a typo in a flag would be a wrong alignment, not a crash."""
    from plastid.external.mafft import MAFFT_STRATEGIES

    assert set(MAFFT_STRATEGIES) == {"auto", "einsi", "linsi", "ginsi", "fftns", "nwns"}
    for name, flag in MAFFT_STRATEGIES.items():
        assert flag.startswith("--"), f"{name}: {flag!r} is not a MAFFT option"
        assert flag == f"--{name}"


def test_an_unknown_strategy_is_refused_with_the_list() -> None:
    """A misspelt algorithm must not silently fall back to MAFFT's default."""
    from plastid.external.mafft import strategy_flag

    with pytest.raises(ValueError, match="auto"):
        strategy_flag("einsi-extra")


def test_the_default_strategy_is_the_one_every_number_was_measured_with() -> None:
    """A change of default is a change of measurement; it has to be deliberate."""
    from plastid.external.mafft import DEFAULT_STRATEGY

    assert DEFAULT_STRATEGY == "auto"
