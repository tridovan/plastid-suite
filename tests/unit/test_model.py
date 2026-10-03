"""Unit tests for the core data model — synthetic data, no file I/O."""

from __future__ import annotations

import pytest

from plastid.core.model import Gene, Genome, Interval


def test_interval_converts_between_conventions() -> None:
    interval = Interval.from_1based(10, 20)
    assert (interval.start, interval.end) == (9, 20)  # 0-based half-open internally
    assert interval.to_1based() == (10, 20)  # 1-based inclusive out
    assert interval.length == 11


def test_interval_rejects_bad_input() -> None:
    with pytest.raises(ValueError, match="end"):
        Interval(5, 4)
    with pytest.raises(ValueError, match="start"):
        Interval(-1, 4)
    with pytest.raises(ValueError, match="strand"):
        Interval(0, 4, strand=2)


def test_subseq_wraps_around_a_circular_origin() -> None:
    genome = Genome(name="synthetic", sequence="ACGTACGTAC")
    # positions 8, 9, then 0, 1 of a circular 10 bp genome
    assert genome.subseq(Interval(8, 12)) == "ACAC"


def test_subseq_rejects_intervals_longer_than_the_genome() -> None:
    genome = Genome(name="synthetic", sequence="ACGT")
    with pytest.raises(ValueError, match="exceeds"):
        genome.subseq(Interval(0, 10))


def test_linear_genome_does_not_wrap() -> None:
    genome = Genome(name="synthetic", sequence="ACGTACGTAC", topology="linear")
    assert genome.subseq(Interval(8, 10)) == "AC"


def test_gene_with_mixed_strands_reports_zero() -> None:
    gene = Gene(
        name="rps12",
        kind="CDS",
        parts=(Interval(0, 10, -1), Interval(20, 30, 1)),
    )
    assert gene.strand == 0  # trans-spliced across the IR: parts disagree
    assert gene.length == 20  # exons only
    assert gene.span == (1, 30)  # 1-based inclusive outer bounds


def test_gene_needs_at_least_one_part() -> None:
    with pytest.raises(ValueError, match="no location parts"):
        Gene(name="broken", kind="CDS", parts=())
