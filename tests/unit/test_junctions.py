"""The four regions and the four junctions, in the record's own frame.

Two of these are the golden contract's own numbers (IR 26,214 · LSC 83,127 · SSC 18,169), which is
the point of reading the record's annotation rather than a rotation of it: a junction is a statement
about the deposited coordinates.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank, read_genbank_records
from plastid.core.junctions import JunctionError, junction_regions
from plastid.core.model import Genome, Interval, RepeatRegion

# The geometry of Chasmanthe floribunda NC_072562, whose IRa runs across the origin: the copy is
# written `join(1..1,124854..151293)`, so it arrives as two parts of **one** copy. Dropping the
# 1 bp piece is what makes a copy's bases stop pairing with the other's — measured 17/10/2026,
# identity 0.2924 instead of 1.0000.
WRAPPED_LENGTH = 151_293
WRAPPED = Genome(
    name="Chasmanthe floribunda",
    accession="NC_072562.1",
    sequence="A" * WRAPPED_LENGTH,
    repeat_regions=(
        RepeatRegion((Interval(81_002, 107_443),)),  # IRb, 26,441 bp
        RepeatRegion(
            (Interval(0, 1), Interval(124_853, WRAPPED_LENGTH))
        ),  # IRa, one copy across the origin
    ),
)


def test_golden_regions_are_the_contract_numbers(golden_genome: Genome) -> None:
    """The golden's four regions, in the frame the record was deposited in."""
    regions = junction_regions(golden_genome)
    assert regions.lsc.length == 83_127
    assert regions.irb.length == 26_214
    assert regions.ssc.length == 18_169
    assert regions.ira.length == 26_214
    assert regions.config_line == (
        "LSC:1-83127,IRb:83128-109341,SSC:109342-127510,IRa:127511-153724"
    )
    assert regions.tiles_genome


def test_junctions_are_the_first_base_of_the_region_they_lead_into(golden_genome: Genome) -> None:
    """JLB is IRb's first base; JSA is the base before IRa's — one off is a different figure."""
    regions = junction_regions(golden_genome)
    assert (regions.jlb, regions.jsb, regions.jsa, regions.jla) == (
        83_128,
        109_342,
        127_510,
        127_511,
    )


def test_a_copy_across_the_origin_keeps_its_piece() -> None:
    """The wrapped piece belongs to IRa, and the region says so rather than being drawn short."""
    regions = junction_regions(WRAPPED)
    assert regions.ira.wraps
    assert regions.ira.text == "124854-151293;1-1"
    assert regions.ira.length == 26_441
    assert regions.irb.length == 26_441
    assert regions.lsc.text == "2-81002"
    assert regions.tiles_genome
    assert regions.config_line == (
        "LSC:2-81002,IRb:81003-107443,SSC:107444-124853,IRa:124854-151293;1-1"
    )


def test_either_copy_may_be_the_one_across_the_origin() -> None:
    """Which copy wraps depends on where the record's origin falls — IRb here, IRa above.

    Measured against the first version of this module, which took "the base after a region" from
    the region's main run: with IRb wrapping, JSB came out at 60,001 — one past the end of a
    60,000 bp record — and SSC started there too.
    """
    wrapped_irb = Genome(
        name="IRb across the origin",
        accession="XX_000005.1",
        sequence="A" * 60_000,
        repeat_regions=(
            RepeatRegion(
                (Interval(0, 1_000), Interval(55_000, 60_000))
            ),  # IRb, one copy across the origin
            RepeatRegion((Interval(5_000, 11_000),)),  # IRa
        ),
    )
    regions = junction_regions(wrapped_irb)
    assert regions.irb.wraps
    assert not regions.ira.wraps
    assert regions.irb.text == "55001-60000;1-1000"
    assert regions.irb.length == 6_000
    assert regions.lsc.text == "11001-55000"
    assert regions.ssc.text == "1001-5000"
    assert (regions.jlb, regions.jsb, regions.jsa, regions.jla) == (55_001, 1_001, 5_000, 5_001)
    assert regions.tiles_genome


def test_a_record_without_ir_annotation_is_refused_by_name() -> None:
    """No annotation, no figure: the message says what was looked for, not that something failed."""
    bare = Genome(name="no IR", accession="XX_000001.1", sequence="ACGT" * 4_000)
    with pytest.raises(JunctionError, match="expected two repeat_region copies"):
        junction_regions(bare)


def test_overlapping_copies_are_refused() -> None:
    overlapping = Genome(
        name="overlap",
        accession="XX_000002.1",
        sequence="A" * 60_000,
        repeat_regions=(Interval(5_000, 35_000), Interval(30_000, 55_000)),
    )
    with pytest.raises(JunctionError, match="overlap"):
        junction_regions(overlapping)


def test_a_join_that_is_not_one_arc_is_refused() -> None:
    """A repeat_region whose parts are not contiguous is malformed, not two convenient pieces.

    A bare short repeat_region is a different feature, not a wrap piece (NC_063921's 1 bp repeat is
    ignored). A long *joined* region with two separate gaps claims to be one copy and is not.
    """
    stranded = Genome(
        name="two gaps in one copy",
        accession="XX_000003.1",
        sequence="A" * 60_000,
        repeat_regions=(
            RepeatRegion((Interval(0, 3_000), Interval(10_000, 13_000))),
            RepeatRegion((Interval(20_000, 26_000),)),
        ),
    )
    with pytest.raises(JunctionError, match="not one circular copy"):
        junction_regions(stranded)


def test_the_four_regions_always_tile_the_circle() -> None:
    """Two disjoint copies cut the circle into four pieces with nothing left over — an invariant.

    It is asserted rather than guarded inside :func:`junction_regions`: the derivation cannot leave
    a gap, so a raise there would be a branch no test could ever reach.
    """
    for genome in (
        WRAPPED,
        Genome(
            name="plain",
            accession="XX_000004.1",
            sequence="A" * 60_001,
            repeat_regions=(Interval(5_000, 30_000), Interval(35_000, 60_000)),
        ),
    ):
        regions = junction_regions(genome)
        assert regions.tiles_genome
        assert regions.total_bp == len(genome)


def test_the_iridaceae_edge_cases_read_or_are_refused_loudly(data_dir: Path) -> None:
    """The vendored hard records: each either cuts into four regions or names why it cannot."""
    records = read_genbank_records(data_dir / "iridaceae_edge_cases.gb")
    assert records, "the edge-case fixture must hold records"
    for record in records:
        try:
            regions = junction_regions(record)
        except JunctionError as error:
            assert record.accession in str(error), "a refusal names the record it refuses"
        else:
            assert regions.tiles_genome


def test_the_golden_file_on_disk_agrees_with_the_fixture(data_dir: Path) -> None:
    """Guard the fixture: the numbers above are read from the same file the toolkit ships."""
    regions = junction_regions(read_genbank(data_dir / "NC_050833.1.gb"))
    assert regions.irb.length == 26_214
