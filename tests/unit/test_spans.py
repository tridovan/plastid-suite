"""Derived spans — introns, trans-splicing junctions and intergenic spacers (P2B-01, D-g).

The definition lives in ``core/spans.py`` because two modules need the same boundaries:
``introns`` publishes the spans, ``diversity`` measures π over them. One derivation, so the two
cannot disagree about where a boundary is — which is the whole reason D-g fixed the IGS
definition instead of leaving it to the reading module.

The synthetic plastome below is in the normalized frame with every case placed on purpose: an
intron inside one region, a junction crossing a boundary, a spacer bounded by a tRNA, a
pseudogene that must not split a spacer, an overlapping pair, and an IR whose second copy has to
come back once (and only once).
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from plastid.core.ir import IRResult, detect_ir, enforce_ir_gate, ir_from_annotation
from plastid.core.model import Gene, Genome, Interval
from plastid.core.normalize import rotate_genome
from plastid.core.spans import (
    IntronSet,
    SpacerSet,
    introns_of,
    spacers_of,
    trna_introns_of,
)

# LSC [0,30) → IRa [30,50) → SSC [50,70) → IRb [70,90). The LSC longer than the SSC is what
# makes this the canonical frame (`lsc_origin` picks the longer single-copy gap), so normalize()
# leaves it alone and every coordinate below is the one the derivation sees.
LENGTH = 90
IRA = Interval(30, 50)
IRB = Interval(70, 90)
MIRROR = IRA.start + IRB.end - 1  # 119


def _ir() -> IRResult:
    return IRResult(
        ira=IRA,
        irb=IRB,
        identity=1.0,
        mismatches=0,
        compared_length=IRA.length,
        exact_core=IRA.length,
        mirror_constant=MIRROR,
        source="manual",
    )


def _gene(
    name: str, kind: str, parts: tuple[tuple[int, int], ...], *, pseudo: bool = False
) -> Gene:
    return Gene(name=name, kind=kind, parts=tuple(Interval(*part) for part in parts), pseudo=pseudo)


def _synthetic() -> Genome:
    return Genome(
        name="synthetic",
        sequence="ACGT" * (LENGTH // 4) + "AC",
        genes=(
            _gene("psaA", "CDS", ((1, 5), (8, 12))),  # intron [5,8), LSC
            _gene("psaB", "CDS", ((14, 18),)),
            _gene("ycf1p", "CDS", ((19, 20),), pseudo=True),  # not a called gene: no anchor
            _gene("rps12", "CDS", ((22, 24), (33, 36))),  # gap [24,33) crosses LSC→IRa
            _gene("rpl2", "CDS", ((36, 39), (41, 44))),  # intron [39,41) in IRa
            _gene("trnI", "tRNA", ((46, 49),)),
            _gene("ndhE", "CDS", ((52, 57),)),
            _gene("ndhF", "CDS", ((56, 61),)),  # overlaps ndhE — no spacer between them
            _gene("ndhG", "CDS", ((64, 68),)),
            _gene("trnI", "tRNA", ((71, 74),)),  # the IRb copy of trnI
            _gene("rpl2", "CDS", ((76, 79), (81, 84))),  # the IRb copy of rpl2
        ),
    )


def test_introns_are_the_gaps_between_parts_within_one_region() -> None:
    genome = _synthetic()
    derived = introns_of(genome, _ir())

    assert [
        (span.intron_id, span.gene, span.start, span.end, span.phase, span.region)
        for span in derived.introns
    ] == [
        ("psaA.1.1", "psaA", 5, 8, 1, "LSC"),  # 4 bases of exon before it, mod three
        ("rpl2.1.1", "rpl2", 39, 41, 0, "IRa"),
        ("rpl2.2.1", "rpl2", 79, 81, 0, "IRb"),
    ]


def test_a_gap_crossing_a_region_boundary_is_a_junction_not_an_intron() -> None:
    derived = introns_of(_synthetic(), _ir())

    assert [
        (span.junction_id, span.gene, span.start, span.end, span.regions)
        for span in derived.junctions
    ] == [("rps12.1.1", "rps12", 24, 33, ("LSC", "IRa"))]
    # and it is not in the intron list: 17 introns on the golden, junctions never counted
    assert "rps12.1.1" not in {span.intron_id for span in derived.introns}


def test_an_intron_inside_the_ir_keeps_its_own_region() -> None:
    # both IR copies are published, as the `introns` module publishes them — the region column
    # is what tells a reader the two rows are the same locus in mirrored copies
    derived = introns_of(_synthetic(), _ir())

    assert [(span.intron_id, span.region) for span in derived.introns] == [
        ("psaA.1.1", "LSC"),
        ("rpl2.1.1", "IRa"),
        ("rpl2.2.1", "IRb"),
    ]


def test_spacers_are_the_gaps_between_called_genes_inside_one_region() -> None:
    derived = spacers_of(_synthetic(), _ir())

    assert [
        (span.igs_id, span.left, span.right, span.start, span.end, span.region)
        for span in derived.igs
    ] == [
        ("psaA-psaB", "psaA", "psaB", 12, 14, "LSC"),
        # ycf1p is a pseudogene: it does not split this spacer, a called gene would
        ("psaB-rps12", "psaB", "rps12", 18, 22, "LSC"),
        ("ndhF-ndhG", "ndhF", "ndhG", 61, 64, "SSC"),
    ]


def test_a_spacer_inside_the_ir_is_named_and_not_published() -> None:
    """D-g excludes the IR: its two copies are one sequence annotated twice.

    Measured on the golden 10/10/2026 — the copies' *spacers* do not always mirror exactly (3 of
    the 18 pairs: ycf1's truncated copy carries no CDS, two rRNA spacers differ by a base), so a
    spacer cannot be deduplicated the way a *gene* can. Excluded, counted, and named instead.
    """
    derived = spacers_of(_synthetic(), _ir())

    assert [(span.igs_id, span.start, span.end, span.region) for span in derived.in_ir] == [
        ("rpl2-trnI", 44, 46, "IRa"),
        ("trnI-rpl2", 74, 76, "IRb"),
    ]
    assert not any(span.region.startswith("IR") for span in derived.igs)


def test_a_spacer_is_bounded_by_exons_not_by_a_gene_span() -> None:
    """A trans-spliced gene occupies its exons: the space between them is not its territory.

    Taking the outer span instead loses real spacers — measured on the golden where rps12's LSC
    exon and its IR exons straddle clpP, and matK sits inside trnK's intron.
    """
    # rps12's two exons straddle clpP, as the golden's exon in the LSC straddles it: an outer
    # span would swallow clpP and lose both spacers, and the gap between the exons is rps12's
    # own junction — not a spacer either
    genome = Genome(
        name="split",
        sequence="ACGT" * (LENGTH // 4) + "AC",
        genes=(
            _gene("rps12", "CDS", ((2, 4), (24, 26))),
            _gene("clpP", "CDS", ((10, 14),)),
        ),
    )
    derived = spacers_of(genome, _ir())
    assert [(span.igs_id, span.start, span.end) for span in derived.igs] == [
        ("rps12-clpP", 4, 10),
        ("clpP-rps12", 14, 24),
    ]


def test_gaps_that_are_not_spacers_are_counted_not_dropped() -> None:
    derived = spacers_of(_synthetic(), _ir())

    # three gaps cross a region boundary: IRa→SSC, SSC→IRb, and the one wrapping the origin
    assert [(gap.left, gap.right) for gap in derived.crossing] == [
        ("trnI", "ndhE"),
        ("ndhG", "trnI"),
        ("rpl2", "psaA"),
    ]
    # rps12 abuts rpl2 (0 bp) and ndhE overlaps ndhF (−1 bp): neither leaves a spacer
    assert [(gap.left, gap.right, gap.length) for gap in derived.no_spacer] == [
        ("rps12", "rpl2", 0),
        ("ndhE", "ndhF", -1),
    ]


def test_a_gap_that_is_not_positive_is_refused_not_published() -> None:
    """Parts that overlap or nest leave a negative "gap" — published as a length until P2B-01.

    The check it replaces ("exons + gaps tile the span") was algebraic and could not fail for any
    pair of parts. Measured 10/10/2026 over 131 available records: 2,850 multi-part
    genes, exactly one with nested parts — OX459128.1, a stand-in mitogenome.
    """
    broken = Genome(
        name="broken",
        sequence="ACGT" * (LENGTH // 4) + "AC",
        genes=(_gene("psaA", "CDS", ((36, 39), (38, 41))),),  # the second part starts inside
    )
    with pytest.raises(ValueError, match=r"psaA.*gap of -1 bp"):
        introns_of(broken, _ir())


def test_region_labels_do_not_move_when_the_origin_rotates(golden_genome: Genome) -> None:
    """The record's features are in its own frame; region names live in the normalized one.

    Measured 10/10/2026: the golden's rotation offset is 0, so a lookup that compares the two
    frames directly is right by accident there — and wrong on the first record that starts
    anywhere else. Rotating the golden moves every coordinate and must move no label.
    """
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    straight = introns_of(golden_genome, ir)

    shifted_genome = rotate_genome(golden_genome, 7_500)
    shifted_ir = detect_ir(shifted_genome)
    enforce_ir_gate(shifted_ir)
    shifted = introns_of(shifted_genome, shifted_ir)

    assert [span.region for span in shifted.introns] == [span.region for span in straight.introns]
    assert [span.intron_id for span in shifted.introns] == [
        span.intron_id for span in straight.introns
    ]
    assert len(shifted.introns) == len(straight.introns) == 17


def test_an_exon_the_rotation_splits_is_not_an_intron(golden_genome: Genome) -> None:
    """A part cut by the record's origin is one exon, and the two pieces are not an intron.

    Measured 28/10/2026: rotating the golden by half puts the origin inside `petD`'s second exon,
    which `rotate_gene` writes as two touching parts (153451..153724 and 0..208). The intron
    caller paired them and published `petD.1.2`, 153,243 bp — twice a real intron's ceiling and
    most of the plastome. The pieces are contiguous on the circle (one ends at the last base, the
    other starts at the first), which is also exactly how a circular record writes an exon that
    crosses position 1, so the rule covers the deposited case too.
    """
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    straight = introns_of(golden_genome, ir)

    shifted_genome = rotate_genome(golden_genome, len(golden_genome) // 2)
    shifted_ir = detect_ir(shifted_genome)
    enforce_ir_gate(shifted_ir)
    shifted = introns_of(shifted_genome, shifted_ir)

    assert [(span.intron_id, span.length) for span in shifted.introns] == [
        (span.intron_id, span.length) for span in straight.introns
    ]
    assert len(shifted.introns) == 17
    assert max(span.length for span in shifted.introns) < len(golden_genome) // 2


def test_the_spacer_set_is_deterministic_and_type_stable() -> None:
    """The two dataclasses are what the modules pass around — frozen, so nobody edits a span."""
    derived: SpacerSet = spacers_of(_synthetic(), _ir())
    introns: IntronSet = introns_of(_synthetic(), _ir())
    assert isinstance(derived.igs, tuple) and isinstance(introns.introns, tuple)
    assert isinstance(derived.crossing, tuple) and isinstance(introns.junctions, tuple)
    # a span knows its own length in bases, half-open [start, end)
    assert introns.introns[0].length == 3
    assert derived.igs[0].length == 2
    assert derived.igs[0].to_1based() == (13, 14)


# ---------------------------------------------------------------- tRNA introns (M11.15)


def _with_trna(gene: Gene) -> Genome:
    return replace(_synthetic(), genes=(*_synthetic().genes, gene))


def test_a_trna_intron_is_the_gap_inside_one_region() -> None:
    """H-23: the tRNA exon structure is an intron and is published with its kind stated."""
    genome = _with_trna(_gene("trnK", "tRNA", ((6, 9), (12, 15))))  # gap [9,12), LSC

    derived = introns_of(genome, _ir())

    assert [span.intron_id for span in derived.trna_introns] == ["trnK.1.1"]
    span = derived.trna_introns[0]
    assert (span.start, span.end, span.length, span.kind, span.region) == (9, 12, 3, "tRNA", "LSC")
    assert derived.unassessed == ()


def test_a_trna_that_wraps_the_origin_is_not_an_intron() -> None:
    """A circular join across the origin has most of the genome between its parts."""
    genome = replace(
        _synthetic(),
        genes=(*_synthetic().genes, _gene("trnW", "tRNA", ((85, 90), (0, 5)))),
    )

    derived = introns_of(genome, _ir())

    assert derived.trna_introns == ()
    assert derived.unassessed == ()


def test_a_gap_past_the_callers_maximum_is_not_assessed_not_absent() -> None:
    """H-22: ARAGORN documents 3,000 bp. A longer gap is named, never counted as an intron."""
    genome = Genome(
        name="long-gap",
        sequence="ACGT" * 5_000,
        # 3,151 - 150 = 3,001 bp — one base past ARAGORN's documented maximum
        genes=(_gene("trnX", "tRNA", ((100, 150), (3_151, 3_201))),),
    )

    introns, unassessed = trna_introns_of(genome, None)

    assert introns == ()
    assert len(unassessed) == 1
    assert "3,001 bp" in unassessed[0] and "3,000 bp maximum" in unassessed[0]


def test_a_trna_gap_across_a_region_boundary_is_named_not_counted() -> None:
    genome = replace(
        _synthetic(),
        genes=(*_synthetic().genes, _gene("trnZ", "tRNA", ((22, 24), (33, 36)))),
    )

    derived = introns_of(genome, _ir())

    assert [span.intron_id for span in derived.trna_introns] == []
    assert any("cross-region gap" in reason for reason in derived.unassessed)


def test_the_cds_set_diversity_reads_is_unchanged_by_the_trna_set() -> None:
    """``diversity`` reads ``introns``; adding tRNA exons there would move its π baseline."""
    genome = _with_trna(_gene("trnK", "tRNA", ((6, 9), (12, 15))))

    derived = introns_of(genome, _ir())

    assert all(span.kind == "CDS" for span in derived.introns)
    assert [span.intron_id for span in derived.introns] == ["psaA.1.1", "rpl2.1.1", "rpl2.2.1"]
