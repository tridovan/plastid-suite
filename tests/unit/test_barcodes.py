"""Hypervariable-region calling (P4-04) — the rule, on hand-built window series.

The π values themselves come from ``core.align.window_diversity`` (tested in ``test_align.py``
and cross-checked module-to-module in ``tests/golden/test_tier_b.py``); these tests are about the
calling rule: the threshold, the merging, and what a hole in the series does.
"""

from __future__ import annotations

import pytest

from plastid.core.align import WindowDiversity
from plastid.core.barcodes import hypervariable_regions, threshold_for


def _windows(values: list[float | None], *, window: int = 100, step: int = 100):
    return tuple(
        WindowDiversity(
            start_1b=index * step + 1,
            end_1b=index * step + window,
            usable_columns=window if value is not None else 0,
            pi=value,
        )
        for index, value in enumerate(values)
    )


def test_the_threshold_is_mean_plus_k_standard_deviations() -> None:
    windows = _windows([1.0, 2.0, 3.0, 4.0])
    # mean 2.5, sample sd = sqrt(5/3) = 1.29099… → threshold at k = 1 is 3.79099…
    assert threshold_for(windows) == pytest.approx(2.5 + (5 / 3) ** 0.5)


def test_a_series_with_fewer_than_two_values_has_no_threshold() -> None:
    assert threshold_for(_windows([None, None])) is None
    assert threshold_for(_windows([0.5, None])) is None


def test_consecutive_windows_above_the_threshold_merge_into_one_region() -> None:
    windows = _windows([0.0, 5.0, 6.0, 7.0, 0.0])
    (region,) = hypervariable_regions(windows, threshold=4.0)
    assert (region.start_1b, region.end_1b) == (101, 400)
    assert region.n_windows == 3
    assert region.mean_pi == pytest.approx(6.0)
    assert region.max_pi == pytest.approx(7.0)
    assert region.span == 300


def test_a_single_window_is_not_a_region_by_default() -> None:
    windows = _windows([0.0, 9.0, 0.0])
    assert hypervariable_regions(windows, threshold=4.0) == ()
    (region,) = hypervariable_regions(windows, threshold=4.0, min_windows=1)
    assert region.n_windows == 1


def test_a_hole_is_not_a_peak() -> None:
    """A window with no usable column has no value — it cannot join a region, and it breaks one."""
    windows = _windows([5.0, 5.0, None, 5.0, 5.0])
    regions = hypervariable_regions(windows, threshold=4.0)
    assert [(region.start_1b, region.end_1b) for region in regions] == [(1, 200), (301, 500)]


def test_regions_come_back_in_coordinate_order() -> None:
    # three runs exceed the threshold, but the leading one is a single window and is dropped
    windows = _windows([9.0, 0.0, 5.0, 5.0, 0.0, 7.0, 7.0])
    regions = hypervariable_regions(windows, threshold=4.0)
    starts = [region.start_1b for region in regions]
    assert starts == sorted(starts)
    assert [(region.start_1b, region.end_1b) for region in regions] == [(201, 400), (501, 700)]


def test_every_shipped_locus_carries_its_source() -> None:
    """The list is data with provenance, like the alias table — not names somebody remembers."""
    from plastid.core.barcodes import KINDS, load_barcode_loci

    loci = load_barcode_loci()
    assert {locus.name for locus in loci} >= {"rbcL", "matK", "trnH-psbA", "trnL-trnF", "ycf1"}
    assert all(locus.source for locus in loci)
    assert all(locus.kind in KINDS for locus in loci)
    assert all(len(locus.members) == 1 for locus in loci if locus.kind == "gene")
    assert all(len(locus.members) == 2 for locus in loci if locus.kind == "intergenic")


def test_a_locus_file_without_a_source_is_refused(tmp_path) -> None:
    from plastid.core.barcodes import BarcodeLocusError, load_barcode_loci

    bad = tmp_path / "loci.yaml"
    bad.write_text("loci:\n  - name: x\n    kind: gene\n    members: [x]\n", encoding="utf-8")
    with pytest.raises(BarcodeLocusError, match="no source"):
        load_barcode_loci(bad)
    bad.write_text(
        "loci:\n  - name: x\n    kind: strange\n    members: [x]\n    source: s\n",
        encoding="utf-8",
    )
    with pytest.raises(BarcodeLocusError, match="kind"):
        load_barcode_loci(bad)


def _circular(genes: dict[str, tuple[int, int]], length: int = 1000):
    """A synthetic circular genome whose gene coordinates are given 1-based inclusive."""
    from plastid.core.model import Gene, Genome, Interval

    features = tuple(
        Gene(name, "CDS", (Interval(start - 1, end),)) for name, (start, end) in genes.items()
    )
    return Genome(
        name="synthetic circle",
        accession="CIRC.1",
        sequence="ACGT" * (length // 4),
        genes=features,
        topology="circular",
    )


def test_an_intergenic_locus_that_crosses_the_origin_is_found_and_reported_as_two_parts() -> None:
    """Measured on the golden genome: a plastome's record starts wherever the submitter started
    it, so a spacer can run through position 1. The module must find it — reporting "not carried"
    for a locus that is there would be a wrong claim about the data, not a missing value."""
    from plastid.core.barcodes import BarcodeLocus
    from plastid.modules.barcodes import _locate

    genome = _circular({"late": (900, 950), "early": (10, 60), "middle": (500, 550)})
    copies = {"late": [genome.genes[0]], "early": [genome.genes[1]]}
    locus = BarcodeLocus(
        name="late-early", kind="intergenic", members=("late", "early"), source="t"
    )
    span = _locate(locus, copies, genome.genes, genome)
    assert span is not None
    assert span.wrapped is True
    assert span.parts == ((951, 1000), (1, 9))
    assert span.span == 59
    assert (span.start, span.end) == (951, 9)  # start > end is how a wrapped span reads


def test_an_intergenic_locus_is_not_found_when_no_route_between_its_members_is_empty() -> None:
    """On a circle a named pair has two routes; the module takes the direct one, and the other
    only when it is empty too. Here both are blocked, so the locus is genuinely not a spacer."""
    from plastid.core.barcodes import BarcodeLocus
    from plastid.modules.barcodes import _locate

    genome = _circular(
        {"first": (10, 60), "blocker": (100, 150), "second": (200, 250), "elsewhere": (700, 750)}
    )
    copies = {"first": [genome.genes[0]], "second": [genome.genes[2]]}
    locus = BarcodeLocus(
        name="first-second", kind="intergenic", members=("first", "second"), source="t"
    )
    assert _locate(locus, copies, genome.genes, genome) is None


def test_the_direct_route_wins_when_both_routes_are_empty() -> None:
    from plastid.core.barcodes import BarcodeLocus
    from plastid.modules.barcodes import _locate

    genome = _circular({"first": (10, 60), "second": (200, 250)})
    copies = {"first": [genome.genes[0]], "second": [genome.genes[1]]}
    locus = BarcodeLocus(
        name="first-second", kind="intergenic", members=("first", "second"), source="t"
    )
    span = _locate(locus, copies, genome.genes, genome)
    assert span is not None and span.wrapped is False
    assert span.parts == ((61, 199),)


def test_a_rotated_reference_gives_the_same_windows_and_locus_pi(
    tmp_path, monkeypatch, golden_genome
) -> None:
    """Measured 27/10/2026: the module aligned the **raw** reference but located windows and
    loci on the **normalized** one, so on a record that does not start at the LSC every window
    span and candidate π was off by the rotation. The fix puts the normalized record in row 0
    of the alignment, which is the frame the tables claim."""
    from dataclasses import replace

    import plastid.modules.barcodes as barcodes
    from plastid.core.align import Alignment
    from plastid.core.ir import detect_ir, enforce_ir_gate
    from plastid.core.normalize import rotate_genome
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    def fake_align(rows, ctx, *, strategy="auto"):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    monkeypatch.setattr(barcodes, "align", fake_align)

    changes = {"A": "G", "G": "A", "C": "T", "T": "C"}
    sequence = list(golden_genome.sequence)
    for position in range(
        1_000, 1_100
    ):  # one hypervariable window, so the region table is non-empty
        sequence[position] = changes[sequence[position]]
    mutated = replace(golden_genome, sequence="".join(sequence))

    def tables(reference):
        layout = RunLayout.create(tmp_path / "runs", "barcodes", fingerprint="test")
        layout.make_dirs()
        ctx = RunContext(layout=layout, params={}, module_id="barcodes", module_order=1)
        ctx.ir = enforce_ir_gate(detect_ir(reference))
        result = barcodes.run(
            (reference, mutated), {"window": 1000, "step": 1000, "min_windows": 1}, ctx
        )
        windows = [
            (
                row["alignment_start"],
                row["reference_start"],
                row["reference_end"],
                row["pi"],
            )
            for row in result.tables["barcodes_windows"].rows
        ]
        candidates = [
            (row["locus"], row["reference_start"], row["reference_end"], row["pi"])
            for row in result.tables["barcodes_candidates"].rows
        ]
        return windows, candidates

    assert tables(rotate_genome(golden_genome, 7500)) == tables(golden_genome)


def test_a_pseudogene_is_not_a_barcode_locus() -> None:
    """A pseudogene is not evidence that the locus is present (measured 27/10/2026)."""
    from plastid.core.align import Alignment
    from plastid.core.model import Gene, Genome, Interval
    from plastid.modules.barcodes import _candidates

    sequence = "ACGT" * 50
    genome = Genome(
        name="Testus",
        accession="T.1",
        sequence=sequence,
        genes=(Gene(name="rbcL", kind="CDS", parts=(Interval(0, 24),), pseudo=True),),
    )
    rows = _candidates(genome, Alignment(names=("T.1",), sequences=(sequence,)))
    rbc = next(row for row in rows if row["locus"] == "rbcL")
    assert rbc["found"] == "no"
    assert "rbcL" in rbc["note"]


def test_identical_genomes_report_no_regions_instead_of_an_empty_table(
    tmp_path, monkeypatch, golden_genome
) -> None:
    """A legitimate "no hypervariable region" result must not be a schema exception.

    ``Table`` refuses zero rows, so the first version crashed on two identical genomes — the
    commonest negative case there is. The summary still reports the zero and a note says the
    table is absent rather than empty, so "none found" and "did not run" stay distinguishable.
    """
    import plastid.modules.barcodes as barcodes
    from plastid.core.align import Alignment
    from plastid.core.ir import detect_ir, enforce_ir_gate
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    def fake_align(rows, ctx, *, strategy="auto"):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    monkeypatch.setattr(barcodes, "align", fake_align)
    layout = RunLayout.create(tmp_path / "runs", "barcodes", fingerprint="test")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="barcodes", module_order=1)
    ctx.ir = enforce_ir_gate(detect_ir(golden_genome))

    result = barcodes.run((golden_genome, golden_genome), {"window": 1000, "step": 1000}, ctx)

    assert "barcodes_regions" not in result.tables
    assert result.tables["barcodes_summary"].rows[0]["n_regions"] == 0
    assert any("regions table is absent" in note for note in result.notes)
