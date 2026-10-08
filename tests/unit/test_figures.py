"""The genome-scale figures: circular map and IR junctions, and the frame they must cross.

D5-03 is the rule under test throughout: a figure reads table cells and recomputes nothing, so a
test can compare the series with the table row by row — and it runs without matplotlib, because
the series is plain data.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.result.export import export_json
from plastid.result.figures import all_figures, circular_map, ir_junction
from plastid.result.figures.frames import Frame
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.figures import render_run_figures

LENGTH = 20_000
ROTATION = 0


def _stats_result(rotation: int = 0) -> tuple[Result, Result]:
    """A small plastome in both frames.

    LSC 1–8000 · IRa 8001–12000 · SSC 12001–16000 · IRb 16001–20000, so a ±3 kb junction
    window is a window and not the whole genome.
    """
    stats = Table(
        id="genome_stats",
        title="Genome statistics",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
            Column("n_genes_total", "Genes", dtype="int"),
        ),
        rows=(
            {
                "organism": "Testus plastoma",
                "accession": "NC_000001.1",
                "length_bp": LENGTH,
                "gc_pct": 37.85,
                "n_genes_total": 4,
            },
        ),
        scope=Scope(label="whole plastome; both IR copies"),
    )
    genes = Table(
        id="genes",
        title="Genes",
        columns=(
            Column("name", "Name"),
            Column("kind", "Kind"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("end_1b", "End", dtype="int", unit="bp"),
            Column("strand", "Strand"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("in_ir", "In IR"),
            Column("pseudo", "Pseudo"),
        ),
        rows=(
            {
                "name": "psbA",
                "kind": "CDS",
                "start_1b": 1,
                "end_1b": 1200,
                "strand": "+",
                "length_bp": 1200,
                "in_ir": "no",
                "pseudo": "no",
            },
            {
                "name": "rps3",
                "kind": "CDS",
                "start_1b": 7800,
                "end_1b": 8200,
                "strand": "-",
                "length_bp": 401,
                "in_ir": "yes",
                "pseudo": "no",
            },
            {
                "name": "trnH-GUG",
                "kind": "tRNA",
                "start_1b": 14_000,
                "end_1b": 14_060,
                "strand": "+",
                "length_bp": 61,
                "in_ir": "no",
                "pseudo": "no",
            },
            {
                "name": "rrn16",
                "kind": "rRNA",
                "start_1b": 16_200,
                "end_1b": 16_300,
                "strand": "-",
                "length_bp": 101,
                "in_ir": "yes",
                "pseudo": "no",
            },
            # A `gene` feature row — the diagnostic row the map must not draw as a second gene.
            {
                "name": "rrn16",
                "kind": "gene",
                "start_1b": 16_200,
                "end_1b": 16_300,
                "strand": "-",
                "length_bp": 101,
                "in_ir": "yes",
                "pseudo": "no",
            },
        ),
        scope=Scope(label="GenBank coordinates, every copy"),
    )
    structure_summary = Table(
        id="structure_summary",
        title="Structure",
        columns=(
            Column("accession", "Accession"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("n_windows", "Windows", dtype="int"),
            Column("n_undefined_windows", "Undefined windows", dtype="int"),
            Column("rotation_offset", "Rotation offset", dtype="int", unit="bp"),
            Column("coord_system", "Coordinates"),
        ),
        rows=(
            {
                "accession": "NC_000001.1",
                "window_bp": 1000,
                "step_bp": 1000,
                "n_windows": 20,
                "n_undefined_windows": 0,
                "rotation_offset": rotation,
                "coord_system": "normalized",
            },
        ),
        scope=Scope(label="normalized frame", coord_system="normalized"),
    )
    regions = Table(
        id="regions",
        title="Regions",
        columns=(
            Column("region", "Region"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("end_1b", "End", dtype="int", unit="bp"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
            Column("gc_skew", "GC skew", dtype="float", precision=4),
        ),
        rows=(
            {
                "region": "LSC",
                "start_1b": 1,
                "end_1b": 8000,
                "length_bp": 8000,
                "gc_pct": 36.0,
                "gc_skew": -0.02,
            },
            {
                "region": "IRa",
                "start_1b": 8001,
                "end_1b": 12_000,
                "length_bp": 4000,
                "gc_pct": 40.0,
                "gc_skew": 0.01,
            },
            {
                "region": "SSC",
                "start_1b": 12_001,
                "end_1b": 16_000,
                "length_bp": 4000,
                "gc_pct": 34.0,
                "gc_skew": -0.03,
            },
            {
                "region": "IRb",
                "start_1b": 16_001,
                "end_1b": 20_000,
                "length_bp": 4000,
                "gc_pct": 40.0,
                "gc_skew": 0.01,
            },
        ),
        scope=Scope(label="normalized frame", coord_system="normalized"),
        key="region",
    )
    skew = Table(
        id="gc_skew_windows",
        title="GC skew",
        columns=(
            Column("region", "Region"),
            Column("window_start_1b", "Start", dtype="int", unit="bp"),
            Column("window_end_1b", "End", dtype="int", unit="bp"),
            Column("window_length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_skew", "GC skew", dtype="float", precision=4),
            Column("cumulative_gc_skew", "Cumulative", dtype="float", precision=4),
        ),
        rows=(
            {
                "region": "LSC",
                "window_start_1b": 1,
                "window_end_1b": 1000,
                "window_length_bp": 1000,
                "gc_skew": -0.05,
                "cumulative_gc_skew": -0.05,
            },
            {
                "region": "LSC",
                "window_start_1b": 1001,
                "window_end_1b": 2000,
                "window_length_bp": 1000,
                "gc_skew": None,
                "cumulative_gc_skew": None,
            },
            {
                "region": "SSC",
                "window_start_1b": 12_001,
                "window_end_1b": 13_000,
                "window_length_bp": 1000,
                "gc_skew": 0.04,
                "cumulative_gc_skew": 0.04,
            },
        ),
        scope=Scope(label="normalized frame", coord_system="normalized"),
    )
    return (
        Result(
            module="genome_stats",
            module_version="0.1.0",
            tables={"genome_stats": stats, "genes": genes},
        ),
        Result(
            module="structure",
            module_version="0.2.0",
            tables={
                "structure_summary": structure_summary,
                "regions": regions,
                "gc_skew_windows": skew,
            },
        ),
    )


def _results(rotation: int = 0) -> dict[str, Result]:
    stats, structure = _stats_result(rotation)
    return {"genome_stats": stats, "structure": structure}


def test_the_frame_maps_a_position_both_ways() -> None:
    frame = Frame(length_bp=LENGTH, rotation_offset=0)
    assert frame.to_normalized(1) == 1 and frame.to_genbank(1) == 1
    shifted = Frame(length_bp=LENGTH, rotation_offset=250)
    # The base that became the normalized position 1 is GenBank position 251.
    assert shifted.to_normalized(251) == 1
    assert shifted.to_genbank(1) == 251
    assert shifted.to_normalized(1) == 19_751


def test_a_span_that_crosses_the_origin_becomes_two_spans() -> None:
    frame = Frame(length_bp=LENGTH, rotation_offset=250)
    assert frame.span_to_genbank(1, 8000) == ((251, 8250),)
    # The GenBank origin (position 1) sits inside this span, so it splits into two arcs.
    assert frame.span_to_genbank(19_700, 19_800) == ((19_950, 20_000), (1, 50))


def test_the_map_series_equals_the_cells_of_the_tables_it_read() -> None:
    """D5-03: every arc the map draws is a cell of `genes`, `regions` or `gc_skew_windows`."""
    results = _results()
    data = circular_map.series(results)
    genes = [row for row in results["genome_stats"].tables["genes"].rows if row["kind"] != "gene"]
    assert [(g.name, g.kind, g.start_1b, g.end_1b, g.strand) for g in data.genes] == [
        (row["name"], row["kind"], row["start_1b"], row["end_1b"], row["strand"]) for row in genes
    ]
    assert data.n_gene_features_skipped == 1
    assert data.length_bp == results["genome_stats"].tables["genome_stats"].value("length_bp")
    assert data.gc_pct == results["genome_stats"].tables["genome_stats"].format("gc_pct")
    assert [(r.region, r.start_1b, r.end_1b) for r in data.regions] == [
        ("LSC", 1, 8000),
        ("IRa", 8001, 12_000),
        ("SSC", 12_001, 16_000),
        ("IRb", 16_001, 20_000),
    ]
    skew = results["structure"].tables["gc_skew_windows"]
    assert [(w.start_1b, w.gc_skew) for w in data.skew] == [
        (row["window_start_1b"], row["gc_skew"]) for row in skew.rows
    ]


def test_the_map_carries_an_undefined_skew_as_a_hole_not_a_zero() -> None:
    data = circular_map.series(_results())
    assert any(window.gc_skew is None for window in data.skew)


def test_a_rotated_record_moves_the_region_arcs_and_says_so() -> None:
    """The golden genome has rotation 0; a figure that skipped the conversion would still pass."""
    zero = circular_map.series(_results(rotation=0))
    shifted = circular_map.series(_results(rotation=250))
    assert [(r.start_1b, r.end_1b) for r in zero.regions][0] == (1, 8000)
    assert [(r.start_1b, r.end_1b) for r in shifted.regions][0] == (251, 8250)
    assert [(r.start_1b, r.end_1b) for r in shifted.regions][1] == (8251, 12_250)
    assert shifted.frame.rotation_offset == 250
    assert zero.frame.rotation_offset == 0


def test_the_junction_series_holds_the_genes_of_each_window() -> None:
    data = ir_junction.series(_results())
    assert [junction.name for junction in data.junctions] == ["LSC / IRa", "SSC / IRb"]
    first, second = data.junctions
    assert first.boundary_1b == 8001 and second.boundary_1b == 16_001
    assert (first.window_start_1b, first.window_end_1b) == (5001, 12_000)
    assert "rps3" in {gene.name for gene in first.genes}
    assert "psbA" not in {gene.name for gene in first.genes}  # 1–1200 is outside the ±3 kb window
    assert {gene.name for gene in second.genes} == {"trnH-GUG", "rrn16"}
    assert data.n_genes == 3


def test_a_gene_that_straddles_the_boundary_is_kept_whole() -> None:
    """The boundary gene is the reason the figure exists — clipping it would hide the finding."""
    data = ir_junction.series(_results())
    rps3 = next(gene for gene in data.junctions[0].genes if gene.name == "rps3")
    assert (rps3.start_1b, rps3.end_1b) == (7800, 8200)


def test_a_figure_that_cannot_find_its_module_says_which_one() -> None:
    with pytest.raises(ValueError, match="structure"):
        circular_map.series({"genome_stats": _results()["genome_stats"]})
    with pytest.raises(ValueError, match="genome_stats"):
        ir_junction.series({"structure": _results()["structure"]})


def test_the_catalogue_holds_every_figure_and_no_untyped_source() -> None:
    """M5.1 lists eleven figure types; eight have a source table and the page is P5-04's ninth.

    A tenth joins in Phase 3B — ``dotplot``, the picture a person reads before confirming an
    inverted SSC (D5-08's confirmation layer, P3B-03). An entry appears in this list only when
    its ``requires`` names real tables, which the loop below checks.
    """
    specs = {spec.id: spec for spec in all_figures()}
    assert sorted(specs) == [
        "barcode_assessment",
        "circular_map",
        "dotplot",
        "genome_page",
        "identity_synteny",
        "ir_junction",
        "pi_window",
        "rscu_heatmap",
        "ssr_distribution",
    ]
    assert specs["genome_page"].formats == ("pdf",), "a page is a document, not a vector asset"
    for spec in specs.values():
        assert spec.requires, f"{spec.id} declares no source table"
        for module, table in spec.requires:
            assert "." not in module and "." not in table, f"{spec.id}: {module}.{table}"


def test_the_run_level_pass_renders_every_figure_the_tables_feed(tmp_path: Path) -> None:
    """Seven figure types are catalogued; this run's tables feed three of them."""
    pytest.importorskip("matplotlib", reason="the 'figures' extra is not installed")
    run = tmp_path / "run"
    for order, result in enumerate(_results().values(), start=1):
        module_dir = run / f"{order:02d}_{result.module}"
        module_dir.mkdir(parents=True)
        (module_dir / "result.json").write_text(export_json(result), encoding="utf-8")

    outcome = render_run_figures(run, formats=("svg",))
    assert sorted(path.name for path in outcome.files) == [
        "circular_map.svg",
        "ir_junction.svg",
    ]
    skipped = dict(outcome.skipped)
    assert "pi_window" in skipped and "diversity" in skipped["pi_window"]
    for path in outcome.files:
        assert path.stat().st_size > 1000
