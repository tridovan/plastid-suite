"""D5-03 for the four table-driven figures: the series equals the cells it read.

Each fixture is the smallest table set the figure's own module emits — the columns are the real
ones, taken from the module sources, so a figure that read the wrong column fails here rather
than in a picture nobody checked.
"""

from __future__ import annotations

import pytest

from plastid.result.figures import (
    barcode_assessment,
    identity_synteny,
    rscu_heatmap,
    ssr_distribution,
)
from plastid.result.result import Column, Result, Scope, Table


def _scope(label: str = "whole plastome") -> Scope:
    return Scope(label=label)


def _rscu_result() -> Result:
    rscu = Table(
        id="rscu",
        title="RSCU",
        columns=(
            Column("codon", "Codon"),
            Column("amino_acid", "AA"),
            Column("count", "Count", dtype="int"),
            Column("rscu", "RSCU", dtype="float", precision=4),
        ),
        rows=(
            {"codon": "GCA", "amino_acid": "A", "count": 370, "rscu": 1.0655},
            {"codon": "GCC", "amino_acid": "A", "count": 236, "rscu": 0.6796},
            {"codon": "TGG", "amino_acid": "W", "count": 100, "rscu": 1.0},
        ),
        scope=_scope(),
        key="codon",
    )
    summary = Table(
        id="codon_usage",
        title="Codon usage",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS", dtype="int"),
            Column("enc", "ENC", dtype="float", precision=2),
        ),
        rows=(
            {"accession": "NC_050833.1", "organism": "Iris domestica", "n_cds": 86, "enc": 49.69},
        ),
        scope=_scope(),
    )
    return Result(
        module="codon_usage",
        module_version="0.1.0",
        tables={"rscu": rscu, "codon_usage": summary},
    )


def _ssr_result() -> Result:
    by_region = Table(
        id="ssr_by_region",
        title="SSR by region",
        columns=(
            Column("region", "Region"),
            Column("total", "Total", dtype="int"),
            Column("mono", "Mono", dtype="int"),
            Column("di", "Di", dtype="int"),
            Column("tri", "Tri", dtype="int"),
            Column("other", "Other", dtype="int"),
        ),
        rows=(
            {"region": "LSC", "total": 24, "mono": 23, "di": 1, "tri": 0, "other": 0},
            {"region": "IRa", "total": 2, "mono": 0, "di": 2, "tri": 0, "other": 0},
        ),
        scope=_scope(),
        key="region",
    )
    ssrs = Table(
        id="ssrs",
        title="SSRs",
        columns=(
            Column("motif", "Motif"),
            Column("unit", "Unit"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("end_1b", "End", dtype="int", unit="bp"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("region", "Region"),
            Column("context", "Context"),
            Column("compound", "Compound"),
        ),
        rows=(
            {
                "motif": "(T)10",
                "unit": "mono",
                "start_1b": 11630,
                "end_1b": 11639,
                "length_bp": 10,
                "region": "LSC",
                "context": "IGS",
                "compound": "no",
            },
            {
                "motif": "(AT)6",
                "unit": "di",
                "start_1b": 40000,
                "end_1b": 40011,
                "length_bp": 12,
                "region": "LSC",
                "context": "Gene",
                "compound": "no",
            },
            {
                "motif": "(T)10",
                "unit": "mono",
                "start_1b": 90000,
                "end_1b": 90009,
                "length_bp": 10,
                "region": "IRa",
                "context": "Intron",
                "compound": "yes",
            },
        ),
        scope=_scope(),
    )
    summary = Table(
        id="ssr_summary",
        title="SSR",
        columns=(
            Column("accession", "Accession"),
            Column("total", "Total", dtype="int"),
            Column("mono", "Mono", dtype="int"),
            Column("di", "Di", dtype="int"),
            Column("tri", "Tri", dtype="int"),
            Column("compound", "Compound", dtype="int"),
            Column("tool", "Tool"),
            Column("thresholds", "Thresholds"),
            Column("reconciliation", "Reconciliation"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "total": 3,
                "mono": 2,
                "di": 1,
                "tri": 0,
                "compound": 1,
                "tool": "MISA 2.1",
                "thresholds": "(1/10) (2/6) (3/5)",
                "reconciliation": "MISA and cpstools agree",
            },
        ),
        scope=_scope(),
    )
    return Result(
        module="ssr",
        module_version="0.1.0",
        tables={"ssr_by_region": by_region, "ssrs": ssrs, "ssr_summary": summary},
    )


def _synteny_result() -> Result:
    identity = Table(
        id="synteny_identity",
        title="Identity series",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("window_start", "Start", dtype="int", unit="bp"),
            Column("window_end", "End", dtype="int", unit="bp"),
            Column("compared_bp", "Compared", dtype="int", unit="bp"),
            Column("identity", "% identity", dtype="float", unit="%", precision=2),
        ),
        rows=(
            {
                "genome": "NC_024936.1",
                "group": "Iris",
                "window_start": 1,
                "window_end": 600,
                "compared_bp": 519,
                "identity": 99.42,
            },
            {
                "genome": "NC_024936.1",
                "group": "Iris",
                "window_start": 601,
                "window_end": 1200,
                "compared_bp": 0,
                "identity": None,
            },
            {
                "genome": "NC_050833.1",
                "group": "Iris",
                "window_start": 1,
                "window_end": 600,
                "compared_bp": 600,
                "identity": 97.0,
            },
        ),
        scope=_scope("normalized frame"),
    )
    blocks = Table(
        id="synteny_blocks",
        title="Blocks",
        columns=(
            Column("block", "Block"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("target_start", "Target start", dtype="int", unit="bp"),
            Column("target_end", "Target end", dtype="int", unit="bp"),
            Column("identity", "%IDY", dtype="float", unit="%", precision=2),
            Column("orientation", "Orientation"),
        ),
        rows=(
            {
                "block": "NC_024936.1#1",
                "genome": "NC_024936.1",
                "group": "Iris",
                "target_start": 1,
                "target_end": 83127,
                "identity": 98.39,
                "orientation": "forward",
            },
            {
                "block": "NC_050833.1#1",
                "genome": "NC_050833.1",
                "group": "Iris",
                "target_start": 100,
                "target_end": 900,
                "identity": 95.0,
                "orientation": "reverse",
            },
        ),
        scope=_scope("normalized frame"),
        key="block",
    )
    summary = Table(
        id="synteny_summary",
        title="Synteny",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("reference_length", "Reference", dtype="int", unit="bp"),
            Column("n_genomes", "Genomes", dtype="int"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("tool", "Tool"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "reference_length": 153724,
                "n_genomes": 2,
                "window_bp": 600,
                "step_bp": 200,
                "tool": "MUMmer nucmer",
            },
        ),
        scope=_scope("normalized frame"),
    )
    return Result(
        module="synteny",
        module_version="0.1.0",
        tables={
            "synteny_identity": identity,
            "synteny_blocks": blocks,
            "synteny_summary": summary,
        },
    )


def _barcodes_result() -> Result:
    windows = Table(
        id="barcodes_windows",
        title="Windows",
        columns=(
            Column("window", "Window", dtype="int"),
            Column("alignment_start", "Start", dtype="int"),
            Column("alignment_end", "End", dtype="int"),
            Column("reference_start", "Reference start", dtype="int", unit="bp"),
            Column("reference_end", "Reference end", dtype="int", unit="bp"),
            Column("usable_columns", "Usable", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
            Column("rank", "Rank", dtype="int"),
            Column("region", "Region"),
        ),
        rows=(
            {
                "window": 1,
                "alignment_start": 1,
                "alignment_end": 600,
                "reference_start": 1,
                "reference_end": 429,
                "usable_columns": 243,
                "pi": 0.198129,
                "rank": 2,
                "region": "",
            },
            {
                "window": 2,
                "alignment_start": 201,
                "alignment_end": 800,
                "reference_start": 176,
                "reference_end": 624,
                "usable_columns": 298,
                "pi": None,
                "rank": None,
                "region": "R1",
            },
        ),
        scope=_scope("normalized frame"),
    )
    candidates = Table(
        id="barcodes_candidates",
        title="Candidates",
        columns=(
            Column("locus", "Locus"),
            Column("kind", "Kind"),
            Column("members", "Members"),
            Column("found", "Found"),
            Column("reference_start", "Reference start", dtype="int", unit="bp"),
            Column("reference_end", "Reference end", dtype="int", unit="bp"),
            Column("span_bp", "Span", dtype="int", unit="bp"),
            Column("pi", "π", dtype="float", precision=6),
            Column("rank", "Rank", dtype="int"),
            Column("note", "Note"),
        ),
        rows=(
            {
                "locus": "rbcL",
                "kind": "gene",
                "members": "rbcL",
                "found": "yes",
                "reference_start": 54991,
                "reference_end": 56436,
                "span_bp": 1446,
                "pi": 0.13755,
                "rank": 8,
                "note": "",
            },
            {
                "locus": "rps16",
                "kind": "gene",
                "members": "rps16",
                "found": "yes",
                "reference_start": 6069,
                "reference_end": 5189,
                "span_bp": 258,
                "pi": 0.1,
                "rank": 7,
                "note": "",
            },
            {
                "locus": "trnH-psbA",
                "kind": "intergenic",
                "members": "trnH-GUG, psbA",
                "found": "yes",
                "reference_start": 153273,
                "reference_end": 71,
                "span_bp": 335,
                "pi": 0.2,
                "rank": 3,
                "note": "spans the origin: 153273–153724 and 1–71",
            },
            {
                "locus": "ITS",
                "kind": "nuclear",
                "members": "ITS",
                "found": "no",
                "reference_start": None,
                "reference_end": None,
                "span_bp": None,
                "pi": None,
                "rank": None,
                "note": "not a plastid locus",
            },
        ),
        scope=_scope("normalized frame"),
        key="locus",
    )
    summary = Table(
        id="barcodes_summary",
        title="Barcodes",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_genomes", "Genomes", dtype="int"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("n_windows", "Windows", dtype="int"),
            Column("mean_pi", "Mean π", dtype="float", precision=6),
            Column("threshold_pi", "Threshold π", dtype="float", precision=6),
            Column("threshold_rule", "Threshold rule"),
            Column("n_regions", "Regions", dtype="int"),
            Column("n_candidates", "Candidates", dtype="int"),
            Column("n_found", "Found", dtype="int"),
            Column("tool", "Tool"),
        ),
        rows=(
            {
                "accession": "NC_050833.1",
                "organism": "Iris domestica",
                "n_genomes": 2,
                "window_bp": 600,
                "step_bp": 200,
                "n_windows": 773,
                "mean_pi": 0.007332,
                "threshold_pi": 0.014157,
                "threshold_rule": "mean + 1 sd of this run's windows",
                "n_regions": 26,
                "n_candidates": 10,
                "n_found": 8,
                "tool": "MAFFT",
            },
        ),
        scope=_scope("normalized frame"),
    )
    return Result(
        module="barcodes",
        module_version="0.2.0",
        tables={
            "barcodes_windows": windows,
            "barcodes_candidates": candidates,
            "barcodes_summary": summary,
        },
    )


TABLES_READ = {
    "rscu_heatmap": lambda: rscu_heatmap.series({"codon_usage": _rscu_result()}),
    "ssr_distribution": lambda: ssr_distribution.series({"ssr": _ssr_result()}),
    "identity_synteny": lambda: identity_synteny.series({"synteny": _synteny_result()}),
    "barcode_assessment": lambda: barcode_assessment.series({"barcodes": _barcodes_result()}),
}


def test_the_rscu_series_is_the_rscu_table_row_for_row() -> None:
    result = _rscu_result()
    data = rscu_heatmap.series({"codon_usage": result})
    table = result.tables["rscu"]
    assert data.codons == tuple(
        (row["amino_acid"], row["codon"], row["count"], row["rscu"]) for row in table.rows
    )
    assert data.families[0][0] == "A" and len(data.families[0][1]) == 2
    assert data.enc == result.tables["codon_usage"].format("enc")
    assert data.n_cds == result.tables["codon_usage"].format("n_cds")


def test_the_ssr_series_is_the_ssr_tables_and_keeps_a_blank_skew_blank() -> None:
    result = _ssr_result()
    data = ssr_distribution.series({"ssr": result})
    by_region = result.tables["ssr_by_region"]
    assert [(entry.region, entry.total) for entry in data.regions] == [
        (row["region"], row["total"]) for row in by_region.rows
    ]
    # (unit column, count, the cell's own rendering) — the count is the cell, not a recount
    assert data.regions[0].segments[0] == ("mono", 23, "23")
    ssrs = result.tables["ssrs"]
    assert [(entry.unit, entry.start, entry.end) for entry in data.spans] == [
        (row["unit"], row["start_1b"], row["end_1b"]) for row in ssrs.rows
    ]
    assert data.units == ("mono", "di")
    assert data.thresholds == result.tables["ssr_summary"].value("thresholds")


def test_the_identity_series_keeps_one_entry_per_genome_and_a_blank_window_blank() -> None:
    result = _synteny_result()
    data = identity_synteny.series({"synteny": result})
    table = result.tables["synteny_identity"]
    flattened = [
        (genome, group, window) for genome, group, windows in data.per_genome for window in windows
    ]
    assert flattened == [
        (row["genome"], row["group"], (row["window_start"], row["window_end"], row["identity"]))
        for row in table.rows
    ]
    assert any(window[2] is None for _, _, windows in data.per_genome for window in windows)
    blocks = result.tables["synteny_blocks"]
    assert [(block[0], block[1], block[2], block[3]) for block in data.blocks] == [
        (row["genome"], row["target_start"], row["target_end"], row["orientation"])
        for row in blocks.rows
    ]


def test_the_barcode_series_is_the_barcode_tables() -> None:
    result = _barcodes_result()
    data = barcode_assessment.series({"barcodes": result})
    windows = result.tables["barcodes_windows"]
    assert list(data.windows) == [
        (row["reference_start"], row["reference_end"], row["pi"]) for row in windows.rows
    ]
    candidates = result.tables["barcodes_candidates"]
    assert [(c[0], c[1], c[2], c[3], c[4]) for c in data.candidates] == [
        (
            row["locus"],
            row["found"],
            row["reference_start"],
            row["reference_end"],
            row["note"],
        )
        for row in candidates.rows
    ]
    assert data.threshold_pi == result.tables["barcodes_summary"].value("threshold_pi")
    assert data.threshold_label == result.tables["barcodes_summary"].format("threshold_pi")


def test_a_span_running_backwards_is_the_locus_footprint_and_only_a_note_makes_it_a_wrap() -> None:
    """Measured: `rps16` runs backwards with an empty note; only the module's wording wraps."""
    limits = (1, 153_724)
    assert barcode_assessment._bands(6069, 5189, "", limits) == ((5189, 6069),)
    assert barcode_assessment._bands(
        153_273, 71, "spans the origin: 153273–153724 and 1–71", limits
    ) == (
        (153_273, 153_724),
        (1, 71),
    )
    assert barcode_assessment._bands(54991, 56436, "", limits) == ((54991, 56436),)


@pytest.mark.parametrize(
    "module",
    [rscu_heatmap, ssr_distribution, identity_synteny, barcode_assessment],
    ids=lambda module: module.SPEC.id,
)
def test_a_missing_module_is_named_rather_than_drawn_empty(module) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError, match=module.SPEC.requires[0][0]):
        module.series({})
