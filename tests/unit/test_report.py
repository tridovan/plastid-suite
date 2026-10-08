"""Report engine: cell references, the draft banner, citations, and the template contract."""

from __future__ import annotations

from dataclasses import replace

import pytest

from plastid.report.citations import Citation, format_citation, load_citations
from plastid.report.engine import (
    DRAFT_BANNER,
    CitationError,
    TemplateRefError,
    find_bare_digits,
    find_refs,
    list_templates,
    load_template,
    render,
)
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.defaults import load_defaults


def _result() -> Result:
    return Result(
        module="genome_stats",
        module_version="0.1.0",
        tables={
            "genome_stats": Table(
                id="genome_stats",
                title="Genome statistics",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
                    Column("n_genes_total", "Genes", dtype="int"),
                    Column("n_cds", "CDS", dtype="int"),
                    Column("n_trna", "tRNA", dtype="int"),
                    Column("n_rrna", "rRNA", dtype="int"),
                    Column("ir_length_bp", "IR", dtype="int", unit="bp"),
                    Column("ir_identity_pct", "IR identity", dtype="float", unit="%", precision=2),
                    Column("ir_source", "IR source"),
                    Column("lsc_bp", "LSC", dtype="int", unit="bp"),
                    Column("ssc_bp", "SSC", dtype="int", unit="bp"),
                    Column("coord_system", "Coordinates"),
                ),
                rows=(
                    {
                        "organism": "Iris domestica",
                        "accession": "NC_050833.1",
                        "length_bp": 153_724,
                        "gc_pct": 37.8536,
                        "n_genes_total": 132,
                        "n_cds": 86,
                        "n_trna": 38,
                        "n_rrna": 8,
                        "ir_length_bp": 26_214,
                        "ir_identity_pct": 100.0,
                        "ir_source": "annotation",
                        "lsc_bp": 83_127,
                        "ssc_bp": 18_169,
                        "coord_system": "genbank",
                    },
                ),
                scope=Scope(label="whole plastome; both IR copies"),
            ),
            "ssr_by_unit": Table(
                id="ssr_by_unit",
                title="SSR by unit",
                columns=(Column("unit", "Unit"), Column("n", "Count", dtype="int")),
                rows=({"unit": "mono", "n": 31}, {"unit": "di", "n": 4}, {"unit": "tri", "n": 2}),
                scope=Scope(label="whole plastome, MISA defaults"),
                key="unit",
            ),
        },
    )


def test_find_refs_and_render_numbers() -> None:
    template = "Length **{{genome_stats.length_bp}} bp**, mono SSRs {{ssr_by_unit.mono.n}}."
    assert find_refs(template) == ["genome_stats.length_bp", "ssr_by_unit.mono.n"]
    rendered = render(template, _result())
    assert "**153724 bp**" in rendered
    assert "mono SSRs 31." in rendered


def test_render_formats_floats_with_the_column_precision() -> None:
    rendered = render("GC {{genome_stats.gc_pct}} %", _result())
    assert "GC 37.85 %" in rendered  # precision=2 from the column, not the raw float


def test_unknown_reference_raises_and_names_the_ref() -> None:
    with pytest.raises(TemplateRefError, match="genome_stats.nope"):
        render("{{genome_stats.nope}}", _result())


def test_render_appends_the_draft_banner() -> None:
    rendered = render("Body.", _result())
    assert rendered.rstrip().endswith(DRAFT_BANNER)


def test_render_appends_citations_from_the_registry() -> None:
    result = replace(_result(), citations=("misa",))
    rendered = render("Body.", result, citations=load_citations())
    assert "### References" in rendered
    assert "MISA-web" in rendered
    assert DRAFT_BANNER in rendered


def test_unknown_citation_key_raises() -> None:
    result = replace(_result(), citations=("nope",))
    with pytest.raises(CitationError, match="nope"):
        render("Body.", result, citations=load_citations())


def test_bare_digit_scan() -> None:
    assert find_bare_digits("The IR is {{ir.length}} bp long.") == []
    assert find_bare_digits("The IR is 26214 bp.") == ["26214"]


def _ssr_result() -> Result:
    return Result(
        module="ssr",
        module_version="0.1.0",
        tables={
            "ssr_summary": Table(
                id="ssr_summary",
                title="SSR summary",
                columns=(
                    Column("accession", "Accession"),
                    Column("tool", "Tool"),
                    Column("total", "Total SSRs", dtype="int"),
                    Column("mono", "Monomeric", dtype="int"),
                    Column("di", "Dinucleotide", dtype="int"),
                    Column("tri", "Trinucleotide", dtype="int"),
                    Column("compound", "In compound formation", dtype="int"),
                    Column("thresholds", "Unit/repeat thresholds"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "tool": "MISA 2.1",
                        "total": 37,
                        "mono": 31,
                        "di": 4,
                        "tri": 2,
                        "compound": 3,
                        "thresholds": "(1/10) (2/6) (3/5) (4/5) (5/5) (6/5)",
                    },
                ),
                scope=Scope(label="whole plastome; MISA defaults"),
            ),
            "ssr_by_region": Table(
                id="ssr_by_region",
                title="SSRs by plastid region",
                columns=(
                    Column("region", "Region"),
                    Column("total", "Total", dtype="int"),
                ),
                rows=tuple(
                    {"region": name, "total": total}
                    for name, total in (("LSC", 27), ("IRa", 0), ("SSC", 10), ("IRb", 0))
                ),
                scope=Scope(label="whole plastome"),
                key="region",
            ),
        },
    )


def _structure_result() -> Result:
    return Result(
        module="structure",
        module_version="0.1.0",
        tables={
            "structure_summary": Table(
                id="structure_summary",
                title="Structure summary",
                columns=(
                    Column("accession", "Accession"),
                    Column("organism", "Organism"),
                    Column("window_bp", "Window", dtype="int", unit="bp"),
                    Column("step_bp", "Step", dtype="int", unit="bp"),
                    Column("n_windows", "Windows", dtype="int"),
                    Column("rotation_offset", "Rotation offset", dtype="int", unit="bp"),
                    Column("coord_system", "Coordinates"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "window_bp": 1000,
                        "step_bp": 1000,
                        "n_windows": 154,
                        "rotation_offset": 0,
                        "coord_system": "normalized",
                    },
                ),
                scope=Scope(label="one row", coord_system="normalized"),
            ),
            "regions": Table(
                id="regions",
                title="Regions",
                columns=(
                    Column("region", "Region"),
                    Column("start_1b", "Start", dtype="int"),
                    Column("end_1b", "End", dtype="int"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
                    Column("gc_skew", "GC-skew (G−C)/(G+C)", dtype="float", precision=4),
                ),
                rows=tuple(
                    {
                        "region": name,
                        "start_1b": 1,
                        "end_1b": 83_127,
                        "length_bp": 83_127,
                        "gc_pct": 36.0,
                        "gc_skew": -0.01,
                    }
                    for name in ("LSC", "IRa", "SSC", "IRb")
                ),
                scope=Scope(label="normalized frame", coord_system="normalized"),
                key="region",
            ),
        },
    )


def _codon_usage_result() -> Result:
    return Result(
        module="codon_usage",
        module_version="0.1.0",
        tables={
            "codon_usage": Table(
                id="codon_usage",
                title="Codon usage summary",
                columns=(
                    Column("accession", "Accession"),
                    Column("organism", "Organism"),
                    Column("n_cds", "CDS in the set", dtype="int"),
                    Column("n_codons", "Codons counted", dtype="int"),
                    Column("enc", "ENC", dtype="float", precision=2),
                    Column("cai", "CAI", dtype="float", precision=4),
                    Column("cbi", "CBI", dtype="float", precision=4),
                    Column("fop", "Fop", dtype="float", precision=4),
                    Column("internal_stops", "Internal stop codons", dtype="int"),
                    Column("translation_table", "Translation table"),
                    Column("reference", "Reference for CAI/CBI/Fop"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_cds": 86,
                        "n_codons": 24_000,
                        "enc": 51.23,
                        "cai": 0.17,
                        "cbi": 0.02,
                        "fop": 0.45,
                        "internal_stops": 0,
                        "translation_table": "NCBI table 1 (standard)",
                        "reference": "the analysed CDS set itself (self-reference)",
                    },
                ),
                scope=Scope(label="CDS only"),
            ),
        },
    )


def _amino_acid_result() -> Result:
    return Result(
        module="amino_acid",
        module_version="0.1.0",
        tables={
            "amino_acid": Table(
                id="amino_acid",
                title="Amino-acid summary",
                columns=(
                    Column("accession", "Accession"),
                    Column("organism", "Organism"),
                    Column("n_cds", "CDS in the set", dtype="int"),
                    Column("n_residues", "Residues counted", dtype="int"),
                    Column("ignored_symbols", "Stops and ambiguous symbols", dtype="int"),
                    Column("gravy", "GRAVY (Kyte–Doolittle)", dtype="float", precision=4),
                    Column("aromaticity", "Aromaticity (F+W+Y)", dtype="float", precision=4),
                    Column("translation_table", "Translation table"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_cds": 86,
                        "n_residues": 24_000,
                        "ignored_symbols": 86,
                        "gravy": -0.1234,
                        "aromaticity": 0.0987,
                        "translation_table": "NCBI table 1 (standard)",
                    },
                ),
                scope=Scope(label="CDS only"),
            ),
        },
    )


def _introns_result() -> Result:
    return Result(
        module="introns",
        module_version="0.1.0",
        tables={
            "intron_summary": Table(
                id="intron_summary",
                title="Intron summary",
                columns=(
                    Column("accession", "Accession"),
                    Column("organism", "Organism"),
                    Column("n_genes_with_introns", "Genes with introns", dtype="int"),
                    Column("n_introns", "Introns", dtype="int"),
                    Column("total_length_bp", "Total intron length", dtype="int", unit="bp"),
                    Column("longest_bp", "Longest intron", dtype="int", unit="bp"),
                    Column("n_trans_splicing", "Trans-splicing junctions", dtype="int"),
                    Column("annotation_check", "Annotation cross-check"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_genes_with_introns": 12,
                        "n_introns": 17,
                        "total_length_bp": 12_547,
                        "longest_bp": 1_055,
                        "n_trans_splicing": 2,
                        "annotation_check": "agrees — 1 gene marked trans-splicing",
                    },
                ),
                scope=Scope(label="CDS only"),
            ),
        },
    )


def _diversity_result() -> Result:
    return Result(
        module="diversity",
        module_version="0.1.0",
        tables={
            "diversity_summary": Table(
                id="diversity_summary",
                title="Nucleotide diversity summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("n_genomes", "Genomes aligned", dtype="int"),
                    Column("alignment_length", "Alignment length", dtype="int", unit="bp"),
                    Column("usable_columns", "Usable columns", dtype="int"),
                    Column("pi", "π (whole alignment)", dtype="float", precision=6),
                    Column("window_bp", "Window", dtype="int", unit="bp"),
                    Column("step_bp", "Step", dtype="int", unit="bp"),
                    Column("tool", "Aligner"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_genomes": 2,
                        "alignment_length": 154_000,
                        "usable_columns": 153_500,
                        "pi": 0.004321,
                        "window_bp": 600,
                        "step_bp": 200,
                        "tool": "MAFFT",
                    },
                ),
                scope=Scope(label="whole aligned plastome", coord_system="normalized"),
            ),
        },
    )


def _substitutions_result() -> Result:
    return Result(
        module="substitutions",
        module_version="0.1.0",
        tables={
            "substitutions_summary": Table(
                id="substitutions_summary",
                title="Substitutions summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("n_genomes", "Genomes compared", dtype="int"),
                    Column("alignment_length", "Alignment length", dtype="int", unit="bp"),
                    Column("compared_columns", "Comparable columns (all comparisons)", dtype="int"),
                    Column("substitutions", "Substitutions", dtype="int"),
                    Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
                    Column("tool", "Aligner"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_genomes": 2,
                        "alignment_length": 154_000,
                        "compared_columns": 153_500,
                        "substitutions": 1_234,
                        "ts_tv": 2.6917,
                        "tool": "MAFFT",
                    },
                ),
                scope=Scope(label="whole aligned plastome", coord_system="normalized"),
            ),
        },
    )


def _barcodes_result() -> Result:
    return Result(
        module="barcodes",
        module_version="0.1.0",
        tables={
            "barcodes_summary": Table(
                id="barcodes_summary",
                title="Barcode and variability summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("n_genomes", "Genomes aligned", dtype="int"),
                    Column("window_bp", "Window", dtype="int", unit="bp"),
                    Column("step_bp", "Step", dtype="int", unit="bp"),
                    Column("n_windows", "Windows with a value", dtype="int"),
                    Column("mean_pi", "Mean π", dtype="float", precision=6),
                    Column("threshold_pi", "Region threshold", dtype="float", precision=6),
                    Column("threshold_rule", "Threshold rule"),
                    Column("n_regions", "Hypervariable regions", dtype="int"),
                    Column("n_candidates", "Candidate loci", dtype="int"),
                    Column("n_found", "Candidates found", dtype="int"),
                    Column("tool", "Aligner"),
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
                scope=Scope(label="the reference genome's normalized frame"),
            ),
            "barcodes_windows": Table(
                id="barcodes_windows",
                title="π by window, ranked",
                columns=(
                    Column("window", "Window", dtype="int"),
                    Column("alignment_start", "Alignment start", dtype="int"),
                    Column("alignment_end", "Alignment end", dtype="int"),
                    Column("reference_start", "Reference start", dtype="int", unit="bp"),
                    Column("reference_end", "Reference end", dtype="int", unit="bp"),
                    Column("usable_columns", "Usable columns", dtype="int"),
                    Column("pi", "π", dtype="float", precision=6),
                    Column("rank", "Rank by π", dtype="int"),
                    Column("region", "Region"),
                ),
                rows=(
                    {
                        "window": 1,
                        "alignment_start": 1,
                        "alignment_end": 600,
                        "reference_start": 1,
                        "reference_end": 600,
                        "usable_columns": 600,
                        "pi": 0.0099,
                        "rank": 411,
                        "region": "",
                    },
                ),
                scope=Scope(label="600 bp window / 200 bp step"),
            ),
            "barcodes_regions": Table(
                id="barcodes_regions",
                title="Hypervariable regions",
                columns=(
                    Column("region", "Region"),
                    Column("rank", "Rank by mean π", dtype="int"),
                    Column("reference_start", "Reference start", dtype="int", unit="bp"),
                    Column("reference_end", "Reference end", dtype="int", unit="bp"),
                    Column("span_bp", "Span", dtype="int", unit="bp"),
                    Column("n_windows", "Windows", dtype="int"),
                    Column("mean_pi", "Mean π", dtype="float", precision=6),
                    Column("max_pi", "Max π", dtype="float", precision=6),
                ),
                rows=(
                    {
                        "region": "R1",
                        "rank": 1,
                        "reference_start": 62234,
                        "reference_end": 63828,
                        "span_bp": 1595,
                        "n_windows": 6,
                        "mean_pi": 0.030225,
                        "max_pi": 0.038721,
                    },
                ),
                scope=Scope(label="maximal runs above the threshold"),
                key="region",
            ),
            "barcodes_candidates": Table(
                id="barcodes_candidates",
                title="Candidate barcode loci",
                columns=(
                    Column("locus", "Locus"),
                    Column("kind", "Kind"),
                    Column("members", "Members"),
                    Column("found", "Found"),
                    Column("reference_start", "Reference start", dtype="int", unit="bp"),
                    Column("reference_end", "Reference end", dtype="int", unit="bp"),
                    Column("span_bp", "Span", dtype="int", unit="bp"),
                    Column("pi", "π", dtype="float", precision=6),
                    Column("rank", "Rank by π", dtype="int"),
                    Column("note", "Note"),
                ),
                rows=(
                    {
                        "locus": "rbcL",
                        "kind": "gene",
                        "members": "rbcL",
                        "found": "yes",
                        "reference_start": 55212,
                        "reference_end": 56672,
                        "span_bp": 1461,
                        "pi": 0.008898,
                        "rank": 7,
                        "note": "",
                    },
                ),
                scope=Scope(label="the loci in data/barcode_loci.yaml"),
                key="locus",
            ),
        },
    )


def _popgen_result() -> Result:
    return Result(
        module="popgen",
        module_version="0.1.0",
        tables={
            "popgen_summary": Table(
                id="popgen_summary",
                title="Population-genetics summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("n_genomes", "Genomes", dtype="int"),
                    Column("n_groups", "Groups", dtype="int"),
                    Column("group_sizes", "Group sizes"),
                    Column("n_genes", "Genes analysed", dtype="int"),
                    Column("tool", "Aligner"),
                ),
                rows=(
                    {
                        "accession": "NC_057985.1",
                        "organism": "Scleromitrion diffusum",
                        "n_genomes": 9,
                        "n_groups": 2,
                        "group_sizes": "Coleus: 4, Oldenlandia: 5",
                        "n_genes": 80,
                        "tool": "MAFFT",
                    },
                ),
                scope=Scope(label="per-gene statistics by sample group"),
            ),
            "popgen_groups": Table(
                id="popgen_groups",
                title="Statistics by gene and group",
                columns=(
                    Column("row_id", "Row"),
                    Column("gene", "Gene"),
                    Column("group", "Group"),
                    Column("n", "Sequences", dtype="int"),
                    Column("sites", "Usable sites", dtype="int", unit="bp"),
                    Column("segregating", "Segregating sites", dtype="int"),
                    Column("singletons", "Singletons", dtype="int"),
                    Column("haplotypes", "Haplotypes", dtype="int"),
                    Column("hd", "Haplotype diversity", dtype="float", precision=4),
                    Column("pi", "π per site", dtype="float", precision=8),
                    Column("theta_w", "θ(S) per site", dtype="float", precision=8),
                    Column("tajimas_d", "Tajima's D", dtype="float", precision=4),
                    Column("fu_li_d_star", "Fu & Li's D*", dtype="float", precision=4),
                    Column("fu_li_f_star", "Fu & Li's F*", dtype="float", precision=4),
                    Column("note", "Note"),
                ),
                rows=(
                    {
                        "row_id": "rbcL#Oldenlandia",
                        "gene": "rbcL",
                        "group": "Oldenlandia",
                        "n": 4,
                        "sites": 1428,
                        "segregating": 41,
                        "singletons": 12,
                        "haplotypes": 3,
                        "hd": 0.8333,
                        "pi": 0.01844071,
                        "theta_w": 0.01985102,
                        "tajimas_d": 1.8434,
                        "fu_li_d_star": 1.8434,
                        "fu_li_f_star": 1.5432,
                        "note": "",
                    },
                ),
                scope=Scope(label="one row per gene and group"),
                key="row_id",
            ),
            "popgen_fst": Table(
                id="popgen_fst",
                title="Fst between groups",
                columns=(
                    Column("row_id", "Row"),
                    Column("gene", "Gene"),
                    Column("group_a", "Group A"),
                    Column("group_b", "Group B"),
                    Column("n_a", "n (A)", dtype="int"),
                    Column("n_b", "n (B)", dtype="int"),
                    Column("fst", "Fst", dtype="float", precision=4),
                    Column("note", "Note"),
                ),
                rows=(
                    {
                        "row_id": "rbcL#Coleus#Oldenlandia",
                        "gene": "rbcL",
                        "group_a": "Coleus",
                        "group_b": "Oldenlandia",
                        "n_a": 4,
                        "n_b": 4,
                        "fst": 0.1852,
                        "note": "",
                    },
                ),
                scope=Scope(label="Nei's Gst between groups"),
                key="row_id",
            ),
        },
    )


def _kaks_result() -> Result:
    return Result(
        module="kaks",
        module_version="0.1.0",
        tables={
            "kaks_summary": Table(
                id="kaks_summary",
                title="Ka/Ks summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("n_genomes", "Genomes compared", dtype="int"),
                    Column("n_genes", "Genes considered", dtype="int"),
                    Column("n_shared", "Genes shared with every genome", dtype="int"),
                    Column("estimator", "Estimator"),
                    Column("aligner", "Aligner"),
                    Column("tool", "Tool"),
                    Column("tool_version", "Tool version"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "n_genomes": 1,
                        "n_genes": 79,
                        "n_shared": 79,
                        "estimator": "Yang & Nielsen 2000 (yn00)",
                        "aligner": "MAFFT",
                        "tool": "PAML yn00",
                        "tool_version": "4.10.10",
                    },
                ),
                scope=Scope(label="the reference's CDS names"),
            ),
            "kaks_pooled": Table(
                id="kaks_pooled",
                title="Ka/Ks pooled over genes",
                columns=(
                    Column("genome", "Genome"),
                    Column("group", "Group"),
                    Column("n_shared", "Shared genes", dtype="int"),
                    Column("n_computed", "Genes scored", dtype="int"),
                    Column("n_not_computed", "Genes without a value", dtype="int"),
                    Column("pooled_ka", "Pooled Ka", dtype="float", precision=6),
                    Column("pooled_ks", "Pooled Ks", dtype="float", precision=6),
                    Column("pooled_ka_ks", "Pooled Ka/Ks", dtype="float", precision=4),
                ),
                rows=(
                    {
                        "genome": "NC_024936.1",
                        "group": "Iris",
                        "n_shared": 79,
                        "n_computed": 75,
                        "n_not_computed": 4,
                        "pooled_ka": 0.003701,
                        "pooled_ks": 0.016742,
                        "pooled_ka_ks": 0.2210,
                    },
                ),
                scope=Scope(label="pooled, site-weighted"),
                key="genome",
            ),
            "kaks_genes": Table(
                id="kaks_genes",
                title="Ka/Ks by gene",
                columns=(
                    Column("row_id", "Row"),
                    Column("gene", "Gene"),
                    Column("genome", "Genome"),
                    Column("group", "Group"),
                    Column("codons", "Codons", dtype="int"),
                    Column("sites_s", "Synonymous sites", dtype="float", precision=1),
                    Column("sites_n", "Nonsynonymous sites", dtype="float", precision=1),
                    Column("ka", "Ka (dN)", dtype="float", precision=6),
                    Column("ks", "Ks (dS)", dtype="float", precision=6),
                    Column("ka_ks", "Ka/Ks", dtype="float", precision=4),
                    Column("status", "Status"),
                    Column("detail", "Detail"),
                ),
                rows=(
                    {
                        "row_id": "rbcL#NC_024936.1",
                        "gene": "rbcL",
                        "genome": "NC_024936.1",
                        "group": "Iris",
                        "codons": 486,
                        "sites_s": 378.2,
                        "sites_n": 1079.8,
                        "ka": 0.0019,
                        "ks": 0.0299,
                        "ka_ks": 0.0635,
                        "status": "ok",
                        "detail": "",
                    },
                ),
                scope=Scope(label="one row per gene and genome"),
                key="row_id",
            ),
        },
    )


def _synteny_result() -> Result:
    return Result(
        module="synteny",
        module_version="0.1.0",
        tables={
            "synteny_summary": Table(
                id="synteny_summary",
                title="Synteny summary",
                columns=(
                    Column("accession", "Reference accession"),
                    Column("organism", "Reference organism"),
                    Column("reference_length", "Reference length", dtype="int", unit="bp"),
                    Column("n_genomes", "Genomes aligned", dtype="int"),
                    Column("window_bp", "Window", dtype="int", unit="bp"),
                    Column("step_bp", "Step", dtype="int", unit="bp"),
                    Column("tool", "Aligner"),
                ),
                rows=(
                    {
                        "accession": "NC_050833.1",
                        "organism": "Iris domestica",
                        "reference_length": 153_724,
                        "n_genomes": 1,
                        "window_bp": 600,
                        "step_bp": 200,
                        "tool": "MUMmer nucmer",
                    },
                ),
                scope=Scope(label="the reference and this comparison's parameters"),
            ),
            "synteny_genomes": Table(
                id="synteny_genomes",
                title="Synteny by genome",
                columns=(
                    Column("genome", "Genome"),
                    Column("group", "Group"),
                    Column("n_blocks", "Alignment blocks", dtype="int"),
                    Column("aligned_bp", "Aligned base pairs", dtype="int", unit="bp"),
                    Column("coverage_pct", "Coverage", dtype="float", unit="%", precision=2),
                    Column("min_identity", "Min %IDY", dtype="float", unit="%", precision=2),
                    Column("max_identity", "Max %IDY", dtype="float", unit="%", precision=2),
                    Column("mean_identity", "Mean %IDY", dtype="float", unit="%", precision=2),
                    Column(
                        "aligned_identity",
                        "Identity over base pairs",
                        dtype="float",
                        unit="%",
                        precision=2,
                    ),
                ),
                rows=(
                    {
                        "genome": "NC_024936.1",
                        "group": "Iris",
                        "n_blocks": 7,
                        "aligned_bp": 152_332,
                        "coverage_pct": 99.09,
                        "min_identity": 89.66,
                        "max_identity": 99.11,
                        "mean_identity": 97.37,
                        "aligned_identity": 99.31,
                    },
                ),
                scope=Scope(label="each genome against the reference"),
                key="genome",
            ),
            "synteny_blocks": Table(
                id="synteny_blocks",
                title="Alignment blocks against the reference",
                columns=(
                    Column("block", "Block"),
                    Column("genome", "Genome"),
                    Column("group", "Group"),
                    Column("target_start", "Target start", dtype="int", unit="bp"),
                    Column("target_end", "Target end", dtype="int", unit="bp"),
                    Column("query_start", "Query start", dtype="int", unit="bp"),
                    Column("query_end", "Query end", dtype="int", unit="bp"),
                    Column("length_target", "Target length", dtype="int", unit="bp"),
                    Column("length_query", "Query length", dtype="int", unit="bp"),
                    Column("identity", "%IDY", dtype="float", unit="%", precision=2),
                    Column("orientation", "Orientation"),
                ),
                rows=(
                    {
                        "block": "NC_024936.1#1",
                        "genome": "NC_024936.1",
                        "group": "Iris",
                        "target_start": 1,
                        "target_end": 13_576,
                        "query_start": 95,
                        "query_end": 13_614,
                        "length_target": 13_576,
                        "length_query": 13_520,
                        "identity": 98.0,
                        "orientation": "forward",
                    },
                ),
                scope=Scope(label="one row per alignment"),
                key="block",
            ),
            "synteny_identity": Table(
                id="synteny_identity",
                title="% identity along the reference",
                columns=(
                    Column("genome", "Genome"),
                    Column("group", "Group"),
                    Column("window_start", "Window start", dtype="int", unit="bp"),
                    Column("window_end", "Window end", dtype="int", unit="bp"),
                    Column("compared_bp", "Aligned base pairs", dtype="int", unit="bp"),
                    Column("identity", "% identity", dtype="float", unit="%", precision=2),
                ),
                rows=(
                    {
                        "genome": "NC_024936.1",
                        "group": "Iris",
                        "window_start": 1,
                        "window_end": 600,
                        "compared_bp": 600,
                        "identity": 99.0,
                    },
                ),
                scope=Scope(label="600 bp window / 200 bp step"),
            ),
        },
    )


def _gene_table_result() -> Result:
    return Result(
        module="gene_table",
        module_version="0.1.0",
        tables={
            "gene_table_summary": Table(
                id="gene_table_summary",
                title="Publication gene table — summary",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("n_genes_total", "Genes", dtype="int"),
                    Column("n_groups", "Functional groups", dtype="int"),
                    Column("n_unique_names", "Distinct gene names", dtype="int"),
                    Column("n_unclassified", "Unclassified", dtype="int"),
                    Column("coord_system", "Coordinates"),
                ),
                rows=(
                    {
                        "organism": "Iris domestica",
                        "accession": "NC_050833.1",
                        "n_genes_total": 132,
                        "n_groups": 18,
                        "n_unique_names": 112,
                        "n_unclassified": 0,
                        "coord_system": "genbank",
                    },
                ),
                scope=Scope(label="every counted gene, both IR copies; grouped by function"),
            ),
            "gene_reconciliation": Table(
                id="gene_reconciliation",
                title="Gene table reconciled against the genome-wide count",
                columns=(
                    Column("kind", "Kind"),
                    Column("in_table", "In the gene table", dtype="int"),
                    Column("counted", "count_genes()", dtype="int"),
                    Column("agreement", "Agreement"),
                ),
                rows=tuple(
                    {
                        "kind": kind,
                        "in_table": copies,
                        "counted": copies,
                        "agreement": "agrees",
                    }
                    for kind, copies in (
                        ("CDS", 86),
                        ("tRNA", 38),
                        ("rRNA", 8),
                        ("Total", 132),
                    )
                ),
                scope=Scope(label="every counted gene, both IR copies; grouped by function"),
                key="kind",
            ),
        },
    )


def _orf_result() -> Result:
    return Result(
        module="orf",
        module_version="0.1.0",
        tables={
            "orf_summary": Table(
                id="orf_summary",
                title="ORF check summary",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("n_cds", "CDS checked", dtype="int"),
                    Column("n_failed", "Failed", dtype="int"),
                    Column("n_pseudogene", "Pseudogenes (excluded from failures)", dtype="int"),
                    Column("n_annotated_not_coding", "Annotated, not coding", dtype="int"),
                    Column("translation_tables", "Translation table(s) used"),
                    Column("coord_system", "Coordinates"),
                ),
                rows=(
                    {
                        "organism": "Iris domestica",
                        "accession": "NC_050833.1",
                        "n_cds": 86,
                        "n_failed": 0,
                        "n_pseudogene": 0,
                        "n_annotated_not_coding": 1,
                        "translation_tables": "11",
                        "coord_system": "genbank",
                    },
                ),
                scope=Scope(label="every CDS feature, both IR copies"),
            ),
            "orf_status": Table(
                id="orf_status",
                title="CDS by ORF status",
                columns=(
                    Column("status", "Status"),
                    Column("meaning", "Meaning"),
                    Column("n", "CDS", dtype="int"),
                ),
                rows=tuple(
                    {"status": status, "meaning": "…", "n": n}
                    for status, n in (
                        ("not_multiple_of_three", 0),
                        ("no_start", 0),
                        ("no_stop", 0),
                        ("internal_stop", 0),
                        ("alternative_start", 5),
                        ("pseudogene", 0),
                        ("ok", 81),
                    )
                ),
                scope=Scope(label="every CDS feature, both IR copies"),
                key="status",
            ),
            "orf_checks": Table(
                id="orf_checks",
                title="ORF check per CDS",
                columns=(
                    Column("cds_id", "CDS"),
                    Column("gene", "Gene"),
                    Column("start_1b", "Start", dtype="int"),
                    Column("end_1b", "End", dtype="int"),
                    Column("strand", "Strand"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("transl_table", "Table", dtype="int"),
                    Column("codon_start", "Codon start", dtype="int"),
                    Column("start_codon", "Start codon"),
                    Column("start_class", "Start class"),
                    Column("stop_codon", "Stop codon"),
                    Column("internal_stops", "Internal stops", dtype="int"),
                    Column("internal_stop_codons", "Internal stop positions (codon)"),
                    Column("status", "Status"),
                    Column("pseudo", "Pseudogene"),
                ),
                rows=(
                    {
                        "cds_id": "rpl2:1-100",
                        "gene": "rpl2",
                        "start_1b": 1,
                        "end_1b": 100,
                        "strand": "+",
                        "length_bp": 100,
                        "transl_table": 11,
                        "codon_start": 1,
                        "start_codon": "ACG",
                        "start_class": "non-canonical (not in the table)",
                        "stop_codon": "TAA",
                        "internal_stops": 0,
                        "internal_stop_codons": "",
                        "status": "alternative_start",
                        "pseudo": "no",
                    },
                ),
                scope=Scope(label="every CDS feature, both IR copies"),
                key="cds_id",
            ),
            "annotated_not_coding": Table(
                id="annotated_not_coding",
                title="Annotated as a gene, not coding",
                columns=(
                    Column("name", "Name"),
                    Column("start_1b", "Start", dtype="int"),
                    Column("end_1b", "End", dtype="int"),
                    Column("strand", "Strand"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("product", "Product"),
                    Column("pseudo", "Pseudogene"),
                ),
                rows=(
                    {
                        "name": "ycf1",
                        "start_1b": 108_448,
                        "end_1b": 109_347,
                        "strand": "+",
                        "length_bp": 900,
                        "product": "",
                        "pseudo": "yes",
                    },
                ),
                scope=Scope(label="every CDS feature, both IR copies"),
            ),
            "orf_reconciliation": Table(
                id="orf_reconciliation",
                title="ORF checks reconciled against the genome-wide CDS count",
                columns=(
                    Column("kind", "Kind"),
                    Column("in_table", "Checked here", dtype="int"),
                    Column("counted", "count_genes()", dtype="int"),
                    Column("agreement", "Agreement"),
                ),
                rows=({"kind": "CDS", "in_table": 86, "counted": 86, "agreement": "agrees"},),
                scope=Scope(label="every CDS feature, both IR copies"),
                key="kind",
            ),
        },
    )


def _gene_names_result() -> Result:
    return Result(
        module="gene_names",
        module_version="0.1.0",
        tables={
            "panel_summary": Table(
                id="panel_summary",
                title="Gene names across the run",
                columns=(
                    Column("n_records", "Records", dtype="int"),
                    Column("n_genes", "Counted genes", dtype="int"),
                    Column("n_distinct_names", "Distinct names as they arrived", dtype="int"),
                    Column("n_distinct_canonical", "Distinct canonical names", dtype="int"),
                    Column("n_renamed", "Written differently", dtype="int"),
                    Column("n_aliased", "One locus, two nomenclatures", dtype="int"),
                    Column("n_resolved_by_sequence", "Resolved by sequence", dtype="int"),
                    Column("n_unresolved", "Unresolved", dtype="int"),
                ),
                rows=(
                    {
                        "n_records": 2,
                        "n_genes": 264,
                        "n_distinct_names": 118,
                        "n_distinct_canonical": 112,
                        "n_renamed": 6,
                        "n_aliased": 8,
                        "n_resolved_by_sequence": 1,
                        "n_unresolved": 0,
                    },
                ),
                scope=Scope(label="every counted gene of every record in the run"),
            ),
            "name_summary": Table(
                id="name_summary",
                title="Gene names by record",
                columns=(
                    Column("genome", "Record"),
                    Column("organism", "Organism"),
                    Column("n_genes", "Counted genes", dtype="int"),
                    Column("n_distinct_names", "Distinct names as they arrived", dtype="int"),
                    Column("n_distinct_canonical", "Distinct canonical names", dtype="int"),
                    Column("n_renamed", "Written differently", dtype="int"),
                    Column("n_aliased", "Two nomenclatures", dtype="int"),
                    Column("n_resolved_by_sequence", "Resolved by sequence", dtype="int"),
                    Column("n_unresolved", "Unresolved", dtype="int"),
                ),
                rows=tuple(
                    {
                        "genome": accession,
                        "organism": organism,
                        "n_genes": 132,
                        "n_distinct_names": 112,
                        "n_distinct_canonical": 112,
                        "n_renamed": 0,
                        "n_aliased": 4,
                        "n_resolved_by_sequence": 0,
                        "n_unresolved": 0,
                    }
                    for accession, organism in (
                        ("NC_050833.1", "Iris domestica"),
                        ("NC_024936.1", "Iris gatesii"),
                    )
                ),
                scope=Scope(label="every counted gene of every record in the run"),
                key="genome",
            ),
            "name_mappings": Table(
                id="name_mappings",
                title="Names that were normalized (or carry a synonym)",
                columns=(
                    Column("name_id", "Name id"),
                    Column("genome", "Record"),
                    Column("name", "Name as it arrived"),
                    Column("canonical", "Canonical name"),
                    Column("printed", "Printed as"),
                    Column("rule", "Why"),
                    Column("product", "Product in the record"),
                ),
                rows=(
                    {
                        "name_id": "NC_050833.1:ycf3",
                        "genome": "NC_050833.1",
                        "name": "ycf3",
                        "canonical": "pafI",
                        "printed": "pafI (ycf3)",
                        "rule": "semantic",
                        "product": "photosystem I assembly protein Ycf3",
                    },
                ),
                scope=Scope(label="every counted gene of every record in the run"),
                key="name_id",
            ),
        },
    )


def _submission_result() -> Result:
    return Result(
        module="submission",
        module_version="0.1.0",
        tables={
            "submission_summary": Table(
                id="submission_summary",
                title="Submission files",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("sequence_id", "Sequence id in the files"),
                    Column("n_features", "Features written", dtype="int"),
                    Column("n_cds", "CDS", dtype="int"),
                    Column("n_trna", "tRNA", dtype="int"),
                    Column("n_rrna", "rRNA", dtype="int"),
                    Column("n_partial_starts", "CDS written with a partial start", dtype="int"),
                    Column("n_problems", "First-pass problems", dtype="int"),
                    Column("n_validator_errors", "NCBI errors", dtype="int"),
                    Column("n_validator_warnings", "NCBI warnings", dtype="int"),
                    Column("sqn", "ASN.1 (.sqn)"),
                ),
                rows=(
                    {
                        "organism": "Iris domestica",
                        "accession": "NC_050833.1",
                        "sequence_id": "NC_050833.1",
                        "n_features": 265,
                        "n_cds": 86,
                        "n_trna": 38,
                        "n_rrna": 8,
                        "n_partial_starts": 2,
                        "n_problems": 0,
                        "n_validator_errors": 7,
                        "n_validator_warnings": 98,
                        "sqn": "produced by table2asn; reports kept: val, stats",
                    },
                ),
                scope=Scope(label="every counted feature, both IR copies"),
            ),
        },
    )


def _annotation_result() -> Result:
    return Result(
        module="annotation",
        module_version="0.2.0",
        tables={
            "annotation_summary": Table(
                id="annotation_summary",
                title="Annotation summary",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("n_genes", "Genes", dtype="int"),
                    Column("n_cds", "CDS", dtype="int"),
                    Column("n_trna", "tRNA", dtype="int"),
                    Column("n_rrna", "rRNA", dtype="int"),
                    Column("n_orf_failed", "CDS failing the ORF checks", dtype="int"),
                    Column("n_aliased_names", "Names with a synonym", dtype="int"),
                    Column("n_unclassified", "Unclassified", dtype="int"),
                    Column("annotator", "Annotator"),
                    Column("reference", "References"),
                ),
                rows=(
                    {
                        "organism": "Iris domestica",
                        "accession": "NC_050833.1",
                        "length_bp": 153_724,
                        "n_genes": 132,
                        "n_cds": 86,
                        "n_trna": 38,
                        "n_rrna": 8,
                        "n_orf_failed": 0,
                        "n_aliased_names": 4,
                        "n_unclassified": 0,
                        "annotator": "Chloe 0.1.18",
                        "reference": "its own curated set (Chloe_refs 9163915d99fd)",
                    },
                ),
                scope=Scope(label="the annotator's gene set, read back with our reader"),
            ),
            "orf_status": Table(
                id="orf_status",
                title="ORF status of the annotated CDS",
                columns=(Column("status", "Status"), Column("n", "CDS", dtype="int")),
                rows=tuple(
                    {"status": status, "n": n}
                    for status, n in (
                        ("not_multiple_of_three", 0),
                        ("no_start", 0),
                        ("no_stop", 0),
                        ("internal_stop", 0),
                        ("alternative_start", 5),
                        ("pseudogene", 0),
                        ("ok", 81),
                    )
                ),
                scope=Scope(label="the annotator's gene set, read back with our reader"),
                key="status",
            ),
        },
    )


def _mito_genes_result() -> Result:
    return Result(
        module="mito_genes",
        module_version="0.1.0",
        tables={
            "mito_summary": Table(
                id="mito_summary",
                title="Mitochondrial gene set",
                columns=(
                    Column("organism", "Organism"),
                    Column("accession", "Accession"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                    Column("topology", "Topology"),
                    Column("n_genes", "Genes", dtype="int"),
                    Column("n_cds", "CDS", dtype="int"),
                    Column("n_trna", "tRNA", dtype="int"),
                    Column("n_rrna", "rRNA", dtype="int"),
                    Column("n_named_by_symbol", "Named by a standard symbol", dtype="int"),
                    Column("ir_gate", "Inverted-repeat gate"),
                ),
                rows=(
                    {
                        "organism": "Oldenlandia corymbosa var. corymbosa",
                        "accession": "OX459128.1",
                        "length_bp": 258_274,
                        "topology": "circular",
                        "n_genes": 51,
                        "n_cds": 31,
                        "n_trna": 17,
                        "n_rrna": 3,
                        "n_named_by_symbol": 31,
                        "ir_gate": "not applicable",
                    },
                ),
                scope=Scope(label="every counted mitochondrial gene"),
            ),
        },
    )


def _mito_annotate_result() -> Result:
    """The P3A-03 module's tables, in the smallest shape its template resolves against."""
    scope = Scope(
        label="CDS and rRNA called by transfer from the reference set", coord_system="genbank"
    )
    summary = Table(
        id="mito_annotation_summary",
        title="Mitochondrial annotation",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("n_calls", "Genes called", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_not_called", "Vocabularies not called", dtype="int"),
            Column("n_cds_failing_orf", "Called CDS failing the ORF checks", dtype="int"),
            Column("n_dropped", "Hits not used", dtype="int"),
            Column("n_unnameable", "Hits this vocabulary cannot name", dtype="int"),
            Column("min_identity_pct", "Identity floor", dtype="float", unit="%"),
            Column("min_coverage_pct", "Coverage floor", dtype="float", unit="%"),
            Column("evalue", "E-value ceiling"),
            Column("transl_table", "Translation table", dtype="int"),
            Column("blast", "BLAST"),
            Column("reference", "Reference set"),
        ),
        rows=(
            {
                "organism": "Oldenlandia corymbosa var. corymbosa",
                "accession": "OX459128.1",
                "length_bp": 258_274,
                "n_calls": 34,
                "n_cds": 31,
                "n_rrna": 3,
                "n_not_called": 0,
                "n_cds_failing_orf": 10,
                "n_dropped": 68,
                "n_unnameable": 15,
                "min_identity_pct": 50.0,
                "min_coverage_pct": 50.0,
                "evalue": "1e-05",
                "transl_table": 1,
                "blast": "BLAST+ 2.17.0",
                "reference": "the packaged reference set",
            },
        ),
        scope=scope,
    )
    return Result(
        module="mito_annotate",
        module_version="0.1.0",
        tables={summary.id: summary},
    )


def _stats_result() -> Result:
    """The P5-03 module's summary, in the shape its template resolves against."""
    columns = (
        Column("accession", "Reference"),
        Column("organism", "Organism"),
        Column("n_genomes", "Genomes", dtype="int"),
        Column("n_groups", "Groups", dtype="int"),
        Column("groups", "Group labels"),
        Column("n_loci", "Loci", dtype="int"),
        Column("n_measures", "Measures", dtype="int"),
        Column("measures", "Measure set"),
        Column("test", "Test"),
        Column("correction", "Correction"),
        Column("alpha", "α", dtype="float", precision=4),
        Column("n_tested", "Tests with a p-value", dtype="int"),
        Column("n_not_computable", "Not computable", dtype="int"),
        Column("n_significant", "Significant after correction", dtype="int"),
        Column("n_observations", "Observations", dtype="int"),
        Column("n_collapsed_cds", "CDS collapsed to one per locus", dtype="int"),
        Column("n_genomes_with_cds", "Genomes carrying CDS", dtype="int"),
        Column("n_genomes_without_cds", "Genomes contributing no CDS", dtype="int"),
        Column("genomes_without_cds", "Genomes with no CDS"),
    )
    row = {
        "accession": "AP045470.1",
        "organism": "Coleus hadiensis",
        "n_genomes": 11,
        "n_groups": 2,
        "groups": "Coleus, Oldenlandia",
        "n_loci": 80,
        "n_measures": 4,
        "measures": "length_bp, gc_pct, gc3_pct, enc",
        "test": "Mann–Whitney U",
        "correction": "Benjamini–Hochberg",
        "alpha": 0.05,
        "n_tested": 303,
        "n_not_computable": 17,
        "n_significant": 0,
        "n_observations": 632,
        "n_collapsed_cds": 57,
        "n_genomes_with_cds": 8,
        "n_genomes_without_cds": 3,
        "genomes_without_cds": "AP045470.1, PQ325515.1, PZ428230.1",
    }
    summary = Table(
        id="stats_summary",
        title="Group comparison",
        columns=columns,
        rows=(row,),
        scope=Scope(label="per-locus measurements"),
    )
    return Result(
        module="stats",
        module_version="0.1.0",
        tables={"stats_summary": summary},
        citations=("biopython",),
    )


def _reference_audit_result() -> Result:
    """Built by the module's own table builder, so the template's contract is with the real one."""
    from plastid.core.reference_audit import AuditReport, Confirmation, Finding, RecordFacts
    from plastid.modules import reference_audit as audit_module

    record = RecordFacts(
        file="/references/NC_050833.1.gb",
        sha256="a" * 64,
        accession="NC_050833.1",
        length_bp=153_724,
        cds=86,
        trna=38,
        rrna=8,
        gene_names=("rbcL", "trnG-UCC"),
        ssc_orientation="opposite",
        blocked=False,
    )
    report = AuditReport(
        findings=(
            Finding(
                file=record.file,
                check="orientation",
                verdict="REPORT",
                reason="SSC orientation differs from the majority (2 of 3 read the other way)",
            ),
        ),
        records=(record,),
        confirmations=(
            Confirmation(
                region="SSC",
                technique="mafft --adjustdirection",
                measured=False,
                note="not measured (mafft 7.526 not found)",
            ),
        ),
    )
    return audit_module.audit_result(report)


def _gene_normalize_result() -> Result:
    """The summary a naming run writes — the cells its template reads.

    Hand-built like the fixtures around it: what the guard checks is that the template's
    references resolve against the schema the module documents, and a fixture that called the
    module's own table builders would only be checking that the module agrees with itself.
    """
    summary = Table(
        id="normalize_summary",
        title="Reference-driven naming",
        columns=(
            Column("mode", "Mode"),
            Column("reference", "Naming reference (A)"),
            Column("records", "Records normalized", dtype="int"),
            Column("genes", "Genes read", dtype="int"),
            Column("renamed", "Renamed", dtype="int"),
            Column("confirmed", "Confirmed", dtype="int"),
            Column("conflicts", "Conflicts", dtype="int"),
            Column("measured", "Status: measured", dtype="int"),
            Column("not_measured", "Status: not measured", dtype="int"),
            Column("outside", "Status: outside vocabulary", dtype="int"),
            Column("artifacts", "Normalized records written", dtype="int"),
        ),
        rows=(
            {
                "mode": "normalize",
                "reference": "NC_050833.1",
                "records": 1,
                "genes": 132,
                "renamed": 5,
                "confirmed": 120,
                "conflicts": 0,
                "measured": 120,
                "not_measured": 11,
                "outside": 1,
                "artifacts": 1,
            },
        ),
        scope=Scope(
            label="one record against the project's naming reference (D5-09, vocabulary A)"
        ),
    )
    return Result(
        module="gene_normalize", module_version="0.1.0", tables={"normalize_summary": summary}
    )


TEMPLATE_FIXTURES = {
    "gene_normalize": _gene_normalize_result,
    "reference_audit": _reference_audit_result,
    "mito_annotate": _mito_annotate_result,
    "mito_genes": _mito_genes_result,
    "annotation": _annotation_result,
    "genome_stats": _result,
    "gene_table": _gene_table_result,
    "gene_names": _gene_names_result,
    "orf": _orf_result,
    "submission": _submission_result,
    "ssr": _ssr_result,
    "structure": _structure_result,
    "codon_usage": _codon_usage_result,
    "amino_acid": _amino_acid_result,
    "introns": _introns_result,
    "diversity": _diversity_result,
    "substitutions": _substitutions_result,
    "barcodes": _barcodes_result,
    "kaks": _kaks_result,
    "popgen": _popgen_result,
    "synteny": _synteny_result,
    "stats": _stats_result,
}


def test_shipped_templates_have_no_bare_digits_and_all_refs_resolve() -> None:
    templates = list_templates()
    assert templates, "no report templates found"
    for path in templates:
        module_id = path.name.split(".")[0]
        assert module_id in TEMPLATE_FIXTURES, f"no fixture Result for template {path.name}"
        template = load_template(path.name)
        assert find_bare_digits(template) == [], f"{path.name}: literal digits outside references"
        for ref in find_refs(template):
            TEMPLATE_FIXTURES[module_id]().value(ref)  # schema vs template contract


def test_citations_file_loads_and_formats() -> None:
    citations = load_citations()
    assert "misa" in citations
    line = format_citation(citations["misa"])
    assert "Beier" in line and line.endswith("doi:10.1093/bioinformatics/btx198")
    assert format_citation(Citation(key="x", text="Plain reference.")) == "Plain reference."


def test_defaults_load() -> None:
    defaults = load_defaults()
    assert defaults["ir"]["min_identity"] == 0.999
    assert defaults["ssr"]["misa"]["unit_repeats"][1] == 10  # MISA default, not tuned
