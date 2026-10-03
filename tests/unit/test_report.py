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
                        "ts_tv": 2.686,
                        "tool": "MAFFT",
                    },
                ),
                scope=Scope(label="whole aligned plastome", coord_system="normalized"),
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


TEMPLATE_FIXTURES = {
    "genome_stats": _result,
    "gene_table": _gene_table_result,
    "ssr": _ssr_result,
    "structure": _structure_result,
    "codon_usage": _codon_usage_result,
    "amino_acid": _amino_acid_result,
    "introns": _introns_result,
    "diversity": _diversity_result,
    "substitutions": _substitutions_result,
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
