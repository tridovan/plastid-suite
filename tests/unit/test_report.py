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


def test_shipped_templates_have_no_bare_digits_and_all_refs_resolve() -> None:
    templates = list_templates()
    assert templates, "no report templates found"
    for path in templates:
        template = load_template(path.name)
        assert find_bare_digits(template) == [], f"{path.name}: literal digits outside references"
        for ref in find_refs(template):
            _result().value(ref)  # raises when the template and the table schema disagree


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
