"""The queryable-number check over *generated* reports (Phase 5, M5.4).

A report the toolkit produces may carry a number only if that number is a bounded number
inside the rendered value of a cell of a shipped table — and of cells that render the same
value, so the reader can point at exactly one place it came from. The template scan
(`find_bare_digits`) guards the template; this guards the report.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.report.audit import Audit, Finding, audit_report, audit_run
from plastid.report.citations import load_citations
from plastid.report.engine import DRAFT_BANNER, load_template, render
from plastid.result.result import Column, Result, Scope, Table


def _result() -> Result:
    """Genome statistics with a keyed region table — the shape the shipped template uses."""
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
            "regions": Table(
                id="regions",
                title="Regions",
                columns=(
                    Column("region", "Region"),
                    Column("length_bp", "Length", dtype="int", unit="bp"),
                ),
                rows=(
                    {"region": "LSC", "length_bp": 83_127},
                    {"region": "IRa", "length_bp": 26_214},
                    {"region": "SSC", "length_bp": 18_169},
                ),
                scope=Scope(label="normalized frame", coord_system="normalized"),
                key="region",
            ),
        },
        citations=("biopython",),
    )


def _codes(audit: Audit) -> list[tuple[str, str]]:
    return sorted((finding.kind, finding.number) for finding in audit.findings)


def test_a_number_no_cell_carries_is_an_orphan() -> None:
    audit = audit_report("The inverted repeat is 99999 bp long.\n", _result())
    assert _codes(audit) == [("orphan", "99999")]
    assert not audit.ok
    assert "99999" in str(audit.findings[0])


def test_a_number_carried_by_two_different_values_is_ambiguous() -> None:
    """A coordinate span and a length that begin with the same digits cannot be told apart."""
    result = Result(
        module="ssr",
        module_version="0.1.0",
        tables={
            "ssrs": Table(
                id="ssrs",
                title="SSRs",
                columns=(
                    Column("motif", "Motif"),
                    Column("span", "Span"),
                    Column("length_bp", "Length", dtype="int"),
                ),
                rows=(
                    {"motif": "(T)10", "span": "5741-5780", "length_bp": 5741},
                    {"motif": "(A)10", "span": "100-200", "length_bp": 100},
                ),
                scope=Scope(label="whole plastome"),
            )
        },
    )
    audit = audit_report("A repeat of 5741 units.\n", result)
    assert _codes(audit) == [("ambiguous", "5741")]
    finding = audit.findings[0]
    assert finding.values == ("5741", "5741-5780")
    assert finding.cells == (("ssrs[#1].length_bp", "5741"), ("ssrs[#1].span", "5741-5780"))
    assert "ssrs[#1].span" in str(finding)


def test_a_number_printed_inside_a_quoted_text_cell_is_not_a_claimed_quantity() -> None:
    """The measured false positive of the strict rule: a rule sentence, a note, a block id."""
    result = Result(
        module="barcodes",
        module_version="0.2.0",
        tables={
            "barcodes_summary": Table(
                id="barcodes_summary",
                title="Barcodes",
                columns=(
                    Column("n_candidates", "Candidates", dtype="int"),
                    Column("threshold_rule", "Threshold rule"),
                ),
                rows=({"n_candidates": 10, "threshold_rule": "mean + 1 sd of this run's windows"},),
                scope=Scope(label="whole plastome"),
            )
        },
    )
    audit = audit_report(
        "10 candidates stood above the threshold (mean + 1 sd of this run's windows).\n", result
    )
    assert audit.ok, [str(finding) for finding in audit.findings]
    assert any("text column" in reason for reason in audit.exempt)


def test_the_same_value_in_two_rows_is_still_exactly_one_answer() -> None:
    """Two rows that render the same value cannot disagree — the reader still gets one number."""
    result = Result(
        module="ssr",
        module_version="0.1.0",
        tables={
            "ssr_by_region": Table(
                id="ssr_by_region",
                title="SSR by region",
                columns=(Column("region", "Region"), Column("total", "Total", dtype="int")),
                rows=({"region": "IRa", "total": 7}, {"region": "IRb", "total": 7}),
                scope=Scope(label="whole plastome"),
                key="region",
            )
        },
    )
    audit = audit_report("There are 7 repeats.\n", result)
    assert audit.ok, [str(finding) for finding in audit.findings]


def test_a_digit_inside_a_word_is_not_a_number() -> None:
    """`rps12` and `NC_050833.1` contain digits; neither is a quantity the report claims."""
    audit = audit_report("The gene rps12 sits at the IR junction of NC_050833.1.\n", _result())
    assert audit.ok, [str(finding) for finding in audit.findings]


def test_a_negative_and_a_scientific_number_are_read_whole() -> None:
    result = Result(
        module="structure",
        module_version="0.1.0",
        tables={
            "regions": Table(
                id="regions",
                title="Regions",
                columns=(
                    Column("region", "Region"),
                    Column("gc_skew", "GC skew", dtype="float", precision=4),
                    Column("evalue", "E-value"),
                ),
                rows=({"region": "LSC", "gc_skew": -0.0230, "evalue": "1e-05"},),
                scope=Scope(label="normalized frame", coord_system="normalized"),
                key="region",
            )
        },
    )
    audit = audit_report("Skew is -0.0230 and the ceiling is 1e-05.\n", result)
    assert audit.ok, [str(finding) for finding in audit.findings]


def test_the_reference_block_is_exempt_and_the_audit_says_why() -> None:
    audit = audit_report(
        "A plain sentence.\n\n### References\n\n"
        "- Cock PJA et al. Bioinformatics 25:1422–1423 (2009). doi:10.1093/bioinformatics/btp163\n",
        _result(),
    )
    assert audit.ok, [str(finding) for finding in audit.findings]
    assert any("citation" in reason for reason in audit.exempt)


def test_the_run_data_label_is_exempt_and_the_audit_says_why() -> None:
    audit = audit_report(
        "> **Data:** stand-in batch 2026-10-04\n\nThe length is 153724 bp.\n", _result()
    )
    assert audit.ok, [str(finding) for finding in audit.findings]
    assert any("label" in reason for reason in audit.exempt)


def test_the_draft_banner_is_exempt() -> None:
    audit = audit_report(f"A plain sentence.\n\n{DRAFT_BANNER}\n", _result())
    assert audit.ok, [str(finding) for finding in audit.findings]


def test_a_report_the_engine_renders_passes_the_audit() -> None:
    result = _result()
    text = render(load_template("genome_stats.en.md"), result, citations=load_citations())
    audit = audit_report(text, result)
    assert audit.ok, [str(finding) for finding in audit.findings]
    assert audit.n_numbers > 5, "the audit checked no numbers — it would pass on anything"


def test_a_planted_number_in_a_rendered_report_fails() -> None:
    """The guard the phase asks for: a number no table cell carries stops the build."""
    result = _result()
    text = render(load_template("genome_stats.en.md"), result, citations=load_citations())
    planted = text.replace("Coordinates are reported", "About 26000 bp. Coordinates are reported")
    audit = audit_report(planted, result)
    assert _codes(audit) == [("orphan", "26000")]


def test_a_report_with_no_citations_has_no_citation_exemption() -> None:
    """An exemption nobody needed is not printed — the reasons list what was actually skipped."""
    result = Result(
        module="genome_stats",
        module_version="0.1.0",
        tables={
            "genome_stats": Table(
                id="genome_stats",
                title="Genome statistics",
                columns=(Column("n_cds", "CDS", dtype="int"),),
                rows=({"n_cds": 86},),
                scope=Scope(label="whole plastome"),
            )
        },
    )
    audit = audit_report("There are 86 CDS.\n", result)
    assert audit.ok
    assert audit.exempt == ()


REPO = Path(__file__).resolve().parents[2]
PHASE4_RUNS = sorted((REPO / ".runs" / "phase4").glob("*/*/manifest.json"))


@pytest.mark.skipif(
    not PHASE4_RUNS,
    reason="the Phase-4 deliverable run directory is absent (.runs/ is gitignored)",
)
def test_the_whole_phase_4_deliverable_run_passes_the_audit() -> None:
    """M5.4's Done-when: the check runs over the whole run directory of a Phase-4 deliverable."""
    reports = []
    for manifest in PHASE4_RUNS:
        reports.extend(audit_run(manifest.parent))
    assert len(reports) >= 5, f"only {len(reports)} reports found — the audit covered too little"
    failures = [
        f"{report.module} {report.report}: {finding}"
        for report in reports
        for finding in report.audit.findings
    ]
    assert not failures, "orphan or ambiguous numbers in the deliverable reports:\n" + "\n".join(
        failures
    )


def test_using_a_finding_as_a_message_names_the_number_and_where_it_came_from() -> None:
    finding = Finding(kind="orphan", number="9", cells=())
    assert "9" in str(finding) and "orphan" in str(finding)
