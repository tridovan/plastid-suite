"""Every CDS analysis states whether partial genes are included, and counts them (M11.17)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import amino_acid, codon_usage, genome_stats, orf
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout


def _genome() -> Genome:
    """Two CDS: one complete, one with an uncertain 5' coordinate."""
    coding = "ATGGCTTAA"
    genes = (
        Gene(name="rbcL", kind="CDS", parts=(Interval(0, len(coding)),)),
        Gene(
            name="psbA",
            kind="CDS",
            parts=(Interval(len(coding), 2 * len(coding), partial_start=True),),
        ),
    )
    return Genome(name="Testus", accession="T.1", sequence=coding * 2, genes=genes)


def _ctx(tmp_path: Path, module_id: str) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", module_id, fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module_id, module_order=1)


def test_codon_usage_counts_partial_cds_instead_of_merging_them_silently(
    tmp_path: Path,
) -> None:
    summary = (
        codon_usage.run(_genome(), {}, _ctx(tmp_path, "codon_usage")).tables["codon_usage"].rows[0]
    )
    assert summary["n_cds"] == 2
    assert summary["n_partial_cds"] == 1


def test_amino_acid_counts_partial_cds(tmp_path: Path) -> None:
    summary = (
        amino_acid.run(_genome(), {}, _ctx(tmp_path, "amino_acid")).tables["amino_acid"].rows[0]
    )
    assert (summary["n_cds"], summary["n_partial_cds"]) == (2, 1)


def test_orf_marks_the_partial_cds_row(tmp_path: Path) -> None:
    rows = orf.run(_genome(), {}, _ctx(tmp_path, "orf")).tables["orf_checks"].rows
    by_gene = {row["gene"]: row for row in rows}
    assert by_gene["psbA"]["partial"] == "yes"
    assert by_gene["rbcL"]["partial"] == "no"


def test_genome_stats_marks_the_partial_gene(tmp_path: Path, golden_genome: Genome) -> None:
    extra = Gene(name="psbX", kind="CDS", parts=(Interval(200, 260, partial_start=True),))
    genome = replace(golden_genome, genes=(*golden_genome.genes, extra))
    ctx = _ctx(tmp_path, "genome_stats")
    ctx.ir = enforce_ir_gate(detect_ir(genome))

    rows = genome_stats.run(genome, {}, ctx).tables["genes"].rows

    assert next(row for row in rows if row["name"] == "psbX")["partial"] == "yes"
