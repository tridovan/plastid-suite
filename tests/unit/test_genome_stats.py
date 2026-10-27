"""``genome_stats`` on a record that annotates nothing (M11.17)."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import Genome
from plastid.modules import genome_stats as stats_module
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout


def test_a_record_with_no_genes_reports_counts_without_a_gene_table(
    tmp_path: Path, golden_genome: Genome
) -> None:
    """The summary's zeros are a result; a zero-row per-gene table is not.

    ``Table`` refuses zero rows, so an unannotated record made this module raise while building
    the table it had nothing to put in (measured 27/10/2026). The counts still report 0 and a
    note says the table is absent.
    """
    genome = replace(golden_genome, genes=())
    layout = RunLayout.create(tmp_path / "runs", "genome_stats", fingerprint="test")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="genome_stats", module_order=1)
    ctx.ir = enforce_ir_gate(detect_ir(genome))

    result = stats_module.run(genome, {}, ctx)

    assert "genes" not in result.tables
    summary = result.tables["genome_stats"].rows[0]
    assert (summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (0, 0, 0)
    assert any("per-gene table is absent" in note for note in result.notes)
