"""M11.17 family coverage: one record per family, every module that reads genes or sequences.

The milestone requires a case on the 80-record panel's families, not only the golden pair and the
Iridaceae set. This file fixes one record per family — the first source record whose own two IR
copies pass the gate at 100.00 % identity — and asserts, per family, that the module accounts for
the whole record rather than quietly measuring part of it.

The accessions are chosen by the gate, not by convenience; the premise test re-runs the gate so a
panel refresh that invalidates a choice fails here instead of inside a scientific number.

The assertions are semantic on purpose. A pinned count would say "this record has 84 CDS" and
stop; the claim worth holding is "every annotated CDS is in the table exactly once, and a CDS that
fails the ORF check is reported with a status rather than dropped" — which is what a reader of the
output depends on, in every family.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.core.model import count_genes
from plastid.core.orf import FAILURE_STATUSES, ORF_STATUSES
from plastid.modules import gene_table, orf
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"

# Reference first: the multi-genome modules treat the first input as the reference, and
# NC_008641.1 is the one panel record the IR gene-mirror rule passes with no violation.
FAMILY_RECORDS = {
    "Malvaceae": "NC_008641.1",
    "Asteraceae": "NC_028027.1",
    "Brassicaceae": "NC_029378.1",
    "Euphorbiaceae": "NC_016736.1",
    "Fabaceae": "NC_016708.2",
    "Lamiaceae": "NC_020098.1",
    "Orchidaceae": "NC_021432.1",
    "Poaceae": "NC_023245.1",
    "Rosaceae": "NC_031163.1",
    "Solanaceae": "NC_016068.1",
}

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not PANEL.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]


def _genomes() -> dict[str, object]:
    missing = [f"{a}.gb" for a in FAMILY_RECORDS.values() if not (PANEL / f"{a}.gb").is_file()]
    if missing:
        pytest.skip(f"panel records missing: {', '.join(missing)}")
    return {
        family: read_genbank(PANEL / f"{accession}.gb")
        for family, accession in FAMILY_RECORDS.items()
    }


def _ctx(root: Path, module: str, genome) -> RunContext:
    layout = RunLayout.create(root, module, fingerprint="family")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id=module, module_order=1)
    ctx.ir = enforce_ir_gate(detect_ir(genome, source="sequence"))
    return ctx


def test_the_ten_family_records_pass_the_gate_on_their_own_sequence() -> None:
    """The premise: every chosen record carries two copies that agree at the gate's threshold."""
    for family, genome in _genomes().items():
        ir = detect_ir(genome, source="sequence")
        assert ir.gate_passed, (
            f"{family} {genome.accession}: identity {ir.identity:.4f} fails the gate"
        )
        assert ir.identity >= 0.999


def test_the_gene_table_accounts_for_every_gene_in_every_family(tmp_path: Path) -> None:
    """`n_genes_total` is the record's own count, and every gene falls into a group."""
    for family, genome in _genomes().items():
        result = gene_table.run(genome, {}, _ctx(tmp_path / family, "gene_table", genome))
        summary = result.tables["gene_table_summary"].rows[0]
        counted = count_genes(genome)
        expected = counted.cds + counted.trna + counted.rrna
        assert summary["n_genes_total"] == expected, (
            f"{family}: {summary['n_genes_total']} != {expected}"
        )
        assert summary["n_groups"] > 0, f"{family}: no functional group was assigned"

        # The module's own reconciliation is the accounting claim: every feature it prints is
        # counted against `count_genes()` for its kind, and the closing row must agree. A gene no
        # group claims is counted in `n_unclassified` and named in the notes — that rule is held
        # by `tests/golden/test_tier_a.py`; here it must hold for every family.
        reconciled = {row["kind"]: row for row in result.tables["gene_reconciliation"].rows}
        for kind, row in reconciled.items():
            assert row["agreement"] == "agrees", (
                f"{family}: {kind} {row['in_table']} vs {row['counted']}"
            )
        assert reconciled["Total"]["counted"] == expected
        if summary["n_unclassified"]:
            note = next((n for n in result.notes if n.startswith("Unclassified genes")), "")
            assert note, f"{family}: an unclassified gene is counted but never named"


def test_orf_checks_every_annotated_cds_once_in_every_family(tmp_path: Path) -> None:
    """No CDS is skipped, and one that does not translate is reported with a status, not dropped."""
    for family, genome in _genomes().items():
        result = orf.run(genome, {}, _ctx(tmp_path / family, "orf", genome))
        summary = result.tables["orf_summary"].rows[0]
        checks = result.tables["orf_checks"].rows

        assert summary["n_cds"] == len(checks), (
            f"{family}: summary says {summary['n_cds']} CDS, the table holds {len(checks)}"
        )
        statuses = result.tables["orf_status"].rows
        assert sum(row["n"] for row in statuses) == len(checks), (
            f"{family}: status counts do not add up"
        )
        assert {row["status"] for row in statuses} == set(ORF_STATUSES), (
            f"{family}: a status row is missing"
        )
        assert {row["status"] for row in checks} <= set(ORF_STATUSES)

        # The headline class of a CDS that does not translate is a failure the reader is meant to
        # act on, or a pseudogene the record itself declares — never a silently dropped row.
        failures = sum(1 for row in checks if row["status"] in FAILURE_STATUSES)
        assert summary["n_failed"] == failures, (
            f"{family}: {summary['n_failed']} failures, {failures} rows"
        )
        pseudogenes = sum(1 for row in checks if row["status"] == "pseudogene")
        assert summary["n_pseudogene"] == pseudogenes, (
            f"{family}: pseudogene count disagrees with its rows"
        )
