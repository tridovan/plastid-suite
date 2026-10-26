"""The audit against the real corpus — the Done-when of P3B-03 (D5-08).

The pinned reference set is the 34 Iridaceae plastomes of the corpus, and a hand count measured what
an audit of them has to say: **32/34 share one SSC orientation, 2/34 are inverted with all their
SSC CDS translating cleanly**, and exactly one record carries a partial coordinate. Those three
numbers are the contract here — the audit may not block the two inverted records (blocking would
refuse two published plastomes, one of them an *Iris*) and may not let the partial one through.

The set itself is 37 MB and lives outside the repo (`.panel/`), so these tests
skip on a machine that has not pinned it; the negative cases that need no corpus live in
tests/unit/test_reference_audit.py and run everywhere.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core import reference_audit as audit
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

PANEL = Path(__file__).resolve().parents[2] / ".panel" / "plastid-reference"

# Both measured by hand on 05/10/2026 and reproduced by the audit on 08/10/2026.
INVERTED = {"NC_042827.1.gb", "NC_062351.1.gb"}
PARTIAL = "NC_024936.1.gb"

pytestmark = pytest.mark.skipif(
    not PANEL.is_dir(), reason="the pinned reference set lives outside the repo (.panel/)"
)


def _paths() -> list[Path]:
    return sorted(PANEL.glob("*.gb"))


def test_the_audit_reproduces_the_measured_corpus() -> None:
    """32/34 one orientation · 2/34 inverted and reported · one record blocked, on partial."""
    report = audit.audit_files(_paths(), confirm=None)

    assert len(report.records) == 34
    assert {Path(finding.file).name for finding in report.blocked} == {PARTIAL}
    assert {finding.check for finding in report.blocked} == {"partial"}
    assert not report.ok  # a blocked record stops a run — that is the gate

    orientations = {Path(record.file).name: record.ssc_orientation for record in report.records}
    assert sum(1 for value in orientations.values() if value == "same") == 32
    assert {name for name, value in orientations.items() if value == "opposite"} == INVERTED

    reported = {Path(finding.file).name for finding in report.reported}
    assert reported >= INVERTED  # reported, never blocked


def test_every_record_carries_its_own_translation_and_translates() -> None:
    """The two rules that make blastx transfer possible, measured across the whole set."""
    report = audit.audit_files(_paths(), confirm=None)

    assert not [finding for finding in report.findings if finding.check == "translation"]
    assert not [finding for finding in report.findings if finding.check == "orf"]
    assert not [finding for finding in report.findings if finding.check == "ir"]


def test_the_gene_set_report_names_the_records_that_spell_a_locus_differently() -> None:
    """Measured: 34/34 records spell something differently — so the report cannot be empty, and it
    has to be about names rather than loci, because which spelling wins is D5-09's decision."""
    report = audit.audit_files(_paths(), confirm=None)

    gene_sets = [finding for finding in report.reported if finding.check == "gene_set"]
    assert len(gene_sets) == 20  # measured: the majority name set is held by 14 of 34
    assert all(
        finding.reason.startswith("gene set differs from the majority of 34 records")
        for finding in gene_sets
    )


@pytest.mark.tier_b
def test_mafft_adjustdirection_names_the_same_two_records(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    """The second reading, and the measurement that says it is worth having.

    Measured 08/10/2026 on MAFFT 7.526: aligning the 34 SSC regions (614,658 bp in total, 17.6 s)
    and asking it to adjust direction flags **exactly** the two records the gene-strand rule calls
    inverted — so the confirmation layer does not merely produce a picture, it independently
    reaches the same verdict.
    """
    require_tools(tools, "mafft")
    layout = RunLayout.create(
        tmp_path, "reference-audit", fingerprint="0" * 16, now=datetime(2026, 10, 8, 12, 0, 0)
    )
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="reference_audit", module_order=1)

    from plastid.external.mafft import adjust_direction

    report = audit.audit_files(_paths(), ctx=ctx, confirm=adjust_direction)

    assert report.confirmations[0].measured
    reversed_names = {
        name for name, direction in report.confirmations[0].directions if direction == "reverse"
    }
    assert reversed_names == INVERTED
    # agreement is asserted, not assumed: a disagreement would arrive as its own finding
    assert not [finding for finding in report.findings if finding.check == "confirmation"]
    assert (layout.logs_dir / "001__mafft__stdout.log").is_file()  # the call is on the record
