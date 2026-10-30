"""M11.10 acceptance on the validation data: 80 RefSeq plastomes (10 families) and 30 PGA variants.

The golden pair is Iridaceae; this check is the evidence that the flow is not fitted to it.

* **80 RefSeq records, 10 families.** The audit classifies every record and names its reason. The
  set blocks 8 records with 10 findings — 4 records whose IR fails on both paths, and 4 whose CDS
  fail the ORF check (one of those carries two bad CDS) — and the rest are usable references.
* **30 PGA annotation variants.** These are annotation files, not reference sets: measured
  26/10/2026 every one is blocked because its CDS carry no `/translation`, which is the rule the
  annotator needs. The classification is asserted rather than smoothed into "0 blocked".
* **A non-Iridaceae annotation.** One Solanaceae target is annotated from its own family's
  non-blocked references, and the artifact is re-read.

The panel is gitignored; without it the module skips with a named reason.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.genbank import read_genbank
from plastid.core.reference_audit import REPORT, audit_files

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"
VARIANTS = ROOT / ".panel" / "ir-validation" / "annotation-variants"
MANIFEST = ROOT / ".panel" / "ir-validation" / "manifest.tsv"

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.panel,
    pytest.mark.skipif(
        not PANEL.is_dir() or not VARIANTS.is_dir(),
        reason="the IR validation panel is not in this tree (.panel/ is gitignored)",
    ),
]

PANEL_RECORDS = 80
PANEL_BLOCKED_RECORDS = 8
PANEL_BLOCK_CHECKS = {"ir": 4, "orf": 6}
VARIANTS_RECORDS = 30
TARGET = "NC_016068.1"  # Solanaceae; 15 other Solanaceae records remain as references


@pytest.fixture(scope="module")
def panel_report():
    return audit_files(sorted(PANEL.glob("*.gb")), confirm=None)


@pytest.fixture(scope="module")
def variant_report():
    return audit_files(sorted(VARIANTS.glob("*.gb")), confirm=None)


def test_every_panel_record_is_classified_with_its_reason(panel_report) -> None:
    assert len(panel_report.records) == PANEL_RECORDS
    blocked_files = {Path(finding.file).name for finding in panel_report.blocked}
    assert len(blocked_files) == PANEL_BLOCKED_RECORDS
    assert Counter(finding.check for finding in panel_report.blocked) == PANEL_BLOCK_CHECKS
    assert all(finding.reason for finding in panel_report.blocked)
    # nothing is silently dropped: every file has a facts row, blocked or not
    assert {Path(record.file).name for record in panel_report.records} == {
        path.name for path in PANEL.glob("*.gb")
    }


def test_records_whose_ir_genes_do_not_mirror_are_reported_not_blocked(panel_report) -> None:
    """Measured 29/10/2026: 23/80 panel records break the IR gene-mirror rule.

    The rule: a gene inside one IR copy must have its copy inside the other, on the **opposite
    strand**. These are third-party annotation inconsistencies — a same-strand copy
    (`NC_016736.1`'s `trnI-CAU`/`trnL-CAA`) or a gene annotated on one side only (`NC_033499.1`'s
    `rrn5`) — checked against the record's own sequence. They are REPORT findings, never blocks: a
    one-sided ORF or pseudogene does not by itself make the record unusable, and a person confirms
    it (pinned to the file's sha256) or fixes the record.
    """
    mirror = [finding for finding in panel_report.findings if finding.check == "mirror"]

    assert len({finding.file for finding in mirror}) == 23
    assert all(finding.verdict == REPORT for finding in mirror)
    assert not any(finding.blocks for finding in mirror)


def test_the_pga_variants_are_annotations_not_a_reference_set(
    variant_report,
) -> None:
    """All 30 block on `/translation` — the rule the annotator relies on, not a panel fault."""
    assert len(variant_report.records) == VARIANTS_RECORDS
    assert len(variant_report.blocked) > VARIANTS_RECORDS
    assert {finding.check for finding in variant_report.blocked} == {"translation", "orf"}
    assert {Path(finding.file).name for finding in variant_report.blocked} == {
        path.name for path in VARIANTS.glob("*.gb")
    }


def _family_of(accession: str) -> str:
    with MANIFEST.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["accession"] == accession:
                return row["family"]
    raise AssertionError(f"{accession} is not in the panel manifest")


def test_a_non_iridaceae_target_is_annotated_from_its_own_family(
    panel_report, blast_tools: dict[str, str | None], tmp_path: Path
) -> None:
    require_tools(blast_tools, "blastx", "blastn", "makeblastdb", "aragorn")
    blocked = {Path(finding.file).name for finding in panel_report.blocked}
    family = _family_of(TARGET)
    references = [
        path
        for path in sorted(PANEL.glob("*.gb"))
        if path.name not in blocked
        and path.name != f"{TARGET}.gb"
        and _family_of(path.stem) == family
    ]
    assert len(references) >= 5, f"{family} has too few usable references: {len(references)}"

    reference_dir = tmp_path / "references"
    reference_dir.mkdir()
    for path in references:
        (reference_dir / path.name).write_bytes(path.read_bytes())

    out = tmp_path / "annotate"
    assert (
        main(
            [
                "annotate",
                "-i",
                str(PANEL / f"{TARGET}.gb"),
                "--reference",
                str(reference_dir),
                "--set",
                "reannotate=yes",
                "-o",
                str(out),
            ]
        )
        == 0
    )
    artifact = next(out.glob("**/plastid_annotation.gb"))
    annotated = read_genbank(artifact)

    counts = Counter(gene.kind for gene in annotated.genes)
    assert counts["CDS"] >= 80, counts
    assert counts["tRNA"] >= 35, counts
    assert counts["rRNA"] >= 4, counts
    assert annotated.accession == TARGET
