"""The reference-set audit: what blocks, what is reported, and what a pin may waive (P3B-03).

The negative cases are the point. D5-08's whole gate exists because the failures it catches are
silent ones — a record with no ``/translation`` produces a search that finds nothing rather than
an error, a partial coordinate produced **zero genes and a success exit** on the real query
(§8.15 trap 1), and two published records genuinely invert their SSC without anything being wrong
(§8.34). A file's own fixture here is the PGA annotation of *Iris hermona*: 85 CDS with no
``/translation`` and an IR that fails the gate at 26.36 % (§8.31) — a real record, and exactly
what the gate is for.
"""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from plastid.core import reference_audit as audit
from plastid.core.reference_audit import (
    BLOCK,
    REPORT,
    ExceptionPin,
    Finding,
    RecordFacts,
    audit_files,
    confirm_findings,
    dotplot_runs,
    gene_set_report,
    majority_report,
    read_exceptions,
    ssc_orientation,
    write_exception,
)

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "tests" / "data" / "NC_050833.1.gb"
CORPUS = REPO / "tests" / "data" / "corpus"
PANEL = REPO / ".panel" / "plastid-reference"

# The PGA annotation of Iris hermona: a real record this project did not annotate, whose faults
# are measured in §8.31 — 85 CDS without /translation, and an IR pair 26.36 % identical.
PGA_HERMONA = CORPUS / "Iris_hermona_NC_081931.1.pga.gb"


def _checks(report: audit.AuditReport, verdict: str | None = None) -> list[str]:
    return [
        finding.check
        for finding in report.findings
        if verdict is None or finding.verdict == verdict
    ]


def _record(
    name: str, *, orientation: str = "same", names: tuple[str, ...] = ("rbcL",)
) -> RecordFacts:
    """A fact row with no file behind it — the report rules are about facts, not about bytes."""
    return RecordFacts(
        file=f"/references/{name}",
        sha256=hashlib.sha256(name.encode()).hexdigest(),
        accession=name,
        length_bp=150_000,
        cds=80,
        trna=38,
        rrna=8,
        gene_names=names,
        ssc_orientation=orientation,
        blocked=False,
    )


# ---------------------------------------------------------------- per-record blocks


def _edited_golden(tmp_path: Path, name: str, *, drop: str, add: str = "") -> Path:
    """The contract genome with one thing changed — a real record, and one real fault.

    Editing the golden rather than hand-writing a small record is the measurement's own choice:
    a hand-written file fails the IR gate as well, so a test could not tell which rule fired.
    Here every other check passes by construction, so a finding can only be the edit.
    """
    text = GOLDEN.read_text(encoding="utf-8")
    assert drop in text, f"the golden no longer carries {drop!r} — the fixture is stale"
    edited = text.replace(drop, add, 1)
    target = tmp_path / name
    target.write_text(edited, encoding="utf-8")
    return target


def test_a_cds_without_translation_blocks(tmp_path: Path) -> None:
    """The engine searches the record's own protein (D5-08 §8.21) — no protein, no transfer."""
    record = _edited_golden(
        tmp_path,
        "no_translation.gb",
        # one qualifier key, renamed: the ORF is untouched, so only the translation can fire
        drop='                     /translation="MSCREGLMSPQTETKASVGF',
        add='                     /translation_removed="MSCREGLMSPQTETKASVGF',
    )

    report = audit_files([record], confirm=None)

    assert [finding.line for finding in report.blocked] == [f"{record} — CDS rbcL: no /translation"]
    assert not report.ok


def test_a_partial_coordinate_blocks(tmp_path: Path) -> None:
    """The class §8.34 caught: one published record, and PGA returned zero genes from it."""
    record = _edited_golden(
        tmp_path,
        "partial.gb",
        drop="     CDS             55212..56672",
        add="     CDS             <55212..56672",
    )

    report = audit_files([record], confirm=None)

    assert [finding.line for finding in report.blocked] == [
        f"{record} — CDS rbcL: partial coordinates (<55212)"
    ]


def test_a_record_with_no_gene_features_blocks(tmp_path: Path) -> None:
    """A sequence with no annotation is not a reference: there is nothing to transfer."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    record = tmp_path / "no_genes.gb"
    sequence = "ACGT" * 100
    seq_record = SeqRecord(Seq(sequence), id="X3", name="X3", annotations={"molecule_type": "DNA"})
    seq_record.features.append(
        SeqFeature(SimpleLocation(0, len(sequence)), type="source", qualifiers={"organism": "x"})
    )
    SeqIO.write(seq_record, record, "genbank")

    report = audit_files([record], confirm=None)

    assert "genes" in _checks(report, BLOCK)
    assert (
        next(finding for finding in report.blocked if finding.check == "genes").line
        == f"{record} — no gene features (no CDS, tRNA or rRNA)"
    )


def test_an_unreadable_file_blocks(tmp_path: Path) -> None:
    broken = tmp_path / "broken.gb"
    broken.write_text("this is not a GenBank file\n", encoding="utf-8")

    report = audit_files([broken], confirm=None)

    assert _checks(report, BLOCK) == ["readable"]
    assert report.blocked[0].line.startswith(f"{broken} — not a GenBank record (")
    assert report.records == ()


def test_a_missing_file_blocks_without_raising(tmp_path: Path) -> None:
    """A set that names a file nobody has is a set that cannot be audited — and it says which."""
    report = audit_files([tmp_path / "gone.gb"], confirm=None)

    assert _checks(report, BLOCK) == ["readable"]
    assert "gone.gb" in report.blocked[0].line


def test_the_real_pga_record_blocks_on_translation_and_ir(tmp_path: Path) -> None:
    """Measured 08/10/2026 on the vendored fixture: 83 of 84 CDS without /translation, IR failed.

    This is *Iris hermona* as PGA annotated it (§8.31) — the same sequence as the published
    record beside it, which passes every check. Both halves matter: the audit must block the
    annotation a search cannot use, and must not block the one it can.
    """
    report = audit_files([PGA_HERMONA], confirm=None)

    blocked = _checks(report, BLOCK)
    assert report.records[0].cds == 84
    assert blocked.count("translation") == 83
    assert "ir" in blocked
    assert report.records[0].ssc_orientation == "unknown"  # no IR, no vote

    published = audit_files([CORPUS / "Iris_hermona_NC_081931.1.published.gb"], confirm=None)
    assert published.blocked == ()
    assert published.ok


def test_the_golden_record_blocks_nothing() -> None:
    """The contract genome is a usable reference: every CDS carries its protein, IR holds."""
    report = audit_files([GOLDEN], confirm=None)

    assert report.blocked == ()
    assert report.records[0].ssc_orientation == "same"
    assert report.records[0].cds == 86


def test_the_message_is_the_file_and_the_reason_and_nothing_else(tmp_path: Path) -> None:
    """``file — reason``, and no advice: what to do about a bad reference is not our call."""
    record = _edited_golden(
        tmp_path,
        "advice.gb",
        drop='                     /translation="MSCREGLMSPQTETKASVGF',
        add='                     /translation_removed="MSCREGLMSPQTETKASVGF',
    )

    line = audit_files([record], confirm=None).blocked[0].line

    assert line == f"{record} — CDS rbcL: no /translation"
    lowered = line.lower()
    for advice in ("try", "consider", "should", "fix", "replace", "please"):
        assert advice not in lowered


def test_a_blinding_finding_cannot_be_built_as_confirmed() -> None:
    """The gate is not waivable by data: no JSON edit can produce a confirmed block."""
    with pytest.raises(ValueError, match="cannot be confirmed away"):
        Finding(
            file="x.gb",
            check="translation",
            verdict=BLOCK,
            reason="CDS rbcL: no /translation",
            confirmed_by="someone",
        )


def test_a_finding_refuses_a_check_this_audit_does_not_know() -> None:
    with pytest.raises(ValueError, match="unknown check"):
        Finding(file="x.gb", check="vibes", verdict=REPORT, reason="seems off")


# ---------------------------------------------------------------- set-level reports


def test_an_inverted_ssc_is_reported_and_never_blocks() -> None:
    """§8.34 measured why: 2/34 published records are inverted with every SSC CDS translating."""
    records = [_record(f"r{index}") for index in range(3)]
    records[2] = _record("r2", orientation="opposite")

    findings = majority_report(
        records,
        value=lambda record: record.ssc_orientation,
        check="orientation",
        statement="SSC orientation",
        ignored=("unknown",),
    )

    assert [finding.verdict for finding in findings] == [REPORT]
    assert findings[0].file.endswith("r2")
    assert "2 of 3 records read the other way" in findings[0].reason


def test_two_references_that_disagree_say_so_rather_than_calling_a_majority() -> None:
    """D5-08's own wording: with N = 2 there is no majority to be right about."""
    records = [_record("a"), _record("b", orientation="opposite")]

    findings = majority_report(
        records,
        value=lambda record: record.ssc_orientation,
        check="orientation",
        statement="SSC orientation",
        ignored=("unknown",),
    )

    assert len(findings) == 2
    assert all("the two references disagree" in finding.reason for finding in findings)
    assert all("N = 2" in finding.reason for finding in findings)


def test_a_set_that_agrees_with_itself_reports_nothing() -> None:
    records = [_record(f"r{index}") for index in range(5)]

    assert (
        majority_report(
            records,
            value=lambda record: record.ssc_orientation,
            check="orientation",
            statement="SSC orientation",
        )
        == []
    )


def test_an_unknown_vote_is_not_a_verdict() -> None:
    """A record with no IR has no orientation — it is left out, not counted as agreeing."""
    records = [_record("a"), _record("b"), _record("c", orientation="unknown")]

    findings = majority_report(
        records,
        value=lambda record: record.ssc_orientation,
        check="orientation",
        statement="SSC orientation",
        ignored=("unknown",),
    )

    assert findings == []


def test_a_gene_set_that_differs_is_named_name_by_name() -> None:
    """The report is about spelling, because that is what the ★ pick decides (D5-09 §7.1)."""
    records = [
        _record("a", names=("rbcL", "ndhI", "trnL-UAA")),
        _record("b", names=("rbcL", "ndhI", "trnL-UAA")),
        _record("c", names=("rbcL", "ndhL", "trnL")),
    ]

    findings = gene_set_report(records)

    assert len(findings) == 1
    assert "missing 2: ndhI, trnL-UAA" in findings[0].reason
    assert "extra 2: ndhL, trnL" in findings[0].reason


def test_a_tied_gene_set_reports_that_there_is_no_majority() -> None:
    records = [
        _record("a", names=("rbcL",)),
        _record("b", names=("rbcL", "ndhI")),
    ]

    findings = gene_set_report(records)

    assert len(findings) == 2
    assert all("no majority name set" in finding.reason for finding in findings)


# ---------------------------------------------------------------- exceptions


def _pin(record: RecordFacts, check: str = "orientation") -> ExceptionPin:
    return ExceptionPin(
        file=record.file,
        sha256=record.sha256,
        check=check,
        reason="confirmed flip-flop with clean ORFs",
        confirmed_by="the dev, from the dot plot",
        confirmed_on="2026-10-08",
    )


def test_a_live_pin_marks_the_finding_confirmed_and_leaves_it_visible() -> None:
    record = _record("a", orientation="opposite")
    finding = Finding(
        file=record.file,
        check="orientation",
        verdict=REPORT,
        reason="SSC orientation differs from the majority (2 of 3 records read the other way)",
    )

    marked, statuses = confirm_findings([finding], [record], [record.file], [_pin(record)])

    assert [status.status for status in statuses] == ["applied"]
    assert marked[0].confirmed
    assert marked[0].confirmed_by == "the dev, from the dot plot"
    assert marked[0].line == finding.line  # the message never changes


def test_changing_the_file_voids_its_exception() -> None:
    """D5-08: the pin is bound to the file's digest — change one byte and the audit starts over."""
    record = _record("a", orientation="opposite")
    finding = Finding(
        file=record.file, check="orientation", verdict=REPORT, reason="SSC orientation differs"
    )
    edited = replace(record, sha256="f" * 64)

    marked, statuses = confirm_findings([finding], [edited], [record.file], [_pin(record)])

    assert [status.status for status in statuses] == ["void"]
    assert not marked[0].confirmed


def test_a_pin_for_a_file_that_is_not_in_the_set_is_stale() -> None:
    record = _record("a")
    _, statuses = confirm_findings([], [record], [record.file], [_pin(_record("elsewhere"))])

    assert [status.status for status in statuses] == ["stale"]


def test_a_pin_for_a_check_the_audit_did_not_raise_is_unused() -> None:
    record = _record("a")
    _, statuses = confirm_findings([], [record], [record.file], [_pin(record, check="gene_set")])

    assert [status.status for status in statuses] == ["unused"]


def test_an_exception_file_round_trips_and_refuses_a_block_check(tmp_path: Path) -> None:
    path = tmp_path / "AUDIT_EXCEPTIONS.json"
    record_file = tmp_path / "record.gb"
    record_file.write_text("LOCUS       X5                     12 bp    DNA     linear\n//\n")

    pin = write_exception(
        path,
        file=record_file,
        check="orientation",
        reason="flip-flop, ORFs clean",
        confirmed_by="dev",
        confirmed_on="2026-10-08",
    )

    assert pin.sha256 == hashlib.sha256(record_file.read_bytes()).hexdigest()
    assert read_exceptions(path) == (pin,)
    with pytest.raises(ValueError, match="not a reportable check"):
        write_exception(
            path,
            file=record_file,
            check="translation",
            reason="who cares",
            confirmed_by="dev",
            confirmed_on="2026-10-08",
        )
    assert read_exceptions(tmp_path / "nothing-here.json") == ()


def test_a_malformed_exception_is_refused_rather_than_skipped(tmp_path: Path) -> None:
    path = tmp_path / "AUDIT_EXCEPTIONS.json"
    path.write_text('{"exceptions": [{"file": "a.gb", "check": "orientation"}]}', encoding="utf-8")

    with pytest.raises(ValueError, match="missing"):
        read_exceptions(path)


# ---------------------------------------------------------------- the dot plot's data


def _random_sequence(length: int, seed: int = 7) -> str:
    """A deterministic pseudo-random sequence — no repeated motif to match by accident."""
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def test_dotplot_runs_keeps_only_run_starts_and_extends_them() -> None:
    """Two conserved blocks around one mismatch are two rows, not one per matching position."""
    # Two *different* conserved blocks, so each word occurs once on each side: with a repeated
    # motif the cross matches are real and the test would be measuring the repeat, not the run.
    subject = "ACGTTGCA" + "GGGG" + "GGATCCTA"
    query = "ACGTTGCA" + "TTTT" + "GGATCCTA"

    runs = dotplot_runs(query, subject, word=8)

    assert [(match.query_pos, match.subject_pos, match.length) for match in runs] == [
        (1, 1, 8),
        (13, 13, 8),
    ]


def test_dotplot_runs_finds_the_anti_diagonal_of_an_inverted_sequence() -> None:
    """The picture's whole point: an inverted region matches the *reverse complement*.

    On a sequence with no internal repeat, an inverted query has no forward match at all — the
    same statement the audit's gene-strand rule makes, drawn.
    """
    from plastid.core.sequence import reverse_complement

    subject = _random_sequence(200)
    inverted = reverse_complement(subject)

    assert dotplot_runs(inverted, subject, word=15) == ()
    assert [
        (match.query_pos, match.subject_pos, match.length)
        for match in dotplot_runs(inverted, reverse_complement(subject), word=15)
    ] == [(1, 1, 200)]


# ---------------------------------------------------------------- the real corpus


@pytest.mark.skipif(not PANEL.is_dir(), reason="the pinned reference set lives outside the repo")
def test_the_pinned_set_blocks_on_the_partial_record_and_reports_the_two_inversions() -> None:
    """The Done-when, measured 08/10/2026: 32/34 one SSC orientation, 2/34 inverted, reported.

    Blocking the two inverted records would refuse two published plastomes, one of them an
    *Iris*; blocking the partial record is the price of a search that needs coordinates, and
    §8.34 predicted exactly this one.
    """
    paths = sorted(PANEL.glob("*.gb"))
    report = audit_files(paths, confirm=None)

    assert [Path(finding.file).name for finding in report.blocked] == ["NC_024936.1.gb"] * 2
    assert {finding.check for finding in report.blocked} == {"partial"}

    orientations = {Path(record.file).name: record.ssc_orientation for record in report.records}
    assert sum(1 for value in orientations.values() if value == "same") == 32
    assert {name for name, value in orientations.items() if value == "opposite"} == {
        "NC_042827.1.gb",
        "NC_062351.1.gb",
    }
    assert {
        Path(finding.file).name for finding in report.reported if finding.check == "orientation"
    } == {
        "NC_042827.1.gb",
        "NC_062351.1.gb",
    }
    assert len(report.records) == 34


# ---------------------------------------------------------------- the door


def test_the_cli_blocks_a_set_and_the_manifest_carries_the_reason(tmp_path: Path) -> None:
    """The Done-when's negative case: exit 2, nothing annotated, the reason in the manifest."""
    import json

    from plastid.cli import EXIT_USAGE, main

    panel = tmp_path / "refs"
    panel.mkdir()
    (panel / "Iris_hermona_pga.gb").write_bytes(PGA_HERMONA.read_bytes())
    out = tmp_path / "runs"

    code = main(["reference", "audit", "--panel", str(panel), "--out", str(out), "--no-confirm"])

    assert code == EXIT_USAGE
    run_dir = next((out / "refs").iterdir())
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert "no /translation" in manifest["modules"][0]["error"]
    assert manifest["modules"][0]["status"] == "failed"
    assert (
        manifest["inputs"][0]["sha256"]
        == hashlib.sha256((panel / "Iris_hermona_pga.gb").read_bytes()).hexdigest()
    )
    # nothing was annotated: an audit writes tables and a report, never a record
    assert list(run_dir.rglob("*.gb")) == []


def test_the_cli_audits_a_clean_set_and_writes_the_report(tmp_path: Path) -> None:
    """A usable set exits 0 and leaves a report a person can pick the ★ record from."""
    from plastid.cli import EXIT_OK, main

    panel = tmp_path / "refs"
    panel.mkdir()
    (panel / "golden.gb").write_bytes(GOLDEN.read_bytes())
    out = tmp_path / "runs"

    code = main(["reference", "audit", "--panel", str(panel), "--out", str(out), "--no-confirm"])

    assert code == EXIT_OK
    run_dir = next((out / "refs").iterdir())
    report = (run_dir / "01_reference_audit" / "report.md").read_text(encoding="utf-8")
    assert "Reference-set audit" in report
    assert "not run (no finding needed confirming)" in report
    assert (run_dir / "01_reference_audit" / "tables" / "name_sets.csv").is_file()


def test_the_confirm_exception_door_pins_a_digest_the_person_never_types(tmp_path: Path) -> None:
    """The person confirms; the tool writes the digest. A typed digest is a pin that may be void."""
    from plastid.cli import EXIT_OK, main

    panel = tmp_path / "refs"
    panel.mkdir()
    record = panel / "record.gb"
    record.write_bytes(GOLDEN.read_bytes())

    code = main(
        [
            "reference",
            "audit",
            "--panel",
            str(panel),
            "--confirm-exception",
            "record.gb:orientation",
            "--confirmed-by",
            "the dev",
            "--reason",
            "flip-flop, ORFs clean",
        ]
    )

    assert code == EXIT_OK
    pins = read_exceptions(panel / audit.EXCEPTIONS_FILENAME)
    assert len(pins) == 1
    assert pins[0].sha256 == hashlib.sha256(record.read_bytes()).hexdigest()
    assert pins[0].check == "orientation"


def test_the_audit_result_is_tables_a_report_can_read() -> None:
    """The report's cells: the summary counts, the findings verbatim, the names as spelled."""
    from plastid.modules import reference_audit as audit_module

    report = audit_files([GOLDEN], confirm=None)
    result = audit_module.audit_result(report)

    assert result.module == audit_module.MODULE_ID
    summary = result.tables["audit_summary"]
    assert summary.value("records") == 1
    assert summary.value("blocked") == 0
    assert summary.value("reported") == 0
    assert result.tables["audit_records"].value("cds") == 86
    names = {row["name"] for row in result.tables["name_sets"].rows}
    assert "rbcL" in names and "trnG-UCC" in names


def test_the_dotplot_figure_draws_the_matches_it_was_given() -> None:
    """A figure reads cells and recomputes nothing (D5-03) — so it needs the table, not the data."""
    pytest.importorskip("matplotlib", reason="the 'figures' extra is not installed")

    from plastid.modules import reference_audit as audit_module
    from plastid.result.figures import dotplot, get

    subject = _random_sequence(400)
    from plastid.core.sequence import reverse_complement

    report = audit_files([GOLDEN], confirm=None)
    inverted = replace(report.records[0], ssc_orientation="opposite", file="/references/inv.gb")
    report = replace(
        report,
        findings=(
            *report.findings,
            Finding(
                file=inverted.file,
                check="orientation",
                verdict=REPORT,
                reason="SSC orientation differs from the majority",
            ),
        ),
        records=(report.records[0], inverted),
    )
    sequences = {Path(report.records[0].file).name: subject, "inv.gb": reverse_complement(subject)}
    data = audit_module.dotplot_for(report, sequences, "SSC")
    assert data is not None and data.queries == ("inv.gb",)

    figure = dotplot.build({"reference_audit": audit_module.audit_result(report, dotplot=data)})

    assert figure.source_tables == ("reference_audit.region_matches",)
    assert figure.formats == ("svg", "pdf", "eps")
    assert figure.payload("svg").startswith(b"<?xml")
    spec = get("dotplot")
    assert spec.requires == (("reference_audit", "region_matches"),)
    assert spec.missing({"reference_audit": audit_module.audit_result(report)}) == (
        "reference_audit.region_matches",
    )


@pytest.mark.skipif(not PANEL.is_dir(), reason="the pinned reference set lives outside the repo")
def test_the_gene_strand_rule_reproduces_the_measured_verdicts() -> None:
    """The rule as §8.34 ran it — same two records, from the same comparison of region votes."""
    from plastid.core.genbank import read_genbank
    from plastid.core.ir import detect_ir, enforce_ir_gate

    verdicts = {}
    for path in sorted(PANEL.glob("*.gb")):
        genome = read_genbank(path)
        verdicts[path.name] = ssc_orientation(genome, enforce_ir_gate(detect_ir(genome)))

    assert sorted(name for name, value in verdicts.items() if value == "opposite") == [
        "NC_042827.1.gb",
        "NC_062351.1.gb",
    ]
