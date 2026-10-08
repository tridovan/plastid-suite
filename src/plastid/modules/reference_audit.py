"""The reference-set audit as a :class:`Result` — tables, and the report's source cells.

This file sits in ``modules/`` because that is where this project turns an analysis into tables
(seam 2: modules build data, the export layer writes files), but it is **not registered** in the
module registry and does not follow ``run(genome, params, ctx)``. The audit's input is a *set of
files*, not a genome: a record it cannot parse is one of the things it has to report, so it can
never be handed to the reader that would refuse it first. Its driver is ``plastid reference
audit`` (and, from P3B-09, ``plastid annotate --audit-only``), which owns the run directory and
hands the audit a context for the aligner call.

Nothing here decides anything: every verdict is :mod:`plastid.core.reference_audit`'s, and the
tables carry the audit's own reasons verbatim, so the report a person reads and the gate the run
obeys cannot say different things.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from plastid.core.reference_audit import AuditReport, dotplot_runs, majority_names
from plastid.core.sequence import reverse_complement
from plastid.result.result import Column, Result, Scope, Table

__all__ = ["Dotplot", "MODULE_ID", "MODULE_VERSION", "audit_result"]

MODULE_ID = "reference_audit"
MODULE_VERSION = "0.1.0"

_SCOPE = Scope(label="plastid reference set, as published")
# The region tables are in the normalized frame's region-local coordinates, not in any record's
# GenBank frame — the scope says so rather than leaving a reader to guess which zero is meant.
_REGION_SCOPE = Scope(
    label="one region per record, normalized frame (LSC → IRa → SSC → IRb)",
    coord_system="normalized",
)


@dataclass(frozen=True, slots=True)
class Dotplot:
    """The sequences a dot plot compares: one subject every query is read against.

    Held as data so the row building stays pure — the caller reads the regions
    (:func:`plastid.core.reference_audit.region_sequences`), this module only arranges them.
    ``queries`` are the records the audit called inverted; the subject is a record it did not.
    """

    region: str
    subject: str
    sequences: Mapping[str, str]
    queries: tuple[str, ...]


def audit_result(report: AuditReport, *, dotplot: Dotplot | None = None) -> Result:
    """Everything the audit found, as tables — including the picture's own coordinates.

    A table with nothing to say is **absent** rather than empty, because a table with no rows is
    not a table (the schema refuses one, and rightly: no reader can tell it from a table that
    failed to fill). What is always present is ``audit_summary`` and ``audit_records`` — the
    counts and the measurements — so a run that found nothing still says so in a cell.
    """
    candidates = {
        "audit_summary": _summary(report),
        "audit_findings": _findings(report),
        "audit_records": _records(report),
        "name_sets": _name_sets(report),
    }
    if dotplot is not None:
        candidates["region_matches"] = _region_matches(dotplot)
        candidates["region_orientation"] = _orientation(report)
    tables = {name: table for name, table in candidates.items() if table is not None}
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=report.notes
        + tuple(
            f"{status.pin.file}: exception {status.status} — {status.detail}"
            for status in report.pins
        ),
    )


def _summary(report: AuditReport) -> Table:
    """One row: the counts the report's prose reads, plus the state of the second reading."""
    confirmation = "not run (no finding needed confirming)"
    if report.confirmations:
        confirmation = "; ".join(
            f"{item.technique} over the {item.region} — "
            + (
                "measured: "
                + ", ".join(f"{name} {direction}" for name, direction in item.directions)
                if item.measured
                else (item.note or "not measured")
            )
            for item in report.confirmations
        )
    return Table(
        id="audit_summary",
        title="Reference-set audit",
        columns=(
            Column("records", "Records audited", dtype="int"),
            Column("blocked", "Blocked", dtype="int"),
            Column("reported", "Reported", dtype="int"),
            Column("confirmed", "Confirmed exceptions", dtype="int"),
            Column("confirmation", "Second reading"),
        ),
        rows=(
            {
                "records": len(report.records),
                "blocked": len(report.blocked),
                "reported": len(report.reported),
                "confirmed": len(report.confirmed),
                "confirmation": confirmation,
            },
        ),
        scope=_SCOPE,
    )


def _findings(report: AuditReport) -> Table | None:
    """Every finding, blocking ones first — the reasons verbatim, as the gate prints them."""
    order = {"BLOCK": 0, "REPORT": 1}
    rows = sorted(
        report.findings, key=lambda finding: (order[finding.verdict], finding.file, finding.check)
    )
    if not rows:
        return None
    return Table(
        id="audit_findings",
        title="Findings",
        columns=(
            Column("file", "File"),
            Column("verdict", "Verdict"),
            Column("check", "Check"),
            Column("reason", "Reason"),
            Column("confirmed_by", "Confirmed by"),
            Column("confirmed_on", "Confirmed on"),
        ),
        rows=tuple(
            {
                "file": finding.file,
                "verdict": finding.verdict,
                "check": finding.check,
                "reason": finding.reason,
                "confirmed_by": finding.confirmed_by,
                "confirmed_on": finding.confirmed_on,
            }
            for finding in rows
        ),
        scope=_SCOPE,
    )


def _records(report: AuditReport) -> Table:
    """What was measured per record — the counts and the orientation vote, as data."""
    return Table(
        id="audit_records",
        title="Records",
        columns=(
            Column("file", "File"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("cds", "CDS", dtype="int"),
            Column("trna", "tRNA", dtype="int"),
            Column("rrna", "rRNA", dtype="int"),
            Column("ssc_orientation", "SSC orientation"),
            Column("blocked", "Blocked", dtype="str"),
            Column("sha256", "sha256"),
        ),
        rows=tuple(
            {
                "file": record.file,
                "accession": record.accession,
                "length_bp": record.length_bp,
                "cds": record.cds,
                "trna": record.trna,
                "rrna": record.rrna,
                "ssc_orientation": record.ssc_orientation,
                "blocked": "yes" if record.blocked else "no",
                "sha256": record.sha256,
            }
            for record in report.records
        ),
        scope=_SCOPE,
    )


def _name_sets(report: AuditReport) -> Table | None:
    """Every name every record carries, and whether the majority carries it too.

    The ★ pick (D5-09 §7.1) is a choice between **vocabularies**, so the report has to show the
    names as the records spell them — pivot this table by ``name`` and ``in_majority`` and the
    two spellings of one locus stand next to each other, with the records that use each.
    """
    majority = majority_names(report.records)
    rows = []
    for record in report.records:
        for index, name in enumerate(sorted(record.gene_names), start=1):
            rows.append(
                {
                    "file": record.file,
                    "name": name,
                    "in_majority": _in_majority(name, majority),
                    "copy": index,
                }
            )
    if not rows:
        return None
    return Table(
        id="name_sets",
        title="Gene names as each record spells them",
        columns=(
            Column("file", "File"),
            Column("name", "Name"),
            Column("copy", "Name #", dtype="int"),
            Column("in_majority", "In the majority set"),
        ),
        rows=tuple(rows),
        scope=_SCOPE,
    )


def _in_majority(name: str, majority: frozenset[str] | None) -> str:
    """``yes`` / ``no`` against the majority set — or ``n/a`` when the set has no majority.

    Reported rather than guessed: with one record, or with the top name set tied, there is
    nothing to be in or out of, and saying ``no`` would invent a baseline (the same rule the
    finding text follows).
    """
    if majority is None:
        return "n/a"
    return "yes" if name in majority else "no"


def _orientation(report: AuditReport) -> Table | None:
    """The two readings side by side: what the gene strands said, and what the aligner said."""
    if not report.confirmations or not report.confirmations[0].directions:
        return None
    confirmation = report.confirmations[0]
    deviating = {name for name, direction in confirmation.directions if direction == "reverse"}
    rows = []
    for record in report.records:
        name = Path(record.file).name
        if name not in dict(confirmation.directions):
            continue
        rows.append(
            {
                "entry": f"{confirmation.region}:{name}",
                "file": name,
                "gene_strands": record.ssc_orientation,
                "mafft": confirmation.direction(name) or "not read",
                "agrees": "yes"
                if (record.ssc_orientation == "opposite") == (name in deviating)
                else "no",
            }
        )
    return (
        Table(
            id="region_orientation",
            title=f"{confirmation.region} orientation, two readings",
            columns=(
                Column("entry", "Entry"),
                Column("file", "File"),
                Column("gene_strands", "Gene strands"),
                Column("mafft", confirmation.technique),
                Column("agrees", "Agree"),
            ),
            rows=tuple(rows),
            scope=_REGION_SCOPE,
            key="entry",
        )
        if rows
        else None
    )


def _region_matches(dotplot: Dotplot) -> Table | None:
    """Exact-match runs of each query region against the subject, both orientations.

    Both directions are drawn on purpose: a query that runs the same way as the subject lands on
    one diagonal and an inverted one on the other, which is the whole point of the picture. No
    run at all means no picture — a table with no rows is not a table, and a figure with no
    points would be a claim that nothing matched, drawn as though it had been measured.
    """
    subject = dotplot.sequences[dotplot.subject]
    rows = []
    for query in dotplot.queries:
        sequence = dotplot.sequences[query]
        for orientation, target in (("+", subject), ("-", reverse_complement(subject))):
            for match in dotplot_runs(sequence, target):
                rows.append(
                    {
                        "region": dotplot.region,
                        "query": query,
                        "subject": dotplot.subject,
                        "orientation": orientation,
                        "query_pos": match.query_pos,
                        "subject_pos": match.subject_pos,
                        "length": match.length,
                    }
                )
    if not rows:
        return None
    rows.sort(key=lambda row: (row["query"], row["orientation"], row["subject_pos"]))
    return Table(
        id="region_matches",
        title=f"{dotplot.region} exact matches ({dotplot.subject} as the subject)",
        columns=(
            Column("region", "Region"),
            Column("query", "Query"),
            Column("subject", "Subject"),
            Column("orientation", "Orientation"),
            Column("query_pos", "Query position", dtype="int", unit="bp"),
            Column("subject_pos", "Subject position", dtype="int", unit="bp"),
            Column("length", "Run length", dtype="int", unit="bp"),
        ),
        rows=tuple(rows),
        scope=_REGION_SCOPE,
    )


def deviating_files(report: AuditReport) -> tuple[str, ...]:
    """The files the audit's orientation rule calls inverted, in report order."""
    return tuple(
        Path(finding.file).name for finding in report.reported if finding.check == "orientation"
    )


def dotplot_for(report: AuditReport, sequences: Mapping[str, str], region: str) -> Dotplot | None:
    """Pair the region sequences into a :class:`Dotplot`, or None when there is nothing to draw.

    The subject is the first record the rule did **not** call inverted: the picture is only worth
    drawing against a record the set agrees on.
    """
    deviating = deviating_files(report)
    if not deviating or not sequences:
        return None
    others = [name for name in sequences if name not in deviating]
    if not others:
        return None
    return Dotplot(
        region=region,
        subject=sorted(others)[0],
        sequences=sequences,
        queries=tuple(sorted(deviating)),
    )
