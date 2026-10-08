"""The reference-set audit — the gate in front of every annotation run (P3B-03, D5-08).

A plastome is annotated by transferring a curated annotation from relatives
(:mod:`plastid.core.plastid_reference`). That only works if the relatives themselves are usable:
a record whose CDS carry no ``/translation`` cannot be searched (D5-08 §8.21 — the engine never
self-translates), a record with partial coordinates stamps partial flags on whatever it annotates
(§8.15 trap 1, where PGA returned **zero** genes from exactly such a record), and a record whose
IR pair does not hold up has no frame to normalize positions against. So the set is audited
*before* anything is annotated, and this module is that audit.

Two verdicts, and the difference is measured rather than stylistic (D5-08, 05/10/2026):

* **BLOCK** — this record cannot be used. The run stops with ``file — reason``, and nothing is
  annotated. The classes are: unreadable · no gene features · a CDS without ``/translation`` ·
  a gene feature with partial coordinates · the IR gate · a CDS whose own sequence fails the ORF
  check.
* **REPORT** — the set is not uniform and a person has to know. Measured in §8.34: 32/34
  published records share one SSC orientation and **2/34 are inverted with all their SSC CDS
  translating cleanly** — the plastome's known flip-flop, not an error, and blocking would refuse
  two published records, one of them an *Iris*. Gene-set deviations are reported for the same
  reason: presence genuinely varies between annotations (§8.31: ``ycf1`` 24/34). A majority is
  only called with **N ≥ 3**; with N = 2 the report says the two references disagree.

**A confirmed exception is pinned to bytes.** The REPORT classes are the ones no rule settles
alone, so a person confirms them once: the confirmation is written into
``AUDIT_EXCEPTIONS.json`` beside the reference set **with the file's sha256**. Change one byte of
the file and the pin no longer matches — the exception is void and the finding comes back. A pin
may never waive a BLOCK check: that file cannot be annotated, and pretending otherwise would fail
later, deeper, and without a reason.

**The second reading is injected, not imported.** ``core/`` may not reach into ``external/``
(the architecture guard), and the aligner is an optional tool anyway: the caller passes the
confirmation callable in — :func:`plastid.external.mafft.adjust_direction` from the CLI — or
passes nothing, and the report then says the finding is *not measured* rather than agreeing with
itself.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, IRResult, detect_ir, enforce_ir_gate
from plastid.core.model import GENE_KINDS, Gene, Genome
from plastid.core.normalize import normalize, regions_from_ir
from plastid.core.orf import check_genome

__all__ = [
    "BLOCK",
    "BLOCK_CHECKS",
    "CHECKS",
    "EXCEPTIONS_FILENAME",
    "REPORT",
    "REPORT_CHECKS",
    "AuditReport",
    "Confirmation",
    "ExceptionPin",
    "Finding",
    "Match",
    "PinStatus",
    "RecordFacts",
    "audit_files",
    "confirm_findings",
    "dotplot_runs",
    "gene_set_report",
    "majority_names",
    "majority_report",
    "read_exceptions",
    "region_sequences",
    "ssc_orientation",
    "write_exception",
]

BLOCK = "BLOCK"
REPORT = "REPORT"

# The checks this audit knows, in report order. A pin names one of them, and the split is data so
# "a pin may not waive a block" is a lookup rather than a second rule.
BLOCK_CHECKS: tuple[str, ...] = ("readable", "genes", "translation", "partial", "ir", "orf")
REPORT_CHECKS: tuple[str, ...] = ("orientation", "gene_set", "confirmation")
CHECKS: tuple[str, ...] = BLOCK_CHECKS + REPORT_CHECKS

# The file a confirmed exception is pinned in, beside the reference set it is about.
EXCEPTIONS_FILENAME = "AUDIT_EXCEPTIONS.json"

# The dot plot's word length: long enough that a match is not chance (4^15), short enough to
# survive the divergence between two species of one family.
DOTPLOT_WORD = 15

# A second reading of one region, exactly as :func:`plastid.external.mafft.adjust_direction`
# answers: rows in, name → was-reversed out.
ConfirmFn = Callable[[Sequence[tuple[str, str]], object], Mapping[str, bool]]


@dataclass(frozen=True, slots=True)
class Finding:
    """One audit result about one file: the file, the check, the verdict, and the reason.

    ``line`` is the message the design fixes — ``file — reason``, specific to the feature or
    coordinate, and nothing else. No suggested fix is ever appended: the audit says what is wrong
    with the bytes, and what to do about it belongs to whoever holds the reference set (D5-08).
    A ``confirmed_by`` finding was waived by a pinned exception; it stays in the report, so the
    waiver is visible rather than silent.
    """

    file: str
    check: str
    verdict: str
    reason: str
    confirmed_by: str | None = None
    confirmed_on: str | None = None

    def __post_init__(self) -> None:
        if self.verdict not in (BLOCK, REPORT):
            raise ValueError(f"verdict must be {BLOCK} or {REPORT}, got {self.verdict!r}")
        if self.check not in CHECKS:
            raise ValueError(f"unknown check {self.check!r} (have: {', '.join(CHECKS)})")
        if self.verdict == BLOCK and self.confirmed_by is not None:
            raise ValueError(
                f"{self.file}: a {BLOCK} finding cannot be confirmed away — {self.check} is a "
                f"check only a fixed file passes"
            )

    @property
    def line(self) -> str:
        """The message, whole: ``file — reason``. Nothing else is ever appended."""
        return f"{self.file} — {self.reason}"

    @property
    def blocks(self) -> bool:
        return self.verdict == BLOCK

    @property
    def confirmed(self) -> bool:
        return self.confirmed_by is not None


@dataclass(frozen=True, slots=True)
class RecordFacts:
    """What the audit measured about one record, so the set-level reports are readable as data."""

    file: str
    sha256: str
    accession: str | None
    length_bp: int
    cds: int
    trna: int
    rrna: int
    gene_names: tuple[str, ...]
    ssc_orientation: str  # "same" | "opposite" | "unknown"
    blocked: bool


@dataclass(frozen=True, slots=True)
class Confirmation:
    """A second reading of one region: which records came out reverse-oriented, by which method.

    ``directions`` maps file → ``"forward"`` / ``"reverse"`` as the aligner read them, relative to
    the first row it was given; the audit always puts a record the rule calls *forward* first, so
    "reverse" and "the rule says inverted" are the same statement. ``measured`` is False when no
    aligner was configured — the rule's finding then stands alone, which is what the report says.
    """

    region: str
    technique: str
    measured: bool
    directions: tuple[tuple[str, str], ...] = ()
    note: str | None = None

    def direction(self, file: str) -> str | None:
        for name, direction in self.directions:
            if Path(name).name == Path(file).name:
                return direction
        return None


@dataclass(frozen=True, slots=True)
class ExceptionPin:
    """A confirmed exception, bound to the file's sha256 — change the bytes, void the pin."""

    file: str
    sha256: str
    check: str
    reason: str
    confirmed_by: str
    confirmed_on: str


@dataclass(frozen=True, slots=True)
class PinStatus:
    """What became of one pin: applied · void (bytes changed) · stale · unused."""

    pin: ExceptionPin
    status: str
    detail: str


@dataclass(frozen=True, slots=True)
class Match:
    """One exact word match, extended to a maximal run: 1-based positions within the region."""

    query_pos: int
    subject_pos: int
    length: int


@dataclass(frozen=True, slots=True)
class AuditReport:
    """The audit's whole answer: what is wrong, what a person must know, and what was measured."""

    findings: tuple[Finding, ...]
    records: tuple[RecordFacts, ...]
    confirmations: tuple[Confirmation, ...] = ()
    pins: tuple[PinStatus, ...] = ()

    @property
    def blocked(self) -> tuple[Finding, ...]:
        return tuple(finding for finding in self.findings if finding.blocks)

    @property
    def reported(self) -> tuple[Finding, ...]:
        return tuple(
            finding
            for finding in self.findings
            if finding.verdict == REPORT and not finding.confirmed
        )

    @property
    def confirmed(self) -> tuple[Finding, ...]:
        return tuple(finding for finding in self.findings if finding.confirmed)

    @property
    def notes(self) -> tuple[str, ...]:
        """Statements about the measurement itself — what was not measured, and why."""
        return tuple(
            f"{confirmation.technique}: {confirmation.note}"
            for confirmation in self.confirmations
            if confirmation.note
        )

    @property
    def ok(self) -> bool:
        """True when nothing blocks — REPORT findings do not stop a run (measured, §8.34)."""
        return not self.blocked

    def lines(self) -> tuple[str, ...]:
        return tuple(finding.line for finding in self.findings)


def audit_files(
    paths: Iterable[str | Path],
    *,
    ctx: object = None,
    exceptions: Sequence[ExceptionPin] = (),
    confirm: ConfirmFn | None = None,
    confirm_note: str | None = None,
) -> AuditReport:
    """Audit a reference set: every file gets its findings, then the set gets its reports.

    ``ctx`` is a :class:`plastid.run.context.RunContext`, handed straight to ``confirm`` so the
    aligner call is recorded like every other external call. ``confirm`` is
    :func:`plastid.external.mafft.adjust_direction`, or None — and None is not a failure: the
    minimum install has no aligner, so the orientation finding is reported as **not measured**,
    with ``confirm_note`` saying why. It is never silently reported as agreed.
    """
    ordered = [Path(path) for path in paths]
    findings: list[Finding] = []
    facts: list[RecordFacts] = []
    for path in ordered:
        record_findings, recorded = _audit_one(path)
        findings.extend(record_findings)
        if recorded is not None:
            facts.append(recorded)

    # The set-level reports come first, because a pin may name one of *them* — an exception
    # applied before the finding it is about exists would come back "unused" and quietly stop
    # working (measured 08/10/2026: the first end-to-end run of the exception door did exactly
    # that). Then the pins, then the second reading — which has nothing left to settle for a
    # finding a person has already confirmed.
    findings.extend(_set_reports(facts))
    findings, pins = confirm_findings(findings, facts, ordered, exceptions)

    confirmations = _second_reading(
        facts,
        [
            finding
            for finding in findings
            if finding.check == "orientation" and not finding.confirmed
        ],
        ctx=ctx,
        confirm=confirm,
        confirm_note=confirm_note,
    )
    if confirmations:
        findings.extend(_confirmation_findings(facts, confirmations))
        # a person may accept a disagreement too — "the aligner is confused here, the gene
        # strands are right" is a judgement the audit has no business overruling
        findings, pins = confirm_findings(findings, facts, ordered, exceptions)
    return AuditReport(
        findings=tuple(findings),
        records=tuple(facts),
        confirmations=tuple(confirmations),
        pins=tuple(pins),
    )


# ---------------------------------------------------------------- per-record checks


def _audit_one(path: Path) -> tuple[list[Finding], RecordFacts | None]:
    findings: list[Finding] = []
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else ""
    try:
        genome = read_genbank(path)
    except (OSError, ValueError) as error:
        return [_block(path, "readable", f"not a GenBank record ({error})")], None

    if not [gene for gene in genome.genes if gene.kind in GENE_KINDS]:
        findings.append(_block(path, "genes", "no gene features (no CDS, tRNA or rRNA)"))

    findings.extend(_translation_findings(path, genome))
    findings.extend(_partial_findings(path, genome))

    ir, ir_finding = _ir_for(path, genome)
    if ir_finding is not None:
        findings.append(ir_finding)
    findings.extend(_orf_findings(path, genome))

    orientation = ssc_orientation(genome, ir) if ir is not None else "unknown"
    counts = Counter(gene.kind for gene in genome.genes)
    facts = RecordFacts(
        file=str(path),
        sha256=digest,
        accession=genome.accession,
        length_bp=len(genome),
        cds=counts.get("CDS", 0),
        trna=counts.get("tRNA", 0),
        rrna=counts.get("rRNA", 0),
        gene_names=tuple(gene.name for gene in genome.genes),
        ssc_orientation=orientation,
        blocked=any(finding.blocks for finding in findings),
    )
    return findings, facts


def _translation_findings(path: Path, genome: Genome) -> list[Finding]:
    """A CDS with no ``/translation`` cannot be searched — the engine never self-translates."""
    return [
        _block(path, "translation", f"CDS {gene.name}: no /translation")
        for gene in genome.genes
        if gene.kind == "CDS" and not gene.translation
    ]


def _partial_findings(path: Path, genome: Genome) -> list[Finding]:
    """A partial coordinate in a reference stamps partial flags on whatever it annotates."""
    findings: list[Finding] = []
    for gene in genome.genes:
        if gene.kind not in GENE_KINDS:
            continue
        marks: list[str] = []
        for part in gene.parts:
            start_1b, end_1b = part.to_1based()
            if part.partial_start:
                marks.append(f"<{start_1b}")
            if part.partial_end:
                marks.append(f">{end_1b}")
        if marks:
            findings.append(
                _block(
                    path,
                    "partial",
                    f"{gene.kind} {gene.name}: partial coordinates ({', '.join(marks)})",
                )
            )
    return findings


def _ir_for(path: Path, genome: Genome) -> tuple[IRResult | None, Finding | None]:
    try:
        return enforce_ir_gate(detect_ir(genome)), None
    except (IRGateError, IRDetectionError) as error:
        return None, _block(path, "ir", f"IR gate failed: {error}")


def _orf_findings(path: Path, genome: Genome) -> list[Finding]:
    """A CDS whose own sequence does not read as protein is not a model to transfer."""
    return [
        _block(path, "orf", f"CDS {check.gene}: {check.status}")
        for check in check_genome(genome)
        if check.is_failure
    ]


def _block(path: Path, check: str, reason: str) -> Finding:
    return Finding(file=str(path), check=check, verdict=BLOCK, reason=reason)


# ---------------------------------------------------------------- per-set reports


def _set_reports(facts: Sequence[RecordFacts]) -> list[Finding]:
    """The two majority reports — orientation and gene set — and nothing else.

    Both are measured over **every record that could be read**, blocked ones included: a blocked
    record is still a published record, its gene strands still say which way its SSC runs, and
    §8.34's 32/34 is a statement about the corpus rather than about the records that survived
    this audit. Usability is the block list's business.
    """
    findings = majority_report(
        facts,
        value=lambda record: record.ssc_orientation,
        check="orientation",
        statement="SSC orientation",
        ignored=("unknown",),
    )
    return findings + gene_set_report(facts)


def _second_reading(
    facts: Sequence[RecordFacts],
    deviating: Sequence[Finding],
    *,
    ctx: object,
    confirm: ConfirmFn | None,
    confirm_note: str | None,
) -> list[Confirmation]:
    """The confirmations to record — none when the set has nothing left to settle."""
    names = [Path(finding.file).name for finding in deviating]
    if not names:
        return []
    if confirm is None:
        return [
            Confirmation(
                region="SSC",
                technique="mafft --adjustdirection",
                measured=False,
                note=confirm_note or "not measured (no aligner configured)",
            )
        ]
    rows, exemplar = _ssc_rows(facts, names)
    if not exemplar:
        return [
            Confirmation(
                region="SSC",
                technique="mafft --adjustdirection",
                measured=False,
                note="not measured (no record the rule calls forward to align against)",
            )
        ]
    reversed_names = confirm(rows, ctx)
    return [
        Confirmation(
            region="SSC",
            technique="mafft --adjustdirection",
            measured=True,
            directions=tuple(
                (name, "reverse" if reversed_names.get(name) else "forward") for name, _ in rows
            ),
        )
    ]


def _confirmation_findings(
    facts: Sequence[RecordFacts], confirmations: Sequence[Confirmation]
) -> list[Finding]:
    """Where the second reading contradicts the rule — reported, never settled by a score."""
    out: list[Finding] = []
    for confirmation in confirmations:
        directions = dict(confirmation.directions)
        for name, direction in confirmation.directions:
            record = next((item for item in facts if Path(item.file).name == Path(name).name), None)
            if record is None or record.ssc_orientation != "opposite":
                continue
            if direction != "reverse":
                out.append(
                    Finding(
                        file=record.file,
                        check="confirmation",
                        verdict=REPORT,
                        reason="SSC orientation disagrees between the gene strands and "
                        f"mafft --adjustdirection ({directions.get(name, 'not read')})",
                    )
                )
    return out


def majority_report(
    records: Sequence[RecordFacts],
    *,
    value: Callable[[RecordFacts], str],
    check: str,
    statement: str,
    ignored: Sequence[str] = (),
) -> list[Finding]:
    """Findings for the records that differ from the most common value — the N ≥ 3 rule.

    "Majority" is the value most records carry; a tie at the top is reported as *no majority*
    rather than broken arbitrarily, and with two records there is nothing to be a majority of —
    the report says the two references disagree, naming both (D5-08's own wording).
    """
    considered = [record for record in records if value(record) not in ignored]
    if len(considered) < 2:
        return []
    tally = Counter(value(record) for record in considered)
    top, top_count = tally.most_common(1)[0]
    tied = [name for name, count in tally.items() if count == top_count]
    if len(tied) == 1 and top_count == len(considered):
        return []  # the set agrees with itself
    if len(considered) == 2:
        names = ", ".join(Path(record.file).name for record in considered)
        return [
            Finding(
                file=record.file,
                check=check,
                verdict=REPORT,
                reason=f"the two references disagree about the {statement} ({names}) — "
                f"N = 2 is too few to call a majority",
            )
            for record in considered
        ]
    if len(tied) > 1:
        spread = ", ".join(f"{name}: {tally[name]}" for name in sorted(tied))
        return [
            Finding(
                file=record.file,
                check=check,
                verdict=REPORT,
                reason=f"no majority over {len(considered)} records for the {statement} ({spread})",
            )
            for record in considered
        ]
    return [
        Finding(
            file=record.file,
            check=check,
            verdict=REPORT,
            reason=f"{statement} differs from the majority "
            f"({top_count} of {len(considered)} records read the other way)",
        )
        for record in considered
        if value(record) != top
    ]


def gene_set_report(records: Sequence[RecordFacts]) -> list[Finding]:
    """Name-set deviations from the most common set — compared as names, never as loci.

    Deliberately a statement about *spelling*: which names a record uses. D5-09 settles which
    vocabulary the project writes (the ★ record) and forbids mixing the two, so the audit may not
    presuppose that answer by canonicalizing first. It reports what the set says, which is also
    what the ★ pick has to be made from (D5-09 §7.1).
    """
    if len(records) < 2:
        return []
    majority = majority_names(records)
    if majority is None:
        spread = Counter(frozenset(record.gene_names) for record in records)
        top = spread.most_common(1)[0][1]
        return [
            Finding(
                file=record.file,
                check="gene_set",
                verdict=REPORT,
                reason=f"no majority name set over {len(records)} records "
                f"(the largest set is held by {top})",
            )
            for record in records
        ]
    findings: list[Finding] = []
    for record in records:
        names = set(record.gene_names)
        missing = sorted(majority - names)
        extra = sorted(names - majority)
        if not missing and not extra:
            continue
        findings.append(
            Finding(
                file=record.file,
                check="gene_set",
                verdict=REPORT,
                reason=f"gene set differs from the majority of {len(records)} records — "
                f"missing {len(missing)}: {_join(missing)}; extra {len(extra)}: {_join(extra)}",
            )
        )
    return findings


def majority_names(records: Sequence[RecordFacts]) -> frozenset[str] | None:
    """The name set most records carry; None when the top is tied (or only one record)."""
    if len(records) < 2:
        return None
    tally = Counter(frozenset(record.gene_names) for record in records)
    top, top_count = tally.most_common(1)[0]
    if sum(1 for count in tally.values() if count == top_count) > 1:
        return None
    return top


def _join(names: Sequence[str]) -> str:
    return ", ".join(names) if names else "none"


# ---------------------------------------------------------------- the second reading


def _ssc_rows(
    records: Sequence[RecordFacts], deviating: Sequence[str]
) -> tuple[tuple[tuple[str, str], ...], str | None]:
    """The SSC sequences to align: a forward record first, so ``_R_`` means "inverted".

    MAFFT orients every sequence against the **first** row it is given, so a record the rule calls
    forward is put there — that makes the aligner's marker and the rule's verdict the same
    statement, and the comparison between them meaningful.
    """
    by_name = {Path(record.file).name: record for record in records}
    wanted = [name for name in deviating if name in by_name]
    exemplar = next(
        (
            Path(record.file).name
            for record in records
            if record.ssc_orientation == "same" and Path(record.file).name not in wanted
        ),
        None,
    )
    if exemplar is None:
        return (), None
    rows: list[tuple[str, str]] = []
    for name in [exemplar, *wanted]:
        record = by_name[name]
        genome = read_genbank(record.file)
        ir = enforce_ir_gate(detect_ir(genome))
        rows.append((name, _region_sequence(genome, ir, "SSC")))
    return tuple(rows), exemplar


# ---------------------------------------------------------------- exceptions


def read_exceptions(path: str | Path) -> tuple[ExceptionPin, ...]:
    """Read the pinned exceptions; a malformed pin raises rather than being skipped.

    Silently dropping a pin would make the audit report findings the person believes they have
    already dealt with — and make a typo in the file invisible.
    """
    target = Path(path)
    if not target.is_file():
        return ()
    payload = json.loads(target.read_text(encoding="utf-8"))
    entries = payload["exceptions"] if isinstance(payload, dict) else payload
    required = ("file", "sha256", "check", "reason", "confirmed_by", "confirmed_on")
    pins: list[ExceptionPin] = []
    for entry in entries:
        missing = [key for key in required if not entry.get(key)]
        if missing:
            raise ValueError(f"{target}: an exception is missing {', '.join(missing)}")
        if entry["check"] not in REPORT_CHECKS:
            raise ValueError(
                f"{target}: {entry['check']!r} is not a reportable check — an exception can only "
                f"confirm one of {', '.join(REPORT_CHECKS)}"
            )
        pins.append(ExceptionPin(**{key: str(entry[key]) for key in required}))
    return tuple(pins)


def write_exception(
    path: str | Path,
    *,
    file: str | Path,
    check: str,
    reason: str,
    confirmed_by: str,
    confirmed_on: str,
) -> ExceptionPin:
    """Pin a confirmed exception to the file's **current** sha256, appending to the file.

    Writing the digest here rather than accepting one is the point: a hand-typed sha256 is a pin
    that may never have matched the bytes, which is the same as no pin at all until it matters.
    """
    if check not in REPORT_CHECKS:
        raise ValueError(
            f"{check!r} is not a reportable check — an exception can only confirm one of "
            f"{', '.join(REPORT_CHECKS)}"
        )
    target = Path(file)
    if not target.is_file():
        raise ValueError(f"{target} is not a file — there is nothing to pin an exception to")
    pin = ExceptionPin(
        file=str(target),
        sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        check=check,
        reason=reason,
        confirmed_by=confirmed_by,
        confirmed_on=confirmed_on,
    )
    destination = Path(path)
    kept = [
        other
        for other in read_exceptions(destination)
        if not (Path(other.file).name == target.name and other.check == check)
    ]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(
            {"exceptions": [asdict(entry) for entry in [*kept, pin]]},
            indent=1,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return pin


def confirm_findings(
    findings: Sequence[Finding],
    records: Sequence[RecordFacts],
    paths: Sequence[str | Path],
    pins: Sequence[ExceptionPin],
) -> tuple[list[Finding], list[PinStatus]]:
    """Apply the pins: mark what a live pin covers, and say where every other pin stands.

    A pin is live only when the file it names is in this set, its **current** digest is the one
    pinned, and the audit really raised that check for that file. Anything else is reported as
    void (bytes changed), stale (no such file) or unused — never quietly ignored, because a pin
    that stopped applying looks exactly like one that is still doing its job.

    Only a REPORT finding can be marked; a BLOCK finding has no confirmed state at all
    (:class:`Finding` refuses to build one), so the gate cannot be neutered by editing a JSON file.
    """
    known = {Path(path).name for path in paths}
    digests = {Path(record.file).name: record.sha256 for record in records}
    raised = {(Path(finding.file).name, finding.check) for finding in findings}
    statuses: list[PinStatus] = []
    for pin in pins:
        name = Path(pin.file).name
        if name not in known:
            statuses.append(PinStatus(pin, "stale", "no such file in this reference set"))
        elif digests.get(name) != pin.sha256:
            statuses.append(PinStatus(pin, "void", "the file's sha256 no longer matches the pin"))
        elif (name, pin.check) not in raised:
            statuses.append(
                PinStatus(pin, "unused", f"the audit raised no {pin.check} finding for this file")
            )
        else:
            statuses.append(PinStatus(pin, "applied", pin.reason))

    live = {
        (Path(status.pin.file).name, status.pin.check): status.pin
        for status in statuses
        if status.status == "applied"
    }
    marked: list[Finding] = []
    for finding in findings:
        pin = live.get((Path(finding.file).name, finding.check))
        if pin is not None and finding.verdict == REPORT:
            marked.append(
                Finding(
                    file=finding.file,
                    check=finding.check,
                    verdict=finding.verdict,
                    reason=finding.reason,
                    confirmed_by=pin.confirmed_by,
                    confirmed_on=pin.confirmed_on,
                )
            )
        else:
            marked.append(finding)
    return marked, statuses


# ---------------------------------------------------------------- measurements


def ssc_orientation(genome: Genome, ir: IRResult) -> str:
    """Is the SSC's coding strand the same as the LSC's, inside this one record?

    Measured as §8.34 measured it: after :func:`plastid.core.normalize.normalize` (which rotates
    the origin and never flips a strand), each region votes with the strand of its own CDS and the
    two votes are compared. Comparing raw coordinates between records would be meaningless — every
    record has its own origin. ``"unknown"`` when a region holds no CDS, or its vote is tied: an
    empty or split vote is not a verdict.
    """
    frame = normalize(genome, ir)
    regions = {region.name: region for region in regions_from_ir(frame.genome, frame.ir)}
    lsc = _majority_strand(frame.genome, regions["LSC"])
    ssc = _majority_strand(frame.genome, regions["SSC"])
    if lsc == 0 or ssc == 0:
        return "unknown"
    return "same" if lsc == ssc else "opposite"


def _majority_strand(genome: Genome, region) -> int:
    votes = Counter(
        gene.strand
        for gene in genome.genes
        if gene.kind == "CDS" and not gene.pseudo and region.start <= _midpoint(gene) < region.end
    )
    if not votes:
        return 0
    top, count = votes.most_common(1)[0]
    if len(votes) > 1 and count == votes.most_common(2)[1][1]:
        return 0
    return top


def _midpoint(gene: Gene) -> int:
    """0-based midpoint of a gene's outer span — the point that decides which region it is in."""
    start_1b, end_1b = gene.span
    return (start_1b - 1 + end_1b) // 2


def dotplot_runs(query: str, subject: str, *, word: int = DOTPLOT_WORD) -> tuple[Match, ...]:
    """Exact-match runs between two sequences — the dots of a dot plot.

    Every position where the two share a ``word``-long word starts a run, extended as far as they
    stay equal; only run starts are kept, so an uninterrupted match is one row rather than one per
    base. Draw the query against the subject **and** against its reverse complement: a region that
    runs the same way lands on one diagonal, an inverted one on the other. That picture is the
    confirmation the user reads when a gene-strand verdict needs eyes on it.
    """
    if word < 1:
        raise ValueError(f"the word length must be at least 1, got {word}")
    starts: dict[str, list[int]] = {}
    for index in range(len(subject) - word + 1):
        starts.setdefault(subject[index : index + word], []).append(index)
    runs: list[Match] = []
    for index in range(len(query) - word + 1):
        for target in starts.get(query[index : index + word], ()):
            if (
                target > 0
                and index > 0
                and query[index - 1 : index - 1 + word] == subject[target - 1 : target - 1 + word]
            ):
                continue  # not the start of the run
            length = word
            while (
                index + length < len(query)
                and target + length < len(subject)
                and query[index + length] == subject[target + length]
            ):
                length += 1
            runs.append(Match(query_pos=index + 1, subject_pos=target + 1, length=length))
    return tuple(runs)


def region_sequences(paths: Mapping[str, str | Path], region: str) -> dict[str, str]:
    """One region out of every named record, in the normalized frame (LSC → IRa → SSC → IRb)."""
    out: dict[str, str] = {}
    for name, path in paths.items():
        genome = read_genbank(path)
        out[name] = _region_sequence(genome, enforce_ir_gate(detect_ir(genome)), region)
    return out


def _region_sequence(genome: Genome, ir: IRResult, region_name: str) -> str:
    frame = normalize(genome, ir)
    regions = {region.name: region for region in regions_from_ir(frame.genome, frame.ir)}
    region = regions[region_name]
    return frame.genome.sequence[region.start : region.end]
