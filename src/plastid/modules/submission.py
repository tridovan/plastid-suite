"""NCBI submission files — FASTA, the five-column table, and optionally the ASN.1 ``.sqn``.

Decision D3-04: the ``.fsa`` and ``.tbl`` are emitted in pure Python and always — NCBI's web
submission accepts exactly that pair, so the minimum install can prepare a submission with no
external tool at all. The ``.sqn`` is produced only when asked for (``--set sqn=true``) **and**
``table2asn`` is present, because it needs a submitter template (``.sbt``) that only the team
has: a template names people and institutions, and this toolkit does not invent those.

The first-pass checks are the ones that fail at upload time if they are missed. They are
reported as findings rather than repaired: a feature past the end of the sequence is a claim
about a sequence that is not there, and only a human knows which side is wrong.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plastid.core.model import Genome, count_genes
from plastid.core.organelle import GROUP_SHARED
from plastid.core.submission import (
    SubmissionInputError,
    partial_start_genes,
    source_qualifiers_left_out,
    submission_problems,
    submittable_features,
    to_fasta,
    to_tbl,
)
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.layout import slugify
from plastid.run.tools import ToolError, which


@dataclass(frozen=True, slots=True)
class ValidatorCounts:
    """What NCBI's validator said about the table, read from its own ``.stats`` file.

    Kept because a report nobody reads is not a check: the ``.val`` is an artifact, but the
    count belongs in the table a person actually looks at.
    """

    errors: int
    warnings: int
    by_class: tuple[tuple[str, str, int], ...]  # (code, level, messages), errors first


def validator_counts(text: str) -> ValidatorCounts:
    """Read table2asn's ``.stats``: totals plus the per-code breakdown it groups them by."""
    errors = warnings = 0
    findings: list[tuple[str, str, int]] = []
    level = ""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.endswith("ERROR-level messages exist"):
            level, errors = "error", int(stripped.split()[0])
            continue
        if stripped.endswith("WARNING-level messages exist"):
            level, warnings = "warning", int(stripped.split()[0])
            continue
        code, separator, count = stripped.partition(":")
        if separator and level and count.strip().isdigit():
            findings.append((code.strip(), level, int(count)))
    order = {"error": 0, "warning": 1}
    return ValidatorCounts(
        errors=errors,
        warnings=warnings,
        by_class=tuple(sorted(findings, key=lambda finding: (order[finding[1]], -finding[2]))),
    )


MODULE_ID = "submission"
MODULE_VERSION = "0.1.0"

TOOL = "table2asn"
SCOPE_LABEL = "every CDS/tRNA/rRNA feature and its parent gene feature, both IR copies"

# The tool's own report files, kept beside the .sqn when it produces them: NCBI's instructions
# are to read the .val and .stats before submitting.
_TOOL_REPORTS = ("val", "stats", "dr")


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="NCBI submission files",
    requires_ir=False,  # writing features is not a position-based analysis
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    seq_id = genome.accession or slugify(genome.name)
    if not submittable_features(genome):
        raise SubmissionInputError(
            f"{genome.accession or genome.name}: the record carries no gene, CDS, tRNA or rRNA "
            f"feature — a submission is built from a checked GenBank annotation, not from a "
            f"bare sequence. Annotate the sample (`plastid annotate`), check the result "
            f"(`plastid run -m gene_normalize --set mode=verify`), and run submission on that "
            f"GenBank"
        )
    problems = submission_problems(genome)
    partial_starts = partial_start_genes(genome)

    artifacts = {
        "submission.fsa": to_fasta(genome, seq_id),
        "submission.tbl": to_tbl(genome, seq_id),
    }
    sqn_note = "not requested"
    validator_status = "not requested"
    validator: ValidatorCounts | None = None
    if _flag(params.get("sqn")):
        artifacts, sqn_note, validator, validator_status = _build_sqn(
            genome, seq_id, params, ctx, artifacts
        )

    counts = count_genes(genome)
    summary = Table(
        id="submission_summary",
        title="Submission files",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("sequence_id", "Sequence id in the files"),
            Column("n_features", "Features written", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_partial_starts", "CDS written with a partial start", dtype="int"),
            Column("n_problems", "First-pass problems", dtype="int"),
            Column("validator_status", "NCBI validator status"),
            Column("n_validator_errors", "NCBI errors", dtype="int"),
            Column("n_validator_warnings", "NCBI warnings", dtype="int"),
            Column("sqn", "ASN.1 (.sqn)"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "sequence_id": seq_id,
                "n_features": sum(
                    1 for gene in genome.genes if gene.kind in ("gene", "CDS", "tRNA", "rRNA")
                ),
                "n_cds": counts.cds,
                "n_trna": counts.trna,
                "n_rrna": counts.rrna,
                "n_partial_starts": len(partial_starts),
                "n_problems": len(problems),
                # A validator that did not run has no counts. Zero would be a measurement this
                # run never made; the scope rule is "not measured", carried as a blank cell
                # with the reason in validator_status (M11.18).
                "validator_status": validator_status,
                "n_validator_errors": validator.errors if validator else None,
                "n_validator_warnings": validator.warnings if validator else None,
                "sqn": sqn_note,
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    tables: dict[str, Table] = {summary.id: summary}
    if problems:
        tables["submission_problems"] = Table(
            id="submission_problems",
            title="First-pass problems — to fix, not to submit",
            columns=(
                Column("problem_id", "Problem"),
                Column("finding", "Finding"),
            ),
            rows=tuple(
                {"problem_id": f"P{index}", "finding": finding}
                for index, finding in enumerate(problems, start=1)
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="problem_id",
        )

    notes = [
        "The .tbl carries the name each gene arrived with: it describes this record, and "
        "renaming a gene would change a submitted record.",
        "Partial coordinates travel from the record into the table; a partial marker on an "
        "internal interval is refused rather than dropped.",
    ]
    left_out = source_qualifiers_left_out(genome)
    if left_out:
        notes.append(
            f"The record's source carries {', '.join(left_out)}, which the five-column table has "
            f"no place for: NCBI's own converter calls them unrecognized qualifiers and exits "
            f"non-zero (while still writing the .sqn), and it derives both itself. They stay on "
            f"the genome read from the record — only the table leaves them out."
        )
    if partial_starts:
        names = ", ".join(sorted({gene.name for gene in partial_starts}))
        notes.append(
            f"{len(partial_starts)} CDS ({names}) do not begin with a codon their translation "
            f"table initiates, so their start is written partial (a `<` on the 5' coordinate). "
            f"The coordinates in the table are the record's own — only the marker is added, and "
            f"it is what makes the feature acceptable to NCBI, whose validator rejects the "
            f"alternative. The usual cause in a plastome is C-to-U RNA editing creating the start "
            f"codon (rpl2, ndhD); anything else is worth a look before submitting."
        )
    if problems:
        notes.append(
            f"{len(problems)} first-pass problem(s) — the files are written so they can be "
            f"inspected, but they are not ready to submit."
        )
    if _flag(params.get("sqn")):
        notes.append(
            "A .sqn is a conversion, not a validation: it does not mean NCBI has accepted the "
            "record. validator_status says whether table2asn's own validator actually ran; when "
            "it did not, the counts are blank rather than zero, because zero would be a "
            "measurement this run never made."
        )
    if validator and validator.by_class:
        tables["validator_findings"] = Table(
            id="validator_findings",
            title=f"NCBI's validator on the .tbl — {validator.errors} error(s), "
            f"{validator.warnings} warning(s)",
            columns=(
                Column("code", "Validator code"),
                Column("level", "Level"),
                Column("n", "Messages", dtype="int"),
            ),
            rows=tuple(
                {"code": code, "level": level, "n": count}
                for code, level, count in validator.by_class
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="code",
        )
        notes.append(
            "The full report is kept beside the .sqn as submission.val — read it before "
            "submitting, which is what NCBI asks of every submission. A finding here is about "
            "the record, not about the writer: a feature whose gene and CDS disagree about the "
            "strand is one this format produces for every duplicated inverted-repeat gene "
            "(measured against a published plastome submission)."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        artifacts=artifacts,
    )


def _flag(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _build_sqn(
    genome: Genome,
    seq_id: str,
    params: dict[str, Any],
    ctx: RunContext,
    artifacts: dict[str, str],
) -> tuple[dict[str, str], str, ValidatorCounts | None, str]:
    """Run table2asn and keep its outputs. A missing tool or template is a named error."""
    explicit = params.get("table2asn")
    binary = None
    if explicit:
        candidate = Path(str(explicit))
        binary = str(candidate) if candidate.is_file() else which(TOOL, extra=(candidate,))
    else:
        binary = which(TOOL)
    if binary is None:
        raise ToolError(
            f"{TOOL} not found — the .sqn step needs it. It is a single binary from NCBI's "
            f"asn1-converters (see docs/tools/table2asn.md for the exact download line); pass "
            f"--set table2asn=/path/to/table2asn if it is not on PATH. The .fsa and .tbl in "
            f"this run are complete and are what NCBI's web submission accepts, so dropping "
            f"--set sqn=true loses nothing but the conversion"
        )
    template = params.get("template")
    if not template:
        raise ToolError(
            "--set sqn=true needs a submitter template: pass --set template=/path/to/your.sbt. "
            "NCBI generates it from a web form (submit.ncbi.nlm.nih.gov/genbank/template/"
            "submission) because it names people and institutions — this toolkit does not "
            "invent those, and a hand-written .sbt is rejected as an unreadable template"
        )
    if not Path(str(template)).is_file():
        raise ToolError(f"the submitter template {template!r} does not exist")

    # table2asn pairs a table with a sequence by filename prefix, so both files are written
    # under the sequence id in a fresh directory — never over the run's own outputs.
    counts: ValidatorCounts | None = None
    with tempfile.TemporaryDirectory(prefix="table2asn_") as work:
        workdir = Path(work)
        (workdir / f"{seq_id}.fsa").write_text(artifacts["submission.fsa"], encoding="utf-8")
        (workdir / f"{seq_id}.tbl").write_text(artifacts["submission.tbl"], encoding="utf-8")
        ctx.tool(
            TOOL,
            [
                binary,
                "-t",
                str(template),
                "-i",
                f"{seq_id}.fsa",
                "-o",
                f"{seq_id}.sqn",
                "-M",
                "n",
            ],
            cwd=workdir,
        )
        produced = workdir / f"{seq_id}.sqn"
        if not produced.is_file():
            raise ToolError(
                f"{TOOL} exited 0 but produced no {produced.name} — see the captured logs"
            )
        artifacts = {**artifacts, "submission.sqn": produced.read_text(encoding="utf-8")}
        kept: list[str] = []
        for extension in _TOOL_REPORTS:
            report = workdir / f"{seq_id}.{extension}"
            if report.is_file():
                text = report.read_text(encoding="utf-8", errors="replace")
                artifacts[f"submission.{extension}"] = text
                kept.append(extension)
                if extension == "stats":
                    counts = validator_counts(text)
    produced_note = f"produced by {Path(binary).name}"
    if kept:
        produced_note += f"; reports kept: {', '.join(kept)}"
    validator_status = (
        "ran — table2asn wrote .stats"
        if counts is not None
        else "not run — the tool wrote no .stats; a .sqn is not a validation"
    )
    return artifacts, produced_note, counts, validator_status
