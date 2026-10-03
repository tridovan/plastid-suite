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
from pathlib import Path
from typing import Any

from plastid.core.model import Genome, count_genes
from plastid.core.submission import submission_problems, to_fasta, to_tbl
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.layout import slugify
from plastid.run.tools import ToolError, which

MODULE_ID = "submission"
MODULE_VERSION = "0.1.0"

TOOL = "table2asn"
SCOPE_LABEL = "every CDS/tRNA/rRNA feature and its parent gene feature, both IR copies"

# The tool's own report files, kept beside the .sqn when it produces them: NCBI's instructions
# are to read the .val and .stats before submitting.
_TOOL_REPORTS = ("val", "stats", "dr")


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="NCBI submission files",
    requires_ir=False,  # writing features is not a position-based analysis
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    seq_id = genome.accession or slugify(genome.name)
    problems = submission_problems(genome)

    artifacts = {
        "submission.fsa": to_fasta(genome, seq_id),
        "submission.tbl": to_tbl(genome, seq_id),
    }
    sqn_note = "not requested"
    if _flag(params.get("sqn")):
        artifacts, sqn_note = _build_sqn(genome, seq_id, params, ctx, artifacts)

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
            Column("n_problems", "First-pass problems", dtype="int"),
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
                "n_problems": len(problems),
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
        "renaming a gene would change a submitted record (TASKS.md pitfall).",
        "Partial coordinates travel from the record into the table; a partial marker on an "
        "internal interval is refused rather than dropped.",
    ]
    if problems:
        notes.append(
            f"{len(problems)} first-pass problem(s) — the files are written so they can be "
            f"inspected, but they are not ready to submit."
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
) -> tuple[dict[str, str], str]:
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
                artifacts[f"submission.{extension}"] = report.read_text(
                    encoding="utf-8", errors="replace"
                )
                kept.append(extension)
    produced_note = f"produced by {TOOL} {Path(binary).name}"
    if kept:
        produced_note += f"; reports kept: {', '.join(kept)}"
    return artifacts, produced_note
