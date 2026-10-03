"""BLAST+ adapter — the searches the mito caller transfers an annotation from (P3A-03).

Every call goes through :meth:`RunContext.tool`, so the run records name · version · path · argv ·
time · exit code · logs like any other tool (seam 3), and the parsed hits come from the search's
own captured stdout rather than a second run.

Three decisions in here were measured, and each of them is silent when it is wrong:

* **``-seg no``.** SEG masks low-complexity stretches of the query — sensible when hunting an
  unknown protein, wrong when transferring an annotation to a plant mitogenome, which is wall to
  wall hydrophobic membrane proteins. Measured on the answer-key record: with the default mask,
  atp9's best alignment drops from 137 to 67 bits and 21 other loci lose score as well. The
  identity floor, not the mask, is what keeps junk out.
* **no ``-parse_seqids``.** The obvious flag for a FASTA whose ids are ``accession|gene|…`` is the
  one that breaks: measured 03/10/2026, ``makeblastdb -parse_seqids`` fails with ``Error:
  Duplicate seq_ids are found: LCL|NC_006581.1|SDH3|-|NICOTIANA`` (the reference set legitimately
  holds repeat copies), and for the ids it does parse it upper-cases them, which would break every
  lowercase gene symbol. Without it, ``sseqid`` is the defline's first token — which is what the
  subject key wants.
* **``-task blastn`` for rRNA.** Megablast and blastn found the same three rRNAs on the answer-key
  record, but the target of a real run may be a divergent sample and the database is 25 sequences:
  sensitivity costs almost nothing here.

BLAST truncates ``sseqid`` at the first whitespace, so a hit's subject is the **first three** header
fields — the reference headers end with an organism name, which never survives. :func:`parse_hits`
refuses anything it cannot read rather than guessing: a malformed line, a hit for another query, or
a subject id with too few fields all raise with the line number.
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.mito_caller import Hsp
from plastid.core.model import Genome
from plastid.core.submission import to_fasta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "BLAST_VERSION",
    "DEFAULT_EVALUE",
    "MAX_HSPS",
    "MAX_TARGET_SEQS",
    "SEGMENT_OFF",
    "BlastResult",
    "blast_version",
    "makeblastdb_argv",
    "parse_hits",
    "resolve_binary",
    "run_blastn",
    "run_blastx",
    "search_argv",
    "subject_key",
]

# The build the toolkit's numbers were measured with. A different BLAST is a different
# annotation, so a mismatch is refused by name rather than recorded as a footnote.
BLAST_VERSION = "2.17.0"

SEARCH_TOOLS: tuple[str, ...] = ("blastx", "blastn")
DB_TOOL = "makeblastdb"

DEFAULT_EVALUE = 1e-5
MAX_HSPS = 10  # per query-subject pair; a pair that returns exactly this is reported as capped
MAX_TARGET_SEQS = 500  # the reference set is 287 sequences — nothing may be truncated away
SEGMENT_OFF = "-seg"  # see the module docstring; the value is "no"

# The tabular columns, in the order the parser reads them. This tuple and the `-outfmt` string
# are one contract, which is why they live in the same file.
_COLUMNS: tuple[str, ...] = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
_FRAME_COLUMN = "qframe"
_VERSION = re.compile(r"\b(\d+)\.(\d+)\.(\d+)")
_SUBJECT_FIELDS = 3


@dataclass(frozen=True, slots=True)
class BlastResult:
    """What one search produced: the hits, the raw text they were parsed from, and what ran."""

    hits: tuple[Hsp, ...]
    raw: str
    search: str
    version: str
    reference: str
    molecule: str  # "prot" for blastx, "nucl" for blastn


def subject_key(sseqid: str) -> str:
    """``accession|gene|vocabulary`` from a BLAST subject id.

    Two measured reasons this is not simply ``sseqid``: BLAST cuts the id at the first whitespace
    (the fourth header field is an organism name with spaces in it), and the dataset itself may
    hold several sequences under one key. Raises rather than returning something short.
    """
    fields = sseqid.split("|")
    if len(fields) < _SUBJECT_FIELDS:
        raise ValueError(
            f"subject id {sseqid!r} carries fewer than three fields — a subject key is "
            f"'accession|gene|vocabulary'. BLAST truncates the id at the first whitespace, so a "
            f"key with more fields would never match a hit; one with fewer means the FASTA was "
            f"not the one this toolkit wrote (data/mito_reference.md has the contract)"
        )
    return "|".join(fields[:_SUBJECT_FIELDS])


def parse_hits(text: str, *, seq_id: str, frame: bool) -> tuple[Hsp, ...]:
    """Read BLAST's tabular output into :class:`Hsp` records, one line at a time, strictly.

    An empty result is an empty tuple: finding nothing is a fact the caller must decide about,
    not a parse error.
    """
    expected = len(_COLUMNS) + (1 if frame else 0)
    hits: list[Hsp] = []
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != expected:
            raise ValueError(
                f"line {number}: expected {expected} tab-separated fields, found {len(fields)}: "
                f"{line[:120]!r}"
            )
        if fields[0] != seq_id:
            raise ValueError(
                f"line {number}: this is a hit for query {fields[0]!r}, not {seq_id!r} — the file "
                f"is not this run's output"
            )
        try:
            hits.append(
                Hsp(
                    subject_key=subject_key(fields[1]),
                    pident=float(fields[2]),
                    length=int(fields[3]),
                    qstart=int(fields[6]),
                    qend=int(fields[7]),
                    sstart=int(fields[8]),
                    send=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                    frame=int(fields[12]) if frame else None,
                )
            )
        except ValueError as error:
            raise ValueError(f"line {number}: {error}") from error
    return tuple(hits)


def _outfmt(frame: bool) -> str:
    columns = _COLUMNS + ((_FRAME_COLUMN,) if frame else ())
    return "6 " + " ".join(columns)


def search_argv(
    tool: str,
    *,
    binary: str,
    query: str | Path,
    database: str | Path,
    evalue: float,
    max_hsps: int,
) -> list[str]:
    """The exact command line, so a test can hold it still without running BLAST."""
    if tool not in SEARCH_TOOLS:
        raise ValueError(f"{tool!r} is not one of {SEARCH_TOOLS}")
    argv = [
        binary,
        "-query",
        str(query),
        "-db",
        str(database),
        "-outfmt",
        _outfmt(tool == "blastx"),
        "-evalue",
        f"{evalue:g}",
        "-max_hsps",
        str(max_hsps),
        "-max_target_seqs",
        str(MAX_TARGET_SEQS),
    ]
    if tool == "blastx":
        argv += [SEGMENT_OFF, "no"]
    else:
        argv += ["-task", "blastn"]
    return argv


def makeblastdb_argv(
    *, binary: str, reference: str | Path, database: str | Path, molecule: str = "prot"
) -> list[str]:
    """The database-building command line — deliberately without ``-parse_seqids`` (docstring)."""
    return [binary, "-in", str(reference), "-dbtype", molecule, "-out", str(database)]


def resolve_binary(
    tool: str, *, explicit: str | Path | None = None, directory: str | Path | None = None
) -> str:
    """Find a BLAST binary: an explicit path, then a directory, then ``PATH``."""
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise ToolError(f"{tool} not found at {path} — see docs/tools/blast.md")
        return str(path)
    if directory is not None:
        candidate = Path(directory) / tool
        if not candidate.is_file():
            raise ToolError(f"{tool} not found in {directory} — see docs/tools/blast.md")
        return str(candidate)
    found = which(tool)
    if found is None:
        raise ToolError(
            f"{tool} not found — install the pinned BLAST+ {BLAST_VERSION}, or point the run at "
            f"one with --set blast=<directory>. See docs/tools/blast.md."
        )
    return found


def blast_version(binary: str) -> str:
    """The BLAST release this binary is, checked against the pin.

    Refused rather than noted: an annotation's coordinates come out of this program, and two BLAST
    releases are two answers nobody has compared. The message names both numbers so the fix is
    either "install the pinned build" or "re-measure and update the pin".
    """
    line, _ = tool_version(Path(binary).name, "binary", binary)
    match = _VERSION.search(line or "")
    if match is None:
        raise ToolError(
            f"{binary}: cannot read a version from {line!r} — a step whose version nobody can "
            f"name is not part of the toolkit"
        )
    version = match.group(0)
    if version != BLAST_VERSION:
        raise ToolError(
            f"{binary} is BLAST {version}, but the toolkit's numbers were measured with "
            f"{BLAST_VERSION} (docs/tools/blast.md) — a different BLAST is a different "
            f"annotation. Install the pinned build, or re-measure and update BLAST_VERSION."
        )
    return version


def run_blastx(
    genome: Genome,
    reference: str | Path,
    ctx: RunContext,
    *,
    seq_id: str,
    blastx: str | Path | None = None,
    makeblastdb: str | Path | None = None,
    directory: str | Path | None = None,
    evalue: float = DEFAULT_EVALUE,
    max_hsps: int = MAX_HSPS,
) -> BlastResult:
    """Search a genome's six reading frames against a protein reference set."""
    return _run(
        "blastx",
        genome,
        reference,
        ctx,
        seq_id=seq_id,
        search=blastx,
        makeblastdb=makeblastdb,
        directory=directory,
        evalue=evalue,
        max_hsps=max_hsps,
    )


def run_blastn(
    genome: Genome,
    reference: str | Path,
    ctx: RunContext,
    *,
    seq_id: str,
    blastn: str | Path | None = None,
    makeblastdb: str | Path | None = None,
    directory: str | Path | None = None,
    evalue: float = DEFAULT_EVALUE,
    max_hsps: int = MAX_HSPS,
) -> BlastResult:
    """Search a genome's two strands against a nucleotide reference set (the rRNA)."""
    return _run(
        "blastn",
        genome,
        reference,
        ctx,
        seq_id=seq_id,
        search=blastn,
        makeblastdb=makeblastdb,
        directory=directory,
        evalue=evalue,
        max_hsps=max_hsps,
    )


def _run(
    tool: str,
    genome: Genome,
    reference: str | Path,
    ctx: RunContext,
    *,
    seq_id: str,
    search: str | Path | None,
    makeblastdb: str | Path | None,
    directory: str | Path | None,
    evalue: float,
    max_hsps: int,
) -> BlastResult:
    binary = resolve_binary(tool, explicit=search, directory=directory)
    db_binary = resolve_binary(DB_TOOL, explicit=makeblastdb, directory=directory)
    version = blast_version(binary)
    molecule = "prot" if tool == "blastx" else "nucl"
    source = Path(reference)
    if not source.is_file():
        raise ToolError(f"the reference set is not a file: {source}")

    with tempfile.TemporaryDirectory(prefix=f"{tool}_") as work:
        workdir = Path(work)
        query = workdir / "target.fna"
        query.write_text(to_fasta(genome, seq_id), encoding="utf-8")
        database = workdir / "reference"
        ctx.tool(
            DB_TOOL,
            makeblastdb_argv(
                binary=db_binary, reference=source, database=database, molecule=molecule
            ),
            kind="binary",
            tool_path=db_binary,
            cwd=workdir,
            timeout=600,
        )
        record = ctx.tool(
            tool,
            search_argv(
                tool,
                binary=binary,
                query=query,
                database=database,
                evalue=evalue,
                max_hsps=max_hsps,
            ),
            kind="binary",
            tool_path=binary,
            cwd=workdir,
            timeout=3600,
        )
        # the search's own captured stdout is the result; re-running it to get the text would be
        # a second measurement of the same thing (and could disagree with the recorded argv)
        raw = Path(record.stdout_file).read_text(encoding="utf-8")

    return BlastResult(
        hits=parse_hits(raw, seq_id=seq_id, frame=tool == "blastx"),
        raw=raw,
        search=binary,
        version=version,
        reference=str(source),
        molecule=molecule,
    )
