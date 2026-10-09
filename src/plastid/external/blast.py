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
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plastid.core.mito_caller import Hsp
from plastid.core.model import Genome
from plastid.core.submission import to_fasta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "BLAST_VERSION",
    "DEFAULT_EVALUE",
    "MAX_HSPS",
    "NAMING_EVALUE",
    "SEGMENT_OFF",
    "BlastResult",
    "ProteinHits",
    "ProteinHit",
    "TraceHit",
    "TraceResult",
    "blast_version",
    "blastp_argv",
    "makeblastdb_argv",
    "parse_hits",
    "parse_protein_hits",
    "parse_trace_hits",
    "resolve_binary",
    "run_blastn",
    "run_blastn_trace",
    "run_blastp",
    "run_blastx",
    "search_argv",
    "subject_count",
    "subject_key",
]

# The build the toolkit's numbers were measured with. A different BLAST is a different
# annotation, so a mismatch is refused by name rather than recorded as a footnote.
BLAST_VERSION = "2.17.0"

SEARCH_TOOLS: tuple[str, ...] = ("blastx", "blastn")
DB_TOOL = "makeblastdb"

DEFAULT_EVALUE = 1e-5
MAX_HSPS = 10  # per query-subject pair; a pair that returns exactly this is reported as capped
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

# The columns of a search that has to answer *where inside an alignment* a position sits, as
# `run_blastn_trace` does for P3B-10's exon structure. Only its own parser reads them.
_TRACE_COLUMNS: tuple[str, ...] = _COLUMNS + ("qseq", "sseq")


@dataclass(frozen=True, slots=True)
class BlastResult:
    """What one search produced: the hits, the raw text they were parsed from, and what ran."""

    hits: tuple[Hsp, ...]
    raw: str
    search: str
    version: str
    reference: str
    molecule: str  # "prot" for blastx, "nucl" for blastn


@dataclass(frozen=True, slots=True)
class TraceHit:
    """One hit **with its alignment** — the only way to place a position inside an alignment.

    Coordinate arithmetic maps an alignment's two ends and nothing in between: the reference's own
    introns carry indels against the target, so the k-th base of the query is not the k-th base of
    the subject, and an exon boundary sitting between them is off by every indel before it.
    Measured on the golden (10/10/2026): the petB junction's 741 bp intron aligns with one gap,
    and the naive mapping put the first exon 1 bp away from where the published record has it.
    """

    subject_key: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qseq: str
    sseq: str

    def map_subject(self, position: int) -> int | None:
        """The 1-based **query** position aligned to a 1-based subject position, or ``None``.

        ``None`` has two meanings and both are "this base is not on the target": the position is
        outside the alignment, or the alignment wrote a gap across from it (an insertion in the
        subject). Either way a caller must not invent a coordinate for it.
        """
        if not min(self.sstart, self.send) <= position <= max(self.sstart, self.send):
            return None
        query, subject = self.qstart, self.sstart
        step_query = 1 if self.qend >= self.qstart else -1
        step_subject = 1 if self.send >= self.sstart else -1
        for query_base, subject_base in zip(self.qseq, self.sseq, strict=True):
            if subject == position:
                return None if query_base == "-" else query
            if query_base != "-":
                query += step_query
            if subject_base != "-":
                subject += step_subject
        return None


@dataclass(frozen=True, slots=True)
class TraceResult:
    """What one traced search produced: the hits with their alignments, and what ran."""

    hits: tuple[TraceHit, ...]
    raw: str
    search: str
    version: str
    reference: str


def subject_count(reference: str | Path) -> int:
    """How many sequences a reference FASTA holds — the only honest ``-max_target_seqs``.

    Counted from the file the search is about to build its database from, so the cap cannot
    disagree with the set: a header is a line whose first character is ``>``, which is what makes
    a FASTA a FASTA.
    """
    with Path(reference).open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.startswith(">"))


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


def parse_trace_hits(text: str, *, seq_id: str) -> tuple[TraceHit, ...]:
    """Read a traced search's tabular output, alignment strings included.

    Strict in the same way :func:`parse_hits` is, and for the same reason: a malformed line read
    leniently is a boundary placed from the wrong row.
    """
    expected = len(_TRACE_COLUMNS)
    hits: list[TraceHit] = []
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != expected:
            raise ValueError(
                f"line {number}: expected {expected} tab-separated fields, found {len(fields)}: "
                f"{line[:120]!r} — a traced search needs the alignment strings"
            )
        if fields[0] != seq_id:
            raise ValueError(
                f"line {number}: this is a hit for query {fields[0]!r}, not {seq_id!r} — the file "
                f"is not this run's output"
            )
        query, subject = fields[-2], fields[-1]
        if len(query) != len(subject):
            raise ValueError(
                f"line {number}: the alignment's two strings are {len(query)} and {len(subject)} "
                f"long — an alignment whose halves do not line up cannot answer where a base is"
            )
        try:
            hits.append(
                TraceHit(
                    subject_key=subject_key(fields[1]),
                    pident=float(fields[2]),
                    length=int(fields[3]),
                    qstart=int(fields[6]),
                    qend=int(fields[7]),
                    sstart=int(fields[8]),
                    send=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                    qseq=query,
                    sseq=subject,
                )
            )
        except ValueError as error:
            raise ValueError(f"line {number}: {error}") from error
    return tuple(hits)


def _outfmt(frame: bool, trace: bool = False) -> str:
    columns = _COLUMNS + ((_FRAME_COLUMN,) if frame else ())
    if trace:
        columns += ("qseq", "sseq")
    return "6 " + " ".join(columns)


def search_argv(
    tool: str,
    *,
    binary: str,
    query: str | Path,
    database: str | Path,
    evalue: float,
    max_hsps: int,
    max_target_seqs: int,
    trace: bool = False,
) -> list[str]:
    """The exact command line, so a test can hold it still without running BLAST.

    ``max_target_seqs`` is the size of the reference set, counted from the file — never a
    constant. Measured 08/10/2026 (P3B-07): the 500 that sufficed for the mito path's 287
    sequences silently truncated the plastid search, which searches 2,838 — BLAST reported 929
    HSPs over 467 subjects where the same command with the cap lifted reported 4,382 over 2,603,
    and **nothing in the output said so**. A cap that cannot truncate is the only kind worth
    having.
    """
    if tool not in SEARCH_TOOLS:
        raise ValueError(f"{tool!r} is not one of {SEARCH_TOOLS}")
    argv = [
        binary,
        "-query",
        str(query),
        "-db",
        str(database),
        "-outfmt",
        _outfmt(tool == "blastx", trace),
        "-evalue",
        f"{evalue:g}",
        "-max_hsps",
        str(max_hsps),
        "-max_target_seqs",
        str(max_target_seqs),
    ]
    if tool == "blastx":
        argv += [SEGMENT_OFF, "no"]
    else:
        argv += ["-task", "blastn"]
    return argv


# --- the naming search (P3B-04) ---------------------------------------------------------------
#
# `blastp`, protein against protein, is the one search D5-08 keeps for the naming step: both
# sides already carry a curated protein (/translation), so blastx's six-frame translation would
# add cost and no sensitivity, and it would let a frame error hide inside a good score.
#
# The E-value and the flags are the ones §8.41's thresholds were measured with, and that is the
# only reason they are not the module's defaults: identity ≥ 50 %, coverage ≥ 0.7 and a margin
# ≥ 10 bits are numbers *of this command line*, and a different one would make them folklore.

NAMING_EVALUE = 10.0
_PROTEIN_COLUMNS: tuple[str, ...] = ("qseqid", "sseqid", "pident", "length", "bitscore", "qlen")


@dataclass(frozen=True, slots=True)
class ProteinHit:
    """One blastp hit for one query gene: what it matched, and how well."""

    query_id: str
    subject: str  # the subject key's gene field — the naming reference's own spelling
    accession: str  # the subject key's accession field
    identity: float  # % identity over the aligned length
    length: int  # aligned length, in amino acids
    bits: float
    query_length: int  # the query protein's length

    @property
    def coverage(self) -> float:
        return self.length / self.query_length if self.query_length else 0.0


@dataclass(frozen=True, slots=True)
class ProteinHits:
    """One search's hits, grouped by query and ranked best first — plus what produced them."""

    by_query: dict[str, tuple[ProteinHit, ...]]
    raw: str
    search: str
    version: str

    def for_query(self, query_id: str) -> tuple[ProteinHit, ...]:
        return self.by_query.get(query_id, ())


def parse_protein_hits(
    text: str, *, subject_vocabulary: str = "A"
) -> dict[str, tuple[ProteinHit, ...]]:
    """Read ``blastp -outfmt 6 qseqid sseqid pident length bitscore qlen``, ranked per query.

    The subject key is the same three-field contract the mito reference writes
    (:func:`subject_key`) — ``accession|gene|vocabulary`` — so a naming run's own FASTA says which
    vocabulary every subject belongs to (``A`` here) and a hit cannot be mistaken for one from
    another table. A line with the wrong field count or a foreign vocabulary raises: this is a
    parser, and a silent skip is how a missing gene becomes an absence.
    """
    grouped: dict[str, list[ProteinHit]] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(_PROTEIN_COLUMNS):
            raise ValueError(
                f"line {number}: expected {len(_PROTEIN_COLUMNS)} tab-separated fields, found "
                f"{len(fields)}: {line[:120]!r}"
            )
        query_id, subject_id, pident, length, bits, qlen = fields
        key = subject_key(subject_id)
        accession, gene, vocabulary = key.split("|", 2)
        if vocabulary != subject_vocabulary:
            raise ValueError(
                f"line {number}: subject {key!r} is vocabulary {vocabulary!r}, not "
                f"{subject_vocabulary!r} — the reference FASTA is not the naming reference"
            )
        grouped.setdefault(query_id, []).append(
            ProteinHit(
                query_id=query_id,
                subject=gene,
                accession=accession,
                identity=float(pident),
                length=int(length),
                bits=float(bits),
                query_length=int(qlen),
            )
        )
    return {
        query_id: tuple(sorted(hits, key=lambda hit: (-hit.bits, hit.subject)))
        for query_id, hits in grouped.items()
    }


def blastp_argv(
    *,
    binary: str,
    query: str | Path,
    database: str | Path,
    evalue: float,
    max_hsps: int,
    max_target_seqs: int,
) -> list[str]:
    """The naming search's command line — held still by a test, so the thresholds stay attached."""
    return [
        binary,
        "-query",
        str(query),
        "-db",
        str(database),
        "-outfmt",
        "6 " + " ".join(_PROTEIN_COLUMNS),
        "-evalue",
        f"{evalue:g}",
        "-max_hsps",
        str(max_hsps),
        "-max_target_seqs",
        str(max_target_seqs),
        SEGMENT_OFF,
        "no",
    ]


def run_blastp(
    queries: Sequence[tuple[str, str]],
    reference: str | Path,
    ctx: RunContext,
    *,
    search: str | Path | None = None,
    makeblastdb: str | Path | None = None,
    directory: str | Path | None = None,
    evalue: float = NAMING_EVALUE,
    max_hsps: int = MAX_HSPS,
) -> ProteinHits:
    """Search ``(query_id, protein)`` rows against a protein reference FASTA on disk.

    Both FASTAs are built by the caller — the reference is the naming record's own curated
    proteins, written with :func:`plastid.core.gene_names.protein_fasta_rows`, and the queries are
    the record's genes. The reference is a path rather than data because it is shared by every
    query and rebuilding it per query would be a second, silently different database.
    """
    if not queries:
        raise ValueError("a naming search needs at least one query protein")
    binary = resolve_binary("blastp", explicit=search, directory=directory)
    db_binary = resolve_binary(DB_TOOL, explicit=makeblastdb, directory=directory)
    version = blast_version(binary)
    source = Path(reference)
    if not source.is_file():
        raise ToolError(f"the naming reference is not a file: {source}")

    with tempfile.TemporaryDirectory(prefix="blastp_") as work:
        workdir = Path(work)
        query = workdir / "queries.faa"
        query.write_text(
            "".join(f">{name}\n{protein}\n" for name, protein in queries), encoding="utf-8"
        )
        database = workdir / "reference"
        ctx.tool(
            DB_TOOL,
            makeblastdb_argv(
                binary=db_binary, reference=source, database=database, molecule="prot"
            ),
            kind="binary",
            tool_path=db_binary,
            cwd=workdir,
            timeout=600,
        )
        record = ctx.tool(
            "blastp",
            blastp_argv(
                binary=binary,
                query=query,
                database=database,
                evalue=evalue,
                max_hsps=max_hsps,
                max_target_seqs=subject_count(source),
            ),
            kind="binary",
            tool_path=binary,
            cwd=workdir,
            timeout=3600,
        )
        raw = Path(record.stdout_file).read_text(encoding="utf-8")
    return ProteinHits(
        by_query=parse_protein_hits(raw),
        raw=raw,
        search=binary,
        version=version,
    )


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
            f"one with --set tools=<directory>. See docs/tools/blast.md."
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
    return _search(
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
        trace=False,
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
    return _search(
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
        trace=False,
    )


def run_blastn_trace(
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
) -> TraceResult:
    """Search a genome against a nucleotide reference and keep **where** every hit aligns.

    The same search :func:`run_blastn` runs, with the alignment strings added to its output: the
    tool's own record is the evidence, and re-running it to ask a second question would be a
    second measurement of the same thing (and could disagree with the recorded argv).
    """
    return _search(
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
        trace=True,
    )


def _search(
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
    trace: bool,
) -> Any:
    """One search, traced or not — the caller decides which result type it gets."""
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
                max_target_seqs=subject_count(source),
                trace=trace,
            ),
            kind="binary",
            tool_path=binary,
            cwd=workdir,
            timeout=3600,
        )
        # the search's own captured stdout is the result; re-running it to get the text would be
        # a second measurement of the same thing (and could disagree with the recorded argv)
        raw = Path(record.stdout_file).read_text(encoding="utf-8")

    if trace:
        return TraceResult(
            hits=parse_trace_hits(raw, seq_id=seq_id),
            raw=raw,
            search=binary,
            version=version,
            reference=str(source),
        )

    return BlastResult(
        hits=parse_hits(raw, seq_id=seq_id, frame=tool == "blastx"),
        raw=raw,
        search=binary,
        version=version,
        reference=str(source),
        molecule=molecule,
    )
