"""The assembly step: reads in; a sequence, its graph, its log and the evidence for them out.

This file sits in ``modules/`` because that is where this project turns work into tables, but like
the reference audit it is **not registered** and does not follow ``run(genome, params, ctx)``: its
input is a pair of read files and its output is the genome every registered module starts from.
``plastid assemble`` drives it and owns the run directory.

**One step, two backends, one output contract.** The local backend runs fastp and GetOrganelle
through ``run_tool()``, or Oatk on HiFi reads. The import backend reads an assembly made elsewhere,
from reads this toolkit prepared and hashed, and accepts it only through
``validate_external_result``. From there all of them are the same object and go through the same
evaluation, so nothing downstream can tell them apart except the provenance that says where and with
what the assembly ran. The toolkit never moves reads to another machine by itself.

**What "done" means.** The assembler's word that a path closed is a claim. Reads mapped back with
an aligner the assembler does not use test it position by position, and a join the reads do not
cross is reported as undetermined rather than passed. Reads cannot say which arrangement of the
pieces is the molecule, so the selected graph is read as well, and a graph through which several
structures close is reported as unresolved. A sequence becomes an artifact for the next step only
when it closed, has every base known, the reads support it, and its graph is resolved.
"""

from __future__ import annotations

import gzip
import itertools
import json
import os
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from plastid.core.assembly import (
    STATUS_COMPLETE,
    STATUS_FRAGMENTED,
    STATUS_GAPPED,
    STATUS_NO_RESULT,
    SUPPORT_NOT_MEASURED,
    SUPPORT_SUPPORTED,
    SUPPORT_UNDETERMINED,
    AssembledPath,
    AssemblyError,
    BoundaryCall,
    ReferenceComparison,
    SupportProfile,
    SupportRule,
    WeakRegion,
    build_profile,
    call_boundary,
    compare_to_reference,
    ir_junctions,
    parse_sam,
    pool_repeated_joins,
    weak_regions,
)
from plastid.core.assembly_graph import (
    GRAPH_NOT_MEASURED,
    GRAPH_RESOLVED,
    AssemblerOutput,
    GraphReading,
    graph_verdict,
)
from plastid.core.ir import IR_MIN_LENGTH, IRResult, ir_de_novo
from plastid.core.model import Genome, Interval
from plastid.core.sequence import gc_percent
from plastid.external.fastp import FastpResult, fastp_version, resolve_fastp
from plastid.external.fastp import trim as run_fastp
from plastid.external.getorganelle import (
    ASSEMBLER_WORD_LIMIT,
    ENTRY_POINT,
    GETORGANELLE_VERSION,
    ORGANELLE,
    SPADES_VERSION,
    GetOrganelleTools,
    read_output,
    resolve_database,
    resolve_tools,
)
from plastid.external.getorganelle import assemble as run_getorganelle
from plastid.external.minimap2 import map_long, map_pair, minimap2_version, resolve_minimap2
from plastid.external.oatk import MITOCHONDRION_PROFILES, OATK_VERSION, PLASTID_PROFILES, OatkTools
from plastid.external.oatk import assemble as run_oatk
from plastid.external.oatk import read_output as read_oatk_output
from plastid.external.oatk import resolve_database as resolve_oatk_database
from plastid.external.oatk import resolve_tools as resolve_oatk_tools
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.contract import (
    ContractError,
    describe_external_result,
    sha256_file,
    validate_external_result,
)
from plastid.run.defaults import load_defaults
from plastid.run.tools import ToolError

__all__ = [
    "ASSEMBLY_PARAMS",
    "BACKEND_IMPORT",
    "BACKEND_LOCAL",
    "MODULE_ID",
    "MODULE_VERSION",
    "PREPARED_FILENAME",
    "Assembled",
    "BackendUnavailable",
    "Evaluation",
    "LongReads",
    "Machine",
    "Measured",
    "ReadPair",
    "assemble_imported",
    "assemble_local",
    "assembly_result",
    "clean_work",
    "evaluate",
    "load_prepared",
    "measure_machine",
    "preflight",
    "prepare",
    "resolve_params",
    "word_limit",
]

MODULE_ID = "assembly"
MODULE_VERSION = "0.1.0"

BACKEND_LOCAL = "local"
BACKEND_IMPORT = "import"

_HIFI_PRESET = "map-hifi"

PREPARED_FILENAME = "prepared.json"
PREPARED_SCHEMA_VERSION = 1

_BYTES_PER_GB = 1_000_000_000
_FASTA_WIDTH = 70

# Every key `--set` may carry on this step, so a misspelling is refused rather than recorded.
ASSEMBLY_PARAMS: tuple[str, ...] = (
    "rounds",
    "kmers",
    "random_seed",
    "threads",
    "memory_gb",
    "max_words",
    "syncmer",
    "min_coverage",
    "max_bases",
    "oatk_db",
    "long_extension",
    "memory_reserve_gb",
    "word_bytes",
    "word_bytes_per_base",
    "word_size_share",
    "trim",
    "tools",
    "database",
    "keep_work",
    "flank",
    "min_clip",
    "min_spanning",
    "clip_window",
    "extension",
    "min_ram_gb",
    "disk_factor",
    "data_label",
)


class BackendUnavailable(Exception):
    """This machine cannot run the step. The message names the reason; nothing fell back."""


@dataclass(frozen=True, slots=True)
class ReadPair:
    """One sample's paired reads, as given. ``sample`` names the run and every output."""

    sample: str
    reads_1: Path
    reads_2: Path

    @property
    def files(self) -> tuple[Path, ...]:
        """The read files, in the order they are hashed."""
        return self.reads_1, self.reads_2


@dataclass(frozen=True, slots=True)
class LongReads:
    """One sample's PacBio HiFi reads, as given, all of them: nuclear and organellar together."""

    sample: str
    reads: Path

    @property
    def files(self) -> tuple[Path, ...]:
        """The read files, in the order they are hashed."""
        return (self.reads,)


@dataclass(frozen=True, slots=True)
class Machine:
    """What the step can use here: memory in GB, free disk in GB where the run writes, cores."""

    ram_gb: float
    disk_free_gb: float
    cores: int


@dataclass(frozen=True, slots=True)
class ReadQC:
    """Read counts around quality control, both mates together; not measured when it did not run."""

    reads_before: int
    bases_before: int
    reads_after: int
    bases_after: int
    q30_before: float
    q30_after: float
    adapter_trimmed_reads: int
    tool_version: str


@dataclass(frozen=True, slots=True)
class Assembled:
    """An assembly from either backend: the one object everything downstream reads."""

    sample: str
    backend: str
    status: str
    status_text: str
    paths: tuple[AssembledPath, ...]
    graph_texts: Mapping[str, str]
    log_text: str
    tool_version: str
    dependencies: Mapping[str, str]
    seed_db: Mapping[str, str]
    base_coverage: float | None
    reads_used: int | None
    architecture: tuple[str, int, int, int, int] | None
    where: str  # a person's sentence: where this assembly ran
    mapping_reads: tuple[Path, ...]  # always the reads as they arrived, whatever the assembler got
    read_qc: ReadQC | None = None
    external_result: Mapping[str, Any] | None = None
    graph_readings: tuple[GraphReading, ...] = ()
    word_limit: int | None = None  # the limit this run was given; None when it ran elsewhere
    word_limit_hit: bool = False
    assembler_name: str = "GetOrganelle"
    citations: tuple[str, ...] = ("getorganelle", "spades", "fastp", "minimap2")
    mapping_preset: str = "sr"  # how minimap2 reads these reads: its preset for the technology


@dataclass(frozen=True, slots=True)
class Measured:
    """What the mapped reads say about one sequence, whoever wrote that sequence."""

    support: str
    profile: SupportProfile
    calls: tuple[BoundaryCall, ...]
    weak: tuple[WeakRegion, ...]
    ir: IRResult | None


@dataclass(frozen=True, slots=True)
class Evaluation:
    """What this toolkit measured about an assembly, none of it taken from the assembler."""

    support: str
    profile: SupportProfile | None
    calls: tuple[BoundaryCall, ...]
    weak: tuple[WeakRegion, ...]
    ir: IRResult | None
    comparisons: tuple[tuple[str, ReferenceComparison], ...]
    reference_name: str | None
    reference_sha256: str | None
    # The same measurement on the reference, made when no path equals it: a disagreement has two
    # possible causes, and the reads can speak to both sides.
    reference_measured: Measured | None = None
    graph_call: str = GRAPH_NOT_MEASURED
    graph_reasons: tuple[str, ...] = ()

    @property
    def usable(self) -> bool:
        """Closed, every base known, read-supported, one molecule: all a next step may consume."""
        return self.support == SUPPORT_SUPPORTED and self.graph_call == GRAPH_RESOLVED


def resolve_params(params: Mapping[str, Any]) -> dict[str, Any]:
    """The step's parameters: packaged defaults, overridden by ``--set``, typed once here."""
    defaults = dict(load_defaults()["assembly"])
    support = dict(defaults.pop("support"))
    merged: dict[str, Any] = {**defaults, **support, **params}
    kmers = merged["kmers"]
    if isinstance(kmers, str):
        kmers = [int(part) for part in kmers.split(",") if part.strip()]
    elif isinstance(kmers, int):
        kmers = [kmers]
    merged["kmers"] = tuple(int(kmer) for kmer in kmers)
    for key in (
        "rounds",
        "random_seed",
        "threads",
        "flank",
        "min_clip",
        "min_spanning",
        "clip_window",
        "extension",
        "syncmer",
        "min_coverage",
        "long_extension",
    ):
        merged[key] = int(merged[key])
    merged["max_bases"] = int(float(merged["max_bases"]))
    for key in (
        "min_ram_gb",
        "disk_factor",
        "memory_reserve_gb",
        "word_bytes",
        "word_bytes_per_base",
        "word_size_share",
    ):
        merged[key] = float(merged[key])
    if merged["memory_gb"] is not None:
        merged["memory_gb"] = float(merged["memory_gb"])
    # Absolute, because a tool runs in the run directory and these are named from wherever the
    # command was typed.
    for key in ("tools", "database", "oatk_db"):
        if merged[key]:
            merged[key] = str(Path(merged[key]).resolve())
    # The assembler's own notation for this number is 3E7, so that spelling has to work here.
    if merged["max_words"] is not None:
        merged["max_words"] = int(float(merged["max_words"]))
    for key in ("trim", "keep_work"):
        value = merged[key]
        merged[key] = value if isinstance(value, bool) else str(value).lower() in ("yes", "true")
    return merged


def measure_machine(directory: Path) -> Machine:
    """Physical memory and the free disk under ``directory``, which is where the run writes."""
    existing = directory
    while not existing.exists():
        existing = existing.parent
    ram = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    return Machine(
        ram_gb=ram / _BYTES_PER_GB,
        disk_free_gb=shutil.disk_usage(existing).free / _BYTES_PER_GB,
        cores=os.cpu_count() or 1,
    )


def word_limit(memory_gb: float, read_length: int, params: Mapping[str, Any]) -> int:
    """The words the assembler may hold in this much memory.

    Peak memory follows the words read recruitment holds, and neither the input's size nor its
    plastid share predicts how many that will be (``docs/tools/getorganelle.md``), so the words are
    limited rather than guessed. A run that reaches the limit stops recruiting and assembles what it
    has, and says so. The word size is the assembler's choice, made after it starts; the limit plans
    for a share of the read length.
    """
    room = max(memory_gb - params["memory_reserve_gb"], 1.0) * _BYTES_PER_GB
    per_word = (
        params["word_bytes"]
        + params["word_bytes_per_base"] * params["word_size_share"] * read_length
    )
    return min(ASSEMBLER_WORD_LIMIT, int(room / per_word))


def _read_length(path: Path, *, reads: int = 1_000) -> int:
    """The longest of a file's first reads, in bases. The assembler takes both forms of FASTQ."""
    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rt", encoding="ascii") if compressed else path.open() as handle:
        first = itertools.islice(FastqGeneralIterator(handle), reads)
        lengths = [len(sequence) for _, sequence, _ in first]
    if not lengths:
        raise AssemblyError(f"{path}: no reads")
    return max(lengths)


def _input_gb(reads: ReadPair | LongReads) -> float:
    return sum(file.stat().st_size for file in reads.files) / _BYTES_PER_GB


def preflight(
    reads: ReadPair | LongReads,
    params: Mapping[str, Any],
    *,
    machine: Machine,
    check_resources: bool,
) -> tuple[GetOrganelleTools | OatkTools, dict[str, Any]]:
    """Everything the local backend needs, or one named reason it cannot run here.

    Memory is a floor, not a prediction: peak memory follows the share of plastid reads in the
    sample much more than the size of the input, so a file's size cannot say what a run will use
    (``docs/tools/getorganelle.md``). What keeps a run inside the machine is the limit on the words
    the assembler holds, and the limit given to SPAdes, whose refusal the caller reports the same
    way as this one. Oatk's limit is the bases it takes, and it writes megabytes, so HiFi reads ask
    nothing of the disk however large the input is.
    """
    long_reads = isinstance(reads, LongReads)
    try:
        if long_reads:
            tools = resolve_oatk_tools(params.get("tools"))
            database = resolve_oatk_database(params.get("oatk_db"))
        else:
            tools = resolve_tools(params.get("tools"))
            database = resolve_database(params.get("database"))
            if params["trim"]:
                fastp_version(resolve_fastp(params.get("tools")))
        minimap2_version(resolve_minimap2(params.get("tools")))
    except ToolError as error:
        raise BackendUnavailable(str(error)) from error
    if check_resources:
        if machine.ram_gb < params["min_ram_gb"]:
            raise BackendUnavailable(
                f"this step needs at least {params['min_ram_gb']:g} GB of memory, the machine "
                f"has {machine.ram_gb:.1f} GB"
            )
        needed = 0.0 if long_reads else _input_gb(reads) * params["disk_factor"]
        if machine.disk_free_gb < needed:
            raise BackendUnavailable(
                f"the assembly needs ~{needed:.1f} GB of disk for this input, the output "
                f"directory has {machine.disk_free_gb:.1f} GB free"
            )
    return tools, database


def _read_qc(result: FastpResult) -> ReadQC:
    return ReadQC(
        reads_before=result.reads_before,
        bases_before=result.bases_before,
        reads_after=result.reads_after,
        bases_after=result.bases_after,
        q30_before=result.q30_before,
        q30_after=result.q30_after,
        adapter_trimmed_reads=result.adapter_trimmed_reads,
        tool_version=result.version,
    )


def _trim(reads: ReadPair, params: Mapping[str, Any], ctx: RunContext) -> FastpResult:
    return run_fastp(
        reads.reads_1,
        reads.reads_2,
        out_1=ctx.out("work/reads_1", "fastq.gz"),
        out_2=ctx.out("work/reads_2", "fastq.gz"),
        report=ctx.out("fastp", "json"),
        html=ctx.out("fastp", "html"),
        threads=params["threads"],
        ctx=ctx,
        binary=resolve_fastp(params.get("tools")),
    )


def _assembled(
    reads: ReadPair | LongReads,
    output: AssemblerOutput,
    *,
    backend: str,
    where: str,
    mapping_reads: tuple[Path, ...],
    read_qc: ReadQC | None,
    external_result: Mapping[str, Any] | None = None,
    word_limit: int | None = None,
) -> Assembled:
    facts = output.facts
    from_hifi = {
        "assembler_name": "Oatk",
        "citations": ("oatk", "minimap2"),
        "mapping_preset": _HIFI_PRESET,
    }
    return Assembled(
        sample=reads.sample,
        backend=backend,
        status=output.status,
        status_text=facts.status_text or "no result stated in the log",
        paths=output.paths,
        graph_texts={graph.name: graph.read_text(encoding="utf-8") for graph in output.graphs},
        log_text=output.log.read_text(encoding="utf-8", errors="replace"),
        tool_version=facts.version or "not stated in the log",
        dependencies=facts.dependencies,
        seed_db=facts.seed_db,
        base_coverage=facts.base_coverage,
        reads_used=facts.reads_used,
        architecture=facts.architecture,
        where=where,
        mapping_reads=mapping_reads,
        read_qc=read_qc,
        external_result=external_result,
        graph_readings=output.graph_readings,
        word_limit=word_limit,
        word_limit_hit=facts.word_limit_hit,
        **(from_hifi if isinstance(reads, LongReads) else {}),
    )


def _assemble_hifi(
    reads: LongReads,
    params: Mapping[str, Any],
    ctx: RunContext,
    *,
    tools: OatkTools,
    database: Mapping[str, Any],
    machine: Machine,
) -> Assembled:
    """HiFi reads through Oatk. A HiFi read is already a consensus, so no quality control runs."""
    try:
        output = run_oatk(
            reads.reads,
            out_dir=ctx.out("work/oatk", "d"),
            tools=tools,
            database=str(database["path"]),
            syncmer=params["syncmer"],
            min_coverage=params["min_coverage"],
            max_bases=params["max_bases"],
            threads=params["threads"],
            ctx=ctx,
        )
    except ToolError as error:
        raise BackendUnavailable(str(error)) from error
    return _assembled(
        reads,
        output,
        backend=BACKEND_LOCAL,
        where=(
            f"local: {machine.cores} cores, {machine.ram_gb:.1f} GB memory, "
            f"{params['threads']} threads, the first {params['max_bases']} bases"
        ),
        mapping_reads=reads.files,
        read_qc=None,
    )


def assemble_local(
    reads: ReadPair | LongReads,
    params: Mapping[str, Any],
    ctx: RunContext,
    *,
    tools: GetOrganelleTools | OatkTools,
    database: Mapping[str, Any],
    machine: Machine,
) -> Assembled:
    """The local backend: quality control, then the assembler, on this machine."""
    if isinstance(reads, LongReads):
        return _assemble_hifi(reads, params, ctx, tools=tools, database=database, machine=machine)
    read_qc: ReadQC | None = None
    pair = (reads.reads_1, reads.reads_2)
    if params["trim"]:
        trimmed = _trim(reads, params, ctx)
        pair = (trimmed.reads_1, trimmed.reads_2)
        read_qc = _read_qc(trimmed)
        read_length = trimmed.read_length_before
    else:
        read_length = _read_length(reads.reads_1)
    memory_gb = params.get("memory_gb") or machine.ram_gb
    max_words = params.get("max_words") or word_limit(memory_gb, read_length, params)
    try:
        output = run_getorganelle(
            pair[0],
            pair[1],
            out_dir=ctx.out("work/getorganelle", "d"),
            tools=tools,
            database=str(database["path"]),
            rounds=params["rounds"],
            kmers=params["kmers"],
            threads=params["threads"],
            random_seed=params["random_seed"],
            memory_gb=max(1, int(memory_gb)),
            max_words=max_words,
            ctx=ctx,
        )
    except ToolError as error:
        raise BackendUnavailable(
            f"the assembler stopped on this machine ({machine.ram_gb:.1f} GB memory, SPAdes limit "
            f"{memory_gb} GB): {error}"
        ) from error
    return _assembled(
        reads,
        output,
        backend=BACKEND_LOCAL,
        where=(
            f"local: {machine.cores} cores, {machine.ram_gb:.1f} GB memory, "
            f"{params['threads']} threads"
        ),
        mapping_reads=(reads.reads_1, reads.reads_2),
        read_qc=read_qc,
        word_limit=max_words,
    )


def clean_work(ctx: RunContext, *, keep: bool = False) -> None:
    """Drop the working files. They are tens of times the size of what the run keeps."""
    work = ctx.out("work/marker", "tmp").parent
    if not keep:
        shutil.rmtree(work, ignore_errors=True)


def _hashed(path: Path, *, relative_to: Path | None = None) -> dict[str, Any]:
    """A file by digest. Inside the prepared directory it is named relatively, so it can move."""
    inside = relative_to is not None and path.parent == relative_to
    name = path.name if inside else str(path.resolve())
    return {"path": name, "sha256": sha256_file(path), "bytes": path.stat().st_size}


def _prepared_reads(payload: Mapping[str, Any], source: Path, key: str) -> tuple[Path, ...]:
    """The read files a prepared record names, resolved against the record's own directory."""
    files = [Path(entry["path"]) for entry in payload[key]]
    return tuple(file if file.is_absolute() else source.parent / file for file in files)


def _command(prepared_1: str, prepared_2: str, params: Mapping[str, Any]) -> list[str]:
    """The assembler's command for a run elsewhere, as close to the local backend as it can be."""
    return [
        ENTRY_POINT,
        "-1",
        prepared_1,
        "-2",
        prepared_2,
        "-F",
        ORGANELLE,
        "-o",
        "output",
        "-R",
        str(params["rounds"]),
        "-k",
        ",".join(str(kmer) for kmer in params["kmers"]),
        "--random-seed",
        str(params["random_seed"]),
    ]


def _hifi_command(reads: str, params: Mapping[str, Any]) -> list[str]:
    """Oatk's command for a run elsewhere. The stderr of that run is the log the import reads."""
    return [
        "oatk",
        "-k",
        str(params["syncmer"]),
        "-c",
        str(params["min_coverage"]),
        "-D",
        str(params["max_bases"]),
        "-m",
        MITOCHONDRION_PROFILES,
        "-p",
        PLASTID_PROFILES,
        "-o",
        "oatk",
        reads,
    ]


def _prepare_hifi(reads: LongReads, params: Mapping[str, Any]) -> dict[str, Any]:
    entry = _hashed(reads.reads)
    return {
        "prepared_schema_version": PREPARED_SCHEMA_VERSION,
        "sample": reads.sample,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "raw_reads": [entry],
        "input_reads": [entry],
        "read_qc": None,
        "tool": {
            "name": "oatk",
            "version": OATK_VERSION,
            "command": _hifi_command(reads.reads.name, params),
        },
        "params": {
            key: params[key] for key in ("syncmer", "min_coverage", "max_bases", "random_seed")
        },
    }


def prepare(
    reads: ReadPair | LongReads, params: Mapping[str, Any], ctx: RunContext
) -> dict[str, Any]:
    """Quality-control the reads here and record what a run elsewhere must be given.

    The reads an outside service receives are the ones hashed here, so the import can prove the
    result it receives was asked for. Nothing leaves this machine in this step.
    """
    if isinstance(reads, LongReads):
        return _prepare_hifi(reads, params)
    read_qc: dict[str, Any] | None = None
    here = ctx.out("prepared", "json").parent
    if params["trim"]:
        trimmed = _trim(reads, params, ctx)
        out_1 = ctx.out("reads_1", "fastq.gz")
        out_2 = ctx.out("reads_2", "fastq.gz")
        shutil.move(trimmed.reads_1, out_1)
        shutil.move(trimmed.reads_2, out_2)
        read_qc = {
            "reads_before": trimmed.reads_before,
            "bases_before": trimmed.bases_before,
            "reads_after": trimmed.reads_after,
            "bases_after": trimmed.bases_after,
            "q30_before": trimmed.q30_before,
            "q30_after": trimmed.q30_after,
            "adapter_trimmed_reads": trimmed.adapter_trimmed_reads,
            "tool_version": trimmed.version,
        }
    else:
        out_1, out_2 = reads.reads_1, reads.reads_2
    return {
        "prepared_schema_version": PREPARED_SCHEMA_VERSION,
        "sample": reads.sample,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "raw_reads": [_hashed(reads.reads_1), _hashed(reads.reads_2)],
        "input_reads": [_hashed(out_1, relative_to=here), _hashed(out_2, relative_to=here)],
        "read_qc": read_qc,
        "tool": {
            "name": "getorganelle",
            "version": GETORGANELLE_VERSION,
            "command": _command(out_1.name, out_2.name, params),
        },
        "params": {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in params.items()
            if key in ("rounds", "kmers", "random_seed", "trim")
        },
    }


def load_prepared(path: str | Path) -> dict[str, Any]:
    """Read a prepared record and check that the reads it names are still those bytes."""
    source = Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise AssemblyError(f"{source}: cannot be read as a prepared record: {error}") from error
    if payload.get("prepared_schema_version") != PREPARED_SCHEMA_VERSION:
        raise AssemblyError(
            f"{source}: prepared_schema_version "
            f"{payload.get('prepared_schema_version')!r} is not {PREPARED_SCHEMA_VERSION}"
        )
    # The raw reads as well as those handed over: the evaluation of an import reads them.
    for key in ("input_reads", "raw_reads"):
        for entry, file in zip(payload[key], _prepared_reads(payload, source, key), strict=True):
            if not file.is_file():
                raise AssemblyError(f"{source}: the prepared reads {file} are gone")
            digest = sha256_file(file)
            if digest != entry["sha256"]:
                raise AssemblyError(
                    f"{file}: sha256 {digest} is not the {entry['sha256']} that was prepared — "
                    "this is not the input the imported result was asked for"
                )
    return payload


def assemble_imported(
    prepared_path: str | Path,
    result_dir: str | Path,
    *,
    source_name: str,
    retrieved_at: str,
    source_url: str | None = None,
) -> Assembled:
    """The import backend: an assembly made elsewhere, accepted only through the contract."""
    prepared_file = Path(prepared_path)
    prepared = load_prepared(prepared_file)
    tool_name = prepared["tool"]["name"]
    from_hifi = tool_name == "oatk"
    output = read_oatk_output(result_dir) if from_hifi else read_output(result_dir)
    facts = output.facts
    if facts.version is None:
        raise AssemblyError(f"{output.log}: the log does not state the assembler's version")
    source = {"name": source_name, "retrieved_at": retrieved_at}
    if source_url:
        source["url"] = source_url
    raw_reads = _prepared_reads(prepared, prepared_file, "raw_reads")
    read_qc = prepared.get("read_qc")
    expected = read_qc["reads_after"] if read_qc else None
    if expected is not None and facts.reads_used is not None and facts.reads_used > expected:
        raise AssemblyError(
            f"{output.log}: the run used {facts.reads_used} reads, more than the {expected} that "
            "were prepared — it did not run on the prepared input"
        )
    has_paths = bool(output.path_files)
    try:
        description = describe_external_result(
            tool_name=tool_name,
            tool_version=facts.version,
            argv=prepared["tool"]["command"],
            input_path=prepared_file,
            query_map={prepared["sample"]: prepared["sample"]},
            source=source,
            raw_path=output.path_files[0] if has_paths else output.log,
            status="ok" if has_paths else "no_result",
            count=len(output.paths) if has_paths else None,
        )
        validate_external_result(description)
    except ContractError as error:
        raise AssemblyError(f"the imported result was refused: {error}") from error
    reads = (
        LongReads(prepared["sample"], raw_reads[0])
        if from_hifi
        else ReadPair(prepared["sample"], raw_reads[0], raw_reads[1])
    )
    return _assembled(
        reads,
        output,
        backend=BACKEND_IMPORT,
        where=f"imported: {source_name}, retrieved {date.fromisoformat(retrieved_at).isoformat()}",
        mapping_reads=raw_reads,
        read_qc=ReadQC(**read_qc) if read_qc else None,
        external_result=description,
    )


def _target_fasta(sequence: str, extension: int, path: Path) -> None:
    """The mapping target: the sequence and its own first bases, so a read can cross the origin."""
    extended = sequence + sequence[: min(extension, len(sequence))]
    lines = [">target"]
    lines.extend(
        extended[start : start + _FASTA_WIDTH] for start in range(0, len(extended), _FASTA_WIDTH)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _named_boundaries(sequence: str, ir: IRResult | None) -> list[tuple[str, int]]:
    """The origin, and the four junctions when an inverted repeat was found."""
    boundaries = [("origin", 0)]
    if ir is not None:
        first, second = sorted((ir.ira, ir.irb), key=lambda interval: interval.start)
        boundaries.extend(
            ir_junctions(first, second, length=len(sequence), frame_offset=ir.frame_offset)
        )
    return boundaries


def _measure(
    sequence: str,
    *,
    name: str,
    tag: str,
    reads: tuple[Path, ...],
    preset: str,
    params: Mapping[str, Any],
    ctx: RunContext,
) -> Measured:
    """Map the reads onto one sequence as a circle and read every join."""
    long_reads = len(reads) == 1
    reach = params["long_extension"] if long_reads else params["extension"]
    extension = min(reach, len(sequence))
    target = ctx.out(f"work/{tag}_target", "fasta")
    _target_fasta(sequence, extension, target)
    mapped = ctx.out(f"work/{tag}_mapped", "sam")
    binary = resolve_minimap2(params.get("tools"))
    if long_reads:
        sam = map_long(
            target,
            reads[0],
            preset=preset,
            out=mapped,
            threads=params["threads"],
            ctx=ctx,
            binary=binary,
        )
    else:
        sam = map_pair(
            target,
            reads[0],
            reads[1],
            out=mapped,
            threads=params["threads"],
            ctx=ctx,
            binary=binary,
        )
    with sam.open(encoding="utf-8") as handle:
        profile = pool_repeated_joins(
            build_profile(
                parse_sam(handle),
                length=len(sequence),
                extension=extension,
                flank=params["flank"],
                min_clip=params["min_clip"],
            ),
            sequence,
        )
    sam.unlink()
    ir = ir_de_novo(Genome(name=name, sequence=sequence, topology="circular"))
    rule = SupportRule(min_spanning=params["min_spanning"], clip_window=params["clip_window"])
    calls = tuple(
        call_boundary(profile, rule, name=label, boundary=boundary)
        for label, boundary in _named_boundaries(sequence, ir)
    )
    weak = weak_regions(profile, rule)
    all_supported = all(call.call == SUPPORT_SUPPORTED for call in calls)
    return Measured(
        support=SUPPORT_SUPPORTED if all_supported and not weak else SUPPORT_UNDETERMINED,
        profile=profile,
        calls=calls,
        weak=weak,
        ir=ir,
    )


def evaluate(
    assembled: Assembled,
    params: Mapping[str, Any],
    ctx: RunContext,
    *,
    reference: Genome | None = None,
    reference_sha256: str | None = None,
) -> Evaluation:
    """Measure the assembly: map the reads back, read every join, compare with the reference.

    Only a closed path with every base known can be mapped as a circle, so an open assembly
    reports its read support as not measured, not as zero. Mapping runs on the first path:
    two isomers share every junction, because they differ only in the orientation of a region
    between identical repeats, and the graph has shown that of this run when it calls it resolved.
    Reads shorter than the repeat cannot separate them, and cannot separate wider rearrangements
    either, which is why a graph with more structures than that is never passed on the reads.

    The reads mapped are the raw ones, whatever the assembler was given. Quality control serves the
    assembler; a check that read its output would lose the evidence it removed, and the aligner
    clips adapter and poor ends by itself. Trimming took a correct *Selaginella* from five crossing
    reads to one at one spot, with the sequence unchanged.

    When a reference is given and no path equals it, the reference is measured the same way, even
    when the assembly stayed open. An assembly and a published record that disagree are then both
    tested by the same reads, and neither is taken as the truth.
    """
    comparisons: tuple[tuple[str, ReferenceComparison], ...] = ()
    closed = assembled.status == STATUS_COMPLETE and bool(assembled.paths)
    if reference is not None and closed:
        comparisons = tuple(
            (path.name, compare_to_reference(path.contigs[0], reference.sequence))
            for path in assembled.paths
        )
    reference_name = reference.accession or reference.name if reference is not None else None
    graph_call, graph_reasons = graph_verdict(assembled.graph_readings)
    reference_measured = None
    if reference is not None and not any(item.identical for _, item in comparisons):
        reference_measured = _measure(
            reference.sequence,
            name=reference.accession or reference.name,
            tag="reference",
            reads=assembled.mapping_reads,
            preset=assembled.mapping_preset,
            params=params,
            ctx=ctx,
        )
    if not closed:
        return Evaluation(
            support=SUPPORT_NOT_MEASURED,
            profile=None,
            calls=(),
            weak=(),
            ir=None,
            comparisons=comparisons,
            reference_name=reference_name,
            reference_sha256=reference_sha256,
            reference_measured=reference_measured,
            graph_call=graph_call,
            graph_reasons=graph_reasons,
        )

    measured = _measure(
        assembled.paths[0].contigs[0],
        name=assembled.sample,
        tag="assembly",
        reads=assembled.mapping_reads,
        preset=assembled.mapping_preset,
        params=params,
        ctx=ctx,
    )
    return Evaluation(
        support=measured.support,
        profile=measured.profile,
        calls=measured.calls,
        weak=measured.weak,
        ir=measured.ir,
        comparisons=comparisons,
        reference_name=reference_name,
        reference_sha256=reference_sha256,
        reference_measured=reference_measured,
        graph_call=graph_call,
        graph_reasons=graph_reasons,
    )


def _fasta_text(name: str, sequences: Sequence[str]) -> str:
    """FASTA with a bare identifier: the next step reads the description as the organism's name."""
    records: list[str] = []
    identifier = "_".join(name.split())
    for index, sequence in enumerate(sequences, start=1):
        records.append(f">{identifier}" if len(sequences) == 1 else f">{identifier}_{index}")
        records.extend(
            sequence[start : start + _FASTA_WIDTH]
            for start in range(0, len(sequence), _FASTA_WIDTH)
        )
    return "\n".join(records) + "\n"


def _base_1based(boundary: int, length: int) -> int:
    """The base just before a boundary, in GenBank coordinates; boundary 0 follows the last base."""
    before = (boundary - 1) % length
    return Interval(before, before + 1).to_1based()[0]


def _regions(assembled: Assembled, evaluation: Evaluation) -> dict[str, int | str]:
    """The four regions from the sequence's own repeat, or a named absence."""
    ir = evaluation.ir
    if ir is None or not assembled.paths:
        missing = "not measured" if evaluation.profile is None else "no inverted repeat found"
        return {"ir_bp": missing, "lsc_bp": missing, "ssc_bp": missing}
    length = assembled.paths[0].length
    first, second = sorted((ir.ira, ir.irb), key=lambda interval: interval.start)
    arcs = sorted((second.start - first.end, length - second.end + first.start), reverse=True)
    return {"ir_bp": ir.length, "lsc_bp": arcs[0], "ssc_bp": arcs[1]}


def _versions(stated: Mapping[str, str]) -> str:
    """Name and version of each entry, as the assembler's log states them."""
    return "; ".join(f"{name} {version}" for name, version in stated.items()) or "not stated"


def _word_limit_cell(assembled: Assembled) -> int | str:
    if assembled.word_limit is not None:
        return assembled.word_limit
    return "not stated" if assembled.assembler_name == "GetOrganelle" else "not applicable"


def _summary(assembled: Assembled, evaluation: Evaluation, scope: Scope) -> Table:
    path = assembled.paths[0] if assembled.paths else None
    profile = evaluation.profile
    sequence = "".join(path.contigs) if path else ""
    row: dict[str, Any] = {
        "sample": assembled.sample,
        "backend": assembled.backend,
        "where": assembled.where,
        "assembler": f"{assembled.assembler_name} {assembled.tool_version}",
        "assembler_programs": _versions(assembled.dependencies),
        "seed_database": _versions(assembled.seed_db),
        "status": assembled.status,
        "assembler_says": assembled.status_text,
        "read_support": evaluation.support,
        "graph_call": evaluation.graph_call,
        "graph_reasons": "; ".join(evaluation.graph_reasons) or "nothing",
        "n_paths": len(assembled.paths),
        "length_bp": path.length if path else "not measured",
        "gc_percent": round(gc_percent(sequence), 2) if sequence else "not measured",
        "assembler_coverage": (
            assembled.base_coverage if assembled.base_coverage is not None else "not measured"
        ),
        "reads_used": assembled.reads_used if assembled.reads_used is not None else "not measured",
        "word_limit": _word_limit_cell(assembled),
        "word_limit_hit": "yes" if assembled.word_limit_hit else "no",
        "mapped_reads": profile.reads if profile else "not measured",
        "mean_depth": round(profile.mean_depth, 2) if profile else "not measured",
        "sd_depth": round(profile.sd_depth, 2) if profile else "not measured",
        "error_rate": (
            round(profile.error_rate, 5)
            if profile and profile.error_rate is not None
            else "not measured"
        ),
        "fewest_spanning": min(profile.spanning) if profile else "not measured",
        **_regions(assembled, evaluation),
        "assembler_repeat": (
            f"{assembled.architecture[4]} bp, {assembled.architecture[0]}"
            if assembled.architecture
            else "not measured"
        ),
    }
    columns = (
        Column("sample", "Sample"),
        Column("backend", "Backend"),
        Column("where", "Where it ran"),
        Column("assembler", "Assembler"),
        Column("assembler_programs", "Programs it ran"),
        Column("seed_database", "Seed database"),
        Column("status", "Status"),
        Column("assembler_says", "The assembler's own words"),
        Column("read_support", "Read support"),
        Column("graph_call", "Graph"),
        Column("graph_reasons", "What the graph leaves open"),
        Column("n_paths", "Paths", dtype="int"),
        Column("length_bp", "Length", dtype="int", unit="bp"),
        Column("gc_percent", "GC", dtype="float", unit="%", precision=2),
        Column("assembler_coverage", "Coverage, by the assembler", dtype="float", precision=2),
        Column("reads_used", "Reads used by the assembler", dtype="int"),
        Column("word_limit", "Word limit given to the assembler", dtype="int"),
        Column("word_limit_hit", "Reached it"),
        Column("mapped_reads", "Reads mapped back", dtype="int"),
        Column("mean_depth", "Mean depth", dtype="float", precision=2),
        Column("sd_depth", "SD of depth", dtype="float", precision=2),
        Column("error_rate", "Errors per matched base", dtype="float", precision=5),
        Column("fewest_spanning", "Fewest spanning reads at any boundary", dtype="int"),
        Column("ir_bp", "Inverted repeat", dtype="int", unit="bp"),
        Column("lsc_bp", "LSC", dtype="int", unit="bp"),
        Column("ssc_bp", "SSC", dtype="int", unit="bp"),
        Column("assembler_repeat", "The assembler's reading of the largest repeat"),
    )
    return Table(
        id="assembly_summary",
        title="Assembly summary",
        columns=columns,
        rows=(row,),
        scope=scope,
    )


def _paths_table(assembled: Assembled, scope: Scope) -> Table:
    return Table(
        id="assembly_paths",
        title="Assembled paths",
        columns=(
            Column("path", "Path"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("contigs", "Contigs", dtype="int"),
            Column("circular", "Circular"),
        ),
        rows=tuple(
            {
                "path": path.name,
                "length_bp": path.length,
                "contigs": len(path.contigs),
                "circular": "yes" if path.circular else "no",
            }
            for path in assembled.paths
        ),
        key="path",
        scope=scope,
    )


def _isomers_cell(reading: GraphReading) -> str:
    if reading.structures is None:
        return "not measured"
    if reading.isomers is None:
        return "not applicable"  # measured: not two structures, so there is no pair to call
    return "yes" if reading.isomers else "no"


def _repeat_cell(reading: GraphReading) -> str:
    if reading.structures is None:
        return "not measured"
    if reading.repeat is None:
        return "none"
    return f"{reading.repeat[1]} bp, {reading.repeat[0]}"


def _graph_table(readings: Sequence[GraphReading], scope: Scope) -> Table:
    missing = "not measured"
    return Table(
        id="assembly_graph",
        title="The selected graph, read without the assembler's words",
        columns=(
            Column("graph", "Graph"),
            Column("call", "Call"),
            Column("segments", "Segments", dtype="int"),
            Column("links", "Links", dtype="int"),
            Column("components", "Components", dtype="int"),
            Column("free_ends", "Free ends", dtype="int"),
            Column("paths_written", "Paths written", dtype="int"),
            Column("structures", "Structures that close", dtype="int"),
            Column("isomers", "Flip-flop isomers"),
            Column("repeat", "Longest repeat of the first path"),
            Column("reasons", "What it leaves open"),
        ),
        rows=tuple(
            {
                "graph": reading.name,
                "call": reading.call,
                "segments": reading.segments,
                "links": reading.links,
                "components": reading.components,
                "free_ends": reading.free_ends,
                "paths_written": reading.paths_written,
                "structures": reading.structures if reading.structures is not None else missing,
                "isomers": _isomers_cell(reading),
                "repeat": _repeat_cell(reading),
                "reasons": "; ".join(reading.reasons) or "nothing",
            }
            for reading in readings
        ),
        key="graph",
        scope=scope,
    )


def _segments_table(readings: Sequence[GraphReading], scope: Scope) -> Table:
    missing = "not measured"
    return Table(
        id="graph_segments",
        title="Segments of the selected graph, longest first",
        columns=(
            Column("graph", "Graph"),
            Column("segment", "Segment"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("depth", "Depth, as the assembler states it", dtype="float", precision=1),
            Column("copies", "Copies in the first path", dtype="int"),
            Column("depth_ratio", "Depth over one copy's depth", dtype="float", precision=2),
        ),
        rows=tuple(
            {
                "graph": reading.name,
                "segment": row.name,
                "length_bp": row.length,
                "depth": round(row.depth, 1) if row.depth is not None else missing,
                "copies": row.copies if row.copies is not None else missing,
                "depth_ratio": (
                    round(row.depth_ratio, 2) if row.depth_ratio is not None else missing
                ),
            }
            for reading in readings
            for row in reading.rows
        ),
        scope=scope,
    )


def _read_qc_table(read_qc: ReadQC, scope: Scope) -> Table:
    return Table(
        id="read_qc",
        title="Read quality control",
        columns=(
            Column("reads_before", "Reads in", dtype="int"),
            Column("bases_before", "Bases in", dtype="int"),
            Column("q30_before", "Q30 in", dtype="float", precision=4),
            Column("reads_after", "Reads kept", dtype="int"),
            Column("bases_after", "Bases kept", dtype="int"),
            Column("q30_after", "Q30 kept", dtype="float", precision=4),
            Column("adapter_trimmed_reads", "Reads with adapter trimmed", dtype="int"),
            Column("tool", "Tool"),
        ),
        rows=(
            {
                "reads_before": read_qc.reads_before,
                "bases_before": read_qc.bases_before,
                "q30_before": read_qc.q30_before,
                "reads_after": read_qc.reads_after,
                "bases_after": read_qc.bases_after,
                "q30_after": read_qc.q30_after,
                "adapter_trimmed_reads": read_qc.adapter_trimmed_reads,
                "tool": f"fastp {read_qc.tool_version}",
            },
        ),
        scope=scope,
    )


def _junction_table(
    calls: Sequence[BoundaryCall],
    length: int,
    scope: Scope,
    *,
    table_id: str = "junction_support",
    title: str = "Read support at the origin and the inverted-repeat junctions",
) -> Table:
    return Table(
        id=table_id,
        title=title,
        columns=(
            Column("junction", "Junction"),
            Column("after_base", "After base", dtype="int"),
            Column("spanning", "Spanning reads", dtype="int"),
            Column("clipped", "Reads clipped here", dtype="int"),
            Column("fraction", "Spanning fraction", dtype="float", precision=3),
            Column("depth_before", "Depth before", dtype="float", precision=1),
            Column("depth_after", "Depth after", dtype="float", precision=1),
            Column("call", "Call"),
        ),
        rows=tuple(
            {
                "junction": call.name,
                "after_base": _base_1based(call.boundary, length),
                "spanning": call.spanning,
                "clipped": call.clipped,
                "fraction": call.fraction if call.fraction is not None else "not measured",
                "depth_before": call.depth_before,
                "depth_after": call.depth_after,
                "call": call.call,
            }
            for call in calls
        ),
        key="junction",
        scope=scope,
    )


def _weak_table(
    weak: Sequence[WeakRegion],
    length: int,
    scope: Scope,
    *,
    table_id: str = "weak_regions",
    title: str = "Boundaries too few reads cross",
) -> Table:
    return Table(
        id=table_id,
        title=title,
        columns=(
            Column("first_after_base", "First: after base", dtype="int"),
            Column("last_after_base", "Last: after base", dtype="int"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("fewest_spanning", "Fewest spanning reads", dtype="int"),
            Column("most_clipped", "Most reads stopping here", dtype="int"),
        ),
        rows=tuple(
            {
                "first_after_base": _base_1based(region.start, length),
                "last_after_base": _base_1based(region.end, length),
                "length_bp": region.length,
                "fewest_spanning": region.lowest,
                "most_clipped": region.most_clipped,
            }
            for region in weak
        ),
        scope=scope,
    )


def _differences_table(evaluation: Evaluation, reference_length: int, scope: Scope) -> Table:
    """Every place a path and the reference part, so a number can be walked back to a position."""
    rows = [
        {
            "path": name,
            "after_reference_base": _base_1based(item.reference_start, reference_length),
            "reference_bp": item.reference_end - item.reference_start,
            "path_bp": item.query_end - item.query_start,
            "mismatches": item.mismatches,
            "gap_columns": item.gap_columns,
        }
        for name, comparison in evaluation.comparisons
        for item in comparison.differences
    ]
    return Table(
        id="reference_differences",
        title="Where the assembled paths and the reference differ",
        columns=(
            Column("path", "Path"),
            Column("after_reference_base", "After reference base", dtype="int"),
            Column("reference_bp", "Reference", dtype="int", unit="bp"),
            Column("path_bp", "Path", dtype="int", unit="bp"),
            Column("mismatches", "Mismatches", dtype="int"),
            Column("gap_columns", "Gap columns", dtype="int"),
        ),
        rows=tuple(rows),
        scope=scope,
    )


def _reference_support_table(evaluation: Evaluation, scope: Scope) -> Table:
    """The assembly and the reference side by side, under the same reads and the same rule."""
    measured = evaluation.reference_measured
    assert measured is not None
    rows = []
    sides: list[tuple[str | None, str, SupportProfile | None, tuple[WeakRegion, ...]]] = [
        ("assembly", evaluation.support, evaluation.profile, evaluation.weak),
        (evaluation.reference_name, measured.support, measured.profile, measured.weak),
    ]
    missing = "not measured"
    for name, support, profile, weak in sides:
        rows.append(
            {
                "sequence": name,
                "length_bp": profile.length if profile else missing,
                "read_support": support,
                "fewest_spanning": min(profile.spanning) if profile else missing,
                "weak_regions": len(weak) if profile else missing,
                "mean_depth": round(profile.mean_depth, 2) if profile else missing,
                "sd_depth": round(profile.sd_depth, 2) if profile else missing,
                "error_rate": (
                    round(profile.error_rate, 5)
                    if profile and profile.error_rate is not None
                    else missing
                ),
            }
        )
    return Table(
        id="reference_read_support",
        title="The same reads on the assembly and on the reference",
        columns=(
            Column("sequence", "Sequence"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("read_support", "Read support"),
            Column("fewest_spanning", "Fewest spanning reads at any boundary", dtype="int"),
            Column("weak_regions", "Weak regions", dtype="int"),
            Column("mean_depth", "Mean depth", dtype="float", precision=2),
            Column("sd_depth", "SD of depth", dtype="float", precision=2),
            Column("error_rate", "Errors per matched base", dtype="float", precision=5),
        ),
        rows=tuple(rows),
        key="sequence",
        scope=scope,
    )


def _comparison_table(evaluation: Evaluation, scope: Scope) -> Table:
    missing = "not measured"
    rows = []
    for name, comparison in evaluation.comparisons:
        rows.append(
            {
                "path": name,
                "reference": evaluation.reference_name,
                "identical": "yes" if comparison.identical else "no",
                "length_delta": comparison.length_delta,
                "identity": comparison.identity if comparison.identity is not None else missing,
                "mismatches": (
                    comparison.mismatches if comparison.mismatches is not None else missing
                ),
                "gap_columns": (
                    comparison.gap_columns if comparison.gap_columns is not None else missing
                ),
                "origin_gap_columns": (
                    comparison.origin_gap_columns
                    if comparison.origin_gap_columns is not None
                    else missing
                ),
                "differences": len(comparison.differences),
                "method": comparison.method,
            }
        )
    return Table(
        id="reference_comparison",
        title="Assembled paths against the reference",
        columns=(
            Column("path", "Path"),
            Column("reference", "Reference"),
            Column("identical", "Same circular molecule"),
            Column("length_delta", "Path minus reference", dtype="int", unit="bp"),
            Column("identity", "Identity", dtype="float", precision=5),
            Column("mismatches", "Mismatches", dtype="int"),
            Column("gap_columns", "Gap columns", dtype="int"),
            Column("origin_gap_columns", "Of those, at the reference's origin", dtype="int"),
            Column("differences", "Places they differ", dtype="int"),
            Column("method", "Method"),
        ),
        rows=tuple(rows),
        key="path",
        scope=scope,
    )


def _repeat_note(assembled: Assembled, evaluation: Evaluation) -> str | None:
    """The graph's reading of the inverted repeat against the sequence's, when they differ.

    A short inverted repeat is ordinary in a plastome that has no quadripartite form, so only one
    long enough to count as an IR copy is compared.
    """
    if evaluation.profile is None or not assembled.graph_readings:
        return None
    kind, length = assembled.graph_readings[0].repeat or ("", 0)
    in_graph = length if kind == "inverted" and length >= IR_MIN_LENGTH else None
    in_sequence = evaluation.ir.length if evaluation.ir is not None else None
    if in_graph == in_sequence:
        return None
    if in_sequence is None:
        return (
            f"The graph holds an inverted repeat of {in_graph} bp and the sequence's own reading "
            "reports none: it looks only for the lengths it was validated on."
        )
    return (
        f"The inverted repeat reads {in_sequence} bp in the sequence and "
        f"{in_graph if in_graph is not None else 'none'} in the graph; compare the copies before "
        "using a position-based analysis."
    )


def assembly_result(
    assembled: Assembled, evaluation: Evaluation, *, data_label: str = ""
) -> Result:
    """The step's tables and files. Empty tables are left out: an absence is not a table."""
    label = f"plastid assembly of {assembled.sample}; {assembled.where}"
    if data_label:
        label = f"{label}; {data_label}"
    scope = Scope(label=label)
    tables: dict[str, Table] = {"assembly_summary": _summary(assembled, evaluation, scope)}
    if assembled.paths:
        tables["assembly_paths"] = _paths_table(assembled, scope)
    if assembled.graph_readings:
        tables["assembly_graph"] = _graph_table(assembled.graph_readings, scope)
        tables["graph_segments"] = _segments_table(assembled.graph_readings, scope)
    if assembled.read_qc is not None:
        tables["read_qc"] = _read_qc_table(assembled.read_qc, scope)
    if evaluation.calls:
        length = assembled.paths[0].length
        tables["junction_support"] = _junction_table(evaluation.calls, length, scope)
        if evaluation.weak:
            tables["weak_regions"] = _weak_table(evaluation.weak, length, scope)
    if evaluation.comparisons:
        tables["reference_comparison"] = _comparison_table(evaluation, scope)
        if any(comparison.differences for _, comparison in evaluation.comparisons):
            reference_length = evaluation.comparisons[0][1].reference_length
            tables["reference_differences"] = _differences_table(
                evaluation, reference_length, scope
            )
    measured = evaluation.reference_measured
    if measured is not None:
        tables["reference_read_support"] = _reference_support_table(evaluation, scope)
        tables["reference_junctions"] = _junction_table(
            measured.calls,
            measured.profile.length,
            scope,
            table_id="reference_junctions",
            title="Read support at the reference's own junctions",
        )
        if measured.weak:
            tables["reference_weak_regions"] = _weak_table(
                measured.weak,
                measured.profile.length,
                scope,
                table_id="reference_weak_regions",
                title="Boundaries of the reference too few reads cross",
            )

    artifacts: dict[str, str] = {"assembler.log": assembled.log_text}
    for name, text in assembled.graph_texts.items():
        artifacts[name] = text
    closed = assembled.status in (STATUS_COMPLETE, STATUS_GAPPED)
    if assembled.paths and closed:
        artifacts["assembly.fasta"] = _fasta_text(assembled.sample, [assembled.paths[0].contigs[0]])
        if len(assembled.paths) > 1:
            artifacts["assembly.paths.fasta"] = _fasta_text(
                assembled.sample,
                [path.contigs[0] for path in assembled.paths],
            )
    elif assembled.status == STATUS_FRAGMENTED:
        artifacts["scaffolds.fasta"] = _fasta_text(
            assembled.sample,
            [contig for path in assembled.paths for contig in path.contigs],
        )
    if assembled.external_result is not None:
        artifacts["external_result.json"] = (
            json.dumps(assembled.external_result, indent=2, ensure_ascii=False) + "\n"
        )
    notes = (f"The assembler's own words: {assembled.status_text}.",)
    if assembled.status == STATUS_NO_RESULT:
        notes += ("The assembler produced no plastome sequence; this is not a measured zero.",)
    # A local run is refused on any other version, so only an import can differ, and an import is
    # accepted on what it states. The difference has to be read off the table, not out of a log.
    spades = assembled.dependencies.get("SPAdes")
    if spades is not None and spades != SPADES_VERSION:
        notes += (f"This assembly ran on SPAdes {spades}; the local route pins {SPADES_VERSION}.",)
    pinned = {"GetOrganelle": GETORGANELLE_VERSION, "Oatk": OATK_VERSION}[assembled.assembler_name]
    if assembled.tool_version != pinned:
        notes += (
            f"This assembly ran on {assembled.assembler_name} {assembled.tool_version}; the local "
            f"route pins {pinned}.",
        )
    if assembled.word_limit_hit:
        notes += (
            "The assembler stopped recruiting reads at its word limit, and assembled what it had. "
            "More memory, or a run elsewhere, would recruit more.",
        )
    notes += tuple(f"The graph leaves this open: {reason}." for reason in evaluation.graph_reasons)
    repeat_note = _repeat_note(assembled, evaluation)
    if repeat_note is not None:
        notes += (repeat_note,)
    if evaluation.ir is not None and assembled.architecture is not None:
        theirs = assembled.architecture[4]
        if theirs != evaluation.ir.length:
            # Two independent readings of one repeat. Neither is assumed right: the toolkit's
            # detector has a known failure on a single-base indel inside a copy (K-07).
            notes += (
                f"The inverted repeat reads {evaluation.ir.length} bp here and {theirs} bp in the "
                "assembler's own log; compare the copies before using a position-based analysis.",
            )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=notes,
        citations=assembled.citations,
        artifacts=artifacts,
    )
