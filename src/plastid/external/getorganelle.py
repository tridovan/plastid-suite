"""GetOrganelle — de novo plastome assembly from short reads.

**Why this assembler.** It finished the most plastomes in the independent benchmark of seven tools
(Freudenthal et al. 2020), and its authors published the raw comparison this toolkit is judged
against: 50 public read sets with their published plastomes, including the three it could not close
(Jin et al. 2020, Additional file 2). ``docs/tools/getorganelle.md`` has the numbers.

**What is read, and from where.** The run's own log states the result, the dependency versions and
the seed database, and the path files are named after the same result. Both are read, and a
disagreement is refused: the result of an assembly decides what may follow it, so it cannot rest
on one string. The same reader serves a local run and a directory imported from elsewhere, which
is how both backends come to mean the same thing by the same status.

**What a local run is given.** The dependencies are named on the command line (``--which-*``) and
the database by ``--config-dir``, so the run does not depend on the shell it starts from. SPAdes
gets a memory limit, because on a small machine the alternative is the kernel's out-of-memory
killer, which leaves no message at all.
"""

from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO

from plastid.core.assembly import (
    STATUS_COMPLETE,
    STATUS_FRAGMENTED,
    STATUS_GAPPED,
    STATUS_NO_RESULT,
    AssembledPath,
    AssemblyError,
)
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "GETORGANELLE_VERSION",
    "ORGANELLE",
    "SPADES_VERSION",
    "GetOrganelleOutput",
    "GetOrganelleTools",
    "assemble",
    "parse_log",
    "read_output",
    "resolve_database",
    "resolve_tools",
]

GETORGANELLE_VERSION = "1.7.7.1"
SPADES_VERSION = "4.3.0"
GETORGANELLE_DOC = "docs/tools/getorganelle.md"
LOG_NAME = "get_org.log.txt"
ORGANELLE = "embplant_pt"
# A plastome run also needs the mitochondrial database: contigs that hit it are what the assembler
# separates from the plastome.
DATABASE_TYPES = (ORGANELLE, "embplant_mt")
ENTRY_POINT = "get_organelle_from_reads.py"

_STATUS = re.compile(rf"Result status of {ORGANELLE}: (.+?)\s*$", re.M)
_VERSION = re.compile(r"^GetOrganelle v(\S+)", re.M)
_DEPENDENCIES = re.compile(r"^DEPENDENCIES: (.+?)\s*$", re.M)
_DATABASE = re.compile(r"^(SEED|LABEL)\s+DB: (.+?)\s*$", re.M)
_COVERAGE = re.compile(rf"Average {ORGANELLE} base-coverage = ([0-9.]+)")
_READS_USED = re.compile(r"Reads used = (\d+)\+(\d+)")
_ERROR = re.compile(r" - ERROR: (.+?)\s*$", re.M)
_ARCHITECTURE = re.compile(
    r"Detecting large repeats \(>1000 bp\) in PATH1 with (.+?), "
    r"Total:LSC:SSC:Repeat\(bp\) = (\d+):(\d+):(\d+):(\d+)"
)
_PATH_FILE = re.compile(r"\.(complete|nearly-complete|scaffolds)\.graph(\d+)\.(.*)path_sequence$")

_STATUS_OF_FILE = {
    "complete": STATUS_COMPLETE,
    "nearly-complete": STATUS_GAPPED,
    "scaffolds": STATUS_FRAGMENTED,
}


@dataclass(frozen=True, slots=True)
class GetOrganelleTools:
    """The assembler and the three programs it drives, with the directory each was found in."""

    assembler: str
    spades: str
    bowtie2: str
    blastn: str


@dataclass(frozen=True, slots=True)
class LogFacts:
    """What the assembler's own log states. Every field may be absent from an imported log."""

    status_text: str | None
    version: str | None
    dependencies: dict[str, str]
    seed_db: dict[str, str]
    label_db: dict[str, str]
    base_coverage: float | None
    reads_used: int | None
    # The assembler's own reading of PATH1: (name, total, LSC, SSC, repeat), in bp. Name is its
    # label for the repeat, such as "IRs detected".
    architecture: tuple[str, int, int, int, int] | None
    errors: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class GetOrganelleOutput:
    """One finished assembly directory, read."""

    status: str
    log: Path
    paths: tuple[AssembledPath, ...]
    path_files: tuple[Path, ...]
    # The graph the paths were read from, then the wider one before it was disentangled, which is
    # what a person opens in Bandage when a path looks wrong.
    graphs: tuple[Path, ...]
    facts: LogFacts


def _pairs(text: str) -> dict[str, str]:
    """``"Bowtie2 2.5.5; SPAdes 4.3.0"`` as ``{"Bowtie2": "2.5.5", "SPAdes": "4.3.0"}``."""
    pairs: dict[str, str] = {}
    for item in text.split(";"):
        parts = item.split()
        if len(parts) == 2:
            pairs[parts[0]] = parts[1]
    return pairs


def parse_log(text: str) -> LogFacts:
    """The facts a log states. The last status line wins: a run reports once per organelle type."""
    status = _STATUS.findall(text)
    version = _VERSION.search(text)
    dependencies = _DEPENDENCIES.search(text)
    databases = {kind: _pairs(value) for kind, value in _DATABASE.findall(text)}
    coverage = _COVERAGE.findall(text)
    reads = _READS_USED.search(text)
    architecture = _ARCHITECTURE.search(text)
    return LogFacts(
        status_text=status[-1] if status else None,
        version=version.group(1) if version else None,
        dependencies=_pairs(dependencies.group(1)) if dependencies else {},
        seed_db=databases.get("SEED", {}),
        label_db=databases.get("LABEL", {}),
        base_coverage=float(coverage[-1]) if coverage else None,
        reads_used=int(reads.group(1)) + int(reads.group(2)) if reads else None,
        architecture=(
            (architecture.group(1), *(int(architecture.group(index)) for index in range(2, 6)))
            if architecture
            else None
        ),
        errors=tuple(_ERROR.findall(text)),
    )


def _status_from_text(text: str | None) -> str:
    if text is None:
        return STATUS_NO_RESULT
    if text == "circular genome":
        return STATUS_COMPLETE
    if text == "circular genome with gaps":
        return STATUS_GAPPED
    if text.endswith("scaffold(s)") or text.endswith("contig(s)"):
        return STATUS_FRAGMENTED
    raise AssemblyError(f"GetOrganelle reported a result this step does not know: {text!r}")


def _read_path_file(path: Path, *, circular: bool) -> AssembledPath:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise AssemblyError(f"{path}: a GetOrganelle path file with no sequence")
    return AssembledPath(
        name=path.name.removesuffix(".path_sequence.fasta"),
        contigs=tuple(str(record.seq).upper() for record in records),
        circular=circular,
    )


def read_output(directory: str | Path) -> GetOrganelleOutput:
    """Read a finished output directory, and refuse one that contradicts itself."""
    root = Path(directory)
    log = root / LOG_NAME
    if not log.is_file():
        raise AssemblyError(f"{root}: no {LOG_NAME} — this is not a GetOrganelle output directory")
    facts = parse_log(log.read_text(encoding="utf-8", errors="replace"))
    status = _status_from_text(facts.status_text)

    files = sorted(root.glob("*.path_sequence.fasta"))
    kinds = set()
    for file in files:
        match = _PATH_FILE.search(file.name.removesuffix(".fasta"))
        if match is None:
            raise AssemblyError(f"{file}: a path file this step cannot name")
        kinds.add(_STATUS_OF_FILE[match.group(1)])
    if status == STATUS_NO_RESULT:
        if files:
            raise AssemblyError(f"{log}: states no result, but {len(files)} path file(s) exist")
    elif kinds != {status}:
        raise AssemblyError(
            f"{root}: the log says {facts.status_text!r} and the path files say "
            f"{sorted(kinds) or 'nothing'} — the result cannot be read from this directory"
        )
    circular = status in (STATUS_COMPLETE, STATUS_GAPPED)
    return GetOrganelleOutput(
        status=status,
        log=log,
        paths=tuple(_read_path_file(file, circular=circular) for file in files),
        path_files=tuple(files),
        graphs=(
            *sorted(root.glob("*.selected_graph.gfa")),
            *sorted(root.glob("extended_*.extend-*.fastg")),
            *sorted(root.glob("extended_*.extend-*.csv")),
        ),
        facts=facts,
    )


def _find(name: str, directory: str | Path | None) -> str:
    if directory is not None:
        candidate = Path(directory) / name
        if not candidate.is_file():
            raise ToolError(f"{name} not found in {directory} — see {GETORGANELLE_DOC}")
        return str(candidate)
    found = which(name)
    if found is None:
        raise ToolError(
            f"{name} not found — install the assembly environment, see {GETORGANELLE_DOC}"
        )
    return found


def resolve_tools(directory: str | Path | None = None) -> GetOrganelleTools:
    """Find the assembler and its three dependencies, and check the two pinned versions."""
    tools = GetOrganelleTools(
        assembler=_find(ENTRY_POINT, directory),
        spades=_find("spades.py", directory),
        bowtie2=_find("bowtie2", directory),
        blastn=_find("blastn", directory),
    )
    line, _ = tool_version("getorganelle", "binary", tools.assembler)
    if line is None or GETORGANELLE_VERSION not in line.split()[-1]:
        raise ToolError(
            f"{tools.assembler}: expected GetOrganelle {GETORGANELLE_VERSION}, got {line!r} — "
            f"see {GETORGANELLE_DOC}"
        )
    line, _ = tool_version("spades", "binary", tools.spades)
    if line is None or SPADES_VERSION not in line:
        raise ToolError(
            f"{tools.spades}: expected SPAdes {SPADES_VERSION}, got {line!r} — "
            f"see {GETORGANELLE_DOC}"
        )
    return tools


def resolve_database(directory: str | Path | None = None) -> dict[str, Any]:
    """The seed and label databases, verified against their own VERSION files.

    The location comes from the argument, then ``GETORG_PATH``, then ``~/.GetOrganelle``, which
    is the order GetOrganelle itself uses. The result names each database by version and digest
    so the manifest records what the assembly was seeded with.
    """
    configured = directory or os.environ.get("GETORG_PATH") or Path.home() / ".GetOrganelle"
    root = Path(configured)
    fetch = f"get_organelle_config.py -a {','.join(DATABASE_TYPES)} --config-dir {root}"
    entries: list[dict[str, str]] = []
    for kind in ("SeedDatabase", "LabelDatabase"):
        version_file = root / kind / "VERSION"
        if not version_file.is_file():
            raise ToolError(
                f"{root}: no {kind}/VERSION — fetch the database once with `{fetch}` "
                f"(see {GETORGANELLE_DOC})"
            )
        declared = {
            parts[0]: parts
            for parts in (line.split() for line in version_file.read_text().splitlines())
            if len(parts) == 3
        }
        for name in DATABASE_TYPES:
            if name not in declared:
                raise ToolError(
                    f"{version_file}: no {name} entry — fetch it with `{fetch}` "
                    f"(see {GETORGANELLE_DOC})"
                )
            _, version, sha256 = declared[name]
            fasta = root / kind / f"{name}.fasta"
            digest = hashlib.sha256(fasta.read_bytes()).hexdigest()
            if digest != sha256:
                raise ToolError(
                    f"{fasta}: sha256 {digest} is not the {sha256} its VERSION file names"
                )
            entries.append({"database": f"{kind}/{name}", "version": version, "sha256": sha256})
    listing = "".join(f"{e['database']}\t{e['version']}\t{e['sha256']}\n" for e in entries)
    return {
        "role": "assembly_seed",
        "format": "getorganelle_database",
        "path": str(root),
        "organelle": ORGANELLE,
        "records": len(entries),
        "databases": entries,
        "digest": hashlib.sha256(listing.encode("utf-8")).hexdigest(),
    }


def assemble(
    reads_1: Path,
    reads_2: Path,
    *,
    out_dir: Path,
    tools: GetOrganelleTools,
    database: str,
    rounds: int,
    kmers: Sequence[int],
    threads: int,
    random_seed: int,
    memory_gb: int,
    ctx: RunContext,
) -> GetOrganelleOutput:
    """Run one assembly into ``out_dir``, which must not exist, and read what it left.

    GetOrganelle catches every failure, logs it and exits 0, so a SPAdes that ran out of memory
    leaves a clean exit and no result. The log is what tells that apart from an assembly that found
    nothing, and the first is this machine's failure, not a property of the sample.
    """
    ctx.tool(
        "getorganelle",
        [
            tools.assembler,
            "-1",
            str(reads_1),
            "-2",
            str(reads_2),
            "-F",
            ORGANELLE,
            "-o",
            str(out_dir),
            "-R",
            str(rounds),
            "-k",
            ",".join(str(kmer) for kmer in kmers),
            "-t",
            str(threads),
            "--random-seed",
            str(random_seed),
            "--config-dir",
            database,
            "--which-spades",
            str(Path(tools.spades).parent),
            "--which-bowtie2",
            str(Path(tools.bowtie2).parent),
            "--which-blast",
            str(Path(tools.blastn).parent),
            "--spades-options",
            f"-m {memory_gb}",
        ],
        tool_path=tools.assembler,
    )
    output = read_output(out_dir)
    if output.status == STATUS_NO_RESULT and output.facts.errors:
        raise ToolError(f"GetOrganelle stopped: {'; '.join(output.facts.errors)}")
    return output
