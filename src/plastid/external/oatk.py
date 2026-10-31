"""Oatk — de novo plastome assembly from PacBio HiFi reads.

**Why this assembler.** It is built for plant organelles from HiFi reads, and it is the assembler of
the 195 species in Zhou et al. 2025 (*Genome Biology* 26:235). ``docs/tools/oatk.md`` has what was
measured here against an answer from Sanger reads.

**What it writes.** One contig for a closed plastome, whose header states its walk through the graph
(``path=u6+,u5-,u4+,u5+``). That is one isomer of the two a repeat gives, so the other is spelled
here from the same graph, and both backends and both assemblers hand on the same pair. A graph with
more structures than that has no pair to complete, and the graph reading reports it.

**How it fails.** Unlike GetOrganelle it exits non-zero, and its reason is an ``[E::...]`` line at
the end of stderr, past where a tool error's excerpt stops. So the exit is allowed here and the
reason is read from the log.

**What it is given.** ``-D`` bounds the reads it takes, which bounds its memory: its reads are held
in a table that grows with them, and a plastome is reached long before a sample's reads run out.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from Bio import SeqIO

from plastid.core.assembly import (
    STATUS_COMPLETE,
    STATUS_FRAGMENTED,
    STATUS_NO_RESULT,
    AssembledPath,
    AssemblyError,
)
from plastid.core.assembly_graph import (
    AssemblerOutput,
    LogFacts,
    other_isomer,
    parse_gfa,
    parse_walk,
    read_selected_graph,
    spell,
)
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "DATABASE_RELEASE",
    "HMMER_VERSION",
    "LOG_NAME",
    "OATK_VERSION",
    "OatkTools",
    "assemble",
    "hmmer_version",
    "parse_log",
    "read_output",
    "resolve_database",
    "resolve_tools",
]

OATK_VERSION = "1.0"
HMMER_VERSION = "3.4"
OATK_DOC = "docs/tools/oatk.md"
LOG_NAME = "oatk.log"
ENTRY_POINT = "oatk"
# The release of OatkDB the acceptance numbers were measured with, and the lineage its README
# recommends for land plants. The mitochondrial profiles are there so a contig of a mitochondrion
# carrying plastid sequence is classed as what it is.
DATABASE_RELEASE = "v20230921"
PLASTID_PROFILES = "embryophyta_pltd.fam"
MITOCHONDRION_PROFILES = "embryophyta_mito.fam"
_PRESSED = (".h3f", ".h3i", ".h3m", ".h3p")

_VERSION = re.compile(r"^\[M::main\] Version: (\S+)", re.M)
_READS = re.compile(r"collected syncmers from (\d+) target sequence")
_ERROR = re.compile(r"^\[E::[^\]]+\] (.+?)\s*$", re.M)
_HMMER = re.compile(r"HMMER (\S+)")
_FIELD = re.compile(r"\b(length|circular|path)=(\S+)")


@dataclass(frozen=True, slots=True)
class OatkTools:
    """The assembler and the HMMER program it calls."""

    assembler: str
    nhmmscan: str


def _find(name: str, directory: str | Path | None) -> str:
    if directory is not None:
        candidate = Path(directory) / name
        if not candidate.is_file():
            raise ToolError(f"{name} not found in {directory} — see {OATK_DOC}")
        return str(candidate)
    found = which(name)
    if found is None:
        raise ToolError(f"{name} not found — install the assembly environment, see {OATK_DOC}")
    return found


def hmmer_version(binary: str) -> str:
    """The release this nhmmscan is. HMMER has no version flag, so its help header is read."""
    result = subprocess.run([binary, "-h"], capture_output=True, text=True, timeout=10, check=False)
    match = _HMMER.search(result.stdout)
    return match.group(1) if match else "not stated"


def resolve_tools(directory: str | Path | None = None) -> OatkTools:
    """Find Oatk and nhmmscan, and check the two pinned versions."""
    assembler = _find(ENTRY_POINT, directory)
    nhmmscan = _find("nhmmscan", directory)
    line, _ = tool_version("oatk", "binary", assembler)
    if line is None or line.split()[-1] != OATK_VERSION:
        raise ToolError(f"{assembler}: expected Oatk {OATK_VERSION}, got {line!r} — see {OATK_DOC}")
    version = hmmer_version(nhmmscan)
    if version != HMMER_VERSION:
        raise ToolError(
            f"{nhmmscan}: expected HMMER {HMMER_VERSION}, got {version!r} — see {OATK_DOC}"
        )
    return OatkTools(assembler=assembler, nhmmscan=nhmmscan)


def resolve_database(directory: str | Path | None = None) -> dict[str, Any]:
    """The two profile databases, by digest, so the manifest records what classified the contigs.

    The location comes from the argument, then ``OATK_DB``. OatkDB states no checksums of its own,
    so the digests are this run's reading of the files; ``docs/tools/oatk.md`` gives those measured.
    """
    configured = directory or os.environ.get("OATK_DB")
    if not configured:
        raise ToolError(
            f"no Oatk database: pass --set oatk_db=<dir> or set OATK_DB (see {OATK_DOC})"
        )
    root = Path(configured)
    entries: list[dict[str, str]] = []
    for name in (PLASTID_PROFILES, MITOCHONDRION_PROFILES):
        profile = root / name
        missing = [file for file in (name, *(name + ext for ext in _PRESSED))]
        missing = [file for file in missing if not (root / file).is_file()]
        if missing:
            raise ToolError(f"{root}: no {', '.join(missing)} — see {OATK_DOC}")
        digest = hashlib.sha256(profile.read_bytes()).hexdigest()
        entries.append({"database": name, "version": DATABASE_RELEASE, "sha256": digest})
    listing = "".join(f"{e['database']}\t{e['version']}\t{e['sha256']}\n" for e in entries)
    return {
        "role": "assembly_seed",
        "format": "oatk_database",
        "path": str(root),
        "organelle": "plastid",
        "records": len(entries),
        "databases": entries,
        "digest": hashlib.sha256(listing.encode("utf-8")).hexdigest(),
    }


def parse_log(text: str, *, status_text: str | None) -> LogFacts:
    """What a run's stderr states. Oatk's log carries no result, so that comes from its contigs."""
    version = _VERSION.search(text)
    reads = _READS.search(text)
    return LogFacts(
        status_text=status_text,
        version=version.group(1) if version else None,
        dependencies={},
        seed_db={},
        label_db={},
        base_coverage=None,
        reads_used=int(reads.group(1)) if reads else None,
        architecture=None,
        errors=tuple(_ERROR.findall(text)),
    )


def _fields(description: str) -> dict[str, str]:
    return dict(_FIELD.findall(description))


def _read_paths(contigs: Path) -> list[AssembledPath]:
    """A closed contig as a path whose title states its walk, the way GetOrganelle's does."""
    records = list(SeqIO.parse(contigs, "fasta"))
    if not records:
        return []
    fields = _fields(records[0].description)
    if len(records) == 1 and fields.get("circular") == "true":
        return [
            AssembledPath(
                name=records[0].id,
                contigs=(str(records[0].seq).upper(),),
                circular=True,
                title=f"{fields.get('path', '')}(circular)",
            )
        ]
    return [
        AssembledPath(
            name="contigs",
            contigs=tuple(str(record.seq).upper() for record in records),
            circular=False,
            title=records[0].description,
        )
    ]


def read_output(directory: str | Path) -> AssemblerOutput:
    """Read a finished output directory, and refuse one that contradicts itself."""
    root = Path(directory)
    log = root / LOG_NAME
    if not log.is_file():
        raise AssemblyError(
            f"{root}: no {LOG_NAME} — keep Oatk's stderr under that name, it states the version"
        )
    text = log.read_text(encoding="utf-8", errors="replace")
    found = sorted(root.glob("*.pltd.ctg.fasta"))
    if len(found) > 1:
        raise AssemblyError(f"{root}: {len(found)} *.pltd.ctg.fasta — one run per directory")
    paths = _read_paths(found[0]) if found else []
    if not paths:
        return AssemblerOutput(
            status=STATUS_NO_RESULT,
            log=log,
            paths=(),
            path_files=(),
            graphs=(),
            facts=parse_log(text, status_text=None),
        )
    closed = paths[0].circular
    status = STATUS_COMPLETE if closed else STATUS_FRAGMENTED
    status_text = "circular genome" if closed else f"{len(paths[0].contigs)} contig(s)"
    selected = sorted(root.glob("*.pltd.gfa"))
    readings = ()
    if selected:
        graph = parse_gfa(selected[0].read_text(encoding="utf-8").splitlines())
        if closed:
            walk = parse_walk(paths[0].title)
            other = other_isomer(graph, walk) if walk is not None else None
            if other is not None:
                paths.append(
                    AssembledPath(
                        name=f"{paths[0].name}.isomer",
                        contigs=(spell(graph, other),),
                        circular=True,
                        title=",".join(f"{name}{strand}" for name, strand in other) + "(circular)",
                    )
                )
        try:
            readings = (read_selected_graph("pltd", graph, paths, status=status),)
        except AssemblyError as error:
            raise AssemblyError(f"{selected[0]}: {error}") from error
    return AssemblerOutput(
        status=status,
        log=log,
        paths=tuple(paths),
        path_files=tuple(found),
        graphs=(
            *selected,
            *sorted(root.glob("*.utg.final.gfa")),
            *sorted(root.glob("*.pltd.ctg.bed")),
        ),
        facts=parse_log(text, status_text=status_text),
        graph_readings=readings,
    )


def assemble(
    reads: Path,
    *,
    out_dir: Path,
    tools: OatkTools,
    database: str,
    syncmer: int,
    min_coverage: int,
    max_bases: int,
    threads: int,
    ctx: RunContext,
) -> AssemblerOutput:
    """Run one assembly into ``out_dir``, which must not exist, and read what it left."""
    out_dir.mkdir(parents=True)
    root = Path(database)
    record = ctx.tool(
        ENTRY_POINT,
        [
            tools.assembler,
            "-k",
            str(syncmer),
            "-c",
            str(min_coverage),
            "-D",
            str(max_bases),
            "-t",
            str(threads),
            "--nhmmscan",
            tools.nhmmscan,
            "-m",
            str(root / MITOCHONDRION_PROFILES),
            "-p",
            str(root / PLASTID_PROFILES),
            "-o",
            str(out_dir / "oatk"),
            str(reads),
        ],
        tool_path=tools.assembler,
        pinned_version=OATK_VERSION,
        expect=(0, 1),
    )
    text = Path(record.stderr_file).read_text(encoding="utf-8", errors="replace")
    (out_dir / LOG_NAME).write_text(text, encoding="utf-8")
    if record.rc != 0:
        reasons = "; ".join(_ERROR.findall(text)) or f"exit {record.rc}, and no reason in its log"
        raise ToolError(f"Oatk stopped: {reasons}")
    output = read_output(out_dir)
    stated = replace(
        output.facts,
        dependencies={"HMMER": HMMER_VERSION},
        seed_db={"OatkDB": DATABASE_RELEASE},
    )
    return replace(output, facts=stated)
