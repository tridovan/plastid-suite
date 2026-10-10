"""Vmatch 2.3.1 — the independent second engine for the exact-repeat count (M2.12, D-f).

**Not a module's engine.** ``repeats`` runs on ``repeat-match``, which carries no licence
constraint, so the minimum install keeps working without this. Vmatch is **academic-only**
(PROJECT.md §7.8), and it is here for the check this project publishes a repeat count with: two
independent tools agreeing, the way MISA 2.1 and cpstools 3.0 both give the golden's 37 SSRs. A
toolkit that needs Vmatch cannot ship to everyone; a count that two engines agree on can.

The pipeline is the tool's own two steps, and both are recorded in the provenance:

    mkvtree -db <fasta> -dna -pl -allout -indexname <idx>
    vmatch -l <floor> -d -p -nodist <idx>

``-d`` (direct) and ``-p`` (palindromic) are the two classes the module publishes; ``-nodist``
prints one line per match instead of an alignment dump. **Output coordinates are 0-based starts**
— unlike ``repeat-match``'s 1-based ones — and that conversion lives here and nowhere else.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.model import Genome
from plastid.external.repeatmatch import ExactMatch
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["VMATCH_VERSION", "VmatchResult", "parse_vmatch", "run_vmatch"]

# The pinned release (§7.8, environment.yml). A different Vmatch can define the classes
# differently, and the whole point of this adapter is that a *count* can be checked.
VMATCH_VERSION = "2.3.1"

PIPELINE_TOOLS: tuple[str, ...] = ("mkvtree", "vmatch")


@dataclass(frozen=True, slots=True)
class VmatchResult:
    """What one run found: the matches, in the same shape ``repeat-match``'s parser returns."""

    matches: tuple[ExactMatch, ...]
    min_length_bp: int
    version: str
    output: str


def parse_vmatch(text: str) -> tuple[ExactMatch, ...]:
    """Read ``vmatch -nodist`` output; a row that will not parse raises with its line number.

    The columns are: length, a strand offset, start, class, length, offset, start, E-value,
    score, % identity. Only the length, the two starts and the class are read — the rest belongs
    to Vmatch's own statistics, and copying them into this toolkit's tables would give two
    numbers nobody could reconcile with ``repeat-match``'s.
    """
    matches: list[ExactMatch] = []
    for number, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields or line.startswith("#"):
            continue
        if len(fields) < 7:
            raise ValueError(f"vmatch output line {number}: cannot read {line!r} as a match row")
        try:
            length = int(fields[0])
            start1 = int(fields[2])
            klass = fields[3]
            start2 = int(fields[6])
        except ValueError as error:
            raise ValueError(f"vmatch output line {number}: {line!r} ({error})") from None
        if klass not in ("D", "P"):
            raise ValueError(
                f"vmatch output line {number}: unknown class {klass!r} — expected D (direct) or "
                f"P (palindromic)"
            )
        first = (start1, start1 + length)
        second = (start2, start2 + length)
        if second < first:
            first, second = second, first
        matches.append(
            ExactMatch(
                start1=first[0],
                end1=first[1],
                start2=second[0],
                end2=second[1],
                palindromic=klass == "P",
            )
        )
    return tuple(matches)


def run_vmatch(
    genome: Genome,
    ctx: RunContext,
    *,
    min_length: int,
    directory: str | Path | None = None,
    timeout: int | None = None,
) -> VmatchResult:
    """Index one genome and list its direct and palindromic exact matches of ≥ ``min_length``."""
    binaries = _resolve(directory)
    version = _version(binaries["vmatch"])
    name, sequence = genome.accession or genome.name, genome.sequence
    with tempfile.TemporaryDirectory(prefix="vmatch_") as work:
        workdir = Path(work)
        fasta = workdir / "genome.fa"
        fasta.write_text(
            f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)),
            encoding="utf-8",
        )
        ctx.tool(
            "mkvtree",
            [
                binaries["mkvtree"],
                "-db",
                fasta.name,
                "-dna",
                "-pl",
                "-allout",
                "-indexname",
                "idx",
            ],
            tool_path=binaries["mkvtree"],
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        record = ctx.tool(
            "vmatch",
            [binaries["vmatch"], "-l", str(min_length), "-d", "-p", "-nodist", "idx"],
            tool_path=binaries["vmatch"],
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        output = Path(record.stdout_file).read_text(encoding="utf-8")
    return VmatchResult(
        matches=parse_vmatch(output), min_length_bp=min_length, version=version, output=output
    )


def _resolve(directory: str | Path | None) -> dict[str, str]:
    """Both binaries, from ``--tools``/the argument, or from ``PATH``."""
    resolved: dict[str, str] = {}
    for tool in PIPELINE_TOOLS:
        if directory is not None:
            candidate = Path(directory) / tool
            if not candidate.is_file():
                raise ToolError(f"{tool} not found in {directory} — see docs/tools/vmatch.md")
            resolved[tool] = str(candidate)
            continue
        found = which(tool)
        if found is None:
            raise ToolError(
                f"{tool} not found — install the pinned Vmatch {VMATCH_VERSION} (academic "
                f"licence), or point the run at one with --tools <directory>. See "
                f"docs/tools/vmatch.md."
            )
        resolved[tool] = found
    prefixes = {str(Path(path).parent) for path in resolved.values()}
    if len(prefixes) > 1:
        raise ToolError(
            f"mkvtree and vmatch come from different installations: {', '.join(sorted(prefixes))} "
            f"— they share an index format and must ship together (docs/tools/vmatch.md)"
        )
    return resolved


def _version(binary: str) -> str:
    line, _ = tool_version(Path(binary).name, "binary", binary)
    if line is None or VMATCH_VERSION not in line:
        raise ToolError(
            f"{binary}: expected Vmatch {VMATCH_VERSION} (the build this toolkit's cross-check "
            f"was measured with), got {line!r} — install the pinned build, or re-measure and "
            f"update VMATCH_VERSION (docs/tools/vmatch.md)"
        )
    return VMATCH_VERSION
