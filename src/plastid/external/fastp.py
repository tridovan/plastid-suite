"""fastp — quality control of short reads before assembly.

Adapter read-through is the part that matters here: a read carrying adapter sequence looks like a
join to nowhere, in the assembler's graph and in the map-back. Quality trimming is fastp's default
and stays on. The report's numbers go into the run's own tables, read from fastp's JSON rather than
from its text output, because the JSON is the format fastp documents as machine-readable.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["FASTP_VERSION", "FastpResult", "fastp_version", "parse_report", "resolve_fastp", "trim"]

FASTP_VERSION = "1.4.0"
FASTP_DOC = "docs/tools/fastp.md"


@dataclass(frozen=True, slots=True)
class FastpResult:
    """The trimmed pair and the counts fastp reported, in reads and bases of both mates together."""

    reads_1: Path
    reads_2: Path
    report: Path
    reads_before: int
    bases_before: int
    reads_after: int
    bases_after: int
    q30_before: float
    q30_after: float
    adapter_trimmed_reads: int
    read_length_before: int
    read_length_after: int
    version: str


def resolve_fastp(directory: str | Path | None = None) -> str:
    """Find fastp: a given directory, then ``PATH``."""
    if directory is not None:
        candidate = Path(directory) / "fastp"
        if not candidate.is_file():
            raise ToolError(f"fastp not found in {directory} — see {FASTP_DOC}")
        return str(candidate)
    found = which("fastp")
    if found is None:
        raise ToolError(f"fastp not found — install the pinned {FASTP_VERSION}, see {FASTP_DOC}")
    return found


def fastp_version(binary: str) -> str:
    """The release this binary is, checked against the pin."""
    line, _ = tool_version("fastp", "binary", binary)
    if line is None or FASTP_VERSION not in line.split():
        raise ToolError(
            f"{binary}: expected fastp {FASTP_VERSION}, got {line!r} — install the pinned build, "
            f"or re-measure and update FASTP_VERSION ({FASTP_DOC})"
        )
    return FASTP_VERSION


def parse_report(text: str, *, source: str) -> dict[str, int | float]:
    """The counts this step keeps, from fastp's JSON. A report without them is refused by name."""
    try:
        report = json.loads(text)
        before = report["summary"]["before_filtering"]
        after = report["summary"]["after_filtering"]
        return {
            "reads_before": int(before["total_reads"]),
            "bases_before": int(before["total_bases"]),
            "q30_before": float(before["q30_rate"]),
            "read_length_before": int(before["read1_mean_length"]),
            "reads_after": int(after["total_reads"]),
            "bases_after": int(after["total_bases"]),
            "q30_after": float(after["q30_rate"]),
            "read_length_after": int(after["read1_mean_length"]),
            "adapter_trimmed_reads": int(report["adapter_cutting"]["adapter_trimmed_reads"]),
        }
    except (KeyError, TypeError, ValueError) as error:
        raise ToolError(f"{source}: not a fastp JSON report this step can read: {error}") from error


def trim(
    reads_1: Path,
    reads_2: Path,
    *,
    out_1: Path,
    out_2: Path,
    report: Path,
    html: Path,
    threads: int,
    ctx: RunContext,
    binary: str,
) -> FastpResult:
    """Run fastp on one pair. Adapters are detected from the reads: libraries differ in kit."""
    version = fastp_version(binary)
    ctx.tool(
        "fastp",
        [
            binary,
            "--in1",
            str(reads_1),
            "--in2",
            str(reads_2),
            "--out1",
            str(out_1),
            "--out2",
            str(out_2),
            "--detect_adapter_for_pe",
            "--thread",
            str(threads),
            "--json",
            str(report),
            "--html",
            str(html),
        ],
        tool_path=binary,
        pinned_version=version,
    )
    counts = parse_report(report.read_text(encoding="utf-8"), source=str(report))
    return FastpResult(reads_1=out_1, reads_2=out_2, report=report, version=version, **counts)
