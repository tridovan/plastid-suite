"""minimap2 — mapping the reads back onto the assembled sequence.

The aligner is deliberately not the assembler's: GetOrganelle recruits reads with Bowtie2, so a
check built on the same tool and index would share its blind spots. The options each carry a reason:

* ``-ax sr`` is the short-read preset, the same for every sample;
* ``--eqx`` writes =/X operators, so mismatches are counted from the CIGAR with no second lookup;
* ``--secondary=no`` keeps one placement per read. A read inside an inverted repeat has two equal
  placements, and the primary falls on either copy, which spreads it over both;
* ``--sam-hit-only`` leaves unmapped reads out. Most reads of a genome skim are not plastid, and
  writing them would make the SAM as large as the input.

A circular sequence has no end, so the mapping target is the sequence followed by its own first
``extension`` bases (:mod:`plastid.core.assembly` folds the positions back).
"""

from __future__ import annotations

from pathlib import Path

from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["MINIMAP2_VERSION", "map_pair", "minimap2_version", "resolve_minimap2"]

MINIMAP2_VERSION = "2.31"
MINIMAP2_DOC = "docs/tools/minimap2.md"


def resolve_minimap2(directory: str | Path | None = None) -> str:
    """Find minimap2: a given directory, then ``PATH``."""
    if directory is not None:
        candidate = Path(directory) / "minimap2"
        if not candidate.is_file():
            raise ToolError(f"minimap2 not found in {directory} — see {MINIMAP2_DOC}")
        return str(candidate)
    found = which("minimap2")
    if found is None:
        raise ToolError(
            f"minimap2 not found — install the pinned {MINIMAP2_VERSION}, see {MINIMAP2_DOC}"
        )
    return found


def minimap2_version(binary: str) -> str:
    """The release this binary is, checked against the pin. It prints e.g. ``2.31-r1302``."""
    line, _ = tool_version("minimap2", "binary", binary)
    if line is None or line.split("-")[0] != MINIMAP2_VERSION:
        raise ToolError(
            f"{binary}: expected minimap2 {MINIMAP2_VERSION}, got {line!r} — install the pinned "
            f"build, or re-measure and update MINIMAP2_VERSION ({MINIMAP2_DOC})"
        )
    return MINIMAP2_VERSION


def map_pair(
    target: Path,
    reads_1: Path,
    reads_2: Path,
    *,
    out: Path,
    threads: int,
    ctx: RunContext,
    binary: str,
) -> Path:
    """Map one read pair to ``target`` and write the SAM to ``out``."""
    version = minimap2_version(binary)
    ctx.tool(
        "minimap2",
        [
            binary,
            "-ax",
            "sr",
            "--eqx",
            "--secondary=no",
            "--sam-hit-only",
            "-t",
            str(threads),
            "-o",
            str(out),
            str(target),
            str(reads_1),
            str(reads_2),
        ],
        tool_path=binary,
        pinned_version=version,
    )
    return out
