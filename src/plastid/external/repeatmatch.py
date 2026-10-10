"""MUMmer's ``repeat-match`` — every maximal exact repeat, both strands (M2.12, D-f).

``repeat-match`` reports maximal exact matches and marks the copy that is read on the reverse
strand, which is exactly the two classes the ``repeats`` module publishes: **dispersed** (a
direct repeat) and **palindromic** (an inverted one).

Which MUMmer this binary is — measured 10/10/2026, because PROJECT.md §7.8's line said "MUMmer3":
the program the toolkit runs comes from the **mummer4 4.0.1** bioconda package (``conda-meta``
names it, and it sits beside ``nucmer``). It has **no version flag of its own** — ``--version`` is
an invalid option and prints the usage — so the release is read from its sibling ``mummer``, and a
``repeat-match`` whose sibling disagrees with the pin is refused rather than recorded as a guess.

Coordinates are GenBank-style 1-based inclusive on the way in and 0-based half-open on the way
out, converted here and nowhere else. The ``r`` suffix is the tool's own marking of a match read
on the reverse strand, and it sits on the *start* coordinate of that copy: for ``83128 153724r
26214`` the second copy is the forward interval ``[153724 − 26214, 153724)`` — which is the
golden's IRb copy, the check that this conversion is the right one.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.model import Genome
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "MIN_LENGTH_BP",
    "REPEAT_MATCH_VERSION",
    "ExactMatch",
    "RepeatMatchResult",
    "parse_repeat_match",
    "repeat_match_version",
    "run_repeat_match",
]

# The release the toolkit's numbers were measured with (docs/tools/mummer4.md). repeat-match has
# no version flag, so this pin is checked against the sibling `mummer` in the same installation.
REPEAT_MATCH_VERSION = "4.0.1"

# D-f's floor: every match at least this long is published, and the number is part of the count.
MIN_LENGTH_BP = 30

PROGRAM = "repeat-match"
SIBLING = "mummer"  # the same package, and the only one of the two that answers --version


@dataclass(frozen=True, slots=True)
class ExactMatch:
    """One maximal exact match: two copies of the same sequence, 0-based half-open.

    ``start1 <= start2`` always — the tool's own order is not guaranteed, and a table whose rows
    swap ends between runs is a table nobody can diff. ``palindromic`` is True when either copy is
    read on the reverse strand, i.e. the match is an inverted repeat.
    """

    start1: int
    end1: int
    start2: int
    end2: int
    palindromic: bool

    @property
    def length(self) -> int:
        return self.end1 - self.start1

    def to_1based(self) -> tuple[int, int, int, int]:
        return (self.start1 + 1, self.end1, self.start2 + 1, self.end2)


@dataclass(frozen=True, slots=True)
class RepeatMatchResult:
    """What one run found, plus the parameters it was asked with."""

    matches: tuple[ExactMatch, ...]
    min_length_bp: int
    version: str
    output: str  # the tool's raw stdout, kept for the record


def parse_repeat_match(text: str) -> tuple[ExactMatch, ...]:
    """Read ``repeat-match`` output; a row that will not parse raises with its line number.

    Skimming past a row would drop a repeat nobody notices — the failure mode this project keeps
    meeting — so the only lines skipped are the tool's own headings ("Long Exact Matches:") and
    its column header.
    """
    matches: list[ExactMatch] = []
    for number, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 3 or not fields[2].isdigit() or fields[0][:1].lower() == "s":
            if _is_heading(line):
                continue
            raise ValueError(
                f"repeat-match output line {number}: cannot read {line!r} as a match row"
            )
        start1, end1, reverse1 = _copy(fields[0], int(fields[2]), number)
        start2, end2, reverse2 = _copy(fields[1], int(fields[2]), number)
        # The tool prints one copy first, but not always the lower one; sorting here keeps the
        # table diffable and the two copies distinguishable.
        if (start2, end2) < (start1, end1):
            start1, end1, start2, end2 = start2, end2, start1, end1
        matches.append(
            ExactMatch(
                start1=start1,
                end1=end1,
                start2=start2,
                end2=end2,
                palindromic=reverse1 or reverse2,
            )
        )
    return tuple(matches)


def repeat_match_version(binary: str) -> str:
    """The MUMmer release this ``repeat-match`` belongs to, read from its sibling ``mummer``."""
    sibling = Path(binary).parent / SIBLING
    if not sibling.is_file():
        raise ToolError(
            f"{binary}: no {SIBLING} beside it, so its release cannot be read — {PROGRAM} has no "
            f"version flag of its own. Install the pinned MUMmer package whole "
            f"(docs/tools/mummer4.md)."
        )
    line, _ = tool_version(SIBLING, "binary", str(sibling))
    if line is None or REPEAT_MATCH_VERSION not in line:
        raise ToolError(
            f"{sibling}: expected MUMmer {REPEAT_MATCH_VERSION} (the build this toolkit's repeat "
            f"counts were measured with), got {line!r} — install the pinned build, or re-measure "
            f"and update REPEAT_MATCH_VERSION (docs/tools/mummer4.md)"
        )
    return REPEAT_MATCH_VERSION


def resolve_repeat_match(directory: str | Path | None = None) -> str:
    """The binary to run — from ``--set repeats=<directory>``, or from ``PATH``."""
    if directory is not None:
        candidate = Path(directory) / PROGRAM
        if not candidate.is_file():
            raise ToolError(f"{PROGRAM} not found in {directory} — see docs/tools/mummer4.md")
        return str(candidate)
    found = which(PROGRAM)
    if found is None:
        raise ToolError(
            f"{PROGRAM} not found — install the pinned MUMmer {REPEAT_MATCH_VERSION}, or point "
            f"the run at one with --set repeats=<directory>. See docs/tools/mummer4.md."
        )
    return found


def run_repeat_match(
    genome: Genome,
    ctx: RunContext,
    *,
    min_length: int = MIN_LENGTH_BP,
    directory: str | Path | None = None,
    timeout: int | None = None,
) -> RepeatMatchResult:
    """Find every maximal exact match of at least ``min_length`` bases, both strands."""
    if min_length < 1:
        raise ValueError(f"min_length must be at least 1, got {min_length}")
    binary = resolve_repeat_match(directory)
    version = repeat_match_version(binary)
    name, sequence = genome.accession or genome.name, genome.sequence
    with tempfile.TemporaryDirectory(prefix="repeat_match_") as work:
        workdir = Path(work)
        fasta = workdir / "genome.fa"
        fasta.write_text(
            f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)),
            encoding="utf-8",
        )
        record = ctx.tool(
            PROGRAM,
            [binary, "-n", str(min_length), fasta.name],
            tool_path=binary,
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        output = Path(record.stdout_file).read_text(encoding="utf-8")
    matches = parse_repeat_match(output)
    short = [match for match in matches if match.length < min_length]
    if short:  # the floor is passed to the tool, so a shorter match means the two disagree
        raise ToolError(
            f"{PROGRAM} returned {len(short)} match(es) shorter than the {min_length} bp floor it "
            f"was given (shortest {min(match.length for match in short)} bp) — the count would not "
            f"mean what its label says"
        )
    return RepeatMatchResult(
        matches=matches, min_length_bp=min_length, version=version, output=output
    )


def _copy(field: str, length: int, line_number: int) -> tuple[int, int, bool]:
    """One printed coordinate → ``(start, end, reverse)``, 0-based half-open."""
    reverse = field.endswith(("r", "R"))
    text = field.rstrip("rR")
    if not text.isdigit():
        raise ValueError(f"repeat-match output line {line_number}: {field!r} is not a coordinate")
    value = int(text)
    if reverse:  # the coordinate is the start in the reverse-complement frame
        return value - length, value, True
    return value - 1, value - 1 + length, False


def _is_heading(line: str) -> bool:
    stripped = line.strip()
    return stripped.endswith(":") or stripped.split()[0].lower().startswith("start")
