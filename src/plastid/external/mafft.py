"""MAFFT — the aligner behind the comparison modules (tier B; `environment.yml` pins 7.526).

One job: sequences in, an :class:`Alignment` out, through the single tool door so the version
and the command line land in the run's provenance. A missing aligner raises a named ToolError —
the comparison modules cannot quietly fall back to a hand-rolled alignment (measuring that
fallback against MAFFT once gave 2.667 where the true Ts/Tv was 2.686).
"""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from pathlib import Path

from plastid.core.align import Alignment
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "DEFAULT_STRATEGY",
    "MAFFT_STRATEGIES",
    "MAFFT_VERSION",
    "adjust_direction",
    "align",
    "parse_directions",
    "strategy_flag",
]

MAFFT_VERSION = "7.526"

# The aligner's own entry points, by the name this toolkit gives each. `auto` is MAFFT's default
# and the one every published number here was measured with; the rest are offered because a
# divergence-rich panel can need a different algorithm, and a user who runs MAFFT by hand picks
# one. Adding a strategy changes the alignment, so it travels in argv and in the column.
MAFFT_STRATEGIES: dict[str, str] = {
    "auto": "--auto",
    "einsi": "--einsi",
    "linsi": "--linsi",
    "ginsi": "--ginsi",
    "fftns": "--fftns",
    "nwns": "--nwns",
}
DEFAULT_STRATEGY = "auto"


def strategy_flag(strategy: str) -> str:
    """The MAFFT flag for a strategy name, refused by name when it is not one of ours."""
    try:
        return MAFFT_STRATEGIES[strategy]
    except KeyError:
        raise ValueError(
            f"mafft: unknown strategy {strategy!r} — choose one of "
            f"{', '.join(sorted(MAFFT_STRATEGIES))}"
        ) from None


def align(
    rows: Sequence[tuple[str, str]], ctx: RunContext, *, strategy: str = DEFAULT_STRATEGY
) -> Alignment:
    """Align ``(name, sequence)`` rows and return them in the order given.

    ``strategy`` selects MAFFT's algorithm (M11.20). The default is ``auto``, which is what every
    number in this repository was measured with; anything else is a different alignment and must
    say so — the module records it in the table's scope and the argv in the manifest.
    """
    text = _run(rows, ctx, extra=(), strategy=strategy)
    aligned = _parse_fasta(text)
    missing = [name for name, _ in rows if name not in aligned]
    if missing:
        raise ToolError(f"mafft returned no row for {', '.join(missing)}")
    return Alignment(
        names=tuple(name for name, _ in rows),
        sequences=tuple(aligned[name] for name, _ in rows),
    )


def adjust_direction(rows: Sequence[tuple[str, str]], ctx: RunContext) -> dict[str, bool]:
    """Which rows MAFFT read as reverse-oriented — name → True when it flipped the sequence.

    ``--adjustdirection`` makes MAFFT compare every sequence's direction against **the first row**
    and reverse-complement the ones that came out the other way, marking each with an ``_R_``
    prefix on its output name (measured 08/10/2026 on 7.526, the corpus's SSC region: the two
    records the audit's gene-strand rule calls inverted are exactly the two that come back
    ``_R_``). It answers the one question a rule over a single record cannot settle alone — is
    this region oriented the way its relatives are? — and the caller decides what the first row
    means.
    """
    text = _run(rows, ctx, extra=("--adjustdirection",))
    directions = parse_directions(text)
    missing = [name for name, _ in rows if name not in directions]
    if missing:
        raise ToolError(f"mafft returned no row for {', '.join(missing)}")
    return directions


def parse_directions(text: str) -> dict[str, bool]:
    """MAFFT's answer, read off its headers: ``_R_`` prefix means it flipped that sequence.

    Split out from the call because it is the part with a rule in it, and a rule deserves a test
    that does not need the aligner installed: the marker is a prefix on the **name**, and the name
    is what BLAST and every other tool would truncate at the first whitespace.
    """
    directions: dict[str, bool] = {}
    for line in text.splitlines():
        if not line.startswith(">"):
            continue
        name = line[1:].strip().split(" ")[0]
        reversed_name = name[3:] if name.startswith("_R_") else name
        directions[reversed_name] = name.startswith("_R_")
    return directions


def _run(
    rows: Sequence[tuple[str, str]],
    ctx: RunContext,
    *,
    extra: Sequence[str],
    strategy: str = DEFAULT_STRATEGY,
) -> str:
    """One MAFFT call, through the tool door: a temp file in, its stdout out.

    Both callers need the same three things — the binary located and version-checked, the aligned
    fasta written somewhere we can name, and every argument recorded — so the call is written
    once here rather than twice (CLAUDE.md: DRY after the second real copy).
    """
    if len(rows) < 2:
        raise ValueError(f"an alignment needs at least two sequences, got {len(rows)}")
    binary = which("mafft")
    if binary is None:
        raise ToolError(
            f"mafft not found — the comparison modules need MAFFT {MAFFT_VERSION}; install it"
        )
    version, _ = tool_version("mafft", "mafft", binary)
    with tempfile.TemporaryDirectory(prefix="mafft_") as work:
        workdir = Path(work)
        fasta = workdir / "input.fasta"
        fasta.write_text(_to_fasta(rows), encoding="utf-8")
        record = ctx.tool(
            "mafft",
            [binary, strategy_flag(strategy), "--quiet", *extra, fasta.name],
            tool_path=binary,
            pinned_version=MAFFT_VERSION,
            cwd=workdir,
        )
        return Path(record.stdout_file).read_text(encoding="utf-8")


def _to_fasta(rows: Sequence[tuple[str, str]], width: int = 60) -> str:
    blocks = []
    for name, sequence in rows:
        body = "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))
        blocks.append(f">{name}\n{body}\n")
    return "".join(blocks)


def _parse_fasta(text: str) -> dict[str, str]:
    """MAFFT's output as our rows — **uppercased**, because that is what the counters read.

    Measured 03/10/2026 on MAFFT 7.526: it prints the alignment in **lowercase**
    (``----------------------------------cccgatatcttactaaaaaaagatat``). The comparison modules
    classify a column by ``core.align._DNA``, which is ``ACGT``, so every column of a lowercase
    alignment was silently skipped: the first real-data Ts/Tv run reported *zero comparable
    columns and a blank ratio* rather than a wrong number. The unit tests never saw it because
    they build alignments by hand, in uppercase.
    """
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    for line in text.splitlines():
        if line.startswith(">"):
            current = line[1:].strip().split(" ")[0]
            sequences[current] = []
        elif current is not None:
            sequences[current].append(line.strip())
    return {name: "".join(parts).upper() for name, parts in sequences.items()}
