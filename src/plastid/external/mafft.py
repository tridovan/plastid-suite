"""MAFFT — the aligner behind the comparison modules (tier B; PROJECT.md §7.8 pins 7.526).

One job: sequences in, an :class:`Alignment` out, through the single tool door so the version
and the command line land in the run's provenance. A missing aligner raises a named ToolError —
the comparison modules cannot quietly fall back to a hand-rolled alignment (measuring that
fallback against MAFFT once gave 2.667 where the true Ts/Tv was 2.686; PROJECT.md §9.B).
"""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from pathlib import Path

from plastid.core.align import Alignment
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = ["MAFFT_VERSION", "align"]

MAFFT_VERSION = "7.526"


def align(rows: Sequence[tuple[str, str]], ctx: RunContext) -> Alignment:
    """Align ``(name, sequence)`` rows and return them in the order given."""
    if len(rows) < 2:
        raise ValueError(f"an alignment needs at least two sequences, got {len(rows)}")
    binary = which("mafft")
    if binary is None:
        raise ToolError(
            f"mafft not found — the comparison modules need MAFFT {MAFFT_VERSION}; "
            f"install it as recorded in PROJECT.md §7.8"
        )
    version, _ = tool_version("mafft", "mafft", binary)
    with tempfile.TemporaryDirectory(prefix="mafft_") as work:
        workdir = Path(work)
        fasta = workdir / "input.fasta"
        fasta.write_text(_to_fasta(rows), encoding="utf-8")
        record = ctx.tool(
            "mafft",
            [binary, "--auto", "--quiet", fasta.name],
            tool_path=binary,
            pinned_version=MAFFT_VERSION,
            cwd=workdir,
        )
        text = Path(record.stdout_file).read_text(encoding="utf-8")
    aligned = _parse_fasta(text)
    missing = [name for name, _ in rows if name not in aligned]
    if missing:
        raise ToolError(f"mafft returned no row for {', '.join(missing)}")
    return Alignment(
        names=tuple(name for name, _ in rows),
        sequences=tuple(aligned[name] for name, _ in rows),
    )


def _to_fasta(rows: Sequence[tuple[str, str]], width: int = 60) -> str:
    blocks = []
    for name, sequence in rows:
        body = "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))
        blocks.append(f">{name}\n{body}\n")
    return "".join(blocks)


def _parse_fasta(text: str) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    for line in text.splitlines():
        if line.startswith(">"):
            current = line[1:].strip().split(" ")[0]
            sequences[current] = []
        elif current is not None:
            sequences[current].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}
