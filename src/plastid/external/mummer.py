"""MUMmer 4 — the whole-genome alignment behind the synteny module (P4-01; decision D4-02).

The pipeline is the one PROJECT.md §7 records for the %identity-by-coordinate plot, and the
module's numbers are only comparable to that record while every step of it is the same:

    nucmer --maxmatch -l 20 -c 65 -t 1 -p out <target.fa> <query.fa>
    delta-filter -1 out.delta > out.filter.delta
    show-coords -rcl out.filter.delta

Three choices in here are deliberate:

* **the parameters are passed, not defaulted.** ``-l 20 -c 65`` are MUMmer 4's own defaults, but
  a default is a number nobody wrote down — and the phase's check is a block count that moves
  with them (measured 04/10/2026: the golden pair gives 12 raw alignments at the defaults, 22 at
  ``-c 40``). They travel in the argv the manifest records, and they live in
  ``data/defaults.yaml`` (``synteny:``) rather than as constants here, so there is one place to
  change and the adapter cannot quietly disagree with the module.
* **``-t 1``.** Threads must not be able to change an output that a manuscript number rests on.
* **``delta-filter -1``** (1-to-1). The alternative modes answer different questions — ``-q``
  maps every query position, which can put one target region in several rows — and a block table
  that double-counts an inverted repeat is a coverage number that lies.

The delta file is read twice over: ``show-coords`` gives the block table with the tool's own
%IDY, and :mod:`plastid.core.delta` gives the alignment itself for the per-window series. The
two must agree on the block set — :func:`run_nucmer` refuses a run where they do not, because a
series drawn from a different alignment than the table beside it is a figure nobody can check.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

from plastid.core.delta import DeltaBlock, parse_delta
from plastid.run.context import RunContext
from plastid.run.tools import ToolError, tool_version, which

__all__ = [
    "MUMMER_VERSION",
    "PIPELINE_TOOLS",
    "BlockRow",
    "MummerResult",
    "parse_coords",
    "run_nucmer",
]

# The release the toolkit's numbers were measured with (docs/tools/mummer4.md). A different
# MUMmer can cluster alignments differently, and a block count is the whole check, so a mismatch
# is refused by name rather than recorded as a footnote.
MUMMER_VERSION = "4.0.1"

PIPELINE_TOOLS: tuple[str, ...] = ("nucmer", "delta-filter", "show-coords")


@dataclass(frozen=True, slots=True)
class BlockRow:
    """One alignment as ``show-coords -rcl`` prints it — the tool's own numbers, unconverted."""

    target_name: str
    query_name: str
    target_start: int  # 1-based inclusive, the tool's own frame
    target_end: int
    query_start: int
    query_end: int
    target_length: int
    query_length: int
    length_target: int
    length_query: int
    identity: float  # %IDY: 100 · (1 − errors / alignment length), indels counted as errors
    reverse: bool


@dataclass(frozen=True, slots=True)
class MummerResult:
    """Everything the three tools produced: the block table, the alignment, and what ran."""

    blocks: tuple[BlockRow, ...]
    delta_blocks: tuple[DeltaBlock, ...]
    coords: str  # the raw show-coords text (kept as a run artifact)
    delta: str  # the raw filtered delta (kept as a run artifact)
    version: str

    @property
    def query_name(self) -> str:
        """The query this run aligned — a run aligns exactly one pair."""
        if not self.delta_blocks:
            raise ValueError("this run produced no alignments — it has no query name")
        return self.delta_blocks[0].query_name


def resolve_binary(tool: str, *, directory: str | Path | None = None) -> str:
    """Find one MUMmer binary: a given directory, then ``PATH``."""
    if directory is not None:
        candidate = Path(directory) / tool
        if not candidate.is_file():
            raise ToolError(f"{tool} not found in {directory} — see docs/tools/mummer4.md")
        return str(candidate)
    found = which(tool)
    if found is None:
        raise ToolError(
            f"{tool} not found — install the pinned MUMmer {MUMMER_VERSION}, or point the run "
            f"at one with --set mummer=<directory>. See docs/tools/mummer4.md."
        )
    return found


def nucmer_version(binary: str) -> str:
    """The MUMmer release this binary is, checked against the pin — refused, not noted."""
    line, _ = tool_version(Path(binary).name, "binary", binary)
    if line is None or MUMMER_VERSION not in line:
        raise ToolError(
            f"{binary}: expected MUMmer {MUMMER_VERSION} (the build this toolkit's numbers were "
            f"measured with), got {line!r} — install the pinned build, or re-measure and update "
            f"MUMMER_VERSION (docs/tools/mummer4.md)"
        )
    return MUMMER_VERSION


def _resolve_pipeline(directory: str | Path | None) -> dict[str, str]:
    """All three binaries, checked to come from one installation.

    They ship as one package and share a file format; a ``nucmer`` from one build feeding a
    ``delta-filter`` from another would be a silent mismatch nobody could name afterwards.
    """
    binaries = {tool: resolve_binary(tool, directory=directory) for tool in PIPELINE_TOOLS}
    prefixes = {str(Path(path).parent) for path in binaries.values()}
    if len(prefixes) > 1:
        raise ToolError(
            f"the MUMmer tools come from different installations: "
            f"{', '.join(sorted(prefixes))} — install them together, or point --set mummer= at "
            f"one directory (docs/tools/mummer4.md)"
        )
    return binaries


def parse_coords(text: str) -> tuple[BlockRow, ...]:
    """Read ``show-coords -rcl`` output; anything unreadable raises with its line number."""
    rows: list[BlockRow] = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        # The header is two path lines, a ruler, the column names and another ruler; everything
        # else with a bar in it is a row, and a row that will not parse is an error — not
        # something to skip quietly past (a skipped row is a missing block nobody notices).
        if "|" not in line or stripped.startswith("[") or set(stripped) <= {"="}:
            continue
        first = line.split("|", 1)[0].split()
        if not first or not first[0].isdigit():
            raise ValueError(f"show-coords line {number}: cannot read {line!r} as an alignment row")
        fields = [field.strip() for field in line.split("|")]
        if len(fields) < 7:
            continue
        try:
            target_start, target_end = (int(value) for value in fields[0].split())
            query_start, query_end = (int(value) for value in fields[1].split())
            length_target, length_query = (int(value) for value in fields[2].split())
            identity = float(fields[3])
            target_length, query_length = (int(value) for value in fields[4].split())
            tags = fields[6].split()
            target_name, query_name = tags[0], tags[1]
        except (ValueError, IndexError) as error:
            raise ValueError(f"show-coords line {number}: cannot read {line!r} ({error})") from None
        rows.append(
            BlockRow(
                target_name=target_name,
                query_name=query_name,
                target_start=target_start,
                target_end=target_end,
                query_start=query_start,
                query_end=query_end,
                target_length=target_length,
                query_length=query_length,
                length_target=length_target,
                length_query=length_query,
                identity=identity,
                reverse=query_start > query_end,
            )
        )
    return tuple(rows)


def run_nucmer(
    target: tuple[str, str],
    query: tuple[str, str],
    ctx: RunContext,
    *,
    min_match: int,
    min_cluster: int,
    filter_mode: str,
    directory: str | Path | None = None,
    timeout: int | None = None,
) -> MummerResult:
    """Align one pair and return the block table, the alignment, and the raw outputs.

    ``target`` and ``query`` are ``(name, sequence)``. Names go into FASTA headers, and MUMmer
    cuts a tag at the first whitespace, so whitespace is replaced by ``_`` — a tag that came back
    different from the name it was given would silently orphan a block.
    """
    binaries = _resolve_pipeline(directory)
    version = nucmer_version(binaries["nucmer"])
    target_name, query_name = _tag(target[0]), _tag(query[0])

    with tempfile.TemporaryDirectory(prefix="mummer_") as work:
        workdir = Path(work)
        _write_fasta(workdir / "target.fa", target_name, target[1])
        _write_fasta(workdir / "query.fa", query_name, query[1])
        ctx.tool(
            "nucmer",
            [
                binaries["nucmer"],
                "--maxmatch",
                "-l",
                str(min_match),
                "-c",
                str(min_cluster),
                "-t",
                "1",
                "-p",
                "out",
                "target.fa",
                "query.fa",
            ],
            tool_path=binaries["nucmer"],
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        delta_path = workdir / "out.delta"
        if not delta_path.is_file():
            raise ToolError("nucmer exited 0 but wrote no out.delta — nothing to filter")

        filtered = workdir / "out.filter.delta"
        record = ctx.tool(
            "delta-filter",
            [binaries["delta-filter"], filter_mode, "out.delta"],
            tool_path=binaries["delta-filter"],
            pinned_version=version,  # same package release as nucmer, verified above
            cwd=workdir,
            timeout=timeout,
        )
        filtered.write_text(Path(record.stdout_file).read_text(encoding="utf-8"), encoding="utf-8")

        record = ctx.tool(
            "show-coords",
            [binaries["show-coords"], "-rcl", "out.filter.delta"],
            tool_path=binaries["show-coords"],
            pinned_version=version,
            cwd=workdir,
            timeout=timeout,
        )
        coords = Path(record.stdout_file).read_text(encoding="utf-8")
        delta = filtered.read_text(encoding="utf-8")

    blocks = parse_coords(coords)
    delta_blocks = parse_delta(delta)
    _require_same_blocks(blocks, delta_blocks)
    return MummerResult(
        blocks=blocks, delta_blocks=delta_blocks, coords=coords, delta=delta, version=version
    )


def _require_same_blocks(
    blocks: tuple[BlockRow, ...], delta_blocks: tuple[DeltaBlock, ...]
) -> None:
    """The table and the alignment must describe the same alignments — or nothing is reported."""
    if len(blocks) != len(delta_blocks):
        raise ToolError(
            f"show-coords reports {len(blocks)} alignment(s) but the filtered delta holds "
            f"{len(delta_blocks)} — the two readings of one run disagree, so neither is usable"
        )
    table = sorted(
        (
            row.target_start,
            row.target_end,
            min(row.query_start, row.query_end),
            max(row.query_start, row.query_end),
        )
        for row in blocks
    )
    alignment = sorted(
        (block.target_start + 1, block.target_end, block.query_start + 1, block.query_end)
        for block in delta_blocks
    )
    if table != alignment:
        raise ToolError(
            "show-coords and the delta file disagree about which blocks the run contains — "
            "refusing to report a table and a series drawn from different alignments"
        )


def _tag(name: str) -> str:
    """A FASTA tag MUMmer will hand back unchanged."""
    return "".join("_" if character.isspace() else character for character in name)


def _write_fasta(path: Path, name: str, sequence: str, width: int = 60) -> None:
    body = "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))
    path.write_text(f">{name}\n{body}\n", encoding="utf-8")
