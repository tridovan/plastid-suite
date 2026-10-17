"""The junction table and the configuration CPJSdraw reads (M10.1) — the pure half.

No external tool is imported here: the drawing is the CLI's job, the way every other non-module step
that calls a tool is wired (`cli.py` reaches `external/`; the `run` layer may not). What lives here
is everything that can be checked without a tool installed.

**Two readers of the same record, and they may disagree.** The four regions come from
:mod:`plastid.core.junctions`, which takes the record's ``repeat_region`` annotation as written. The
**IR gate** comes from :mod:`plastid.core.ir`, which asks whether the two copies really are the same
inverted repeat. A record can pass one and fail the other — measured 17/10/2026 on three of the
eight Iridaceae records, whose annotated copies differ in length by 1, 2 and 4 bp (H-25) — so the
gate's verdict is written **into the table** rather than left to be rediscovered. The figure is
still drawn: it is a picture of what the record says, and refusing to draw it would hide the very
disagreement a reader needs to see.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from plastid.core.ir import IRDetectionError, detect_ir
from plastid.core.junctions import JunctionRegions, junction_regions
from plastid.core.model import Genome

__all__ = ["JunctionRow", "junction_config", "junction_rows", "junction_table"]

# Which features the flanking columns name, and how far out they look. The manuscript's figure
# annotates the genes either side of each junction; a table naming only the nearest would not let a
# reader check the drawing.
_GENE_KINDS = ("CDS", "tRNA", "rRNA")
_FLANK_BP = 2_000
_FLANK_MAX = 4

_COLUMNS = (
    "species",
    "accession",
    "genome_bp",
    "LSC",
    "IRb",
    "SSC",
    "IRa",
    "LSC_bp",
    "IRb_bp",
    "SSC_bp",
    "IRa_bp",
    "IRb-IRa_bp",
    "JLB",
    "JSB",
    "JSA",
    "JLA",
    "ir_gate",
    "ir_identity",
    "ir_source",
    "genes_at_JLB",
    "genes_at_JSB",
    "genes_at_JSA",
    "genes_at_JLA",
)


@dataclass(frozen=True, slots=True)
class JunctionRow:
    """One record's line of the junction table, and the file CPJSdraw will be handed."""

    path: Path
    species: str
    accession: str
    genome_bp: int
    regions: JunctionRegions
    gate: str  # "passed" | "failed" | "not found"
    identity: float | None
    source: str | None
    flanks: tuple[str, str, str, str]  # JLB, JSB, JSA, JLA

    @property
    def warning(self) -> str | None:
        """What a reader must be told when the gate refuses a record this figure still draws."""
        if self.gate == "passed":
            return None
        measured = (
            f" — the two annotated copies are {self.identity * 100:.2f} % identical"
            if self.identity is not None
            else ""
        )
        return (
            f"{self.species} ({self.accession}): the IR gate {self.gate}{measured}. The figure is "
            f"drawn from the record's own annotation; a position-based analysis of this record "
            f"will refuse it (H-25, ROADMAP.md → Known defects)"
        )


def junction_rows(
    records: Sequence[tuple[Path, Genome]], order: Sequence[str] | None = None
) -> tuple[JunctionRow, ...]:
    """One row per record, optionally ordered — by a tree's distances, or by the caller's list.

    ``records`` pairs each record with the file it was read from: CPJSdraw is handed those paths, so
    a line's coordinates and the file they describe cannot come apart.
    """
    if not records:
        raise ValueError("ir-junction: no records — there is nothing to draw")
    rows = [_row(Path(path), genome) for path, genome in records]
    if order is not None:
        rank = {name: position for position, name in enumerate(order)}
        rows.sort(key=lambda row: rank.get(row.species, len(rank)))
    return tuple(rows)


def junction_config(rows: Sequence[JunctionRow]) -> str:
    """The configuration CPJSdraw reads: one absolute path and one region string per line."""
    return "".join(f"{row.path.resolve()}\t{row.regions.config_line}\n" for row in rows)


def junction_table(rows: Sequence[JunctionRow]) -> str:
    """The table as CSV text — English headers, one row per species, no thousands separators."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(_COLUMNS)
    for row in rows:
        regions = row.regions
        writer.writerow(
            [
                row.species,
                row.accession,
                row.genome_bp,
                regions.lsc.text,
                regions.irb.text,
                regions.ssc.text,
                regions.ira.text,
                regions.lsc.length,
                regions.irb.length,
                regions.ssc.length,
                regions.ira.length,
                regions.irb.length - regions.ira.length,
                regions.jlb,
                regions.jsb,
                regions.jsa,
                regions.jla,
                row.gate,
                "" if row.identity is None else f"{row.identity:.4f}",
                row.source or "",
                *row.flanks,
            ]
        )
    return buffer.getvalue()


def _row(path: Path, genome: Genome) -> JunctionRow:
    regions = junction_regions(genome)
    gate, identity, source = _gate(genome)
    return JunctionRow(
        path=path,
        species=genome.name,
        accession=genome.accession or "",
        genome_bp=len(genome),
        regions=regions,
        gate=gate,
        identity=identity,
        source=source,
        flanks=(
            _flanking(genome, regions.jlb),
            _flanking(genome, regions.jsb),
            _flanking(genome, regions.jsa),
            _flanking(genome, regions.jla),
        ),
    )


def _gate(genome: Genome) -> tuple[str, float | None, str | None]:
    """The IR gate's verdict for this record — reported here, acted on by the modules."""
    try:
        result = detect_ir(genome)
    except IRDetectionError:
        return "not found", None, None
    return ("passed" if result.gate_passed else "failed"), result.identity, result.source


def _flanking(genome: Genome, point: int) -> str:
    """The genes within ``_FLANK_BP`` of a junction, nearest first.

    Distances are in the record's own frame, 1-based inclusive: a gene *contains* the junction
    (0 bp), lies *before* it, or *after* it. Ties break on the text so the column is stable — a
    table whose cells move between runs cannot be diffed.
    """
    found: list[tuple[int, str]] = []
    for gene in genome.genes:
        if gene.kind not in _GENE_KINDS:
            continue
        for part in gene.parts:
            start, end = part.start + 1, part.end
            if start <= point <= end:
                distance, where = 0, "contains"
            elif end < point:
                distance, where = point - end, "before"
            else:
                distance, where = start - point, "after"
            if distance <= _FLANK_BP:
                found.append((distance, f"{gene.name} ({where} {distance} bp)"))
            break
    found.sort(key=lambda item: (item[0], item[1]))
    return " | ".join(dict.fromkeys(text for _, text in found[:_FLANK_MAX]))
