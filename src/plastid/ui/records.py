"""A record as the interface shows it: what it is, its map, and what the IR gate says of it.

Every value comes from the toolkit's own readers and from `core/ir.py`, so what this page says about
a record is what a run would say. The interface never measures anything of its own.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import Any

from plastid.core.digest import sha256_file
from plastid.core.genbank import read_genbank, read_record
from plastid.core.ir import IRDetectionError, IRResult, detect_ir, identity_percent
from plastid.core.model import Genome, count_genes
from plastid.core.normalize import normalize
from plastid.core.sequence import gc_percent

__all__ = ["RECORD_SUFFIXES", "RecordBook", "brief", "inspect"]

RECORD_SUFFIXES = (".gb", ".gbk", ".gbff", ".fasta", ".fa", ".fna", ".fas")
_COUNTED = ("CDS", "tRNA", "rRNA")
# A gene drawn wider than this share of the map is a trans-spliced span, not a place on it.
_WIDEST_TICK_PCT = 10.0


def _ticks(genome: Genome) -> list[list[float]]:
    """Each exon as [left %, width %, strand], the gene map's only drawing data."""
    length = len(genome)
    ticks: set[tuple[float, float, int]] = set()
    for gene in genome.genes:
        if gene.kind not in _COUNTED:
            continue
        for part in gene.parts:
            width = (part.end - part.start) * 100 / length
            if width <= _WIDEST_TICK_PCT:
                strand = 1 if part.strand == 1 else 0
                ticks.add((round(part.start * 100 / length, 2), round(width, 3), strand))
    return [list(tick) for tick in sorted(ticks)]


def _regions(genome: Genome, ir: IRResult) -> tuple[list[dict[str, Any]], list[list[float]]]:
    """The four regions in the toolkit's normalized frame, and the genes drawn in that frame."""
    frame = normalize(genome, ir)
    regions = [
        {"name": region.name, "start": region.start + 1, "end": region.end, "length": region.length}
        for region in frame.regions
    ]
    return regions, _ticks(frame.genome)


def _gate(genome: Genome, source: str) -> dict[str, Any]:
    """What the gate says through one IR source, in the words a run would use."""
    try:
        ir = detect_ir(genome, source=source)
    except IRDetectionError as error:
        return {"source": source, "status": "undetected", "reason": str(error)}
    answer: dict[str, Any] = {
        "source": source,
        "status": "passed" if ir.gate_passed else "refused",
        "identity_pct": identity_percent(ir.identity),
        "length": ir.length,
        "copies": [ir.ira.length, ir.irb.length],
        "annotated_lengths": list(ir.annotated_lengths),
        "denovo_length": ir.denovo_length,
        "matched": ir.matched,
        "mismatched": ir.mismatches,
        "terminal_gaps": ir.terminal_gaps,
        "internal_gaps": ir.internal_gaps,
        "compared": ir.compared_length,
    }
    try:
        answer["regions"], answer["ticks"] = _regions(genome, ir)
    except ValueError as error:
        answer["regions"], answer["ticks"] = [], _ticks(genome)
        answer["map_note"] = str(error)
    return answer


def inspect(path: Path) -> dict[str, Any]:
    """Everything the records page shows about one file."""
    genome = read_record(path)
    counts = count_genes(genome)
    annotated = bool(genome.genes)
    gates: Sequence[str] = ("annotation", "sequence") if annotated else ("sequence",)
    return {
        "organism": genome.source.get("organism") or genome.name,
        "accession": genome.accession,
        "length": len(genome),
        "gc_pct": round(gc_percent(genome.sequence), 2),
        "annotated": annotated,
        "genes": {
            "total": counts.cds + counts.trna + counts.rrna,
            "cds": counts.cds,
            "trna": counts.trna,
            "rrna": counts.rrna,
        },
        "gates": {source: _gate(genome, source) for source in gates},
    }


def brief(path: Path) -> dict[str, Any]:
    """Who a reference record is and how many genes it carries, read as the annotator reads it.

    No IR gate here: a reference record is searched, not measured, and a set of thirty would make
    the list wait on thirty gates.
    """
    genome = read_genbank(path)
    counts = count_genes(genome)
    return {
        "organism": genome.source.get("organism") or genome.name,
        "accession": genome.accession,
        "genes": counts.cds + counts.trna + counts.rrna,
    }


@dataclass(slots=True)
class _Seen:
    digest: str
    summary: dict[str, Any]


class RecordBook:
    """Inspected records by content, so a page reload does not read 150 kb of sequence again."""

    def __init__(self) -> None:
        self._seen: dict[Path, _Seen] = {}
        self._lock = Lock()

    def read(self, path: Path) -> dict[str, Any]:
        digest = sha256_file(path)
        with self._lock:
            kept = self._seen.get(path)
        if kept is not None and kept.digest == digest:
            return kept.summary
        summary = {**inspect(path), "sha256": digest}
        with self._lock:
            self._seen[path] = _Seen(digest, summary)
        return summary
