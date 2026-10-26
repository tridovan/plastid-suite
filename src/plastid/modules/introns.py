"""Introns — count, length, position and phase per gene (ROADMAP M2.7).

Introns are derived from the CDS exon boundaries: the gap between two consecutive parts of one
gene. Not every gap is an intron — a junction whose parts sit in **different plastid regions**
joins separately transcribed RNAs (trans-splicing) and is reported as such, never counted as an
intron. The derivation itself lives in :mod:`plastid.core.spans`, because ``diversity`` measures
π over exactly these spans and two modules must not be able to disagree about a boundary.

Measured on the golden genome (2026-10-03): the two rps12 junctions cross LSC→IRa (28,695 bp) and
LSC→IRb (69,749 bp) while all 17 real introns (541–1,055 bp) sit inside a single region — and the
record's own ``/trans_splicing`` qualifier marks exactly those genes. Geometry and annotation
agree; a disagreement is reported, not smoothed over.
"""

from __future__ import annotations

from typing import Any

from plastid.core.model import Genome
from plastid.core.organelle import GROUP_SHARED
from plastid.core.spans import introns_of
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "introns"
MODULE_VERSION = "0.1.1"  # 0.1.1: the derivation moved to plastid.core.spans (P2B-01)


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Introns",
    requires_ir=True,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    derived = introns_of(genome, ctx.require_ir())

    rows: list[dict[str, Any]] = [
        {
            "intron_id": span.intron_id,
            "gene": span.gene,
            "start_1b": span.start + 1,
            "end_1b": span.end,
            "length_bp": span.length,
            "phase": span.phase,
        }
        for span in derived.introns
    ]
    junctions: list[dict[str, Any]] = [
        {
            "gene": span.gene,
            "start_1b": span.start + 1,
            "end_1b": span.end,
            "length_bp": span.length,
        }
        for span in derived.junctions
    ]

    summary = Table(
        id="intron_summary",
        title="Intron summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_genes_with_introns", "Genes with introns", dtype="int"),
            Column("n_introns", "Introns", dtype="int"),
            Column("total_length_bp", "Total intron length", dtype="int", unit="bp"),
            Column("longest_bp", "Longest intron", dtype="int", unit="bp"),
            Column("n_trans_splicing", "Trans-splicing junctions", dtype="int"),
            Column("annotation_check", "Annotation cross-check"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "n_genes_with_introns": len({row["gene"] for row in rows}),
                "n_introns": len(rows),
                "total_length_bp": sum(row["length_bp"] for row in rows),
                "longest_bp": max((row["length_bp"] for row in rows), default=0),
                "n_trans_splicing": len(junctions),
                "annotation_check": _annotation_check(genome, junctions),
            },
        ),
        scope=Scope(
            label=(
                "CDS only, pseudogenes excluded, both IR copies included; "
                "gaps inside one plastid region count as introns"
            )
        ),
    )

    table = Table(
        id="introns",
        title="Introns",
        columns=(
            Column("intron_id", "Intron (gene.copy.n)"),
            Column("gene", "Gene"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("phase", "Phase", dtype="int"),
        ),
        rows=tuple(rows),
        scope=Scope(
            label="GenBank coordinates; exon gaps from the CDS parts, both IR copies included"
        ),
        key="intron_id",
    )

    notes = [
        f"{len(junctions)} trans-splicing junction(s) found — not counted as introns.",
        "Phase = bases of the preceding exon(s) modulo three, in transcript order.",
    ]
    notes.extend(
        f"trans-splicing junction in {entry['gene']}: {entry['start_1b']}–{entry['end_1b']} "
        f"({entry['length_bp']} bp) — its two parts lie in different plastid regions."
        for entry in junctions
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"intron_summary": summary, "introns": table},
        notes=tuple(notes),
        citations=("biopython",),
    )


def _annotation_check(genome: Genome, junctions: list[dict[str, Any]]) -> str:
    """Compare the geometry with the record's own ``/trans_splicing`` qualifier."""
    marked = {gene.name for gene in genome.genes if gene.trans_splicing}
    crossing = {entry["gene"] for entry in junctions}
    if marked == crossing:
        noun = "gene" if len(marked) == 1 else "genes"
        return f"agrees — {len(marked)} {noun} marked trans-splicing"
    return f"DISAGREES for {', '.join(sorted(marked ^ crossing))}"
