"""Introns — count, length, position and phase per gene (PROJECT.md §13.2, ROADMAP M2.7).

Introns are derived from the CDS exon boundaries: the gap between two consecutive parts of one
gene. Not every gap is an intron — a junction whose parts sit in **different plastid regions**
joins separately transcribed RNAs (trans-splicing) and is reported as such, never counted as an
intron.

Measured on the golden genome (2026-10-03): the two rps12 junctions cross LSC→IRa (28,695 bp) and
LSC→IRb (69,749 bp) while all 17 real introns (541–1,055 bp) sit inside a single region — and the
record's own ``/trans_splicing`` qualifier marks exactly those genes. Geometry and annotation
agree; a disagreement is reported, not smoothed over.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from plastid.core.model import Gene, Genome, cds_features
from plastid.core.normalize import normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "introns"
MODULE_VERSION = "0.1.0"


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Introns",
    requires_ir=True,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    frame = normalize(genome, ir)
    boundaries = tuple((region.start, region.end) for region in frame.regions)

    rows: list[dict[str, Any]] = []
    junctions: list[dict[str, Any]] = []
    copies: Counter[str] = Counter()
    for gene in cds_features(genome):
        copies[gene.name] += 1
        gene_rows: list[dict[str, Any]] = []
        gene_junctions: list[dict[str, Any]] = []
        preceding = 0
        pairs = zip(gene.parts, gene.parts[1:], strict=False)  # one shorter by design
        for number, (first, second) in enumerate(pairs, start=1):
            preceding += first.length
            start = min(first.end, second.end)
            end = max(first.start, second.start)
            span = {
                "gene": gene.name,
                "start_1b": start + 1,
                "end_1b": end,
                "length_bp": end - start,
            }
            if _region_of(start, boundaries) != _region_of(end - 1, boundaries):
                gene_junctions.append(span)
            else:
                # gene.copy.n: the IR carries a second copy of several genes, and a repeated
                # table key would make every reference to it ambiguous
                intron_id = f"{gene.name}.{copies[gene.name]}.{number}"
                gene_rows.append({"intron_id": intron_id, **span, "phase": preceding % 3})
        _check_span(gene, gene_rows, gene_junctions)
        rows.extend(gene_rows)
        junctions.extend(gene_junctions)

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


def _region_of(position: int, boundaries: tuple[tuple[int, int], ...]) -> int | None:
    """Index of the plastid region holding ``position`` — the discriminator for trans-splicing."""
    for index, (start, end) in enumerate(boundaries):
        if start <= position < end:
            return index
    return None


def _annotation_check(genome: Genome, junctions: list[dict[str, Any]]) -> str:
    """Compare the geometry with the record's own ``/trans_splicing`` qualifier."""
    marked = {gene.name for gene in genome.genes if gene.trans_splicing}
    crossing = {entry["gene"] for entry in junctions}
    if marked == crossing:
        noun = "gene" if len(marked) == 1 else "genes"
        return f"agrees — {len(marked)} {noun} marked trans-splicing"
    return f"DISAGREES for {', '.join(sorted(marked ^ crossing))}"


def _check_span(
    gene: Gene, intron_rows: list[dict[str, Any]], junctions: list[dict[str, Any]]
) -> None:
    """A gene's annotated span is exactly its exons plus the gaps between them.

    This is P2-05's cross-check against the annotation's part boundaries: if the derived gaps do
    not tile the span, the parts and the gaps disagree and the run must stop rather than publish
    an intron list that does not add up. Skipped for trans-spliced genes, whose parts are not
    contiguous by design.
    """
    if gene.trans_splicing or junctions:
        return  # trans-spliced parts are not contiguous: their gaps can even overlap
    low, high = gene.span  # 1-based inclusive
    exons = sum(part.length for part in gene.parts)
    gaps = sum(row["length_bp"] for row in intron_rows)
    if exons + gaps != high - low + 1:
        raise ValueError(
            f"{gene.name}: exons + gaps ({exons + gaps}) do not tile the annotated span "
            f"({high - low + 1}) — the CDS parts and their boundaries disagree"
        )
