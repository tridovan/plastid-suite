"""Introns — count, length, position and phase per gene, CDS and tRNA (ROADMAP M2.7, M11.15).

An intron is the gap between two consecutive parts of one gene, both ends in one plastid region.
The classes are counted separately and published together with their ``kind``, because they are
different facts and a single total hides which one moved:

* **CDS introns** — the set ``diversity`` measures π over (17 on the golden, both IR copies);
* **tRNA introns** — a tRNA's parts are its two exons around an anticodon-loop intron (8 on the
  golden, the number the published record annotates; ARAGORN only returns them with ``-i``);
* **trans-splicing junctions** — a CDS gap whose parts sit in **different plastid regions**
  (``rps12`` on the golden) joins separately transcribed RNAs and is reported, never counted as
  an intron;
* **not assessed** — a multi-part tRNA gap longer than the tRNA caller's documented maximum, a
  gap across a region boundary, or a join across the origin. These are named in the notes exactly
  so "no intron" and "could not tell" do not look alike.

A class is **not measured** when the record annotates no feature of that kind. An annotation with
no CDS cannot say "no CDS introns"; the summary says so rather than printing a zero that reads as
a measurement.

Measured on the golden genome (2026-10-03): the two rps12 junctions cross LSC→IRa (28,695 bp) and
LSC→IRb (69,749 bp) while all 17 CDS introns (541–1,055 bp) sit inside a single region — and the
record's own ``/trans_splicing`` qualifier marks exactly those genes. Geometry and annotation
agree; a disagreement is reported, not smoothed over. The 8 tRNA introns are checked against the
record's own multi-part tRNA features, not against ARAGORN's count (H-23).
"""

from __future__ import annotations

from typing import Any

from plastid.core.model import Genome
from plastid.core.organelle import GROUP_SHARED
from plastid.core.spans import IntronSpan, introns_of
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "introns"
MODULE_VERSION = "0.2.0"  # 0.2.0: tRNA introns, separate counts and an assessment (M11.15)


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

    cdss = [gene for gene in genome.genes if gene.kind == "CDS" and not gene.pseudo]
    trnas = [gene for gene in genome.genes if gene.kind == "tRNA" and not gene.pseudo]
    rows = [_row(span) for span in (*derived.introns, *derived.trna_introns)]
    junctions: list[dict[str, Any]] = [
        {
            "gene": span.gene,
            "start_1b": span.start + 1,
            "end_1b": span.end,
            "length_bp": span.length,
        }
        for span in derived.junctions
    ]
    judged = rows
    cds_rows = [row for row in rows if row["kind"] == "CDS"]
    trna_rows = [row for row in rows if row["kind"] == "tRNA"]

    summary = Table(
        id="intron_summary",
        title="Intron summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_genes_with_introns", "Genes with introns", dtype="int"),
            Column("n_introns", "Introns", dtype="int"),
            Column("n_cds_introns", "CDS introns", dtype="int"),
            Column("n_trna_introns", "tRNA introns", dtype="int"),
            Column("total_length_bp", "Total intron length", dtype="int", unit="bp"),
            Column("longest_bp", "Longest intron", dtype="int", unit="bp"),
            Column("n_trans_splicing", "Trans-splicing junctions", dtype="int"),
            Column("assessment", "Assessment"),
            Column("annotation_check", "Annotation cross-check"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "n_genes_with_introns": len({row["gene"] for row in judged}),
                "n_introns": len(judged),
                # A missing feature class is not a zero: an annotation with no CDS cannot say
                # "no CDS introns", and printing 0 would read as a measurement (M11.15).
                "n_cds_introns": len(cds_rows) if cdss else "not measured",
                "n_trna_introns": len(trna_rows) if trnas else "not measured",
                "total_length_bp": sum(row["length_bp"] for row in judged),
                "longest_bp": max((row["length_bp"] for row in judged), default=0),
                "n_trans_splicing": len(junctions),
                "assessment": _assessment(bool(cdss), bool(trnas), derived.unassessed),
                "annotation_check": _annotation_check(genome, junctions),
            },
        ),
        scope=Scope(
            label=(
                "CDS and tRNA introns, pseudogenes excluded, both IR copies included; a gap "
                "inside one plastid region counts as an intron, a cross-region gap as a "
                "trans-splicing junction"
            )
        ),
    )

    notes = [
        f"{len(derived.introns)} CDS intron(s) and {len(derived.trna_introns)} tRNA intron(s) "
        f"counted.",
        f"{len(junctions)} trans-splicing junction(s) found — not counted as introns.",
        "Phase = bases of the preceding exon(s) modulo three, in transcript order.",
    ]
    if not cdss:
        notes.append("CDS introns: not measured — the record annotates no CDS.")
    if not trnas:
        notes.append("tRNA introns: not measured — the record annotates no tRNA.")
    notes.extend(f"not assessed: {reason}" for reason in derived.unassessed)
    notes.extend(
        f"trans-splicing junction in {entry['gene']}: {entry['start_1b']}–{entry['end_1b']} "
        f"({entry['length_bp']} bp) — its two parts lie in different plastid regions."
        for entry in junctions
    )

    # A table with no rows is not a table (the schema refuses one, and a reader cannot tell it
    # from one that failed to fill): the summary and the notes still say what was and was not
    # measured, which is the fact this run has to carry when nothing is published.
    tables = {"intron_summary": summary}
    if rows:
        tables["introns"] = Table(
            id="introns",
            title="Introns",
            columns=(
                Column("intron_id", "Intron (gene.copy.n)"),
                Column("gene", "Gene"),
                Column("kind", "Kind"),
                Column("start_1b", "Start", dtype="int"),
                Column("end_1b", "End", dtype="int"),
                Column("length_bp", "Length", dtype="int", unit="bp"),
                Column("phase", "Phase", dtype="int"),
            ),
            rows=tuple(rows),
            scope=Scope(
                label=(
                    "GenBank coordinates; exon gaps from each feature's own parts, both IR "
                    "copies included; phase is a reading-frame offset for CDS and has no frame "
                    "meaning for tRNA"
                )
            ),
            key="intron_id",
        )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        citations=("biopython",),
    )


def _row(span: IntronSpan) -> dict[str, Any]:
    return {
        "intron_id": span.intron_id,
        "gene": span.gene,
        "kind": span.kind,
        "start_1b": span.start + 1,
        "end_1b": span.end,
        "length_bp": span.length,
        "phase": span.phase,
    }


def _assessment(cds: bool, trna: bool, unassessed: tuple[str, ...]) -> str:
    """What this record could and could not answer — a cell, not something a reader infers."""
    measured = [name for name, present in (("CDS", cds), ("tRNA", trna)) if present]
    missing = [name for name, present in (("CDS", cds), ("tRNA", trna)) if not present]
    parts = []
    parts.append(f"{' and '.join(measured)} introns measured" if measured else "nothing measured")
    if missing:
        parts.append(f"{' and '.join(missing)} introns not measured (no such feature)")
    if unassessed:
        parts.append(f"{len(unassessed)} tRNA gap(s) not assessed")
    return "; ".join(parts)


def _annotation_check(genome: Genome, junctions: list[dict[str, Any]]) -> str:
    """Compare the geometry with the record's own ``/trans_splicing`` qualifier."""
    marked = {gene.name for gene in genome.genes if gene.trans_splicing}
    crossing = {entry["gene"] for entry in junctions}
    if marked == crossing:
        noun = "gene" if len(marked) == 1 else "genes"
        return f"agrees — {len(marked)} {noun} marked trans-splicing"
    return f"DISAGREES for {', '.join(sorted(marked ^ crossing))}"
