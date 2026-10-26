"""Genome statistics — the first end-to-end module.

Golden expectations live in the tests; this file must simply report what
the reader and the gate measured, with every table labelled by its coordinate frame.
"""

from __future__ import annotations

from typing import Any

from plastid.core.gene_names import GenomeNames, load_gene_aliases, normalize_genome
from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, Interval, count_genes
from plastid.core.normalize import normalize
from plastid.core.organelle import GROUP_PLASTID
from plastid.core.sequence import gc_percent
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "genome_stats"
MODULE_VERSION = "0.1.0"


@module(
    MODULE_ID,
    group=GROUP_PLASTID,
    version=MODULE_VERSION,
    title="Genome statistics",
    requires_ir=True,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    counts = count_genes(genome)
    frame = normalize(genome, ir)
    names = normalize_genome(genome, load_gene_aliases())
    regions = {region.name: region for region in frame.regions}

    stats = Table(
        id="genome_stats",
        title="Genome statistics",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
            Column("n_genes_total", "Genes", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_gene_features", "gene features (diagnostic only)", dtype="int"),
            Column("n_unique_gene_names", "Unique gene names", dtype="int"),
            Column("ir_length_bp", "IR", dtype="int", unit="bp"),
            Column("ir_identity_pct", "IR identity", dtype="float", unit="%", precision=2),
            Column("ir_source", "IR source"),
            Column("lsc_bp", "LSC", dtype="int", unit="bp"),
            Column("ssc_bp", "SSC", dtype="int", unit="bp"),
            Column("coord_system", "Coordinates"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "length_bp": len(genome),
                "gc_pct": round(gc_percent(genome.sequence), 4),
                "n_genes_total": counts.total,
                "n_cds": counts.cds,
                "n_trna": counts.trna,
                "n_rrna": counts.rrna,
                "n_gene_features": counts.gene_features,
                "n_unique_gene_names": counts.unique_names,
                "ir_length_bp": ir.length,
                "ir_identity_pct": round(ir.identity * 100, 4),
                "ir_source": ir.source,
                "lsc_bp": regions["LSC"].length,
                "ssc_bp": regions["SSC"].length,
                "coord_system": "genbank",
            },
        ),
        scope=Scope(
            label="whole plastome; every copy of CDS/tRNA/rRNA; region split from the IR pair",
            coord_system="genbank",
        ),
    )

    region_table = Table(
        id="regions",
        title="Regions",
        columns=(
            Column("region", "Region"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
        ),
        rows=tuple({"region": region.name, "length_bp": region.length} for region in frame.regions),
        scope=Scope(
            label="normalized frame (LSC → IRa → SSC → IRb)",
            coord_system="normalized",
        ),
        key="region",
    )

    gene_table = Table(
        id="genes",
        title="Genes (GenBank coordinates, every copy)",
        columns=(
            Column("name", "Name"),
            Column("kind", "Kind"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("strand", "Strand"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("in_ir", "In IR"),
            Column("pseudo", "Pseudo"),
        ),
        rows=tuple(_gene_row(gene, ir, names) for gene in genome.genes),
        scope=Scope(label="GenBank coordinates; both IR copies included", coord_system="genbank"),
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={stats.id: stats, region_table.id: region_table, gene_table.id: gene_table},
        notes=[f"IR source: {ir.source}; two copies {ir.identity * 100:.2f} % identical."],
        citations=("biopython",),
    )


def _gene_row(gene: Gene, ir: IRResult, names: GenomeNames) -> dict[str, Any]:
    start, end = gene.span
    return {
        # the label, not the raw name: `pafI (ycf3)` is this locus whichever dialect the
        # record used (P3-02) — one name in every table of one run
        "name": names.label(gene),
        "kind": gene.kind,
        "start_1b": start,
        "end_1b": end,
        "strand": "+" if gene.strand >= 0 else "-",
        "length_bp": gene.length,
        "in_ir": "yes" if _in_ir(gene, ir) else "no",
        "pseudo": "yes" if gene.pseudo else "no",
    }


def _in_ir(gene: Gene, ir: IRResult) -> bool:
    return any(_overlaps(part, copy) for part in gene.parts for copy in (ir.ira, ir.irb))


def _overlaps(part: Interval, copy: Interval) -> bool:
    return part.start < copy.end and copy.start < part.end
