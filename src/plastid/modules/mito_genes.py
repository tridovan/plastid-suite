"""Plant mitochondrial gene set (P3-04) — the path that must not touch the IR gate.

A plant mitogenome has no canonical inverted repeat, so nothing on this path may call
``ctx.require_ir()``: measured 2026-10-03, both stand-in mitogenomes raise ``IRDetectionError``,
and they are right to. A module that gated here would refuse to describe a genome it can read
perfectly well — the gate exists to stop *position-based* analysis, not to stop counting genes.

The other half of the mito problem is naming. The stand-in mitogenome ``OX459128.1`` carries no
``/gene`` qualifier on any feature, so the name our reader reports is the record's ``/product``
sentence. ``data/mito_gene_names.yaml`` maps those sentences to standard symbols, and is a
separate file from the plastid aliases because the same sentence occurs in both organelles.
A product text with no mapping is reported as it arrived, never guessed.
"""

from __future__ import annotations

from typing import Any

from plastid.core.gene_names import (
    SEMANTIC_RULE,
    GenomeNames,
    load_mito_aliases,
    normalize_genome,
)
from plastid.core.model import GENE_KINDS, Gene, Genome, count_genes
from plastid.core.organelle import GROUP_MITO
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "mito_genes"
MODULE_VERSION = "0.1.0"

SCOPE_LABEL = "every counted mitochondrial gene, in the record's own coordinates"
NO_IR_NOTE = (
    "No inverted-repeat gate applies: a plant mitogenome has no canonical IR, so this module "
    "never asks for one. The gate stays where position-based analysis needs it."
)


@module(
    MODULE_ID,
    group=GROUP_MITO,
    version=MODULE_VERSION,
    title="Mitochondrial gene set",
    requires_ir=False,  # deliberate: see the module docstring and NO_IR_NOTE
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    aliases = load_mito_aliases()
    names = normalize_genome(genome, aliases)
    counts = count_genes(genome)

    counted = [gene for gene in genome.genes if gene.kind in GENE_KINDS]
    # A name with whitespace in it is a product sentence, not a gene symbol — which is what a
    # record without /gene qualifiers leaves us (measured on OX459128.1). One we cannot map is
    # reported as it arrived; a record that already gives symbols is not second-guessed.
    unmapped = sorted(
        {
            gene.name
            for gene in counted
            if " " in gene.name and names.normalized(gene).rule != SEMANTIC_RULE
        }
    )

    summary = Table(
        id="mito_summary",
        title="Mitochondrial gene set",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("topology", "Topology"),
            Column("n_genes", "Genes", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_named_by_symbol", "Named by a standard symbol", dtype="int"),
            Column("ir_gate", "Inverted-repeat gate"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "length_bp": len(genome),
                "topology": genome.topology,
                "n_genes": counts.total,
                "n_cds": counts.cds,
                "n_trna": counts.trna,
                "n_rrna": counts.rrna,
                "n_named_by_symbol": sum(
                    1 for gene in counted if names.normalized(gene).rule == "semantic"
                ),
                "ir_gate": "not applicable",
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    gene_set = Table(
        id="gene_set",
        title="Mitochondrial genes",
        columns=(
            Column("gene_id", "Gene id"),
            Column("gene", "Gene"),
            Column("name_as_arrived", "Name as it arrived"),
            Column("kind", "Kind"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("product", "Product in the record"),
        ),
        rows=tuple(_gene_row(gene, names) for gene in counted),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="gene_id",
    )

    tables: dict[str, Table] = {summary.id: summary, gene_set.id: gene_set}
    if unmapped:
        tables["unmapped_names"] = Table(
            id="unmapped_names",
            title="Names with no standard symbol — reported, not renamed",
            columns=(
                Column("gene_id", "Gene id"),
                Column("name", "Name as it arrived"),
                Column("product", "Product in the record"),
            ),
            rows=tuple(
                {
                    "gene_id": f"{gene.name}:{gene.span[0]}-{gene.span[1]}",
                    "name": gene.name,
                    "product": gene.product or "",
                }
                for gene in counted
                if gene.name in unmapped
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="gene_id",
        )

    notes = [NO_IR_NOTE]
    if unmapped:
        notes.append(
            f"{len(unmapped)} name(s) have no entry in data/mito_gene_names.yaml and are listed "
            f"as they arrived: a name this toolkit cannot map is a name it does not invent."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
    )


def _gene_row(gene: Gene, names: GenomeNames) -> dict[str, Any]:
    entry = names.normalized(gene)
    start, end = gene.span
    return {
        "gene_id": f"{entry.canonical}:{start}-{end}",
        # the symbol alone, not the `canonical (synonym)` label the plastid tables print: there a
        # synonym is another *name*, here it is the product sentence, which already has its own
        # column — printing it twice would be noise, not information
        "gene": entry.canonical,
        "name_as_arrived": gene.name,
        "kind": gene.kind,
        "length_bp": gene.length,
        "product": gene.product or "",
    }
