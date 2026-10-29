"""Publication gene table — the table a plastome paper prints, grouped by function (M2.2).

**Counting convention.** ``Number`` counts *copies*, both inverted-repeat copies included, so
the table totals the fixed ``132 = 86 CDS + 38 tRNA + 8 rRNA`` and can be reconciled against
:func:`~plastid.core.model.count_genes` — which this module does, in its own table, rather
than trusting the grouping to have covered everything (the same reconciliation the SSR module
gained in P2-06, for the same reason: a silent drop looks exactly like a smaller number).

The grouping itself is data (``data/gene_groups.yaml``); this module never hard-codes a gene
name. A gene no group claims is reported as an unclassified row and named in the notes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from plastid.core.gene_groups import GeneGroup, GeneGroups, assign, load_gene_groups
from plastid.core.gene_names import GenomeNames, load_gene_aliases, normalize_genome
from plastid.core.model import GENE_KINDS, Gene, GeneCounts, Genome, count_genes
from plastid.core.organelle import GROUP_PLASTID
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "gene_table"
MODULE_VERSION = "0.2.0"

SCOPE_LABEL = "every counted gene, both IR copies; grouped by function"

_KINDS = ("CDS", "tRNA", "rRNA")


@module(
    MODULE_ID,
    group=GROUP_PLASTID,
    version=MODULE_VERSION,
    title="Publication gene table",
    requires_ir=False,  # grouping a gene by function is not a position-based analysis
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    groups = load_gene_groups()
    names = normalize_genome(genome, load_gene_aliases())
    rows = assign(genome, groups, canonical_of=names.canonical)
    counts = count_genes(genome)
    # Copies, not names: every other number in the summary counts copies, and a reader
    # subtracting this cell from `n_genes_total` must land on the grouped rows. Measured
    # 29/10/2026 on NC_016068.1 (Solanaceae): 24 copies of `hypothetical protein` had no group
    # while the summary read 1.
    unclassified_genes = [row.gene for row in rows if row.group is None]
    unclassified_labels = sorted({names.label(gene) for gene in unclassified_genes})

    by_group: dict[GeneGroup, list[Gene]] = defaultdict(list)
    by_kind: Counter[str] = Counter()
    for row in rows:
        by_kind[row.gene.kind] += 1
        if row.group is not None:
            by_group[row.group].append(row.gene)

    summary = Table(
        id="gene_table_summary",
        title="Publication gene table — summary",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("n_genes_total", "Genes", dtype="int"),
            Column("n_groups", "Functional groups", dtype="int"),
            Column("n_unique_names", "Distinct gene names", dtype="int"),
            Column("n_unclassified", "Unclassified", dtype="int"),
            Column("coord_system", "Coordinates"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "n_genes_total": counts.total,
                "n_groups": sum(1 for genes in by_group.values() if genes),
                # distinct *canonical* names: the count a reader of this table can reproduce by
                # reading it, which is not the same as the record's own name count once two
                # dialects of one locus are merged
                "n_unique_names": len(
                    {names.canonical(gene) for gene in genome.genes if gene.kind in GENE_KINDS}
                ),
                "n_unclassified": len(unclassified_genes),
                "coord_system": "genbank",
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    groups_table = Table(
        id="gene_groups",
        title="Genes by functional group",
        columns=(
            Column("category", "Category"),
            Column("group", "Group of genes"),
            Column("genes", "Name of genes"),
            Column("copies", "Number", dtype="int"),
        ),
        rows=(
            *_group_rows(groups, by_group, names),
            *_unclassified_rows(unclassified_genes, names),
            _total_row(len(rows)),
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="group",
        emphasis=("Total",),
    )

    reconciliation = Table(
        id="gene_reconciliation",
        title="Gene table reconciled against the genome-wide count",
        columns=(
            Column("kind", "Kind"),
            Column("in_table", "In the gene table", dtype="int"),
            Column("counted", "count_genes()", dtype="int"),
            Column("agreement", "Agreement"),
        ),
        rows=tuple(
            _reconciliation_row(kind, _in_table(by_kind, kind), _counted_for(counts, kind))
            for kind in (*_KINDS, "Total")
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="kind",
    )

    notes = [
        "Number counts copies, both inverted-repeat copies included.",
        "A name marked (pseudo) has at least one pseudogene copy among its copies.",
    ]
    aliased = sorted({names.label(row.gene) for row in rows if names.normalized(row.gene).synonyms})
    if aliased:
        notes.append(
            "A name printed with a synonym is one locus under two nomenclatures (D3-02); the "
            f"canonical name leads: {', '.join(aliased)}"
        )
    if unclassified_labels:
        notes.append(
            f"Unclassified genes (no group in gene_groups.yaml): {', '.join(unclassified_labels)}"
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            summary.id: summary,
            groups_table.id: groups_table,
            reconciliation.id: reconciliation,
        },
        notes=tuple(notes),
    )


def _in_table(by_kind: Counter[str], kind: str) -> int:
    """``Total`` is not a key of the counter — it is every kind summed (the first run of this
    module reported 0 there, which is exactly the silent-wrong the reconciliation exists for)."""
    if kind == "Total":
        return sum(by_kind.values())
    return by_kind.get(kind, 0)


def _counted_for(counts: GeneCounts, kind: str) -> int:
    return {"CDS": counts.cds, "tRNA": counts.trna, "rRNA": counts.rrna, "Total": counts.total}[
        kind
    ]


def _group_rows(
    groups: GeneGroups, by_group: dict[GeneGroup, list[Gene]], names: GenomeNames
) -> tuple[dict[str, Any], ...]:
    """One row per non-empty group, in the order the grouping file declares them.

    File order, not genome order: photosystem before ribosome is the order a reader expects,
    and it is the one thing that must not drift between two genomes of the same paper.
    """
    return tuple(
        {
            "category": group.category,
            "group": group.name,
            "genes": _names(by_group[group], names),
            "copies": len(by_group[group]),
        }
        for group in groups.groups
        if by_group.get(group)
    )


def _unclassified_rows(genes: list[Gene], names: GenomeNames) -> tuple[dict[str, Any], ...]:
    """The row the grouping file promises for a gene no group claims.

    Without it the table's rows do not add up to its own total, and the gene is only counted.
    """
    if not genes:
        return ()
    return (
        {
            "category": "",
            "group": "Unclassified",
            "genes": _names(genes, names),
            "copies": len(genes),
        },
    )


def _total_row(total: int) -> dict[str, Any]:
    """The closing row. It carries the key ``Total`` rather than a blank: every row of a keyed
    table stays addressable, so a report can point at the grand total like any other cell."""
    return {"category": "", "group": "Total", "genes": "", "copies": total}


def _names(genes: list[Gene], names: GenomeNames) -> str:
    """Each locus once, in a stable order, with the copy count and pseudogene marker.

    The label carries both dialects (``pafI (ycf3)``): a reader who knows either name finds the
    gene, and never has to know which dialect this record happened to use.
    """
    labels = {gene: names.label(gene) for gene in genes}
    copies = Counter(labels.values())
    pseudo = {labels[gene] for gene in genes if gene.pseudo}
    parts = []
    for label in sorted(copies):
        text = f"{label} (pseudo)" if label in pseudo else label
        if copies[label] > 1:
            text = f"{text} (\N{MULTIPLICATION SIGN}{copies[label]})"
        parts.append(text)
    return ", ".join(parts)


def _reconciliation_row(kind: str, in_table: int, counted: int) -> dict[str, Any]:
    return {
        "kind": kind,
        "in_table": in_table,
        "counted": counted,
        "agreement": "agrees" if in_table == counted else f"DISAGREES — expected {counted}",
    }
