"""Analysis modules — one job per file.

Module files never import each other; this package registers them explicitly.
The architecture test permits sibling imports only from this file.
"""

from __future__ import annotations

from plastid.modules import (  # noqa: F401  (importing registers them)
    amino_acid,
    barcodes,
    codon_usage,
    diversity,
    gene_names,
    gene_normalize,
    gene_table,
    genome_stats,
    introns,
    kaks,
    mito_annotate,
    mito_genes,
    orf,
    phylogeny,
    plastid_annotate,
    popgen,
    repeats,
    rna_editing,
    ssr,
    stats,
    structure,
    submission,
    substitutions,
    synteny,
)

__all__ = [
    "amino_acid",
    "barcodes",
    "codon_usage",
    "diversity",
    "gene_names",
    "gene_normalize",
    "gene_table",
    "genome_stats",
    "introns",
    "kaks",
    "mito_annotate",
    "plastid_annotate",
    "mito_genes",
    "orf",
    "phylogeny",
    "popgen",
    "repeats",
    "rna_editing",
    "ssr",
    "stats",
    "structure",
    "submission",
    "substitutions",
    "synteny",
]
