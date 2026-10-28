"""Phylogeny — per-gene alignment, concatenated supermatrix, maximum-likelihood tree (M11.21).

**The route is CGAS's, named rather than implied.** `cgas_module14.py` aligns each protein-coding
gene with MAFFT, concatenates the alignments into one matrix, and infers a tree with IQ-TREE
(`-m MFP -bb 1000 -alrt 1000 -bnni`). This module does the same through this toolkit's own
pieces: `external.mafft.align()` per gene, a matrix built here, `external.iqtree.run_iqtree()`.
It is **one** engine — no second aligner, no second inference tool, and the tool version and argv
of every call land in the manifest.

**The rule that decides which genes enter the matrix**, written here because it changes the tree
and a reader has to be able to check it:

* **CDS only, one copy per gene name.** `cds_by_canonical_name` is the same helper `popgen`,
  `kaks` and `substitutions` use, so the inverted repeat's second copy is counted once — the
  gene set of a comparison is *names*, not features.
* **Pseudo genes are excluded** (the helper drops them) and **partial CDS are excluded by
  default**, because a partial gene's alignment is missing an end and would silently lengthen the
  matrix's tail of gaps. `--set include_partial=yes` admits them and the matrix table says which
  gene was partial.
* **A gene needs at least two carriers** to be aligned at all. A gene present in one genome only
  is **excluded and named**, never filled with gaps as if it were absent-but-measured.
* **Missing data is allowed once a gene is in.** A genome that lacks a gene gets a run of gaps
  for that block — that is the supermatrix convention, and the `matrix_genes` table reports the
  per-gene coverage so a tree driven by a sparse gene cannot be read as a complete one.

**What this module does not do.** It does not choose an outgroup, does not root the tree, and
does not claim the tree is right. A topology is not a verifiable number, so the tables report what
was built (taxa, blocks, sites, tool, argv, support present) and leave interpretation to the
reader. `plastid figures` draws the `.treefile` it writes.
"""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from plastid.core.align import Alignment
from plastid.core.gene_names import cds_by_canonical_name
from plastid.core.model import Genome, coding_sequence
from plastid.core.organelle import GROUP_SHARED
from plastid.external.iqtree import run_iqtree
from plastid.external.mafft import align
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "phylogeny"
MODULE_VERSION = "0.1.0"

DEFAULT_THREADS = "1"
DEFAULT_SEED = 1


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Phylogeny (ML tree from a per-gene supermatrix)",
    arity=2,
    citations=("mafft", "iqtree"),
    params=(
        ParamSpec(
            "mafft_strategy",
            "str",
            "MAFFT algorithm: auto (default), einsi, linsi, ginsi, fftns, nwns",
        ),
        ParamSpec("iqtree", "str", "directory holding the pinned iqtree2 binary"),
        ParamSpec("model", "str", "substitution model or MFP for ModelFinder Plus"),
        ParamSpec("ultrafast_bootstrap", "int", "ultrafast bootstrap replicates (at least 1000)"),
        ParamSpec("alrt", "int", "SH-aLRT replicates"),
        ParamSpec("bnni", "bool", "optimise the bootstrap trees by NNI"),
        ParamSpec(
            "threads",
            "str",
            "IQ-TREE threads: 1 (default, reproducible) or a number / AUTO "
            "(faster, not bit-reproducible)",
        ),
        ParamSpec("seed", "int", "IQ-TREE random seed, so the same input gives the same tree"),
        ParamSpec(
            "include_partial", "bool", "include CDS marked partial instead of excluding them"
        ),
        ParamSpec("groups", "str", "sample-group label for each input (supplied by --groups)"),
    ),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    if len(genomes) < 3:
        raise ValueError(
            f"phylogeny takes at least 3 genomes, got {len(genomes)} — a tree of two sequences "
            f"is one comparison, not a topology; use `substitutions` or `diversity` for that"
        )
    keys = [genome.accession or genome.name for genome in genomes]
    labels: dict[str, str] = params.get("groups", {}) or {}
    include_partial = _flag(params.get("include_partial"))

    by_genome = {
        key: cds_by_canonical_name(genome) for key, genome in zip(keys, genomes, strict=True)
    }
    names = sorted({name for mapping in by_genome.values() for name in mapping})

    gene_rows: list[dict[str, Any]] = []
    gene_alignments: dict[str, Alignment] = {}
    kept: list[str] = []
    for name in names:
        carriers = [
            (key, genome, by_genome[key][name])
            for key, genome in zip(keys, genomes, strict=True)
            if name in by_genome[key]
        ]
        reason = None
        if len(carriers) < 2:
            reason = f"present in {len(carriers)} genome(s); a block needs at least 2"
        elif not include_partial and any(_partial(genome, name) for _, genome, _ in carriers):
            reason = "partial CDS in at least one genome (set include_partial=yes to admit)"
        if reason is not None:
            gene_rows.append(
                {
                    "gene": name,
                    "included": "no",
                    "carriers": len(carriers),
                    "sites": None,
                    "reason": reason,
                }
            )
            continue
        alignment = align(
            [(key, coding_sequence(genome, gene)) for key, genome, gene in carriers],
            ctx,
            strategy=str(params.get("mafft_strategy", "auto")),
        )
        gene_alignments[name] = alignment
        kept.append(name)
        gene_rows.append(
            {
                "gene": name,
                "included": "yes",
                "carriers": len(carriers),
                "sites": alignment.length,
                "reason": "",
            }
        )

    if len(kept) < 2:
        raise ValueError(
            f"phylogeny: {len(kept)} gene(s) usable across {len(genomes)} genomes — a "
            f"supermatrix needs at least 2 shared non-partial CDS; see the gene table for why "
            f"the others were left out"
        )

    # Build the matrix: block per gene, in name order, with gaps for a genome that lacks it.
    blocks: list[tuple[str, dict[str, str]]] = []
    for name in kept:
        alignment = gene_alignments[name]
        blocks.append((name, dict(zip(alignment.names, alignment.sequences, strict=True))))
    rows: dict[str, list[str]] = {key: [] for key in keys}
    for _name, block in blocks:
        width = len(next(iter(block.values())))
        for key in keys:
            rows[key].append(block.get(key, "-" * width))
    sequences = {key: "".join(chunks) for key, chunks in rows.items()}
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"phylogeny: supermatrix rows differ in length: {sorted(lengths)}")
    matrix = Alignment(names=tuple(keys), sequences=tuple(sequences[key] for key in keys))

    # The matrix is written to a throwaway path for the tool to read; the *artifact* carries the
    # text, so the module never names a path in the run directory (seam 2: only write_payload
    # writes). `supermatrix.fasta` in the run directory is that artifact, not a side effect.
    matrix_text = _to_fasta(matrix)
    with tempfile.TemporaryDirectory(prefix="phylo_matrix_") as work:
        path = Path(work) / "supermatrix.fasta"
        path.write_text(matrix_text, encoding="utf-8")
        result = run_iqtree(
            path,
            ctx,
            model=str(params.get("model", "MFP")),
            ultrafast_bootstrap=int(params.get("ultrafast_bootstrap", 1000)),
            alrt=int(params.get("alrt", 1000)),
            bnni=_flag(params.get("bnni", True)),
            threads=str(params.get("threads", DEFAULT_THREADS)),
            seed=int(params.get("seed", DEFAULT_SEED)),
            directory=params.get("iqtree"),
        )

    scope = Scope(
        label=(
            f"{len(keys)} genomes; {len(kept)} CDS blocks concatenated in name order; "
            f"one copy per gene name; MAFFT per gene; IQ-TREE {result.version}"
        )
    )
    matrix_rows = [
        {
            "accession": key,
            "group": labels.get(key, ""),
            "sites": len(sequences[key]),
            "gapped_sites": sequences[key].count("-"),
            "blocks": len(kept),
        }
        for key in keys
    ]
    summary = Table(
        id="phylogeny_summary",
        title="Phylogeny summary",
        scope=scope,
        columns=(
            Column("genomes", "Genomes", dtype="int"),
            Column("genes_included", "Genes included", dtype="int"),
            Column("genes_excluded", "Genes excluded", dtype="int"),
            Column("alignment_sites", "Supermatrix sites", dtype="int", unit="bp"),
            Column("model", "Model"),
            Column("tool", "Tool"),
            Column("version", "Version"),
            Column("tool_call", "Command line"),
        ),
        rows=(
            {
                "genomes": len(keys),
                "genes_included": len(kept),
                "genes_excluded": len(gene_rows) - len(kept),
                "alignment_sites": matrix.length,
                "model": str(params.get("model", "MFP")),
                "tool": "iqtree2",
                "version": result.version,
                "tool_call": " ".join(result.argv),
            },
        ),
    )
    genes = Table(
        id="matrix_genes",
        title="Genes in the supermatrix",
        scope=scope,
        columns=(
            Column("gene", "Gene"),
            Column("included", "Included"),
            Column("carriers", "Genomes carrying it", dtype="int"),
            Column("sites", "Alignment sites", dtype="int", unit="bp"),
            Column("reason", "Reason when excluded"),
        ),
        rows=tuple(gene_rows),
    )
    samples = Table(
        id="matrix_samples",
        title="Taxa in the supermatrix",
        scope=scope,
        columns=(
            Column("accession", "Accession"),
            Column("group", "Group"),
            Column("sites", "Sites", dtype="int", unit="bp"),
            Column("gapped_sites", "Gapped sites", dtype="int", unit="bp"),
            Column("blocks", "Blocks", dtype="int"),
        ),
        rows=tuple(matrix_rows),
    )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "phylogeny_summary": summary,
            "matrix_genes": genes,
            "matrix_samples": samples,
        },
        # Keys are `<logical>.<extension>`: the payload writer splits on the last dot, so
        # `supermatrix.fasta` lands as the file of that name in the step's directory.
        artifacts={
            "supermatrix.fasta": matrix_text,
            "tree.treefile": result.treefile,
            "tree.iqtree": result.report,
        },
        citations=("mafft", "iqtree"),
    )


def _flag(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _partial(genome: Genome, name: str) -> bool:
    """Whether this genome's copy of ``name`` is marked partial at either end."""
    from plastid.core.gene_names import cds_by_canonical_name

    gene = cds_by_canonical_name(genome).get(name)
    if gene is None:
        return False
    return any(part.partial_start or part.partial_end for part in gene.parts)


def _to_fasta(alignment: Alignment, width: int = 60) -> str:
    lines: list[str] = []
    for name, sequence in zip(alignment.names, alignment.sequences, strict=True):
        lines.append(f">{name}")
        lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"
