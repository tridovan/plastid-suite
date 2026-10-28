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

**Three matrices, one tree at a time.** `matrix=cds` (the default) builds the CDS supermatrix
described above. `matrix=introns` and `matrix=igs` build a matrix from the intron and intergenic-
spacer loci `core.spans` derives, keyed by canonical gene name so two records that spell a locus
two ways still meet. Both are **position-based**, so each genome passes the IR gate before its
spans are read, and the sequence taken is the forward strand (an intron has no reading strand).
The IR's second copy is one locus and is kept once. A locus needs at least two carriers, exactly
as a gene does. `aligner=macse` is refused for these matrices: MACSE aligns codons through a
translation, and a spacer has no frame.

**Partitions are opt-in, and the aligner decides what a codon is.** `partition=none` (the
default) is the unpartitioned supermatrix. `partition=gene` gives IQ-TREE one partition per gene
block, and `partition=codon` gives three per block (alignment columns 1/2/3). The default aligner
is MAFFT on nucleotides; `aligner=macse` uses MACSE 2.07, which aligns through the amino-acid
translation and writes a frameshift as `!` and a stop as `*` (both read as gaps here). It is
opt-in, as it is in CGAS. Under MAFFT a codon partition is a partition of alignment columns, not
of verified codon positions: an indel that shifts the frame inside a gene makes columns 1/2/3 stop
being codon positions for the sequences after it. Under MACSE the alignment is codon-aware, but
the codon partition is still a column partition of the alignment MACSE returned. The scope label
and the summary's `Partitioning`/`Aligner` columns say which combination ran.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from plastid.core.align import Alignment
from plastid.core.gene_names import aliases_for, cds_by_canonical_name, normalize_name
from plastid.core.ir import IRGateError, detect_ir, enforce_ir_gate
from plastid.core.model import Genome, Interval, coding_sequence
from plastid.core.organelle import GROUP_SHARED
from plastid.core.spans import introns_of, spacers_of
from plastid.external.iqtree import run_iqtree
from plastid.external.macse import align_cds, resolve_macse
from plastid.external.mafft import align
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.params import ParamValueError

MODULE_ID = "phylogeny"
# 0.5.0 (28/10/2026): `macse_timeout` — the 34-genome measurement found one gene over the old
# fixed 600 s. 0.4.0 added `matrix=cds|introns|igs`; the default is the CDS supermatrix and every
# earlier number is unchanged.
MODULE_VERSION = "0.5.0"

DEFAULT_THREADS = "1"
DEFAULT_SEED = 1
PARTITION_MODES = ("none", "gene", "codon")
MATRIX_MODES = ("cds", "introns", "igs")
DEFAULT_MATRIX = "cds"
# Seconds allowed for one MACSE call. Measured 28/10/2026 on 34 genomes: 62 genes finished at
# ~63 s each, then one gene exceeded 600 s and the run failed — the default is enough for a panel
# of four, not for a full family, so it is a parameter and the failure names the value.
MACSE_TIMEOUT = 600
ALIGNERS = ("mafft", "macse")


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
            "MAFFT algorithm: auto (default), einsi, linsi, ginsi, fftns",
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
        ParamSpec(
            "matrix",
            "str",
            "which loci the supermatrix holds: cds (default), introns, or igs (intergenic spacers)",
        ),
        ParamSpec(
            "partition",
            "str",
            "IQ-TREE partitioning: none (default), gene, or codon (alignment columns 1/2/3)",
        ),
        ParamSpec(
            "aligner",
            "str",
            "per-gene aligner: mafft (default) or macse (codon-aware, opt-in)",
        ),
        ParamSpec("macse", "str", "path to macse_v2.07.jar or the `macse` wrapper"),
        ParamSpec(
            "macse_timeout",
            "int",
            "seconds allowed per gene for MACSE (default 600; a 34-genome matrix needs more)",
        ),
        ParamSpec(
            "ir_source",
            "str",
            "annotation or sequence: which IR path matrix=introns/igs reads (M11.17)",
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
    partition_mode = _partition_mode(params.get("partition", "none"))
    matrix_mode = _matrix_mode(params.get("matrix", DEFAULT_MATRIX))
    aligner = _aligner(params.get("aligner", "mafft"))
    if matrix_mode != "cds" and aligner == "macse":
        raise ParamValueError(
            "phylogeny: aligner=macse aligns a CDS through its translation; an intron or an "
            "intergenic spacer has no reading frame — drop aligner=macse for matrix="
            f"{matrix_mode}"
        )
    macse_binary = None
    if aligner == "macse":
        if "mafft_strategy" in params:
            raise ParamValueError(
                "phylogeny: mafft_strategy has no effect with aligner=macse — set the strategy "
                "on the MAFFT run, or drop it here"
            )
        # Resolved once, before the gene loop: a missing jar or an unpinned release fails here
        # rather than after the first gene has been aligned.
        macse_timeout = int(params.get("macse_timeout", MACSE_TIMEOUT))
        if macse_timeout < 1:
            raise ParamValueError(
                f"phylogeny: macse_timeout must be at least 1 second, got {macse_timeout}"
            )
        macse_binary = resolve_macse(params.get("macse"))
    elif "macse" in params:
        raise ParamValueError(
            "phylogeny: macse=… has no effect with aligner=mafft — pass aligner=macse to use it"
        )
    elif "macse_timeout" in params:
        raise ParamValueError(
            "phylogeny: macse_timeout=… has no effect with aligner=mafft — pass aligner=macse "
            "to use it"
        )

    units: dict[str, dict[str, str]] = {}
    partial_by: dict[str, set[str]] = {}
    if matrix_mode == "cds":
        by_genome = {
            key: cds_by_canonical_name(genome) for key, genome in zip(keys, genomes, strict=True)
        }
        for key, genome in zip(keys, genomes, strict=True):
            for name, gene in by_genome[key].items():
                units.setdefault(name, {})[key] = coding_sequence(genome, gene)
                if gene.partial:
                    partial_by.setdefault(name, set()).add(key)
    else:
        for key, genome in zip(keys, genomes, strict=True):
            for locus, sequence in _span_loci(genome, matrix_mode, params).items():
                units.setdefault(locus, {})[key] = sequence
    names = sorted(units)

    gene_rows: list[dict[str, Any]] = []
    gene_alignments: dict[str, Alignment] = {}
    kept: list[str] = []
    for name in names:
        carriers = [(key, units[name][key]) for key in keys if key in units[name]]
        reason = None
        if len(carriers) < 2:
            reason = f"present in {len(carriers)} genome(s); a block needs at least 2"
        elif (
            matrix_mode == "cds"
            and not include_partial
            and partial_by.get(name, set()) & {key for key, _ in carriers}
        ):
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
        if macse_binary is not None:
            alignment = align_cds(carriers, ctx, binary=macse_binary, timeout=macse_timeout)
        else:
            alignment = align(carriers, ctx, strategy=str(params.get("mafft_strategy", "auto")))
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
            f"phylogeny: {len(kept)} {_unit_word(matrix_mode)}(s) usable across {len(genomes)} "
            f"genomes — a supermatrix needs at least 2 shared loci; see the locus table for why "
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
    charsets = _block_charsets(blocks, partition_mode)
    partition_text = _partition_nexus(charsets) if charsets else None

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
            partition=partition_text,
        )

    partition_label = (
        "none" if partition_mode == "none" else f"{partition_mode} ({len(charsets)} partitions)"
    )
    aligner_label = (
        f"MACSE {macse_binary.version} per gene (codon-aware; ! and * read as gaps)"
        if macse_binary is not None
        else "MAFFT per gene"
    )
    matrix_label = {
        "cds": "CDS blocks",
        "introns": "intron blocks (forward strand)",
        "igs": "intergenic-spacer blocks (forward strand)",
    }[matrix_mode]
    scope = Scope(
        label=(
            f"{len(keys)} genomes; {len(kept)} {matrix_label} concatenated in name order; "
            f"one copy per gene name; {aligner_label}; partition: {partition_label}; "
            f"IQ-TREE {result.version}"
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
            Column("matrix", "Matrix"),
            Column("genes_included", "Genes included", dtype="int"),
            Column("genes_excluded", "Genes excluded", dtype="int"),
            Column("alignment_sites", "Supermatrix sites", dtype="int", unit="bp"),
            Column("partition", "Partitioning"),
            Column("aligner", "Aligner"),
            Column("model", "Model"),
            Column("tool", "Tool"),
            Column("version", "Version"),
            Column("tool_call", "Command line"),
        ),
        rows=(
            {
                "genomes": len(keys),
                "matrix": matrix_mode,
                "genes_included": len(kept),
                "genes_excluded": len(gene_rows) - len(kept),
                "alignment_sites": matrix.length,
                "partition": partition_label,
                "aligner": "macse " + macse_binary.version if macse_binary else "mafft",
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
            **({"supermatrix.partition.nex": partition_text} if partition_text is not None else {}),
        },
        citations=("mafft", "iqtree"),
    )


def _matrix_mode(value: object) -> str:
    """``cds`` · ``introns`` · ``igs``, refused by name otherwise (a typo is not a default)."""
    mode = str(value).strip().lower()
    if mode not in MATRIX_MODES:
        raise ParamValueError(
            f"phylogeny: matrix={value!r} is not one of "
            f"{', '.join(MATRIX_MODES)} — choose the matrix explicitly"
        )
    return mode


def _unit_word(mode: str) -> str:
    return {"cds": "gene", "introns": "intron", "igs": "spacer"}[mode]


def _span_loci(genome: Genome, mode: str, params: dict[str, Any]) -> dict[str, str]:
    """The intron or intergenic-spacer sequences of one genome, keyed by a comparable locus id.

    The ids come from `core.spans`, the one derivation the `introns`/`diversity` tables use, so a
    locus here is the locus those tables name. Gene names are canonicalised before they enter an
    id because the whole point of a cross-genome matrix is that two records spell one locus two
    ways; the IR's second copy is the same locus and is kept once (first copy wins), as `stats`
    collapses it.

    **Position-based, so the IR gate runs**: an intron boundary is a coordinate, and a record whose
    two copies disagree past the threshold cannot supply one. Each genome is checked on its own —
    `ctx.ir` is the reference's, and a query genome's IR is its own.

    The sequence is the **forward strand**. An intron is not transcribed, so there is no reading
    strand to follow; taking the record's own strand would make the same locus align in two
    orientations across records. The scope label says so.
    """
    source = str(params.get("ir_source", "annotation"))
    label = genome.accession or genome.name
    try:
        ir = enforce_ir_gate(detect_ir(genome, source=source))
    except IRGateError as error:
        raise IRGateError(f"phylogeny: matrix={mode} needs an IR on {label}: {error}") from error
    aliases = aliases_for(genome)
    sequences: dict[str, str] = {}
    if mode == "introns":
        for span in introns_of(genome, ir).introns:
            gene, _copy, number = span.intron_id.rsplit(".", 2)
            locus = f"{normalize_name(gene, aliases).canonical}.{number}"
            sequences.setdefault(locus, genome.subseq(Interval(span.start, span.end)))
    else:
        for spacer in spacers_of(genome, ir).igs:
            left = normalize_name(spacer.left, aliases).canonical
            right = normalize_name(spacer.right, aliases).canonical
            locus = f"{left}-{right}"
            sequences.setdefault(locus, genome.subseq(Interval(spacer.start, spacer.end)))
    return sequences


def _block_charsets(
    blocks: Sequence[tuple[str, dict[str, str]]], mode: str
) -> list[tuple[str, int, int, int | None]]:
    """1-based inclusive column ranges for the partition file, in supermatrix order.

    ``gene`` yields one range per block; ``codon`` yields three (columns 1/2/3 of the block).
    The cursor is the first column of the block, so the ranges tile the matrix without gaps or
    overlaps. ``none`` yields nothing — the adapter then runs IQ-TREE without ``-p``.
    """
    charsets: list[tuple[str, int, int, int | None]] = []
    used_names: set[str] = set()
    cursor = 1
    for block_name, block in blocks:
        width = len(next(iter(block.values())))
        if mode == "gene":
            charsets.append((_nexus_name(block_name, used_names), cursor, cursor + width - 1, None))
        elif mode == "codon":
            if width < 3:
                raise ValueError(
                    f"phylogeny: gene {block_name!r} aligns to {width} column(s); a codon "
                    f"partition needs at least 3 — use partition=gene for this matrix"
                )
            for position in (1, 2, 3):
                charsets.append(
                    (
                        _nexus_name(f"{block_name}_codon{position}", used_names),
                        cursor + position - 1,
                        cursor + width - 1,
                        3,
                    )
                )
        cursor += width
    return charsets


def _aligner(value: object) -> str:
    """``mafft`` (default) or ``macse``, refused by name otherwise."""
    aligner = str(value).strip().lower()
    if aligner not in ALIGNERS:
        raise ParamValueError(
            f"phylogeny: aligner={value!r} is not one of {', '.join(ALIGNERS)} — choose the "
            f"per-gene aligner explicitly"
        )
    return aligner


def _partition_mode(value: object) -> str:
    """``none`` · ``gene`` · ``codon``, refused by name otherwise (a typo is not a default)."""
    mode = str(value).strip().lower()
    if mode not in PARTITION_MODES:
        raise ParamValueError(
            f"phylogeny: partition={value!r} is not one of "
            f"{', '.join(PARTITION_MODES)} — choose the partitioning explicitly"
        )
    return mode


def _nexus_name(name: str, used: set[str]) -> str:
    """A NEXUS-safe, unique charset name for a gene block (``rps12_codon3`` stays readable)."""
    clean = re.sub(r"[^A-Za-z0-9_.]", "_", name) or "gene"
    if clean[0].isdigit():
        clean = f"g{clean}"
    candidate, suffix = clean, 2
    while candidate in used:
        candidate = f"{clean}_{suffix}"
        suffix += 1
    used.add(candidate)
    return candidate


def _partition_nexus(charsets: Sequence[tuple[str, int, int, int | None]]) -> str:
    """A NEXUS ``begin sets`` block for IQ-TREE's ``-p``.

    Ranges are 1-based inclusive supermatrix columns, which is the coordinate system IQ-TREE
    reads. ``step`` is written with NEXUS's ``\\N`` syntax (codon position 1 is ``1-N\\3``), so
    the file stays a range description rather than a list of every column.
    """
    lines = ["#nexus", "begin sets;"]
    for name, start, end, step in charsets:
        step_text = f"\\{step}" if step is not None else ""
        lines.append(f"  charset {name} = {start}-{end}{step_text};")
    lines += ["end;", ""]
    return "\n".join(lines)


def _flag(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _to_fasta(alignment: Alignment, width: int = 60) -> str:
    lines: list[str] = []
    for name, sequence in zip(alignment.names, alignment.sequences, strict=True):
        lines.append(f">{name}")
        lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"
