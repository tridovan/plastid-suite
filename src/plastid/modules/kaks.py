"""Ka/Ks selection pressure — pairwise dN/dS per shared CDS (P4-02, M4.1; decision D4-01).

The estimator is PAML's **yn00**, the Yang & Nielsen (2000) method, run on a **codon alignment**
of each gene: MAFFT aligns the two CDS as nucleotides, the terminal stop codon is trimmed (yn00's
poison), and the alignment must be codon-aligned before yn00 ever sees it — yn00 reads three
bases at a time and does not check the frame, so a one-base gap is scored as if it were real
(measured 04/10/2026; ``core.align.codon_frame_breaks``).

Every gene gets a row, including the ones that could not be scored: a gene one genome does not
carry, a gene whose annotation is not whole codons, an alignment that breaks the frame, a run
yn00 itself refused. A comparison that silently dropped them would report a ratio over an
unstated subset — which is the failure this project exists to avoid. ``Ks = 0`` is *labelled*:
the ratio is undefined there, never zero and never one.

This module deliberately does **not** require the IR gate: it is gene-based, not position-based,
and the mito half of Phase 4 (D4-05) has no inverted repeat to gate on.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from plastid.core.align import codon_frame_breaks
from plastid.core.codon import translate
from plastid.core.gene_names import cds_by_canonical_name
from plastid.core.model import Gene, Genome, coding_sequence
from plastid.core.organelle import GROUP_SHARED
from plastid.external.mafft import align
from plastid.external.paml import resolve_binary, run_yn00, yn00_version
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.tools import ToolError

MODULE_ID = "kaks"
MODULE_VERSION = "0.3.0"  # behaviour changed 2026-10-04 after the first runs

ESTIMATOR = "Yang & Nielsen 2000 (yn00)"
STOPS = ("TAA", "TAG", "TGA")

# Statuses a row can carry. "ok" means a ratio; "ks = 0" means the ratio is undefined (and is
# labelled, never divided by); the rest say why a gene has no value at all.
STATUS_OK = "ok"
STATUS_KS_ZERO = "ks = 0"
STATUS_ABSENT_REFERENCE = "absent from the reference"
STATUS_ABSENT_GENOME = "absent from this genome"
STATUS_NOT_WHOLE_CODONS = "not whole codons"
STATUS_CODON_START = "codon_start is not one"
STATUS_INTERNAL_STOP = "internal stop codon"
STATUS_FRAMESHIFT = "frameshift in the alignment"
STATUS_TOOL_FAILED = "yn00 refused the pair"


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Ka/Ks selection pressure",
    requires_ir=False,  # gene-based; see the module docstring
    arity=2,
    citations=("paml", "yang_nielsen", "mafft"),
    params=(
        ParamSpec("groups", "str", "sample-group label for each input (supplied by --groups)"),
        ParamSpec("paml", "str", "directory holding yn00 when it is not on PATH"),
        ParamSpec(
            "mafft_strategy",
            "str",
            "MAFFT algorithm: auto (default), einsi, linsi, ginsi, fftns",
        ),
    ),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    reference = genomes[0]
    reference_name = reference.accession or reference.name
    labels: dict[str, str] = params.get("groups", {})
    directory = params.get("paml")

    # A missing or unpinned estimator is a run-level failure, raised before any gene is scored;
    # a single gene yn00 cannot score is data, and becomes a row below.
    binary = resolve_binary(directory=directory)
    version = yn00_version(binary)

    # The record's own organelle picks the vocabulary (core.gene_names.aliases_for).
    reference_genes = cds_by_canonical_name(reference)
    per_genome = {
        (genome.accession or genome.name): cds_by_canonical_name(genome) for genome in genomes[1:]
    }
    names = sorted(set(reference_genes) | {name for genes in per_genome.values() for name in genes})
    # "Shared" means carried by the reference AND by the genome it is compared with — not "carried
    # by the reference". Measured 04/10/2026 on the mito sketch: the reference carries 27 named
    # CDS and the query 23, and labelling the reference's own count "shared with every genome"
    # put a wrong claim on the summary row. The true intersection is computed per genome, and the
    # summary counts the stricter set — genes no annotated genome is missing.
    annotated = [(key, genes) for key, genes in per_genome.items() if genes]
    shared_by_genome = {
        key: [name for name in reference_genes if name in genes]
        for key, genes in per_genome.items()
    }
    shared = [name for name in reference_genes if all(name in genes for _key, genes in annotated)]
    if not annotated:
        raise ValueError(
            "kaks has nothing to compare: no genome other than the reference carries an "
            "annotated CDS"
        )

    gene_rows: list[dict[str, Any]] = []
    pooled_rows: list[dict[str, Any]] = []
    for genome in genomes[1:]:
        key = genome.accession or genome.name
        genes = per_genome[key]
        scored: list[tuple[float, float, float | None, float | None]] = []
        for name in names:
            row = _score_gene(
                name,
                key,
                labels.get(key, ""),
                reference,
                reference_genes,
                genome,
                genes,
                ctx,
                binary=binary,
                strategy=str(params.get("mafft_strategy", "auto")),
            )
            gene_rows.append(row)
            if row["ka"] is not None and row["ks"] is not None:
                scored.append((row["ka"], row["ks"], row["sites_n"], row["sites_s"]))
        pooled_rows.append(
            _pool(key, labels.get(key, ""), scored, len(shared_by_genome.get(key, ())))
        )

    summary_scope = Scope(
        label=(
            "the reference's CDS names; per-gene values are in kaks_genes, pooled values in "
            "kaks_pooled"
        ),
        coord_system="genbank",
    )
    genes_scope = Scope(
        label=(
            "one row per gene and genome, both inverted-repeat copies counted once; the "
            "alignment is per gene, on the CDS itself, the terminal stop trimmed; Ka and Ks are "
            "PAML yn00's dN and dS; Ka/Ks is undefined (blank) when Ks is zero"
        ),
        include=("CDS",),
        exclude_regions=("IRa", "IRb"),
        coord_system="genbank",
    )
    pooled_scope = Scope(
        label=(
            "pooled over the genes that produced values, weighted by yn00's own site counts "
            "(Ka by nonsynonymous sites, Ks by synonymous sites); it is not the mean of the "
            "per-gene ratios"
        ),
        include=("CDS",),
        exclude_regions=("IRa", "IRb"),
        coord_system="genbank",
    )

    summary = Table(
        id="kaks_summary",
        title="Ka/Ks summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes compared", dtype="int"),
            Column("n_genes", "Genes considered", dtype="int"),
            Column("n_shared", "Genes shared with every annotated genome", dtype="int"),
            Column("n_partial_genes", "Reference genes partial", dtype="int"),
            Column("estimator", "Estimator"),
            Column("aligner", "Aligner"),
            Column("tool", "Tool"),
            Column("tool_version", "Tool version"),
        ),
        rows=(
            {
                "accession": reference_name,
                "organism": reference.name,
                "n_genomes": len(genomes) - 1,
                "n_genes": len(names),
                "n_shared": len(shared),
                "n_partial_genes": sum(
                    1 for name in names if name in reference_genes and reference_genes[name].partial
                ),
                "estimator": ESTIMATOR,
                "aligner": "MAFFT",
                "tool": "PAML yn00",
                "tool_version": version,
            },
        ),
        scope=summary_scope,
    )

    pooled = Table(
        id="kaks_pooled",
        title="Ka/Ks pooled over genes",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("n_shared", "Shared genes", dtype="int"),
            Column("n_computed", "Genes scored", dtype="int"),
            Column("n_not_computed", "Genes without a value", dtype="int"),
            Column("pooled_ka", "Pooled Ka", dtype="float", precision=6),
            Column("pooled_ks", "Pooled Ks", dtype="float", precision=6),
            Column("pooled_ka_ks", "Pooled Ka/Ks", dtype="float", precision=4),
        ),
        rows=tuple(pooled_rows),
        scope=pooled_scope,
        key="genome",
    )

    genes = Table(
        id="kaks_genes",
        title="Ka/Ks by gene",
        columns=(
            Column("row_id", "Row"),
            Column("gene", "Gene"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("partial", "Partial CDS"),
            Column("codons", "Codons", dtype="int"),
            Column("sites_s", "Synonymous sites", dtype="float", precision=1),
            Column("sites_n", "Nonsynonymous sites", dtype="float", precision=1),
            Column("ka", "Ka (dN)", dtype="float", precision=6),
            Column("ks", "Ks (dS)", dtype="float", precision=6),
            Column("ka_ks", "Ka/Ks", dtype="float", precision=4),
            Column("status", "Status"),
            Column("detail", "Detail"),
        ),
        rows=tuple(gene_rows),
        scope=genes_scope,
        key="row_id",
    )

    unshared = [name for name in names if name not in reference_genes]
    notes = [
        f"Estimator: {ESTIMATOR}, through PAML yn00 {version}; each gene is aligned on its own "
        f"CDS with MAFFT, terminal stop trimmed.",
        "A gene with Ks = 0 has an undefined Ka/Ks — the cell is blank, never zero or one.",
        "Genes that could not be scored keep a row and a status; nothing is dropped silently.",
    ]
    if unshared:
        notes.append(
            f"{len(unshared)} gene(s) are carried by a query genome but not by the reference "
            f"record and cannot be compared: {', '.join(unshared)}."
        )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"kaks_summary": summary, "kaks_pooled": pooled, "kaks_genes": genes},
        notes=tuple(notes),
        citations=("paml", "yang_nielsen", "mafft"),
    )


def _score_gene(
    name: str,
    genome_key: str,
    group: str,
    reference: Genome,
    reference_genes: dict[str, Gene],
    genome: Genome,
    genes: dict[str, Gene],
    ctx: RunContext,
    *,
    binary: str,
    strategy: str = "auto",
) -> dict[str, Any]:
    """One gene's row — a value, or the reason there is none. Never a silently missing row."""
    row: dict[str, Any] = {
        "row_id": f"{name}#{genome_key}",
        "gene": name,
        "genome": genome_key,
        "group": group,
        "partial": "",
        "codons": None,
        "ka": None,
        "ks": None,
        "ka_ks": None,
        "sites_n": None,
        "sites_s": None,
        "status": "",
        "detail": "",
    }
    if name not in reference_genes:
        row["status"] = STATUS_ABSENT_REFERENCE
        row["detail"] = f"{name} is not a CDS of the reference record"
        return row
    row["partial"] = "yes" if reference_genes[name].partial else "no"
    if name not in genes:
        row["status"] = STATUS_ABSENT_GENOME
        row["detail"] = f"{genome_key} does not carry {name}"
        return row

    reference_gene, query_gene = reference_genes[name], genes[name]
    for gene, label in ((reference_gene, "reference"), (query_gene, genome_key)):
        if gene.codon_start != 1:
            row["status"] = STATUS_CODON_START
            row["detail"] = f"{label}: /codon_start={gene.codon_start}"
            return row

    sequence_a = _coding(reference, reference_gene)
    sequence_b = _coding(genome, query_gene)
    for sequence, label in ((sequence_a, "reference"), (sequence_b, genome_key)):
        if len(sequence) % 3:
            row["status"] = STATUS_NOT_WHOLE_CODONS
            row["detail"] = f"{label}: {len(sequence)} bases"
            return row
        protein = translate(sequence)
        if "*" in protein[:-1]:
            row["status"] = STATUS_INTERNAL_STOP
            row["detail"] = f"{label}: stop at codon {protein.index('*') + 1}"
            return row

    alignment = align(
        [("reference", sequence_a), (genome_key, sequence_b)],
        ctx,
        strategy=strategy,
    )
    breaks = codon_frame_breaks(alignment)
    if breaks:
        row["status"] = STATUS_FRAMESHIFT
        row["detail"] = "; ".join(breaks)
        return row

    try:
        result = run_yn00(
            "reference",
            alignment.sequences[0],
            genome_key,
            alignment.sequences[1],
            ctx,
            binary=binary,
        )
    except (ToolError, ValueError) as error:
        # One gene yn00 cannot score is a finding about that gene; a missing or unpinned tool
        # was already raised above and stopped the run.
        row["status"] = STATUS_TOOL_FAILED
        row["detail"] = str(error)[:300]
        return row

    row["codons"] = result.codons
    row["ka"] = result.dN
    row["ks"] = result.dS
    row["sites_n"] = result.sites_nonsynonymous
    row["sites_s"] = result.sites_synonymous
    if result.dS == 0:
        row["status"] = STATUS_KS_ZERO
        row["detail"] = "Ks = 0 — the ratio is undefined"
    else:
        row["ka_ks"] = result.dN / result.dS
        row["status"] = STATUS_OK
    return row


def _coding(genome: Genome, gene: Any) -> str:
    """The CDS without its terminal stop codon — yn00 stops on a stop codon, so it is trimmed
    here rather than repaired after the tool refuses the pair."""
    sequence = coding_sequence(genome, gene)
    if sequence[-3:].upper() in STOPS:
        return sequence[:-3]
    return sequence


def _pool(
    genome_key: str,
    group: str,
    scored: list[tuple[float, float, float | None, float | None]],
    n_shared: int,
) -> dict[str, Any]:
    """Site-weighted pooled Ka and Ks over the genes that produced values."""
    weighted_n = [(ka, n) for ka, _ks, n, _s in scored if n]
    weighted_s = [(ks, s) for _ka, ks, _n, s in scored if s]
    pooled_ka = (
        sum(ka * n for ka, n in weighted_n) / sum(n for _ka, n in weighted_n)
        if weighted_n
        else None
    )
    pooled_ks = (
        sum(ks * s for ks, s in weighted_s) / sum(s for _ks, s in weighted_s)
        if weighted_s
        else None
    )
    ratio = pooled_ka / pooled_ks if pooled_ka is not None and pooled_ks not in (None, 0) else None
    return {
        "genome": genome_key,
        "group": group,
        "n_shared": n_shared,
        "n_computed": len(scored),
        "n_not_computed": n_shared - len(scored),
        "pooled_ka": pooled_ka,
        "pooled_ks": pooled_ks,
        "pooled_ka_ks": ratio,
    }
