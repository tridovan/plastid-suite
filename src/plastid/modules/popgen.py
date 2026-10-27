"""Population genetics per gene and per group — DnaSP-style, pure Python (P4-03, M4.4; D4-03).

One gene at a time: the module aligns the CDS across every input genome that carries it, splits
the alignment by the run's sample groups, and reports each group's statistics beside the pairwise
Fst between groups. Nothing here needs an external estimator — the arithmetic is in
:mod:`plastid.core.popgen`, so the minimum install runs it (CLAUDE.md).

What this dataset *is* is part of the result. The groups a run is given are species or genera,
not populations of one species, so:

* an **Fst between groups** here measures divergence between taxa, not differentiation within a
  species — the table says so rather than letting the number travel under the population
  meaning;
* **Tajima's D needs at least four sequences** (at n ≤ 3 its variance is identically zero), so a
  group of two or three gets a blank with the reason, never a number;
* every statistic is published with its **n** and its usable site count, because a π from three
  sequences and a π from thirty are different measurements.

Gaps and ambiguous bases are dropped column-wise (complete deletion), the same rule the rest of
the toolkit uses.
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations
from typing import Any

from plastid.core.align import Alignment
from plastid.core.gene_names import cds_by_canonical_name
from plastid.core.model import Genome, coding_sequence
from plastid.core.organelle import GROUP_SHARED
from plastid.core.popgen import fst_gst, summarise
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "popgen"
MODULE_VERSION = "0.2.0"  # behaviour changed 2026-10-04 after the first runs


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Population genetics",
    requires_ir=False,  # gene-based; the mito half of the phase has no IR to gate on
    arity=2,
    citations=("mafft",),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    labels: dict[str, str] = dict(params.get("groups") or {})
    keys = [genome.accession or genome.name for genome in genomes]
    if not labels:
        raise ValueError(
            "popgen compares groups and no groups were given — pass --groups with a sample-group "
            "file, or run a module that does not need groups"
        )
    unlabelled = [key for key in keys if key not in labels]
    if unlabelled:
        raise ValueError(
            f"every genome needs a group, {', '.join(unlabelled)} has none — an unlabelled "
            f"sequence would be silently left out of every group statistic"
        )
    groups = sorted({labels[key] for key in keys})

    reference = genomes[0]
    # The record's own organelle picks the vocabulary (core.gene_names.aliases_for).
    by_genome = {
        key: cds_by_canonical_name(genome) for key, genome in zip(keys, genomes, strict=True)
    }
    names = sorted({name for mapping in by_genome.values() for name in mapping})

    group_rows: list[dict[str, Any]] = []
    fst_rows: list[dict[str, Any]] = []
    aligned_genes = 0
    analysed: set[str] = set()
    skipped: list[str] = []
    for name in names:
        carriers = [
            (key, genome, by_genome[key][name])
            for key, genome in zip(keys, genomes, strict=True)
            if name in by_genome[key]
        ]
        if len(carriers) < 2:
            skipped.append(name)
            continue
        alignment = align(
            [(key, coding_sequence(genome, gene)) for key, genome, gene in carriers], ctx
        )
        aligned_genes += 1
        analysed.add(name)
        subsets: dict[str, Alignment] = {}
        for group in groups:
            indices = [
                index
                for index, (key, _genome, _gene) in enumerate(carriers)
                if labels[key] == group
            ]
            if not indices:
                continue
            subsets[group] = Alignment(
                names=tuple(carriers[index][0] for index in indices),
                sequences=tuple(alignment.sequences[index] for index in indices),
            )
            summary = summarise(subsets[group])
            group_rows.append(
                {
                    "row_id": f"{name}#{group}",
                    "gene": name,
                    "group": group,
                    "n": summary.n,
                    "sites": summary.sites,
                    "segregating": summary.segregating,
                    "singletons": summary.singletons,
                    "haplotypes": summary.haplotypes,
                    "hd": summary.haplotype_diversity,
                    "pi": summary.pi_per_site,
                    "theta_w": summary.theta_w_per_site,
                    "tajimas_d": summary.tajimas_d,
                    "fu_li_d_star": summary.fu_li_d_star,
                    "fu_li_f_star": summary.fu_li_f_star,
                    "note": "; ".join(summary.notes),
                }
            )
        for first, second in combinations(sorted(subsets), 2):
            value = fst_gst((subsets[first], subsets[second]))
            fst_rows.append(
                {
                    "row_id": f"{name}#{first}#{second}",
                    "gene": name,
                    "group_a": first,
                    "group_b": second,
                    "n_a": subsets[first].n_rows,
                    "n_b": subsets[second].n_rows,
                    "fst": value,
                    "estimator": "Nei's Gst",
                    "note": (
                        ""
                        if value is not None
                        else "not computable — either a group has fewer than two sequences "
                        "(within-group diversity is not estimable) or there is no diversity to "
                        "partition between the groups"
                    ),
                }
            )

    if not group_rows:
        raise ValueError(
            "popgen found no gene carried by at least two genomes — nothing is estimable; "
            "an unannotated record cannot contribute to a gene-by-gene comparison (each input "
            "needs CDS, or the module has no unit to analyse)"
        )
    sizes = ", ".join(
        f"{group}: {sum(1 for key in keys if labels[key] == group)}" for group in groups
    )
    summary_scope = Scope(
        label=(
            "per-gene statistics by sample group; sites are columns every sequence in the group "
            "carries unambiguously (complete deletion)"
        ),
        coord_system="genbank",
    )
    fst_scope = Scope(
        label=(
            "Nei's Gst between groups, per gene; the groups of this dataset are taxa, not "
            "populations, so this is a divergence statistic"
        ),
        coord_system="genbank",
    )

    summary = Table(
        id="popgen_summary",
        title="Population-genetics summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes", dtype="int"),
            Column("n_groups", "Groups", dtype="int"),
            Column("group_sizes", "Group sizes"),
            Column("n_genes", "Genes analysed", dtype="int"),
            Column("n_partial_genes", "Analysed genes partial in the reference", dtype="int"),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference.accession or "",
                "organism": reference.name,
                "n_genomes": len(genomes),
                "n_groups": len(groups),
                "group_sizes": sizes,
                "n_genes": aligned_genes,
                "n_partial_genes": sum(
                    1
                    for name in analysed
                    if name in by_genome[keys[0]] and by_genome[keys[0]][name].partial
                ),
                "tool": "MAFFT",
            },
        ),
        scope=summary_scope,
    )

    group_table = Table(
        id="popgen_groups",
        title="Statistics by gene and group",
        columns=(
            Column("row_id", "Row"),
            Column("gene", "Gene"),
            Column("group", "Group"),
            Column("n", "Sequences", dtype="int"),
            Column("sites", "Usable sites", dtype="int", unit="bp"),
            Column("segregating", "Segregating sites", dtype="int"),
            Column("singletons", "Singletons", dtype="int"),
            Column("haplotypes", "Haplotypes", dtype="int"),
            Column("hd", "Haplotype diversity", dtype="float", precision=4),
            Column("pi", "π per site", dtype="float", precision=8),
            Column("theta_w", "θ(S) per site", dtype="float", precision=8),
            Column("tajimas_d", "Tajima's D", dtype="float", precision=4),
            Column("fu_li_d_star", "Fu & Li's D*", dtype="float", precision=4),
            Column("fu_li_f_star", "Fu & Li's F*", dtype="float", precision=4),
            Column("note", "Note"),
        ),
        rows=tuple(group_rows),
        scope=summary_scope,
        key="row_id",
    )

    fst_table = Table(
        id="popgen_fst",
        title="Fst between groups",
        columns=(
            Column("row_id", "Row"),
            Column("gene", "Gene"),
            Column("group_a", "Group A"),
            Column("group_b", "Group B"),
            Column("n_a", "n (A)", dtype="int"),
            Column("n_b", "n (B)", dtype="int"),
            Column("fst", "Fst", dtype="float", precision=4),
            Column("estimator", "Estimator"),
            Column("note", "Note"),
        ),
        rows=tuple(fst_rows),
        scope=fst_scope,
        key="row_id",
    )

    notes = [
        "Estimators are written here from their published formulas; no "
        "external tool is needed, so the minimum install runs this module.",
        "Tajima's D and Fu & Li's D*/F* need at least four sequences — at n ≤ 3 their variance "
        "is identically zero — and a blank cell carries that reason in the Note column.",
        "The groups of this run are taxa, not populations of one species: an Fst between them "
        "measures divergence, and a neutrality test over a handful of sequences is thin — both "
        "are reported with their n.",
        "Fst here is Nei's Gst, which reads haplotype frequencies: genes whose sequences are all "
        "distinct give the same value however far apart the sequences are — how much divergence "
        "there is lives in π, not in this column.",
    ]
    if skipped:
        notes.append(
            f"{len(skipped)} gene(s) are carried by fewer than two genomes and cannot be "
            f"compared: {', '.join(skipped)}."
        )
    empty = [key for key in keys if not by_genome[key]]
    if empty:
        notes.append(
            f"{len(empty)} record(s) carry no CDS and contribute nothing to any gene: "
            f"{', '.join(empty)}."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "popgen_summary": summary,
            "popgen_groups": group_table,
            "popgen_fst": fst_table,
        },
        notes=tuple(notes),
        citations=("mafft",),
    )
