"""Substitutions — pairwise counts against a reference genome (PROJECT.md §13.2, ROADMAP M2.6).

Every count is stated against the genome it was measured from: row 0 of the MAFFT alignment is
the reference, and each other row is compared to it. Columns where either row is a gap or an
ambiguous base are skipped and reported, never scored as a match or a substitution.

Scope matters more here than anywhere: §6's golden Ts/Tv of 2.686 was measured over **CDS only,
76 shared genes, IR excluded**, which is narrower than this module's default whole-alignment
comparison. The two must never be quoted as if they were the same number. Ts/Tv is left blank
when there are no transversions — it is undefined, not infinite.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from plastid.core.align import Alignment, reference_positions, substitution_counts
from plastid.core.gene_names import load_gene_aliases, normalize_genome
from plastid.core.ir import detect_ir
from plastid.core.model import Genome, cds_features, coding_sequence
from plastid.core.normalize import normalize
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "substitutions"
MODULE_VERSION = "0.1.0"


# The two scopes this module can report, and why the difference is not cosmetic: §6's golden
# Ts/Tv (2.686) was measured over CDS only, the 76 genes the two genomes share, IR excluded —
# 66,066 alignment columns. The whole-alignment default is a different number over different
# bases, and the two must never be quoted as one another (PROJECT.md §6, §9.B-13).
SCOPES: tuple[str, ...] = ("whole", "cds")
DEFAULT_SCOPE = "whole"


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="Substitutions",
    requires_ir=True,
    arity=2,
    citations=("mafft",),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    reference = genomes[0]
    scope = str(params.get("scope", DEFAULT_SCOPE))
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}, got {scope!r}")
    if scope == "cds":
        return _cds_scope(genomes, params, ctx)

    alignment = align(
        [(genome.accession or genome.name, genome.sequence) for genome in genomes], ctx
    )
    frame = normalize(reference, detect_ir(reference))
    positions = reference_positions(alignment)

    labels: dict[str, str] = params.get("groups", {})
    rows = []
    for index in range(1, alignment.n_rows):
        counts = substitution_counts(alignment.sequences[0], alignment.sequences[index])
        rows.append(
            {
                "genome": alignment.names[index],
                "group": labels.get(alignment.names[index], ""),
                **_ratio(counts),
            }
        )

    regions = []
    for region in frame.regions:
        columns = tuple(
            column
            for column, position in enumerate(positions)
            if position is not None and region.start <= position < region.end
        )
        totals = {"compared": 0, "substitutions": 0, "transitions": 0, "transversions": 0}
        for index in range(1, alignment.n_rows):
            counts = _counts_in_columns(alignment, index, columns)
            for key in totals:
                totals[key] += counts[key]
        regions.append({"region": region.name, **_ratio(totals)})

    whole = {
        key: sum(row[key] for row in rows if isinstance(row.get(key), int))
        for key in ("compared", "substitutions", "transitions", "transversions")
    }

    summary = Table(
        id="substitutions_summary",
        title="Substitutions summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes compared", dtype="int"),
            Column("alignment_length", "Alignment length", dtype="int", unit="bp"),
            Column("compared_columns", "Comparable columns (all comparisons)", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference.accession or "",
                "organism": reference.name,
                "n_genomes": alignment.n_rows,
                "alignment_length": alignment.length,
                "compared_columns": whole["compared"],
                "substitutions": whole["substitutions"],
                "ts_tv": _ratio(whole)["ts_tv"],
                "tool": "MAFFT",
            },
        ),
        scope=Scope(
            label=(
                "whole aligned plastome against the reference (row one of the alignment); "
                "gap and ambiguous columns skipped; not the CDS-only scope of PROJECT.md §6"
            ),
            coord_system="normalized",
        ),
    )

    table = Table(
        id="substitutions",
        title="Substitutions against the reference",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("compared", "Comparable columns", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=tuple(rows),
        scope=Scope(
            label=(
                "each genome against the reference; gap and ambiguous columns skipped; "
                "the group label comes from the sample-group file given to the run"
            ),
            coord_system="normalized",
        ),
        key="genome",
    )

    region_table = Table(
        id="substitutions_by_region",
        title="Substitutions by region",
        columns=(
            Column("region", "Region"),
            Column("compared", "Comparable columns", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=tuple(regions),
        scope=Scope(
            label=(
                "the reference genome's normalized regions mapped onto alignment columns; "
                "summed over every comparison"
            ),
            coord_system="normalized",
        ),
        key="region",
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "substitutions_summary": summary,
            "substitutions": table,
            "substitutions_by_region": region_table,
        },
        notes=(
            "Ts/Tv is blank when no transversion was found: it is undefined, not infinite.",
            "§6's golden Ts/Tv has the narrower scope CDS only · 76 shared genes · IR excluded.",
        ),
        citations=("mafft",),
    )


def _cds_scope(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    """§6's method: the shared CDS, one copy per gene, IR excluded, aligned gene by gene.

    The numbers this reproduces are the project's own golden ones — 317 transitions, 118
    transversions, 2.686 over 76 genes and 66,066 columns — so this is the scope in which the
    toolkit's Ts/Tv is comparable to the published table at all.
    """
    reference = genomes[0]
    aliases = load_gene_aliases()
    ref_genes = _by_name(reference, aliases)
    shared = [
        name
        for name in sorted(ref_genes)
        if all(_by_name(g, aliases).get(name) for g in genomes[1:])
    ]

    totals = {"compared": 0, "substitutions": 0, "transitions": 0, "transversions": 0}
    per_gene: list[dict[str, Any]] = []
    for name in shared:
        rows = [
            (
                genome.accession or genome.name,
                coding_sequence(genome, _by_name(genome, aliases)[name]),
            )
            for genome in genomes
        ]
        alignment = align(rows, ctx)
        counts = {"compared": 0, "substitutions": 0, "transitions": 0, "transversions": 0}
        for index in range(1, alignment.n_rows):
            for key, value in substitution_counts(
                alignment.sequences[0], alignment.sequences[index]
            ).items():
                if key in counts:
                    counts[key] += value
        for key in totals:
            totals[key] += counts[key]
        per_gene.append({"gene": name, "columns": alignment.length, **_ratio(counts)})

    summary_scope = Scope(
        label=(
            "CDS only; the genes both genomes share, one copy each (IR excluded); each gene "
            "aligned on its own — the scope PROJECT.md §6 measured 2.686 in"
        ),
        include=("CDS",),
        exclude_regions=("IRa", "IRb"),
        coord_system="normalized",
    )
    summary = Table(
        id="substitutions_summary",
        title="Substitutions summary (CDS scope)",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes compared", dtype="int"),
            Column("n_genes", "Shared genes", dtype="int"),
            Column("alignment_length", "Alignment columns", dtype="int", unit="bp"),
            Column("compared_columns", "Comparable columns", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference.accession or "",
                "organism": reference.name,
                "n_genomes": len(genomes),
                "n_genes": len(shared),
                "alignment_length": sum(row["columns"] for row in per_gene),
                "compared_columns": totals["compared"],
                "substitutions": totals["substitutions"],
                "transitions": totals["transitions"],
                "transversions": totals["transversions"],
                "ts_tv": _ratio(totals)["ts_tv"],
                "tool": "MAFFT",
            },
        ),
        scope=summary_scope,
    )
    genes = Table(
        id="substitutions_by_gene",
        title="Substitutions by gene (CDS scope)",
        columns=(
            Column("gene", "Gene"),
            Column("columns", "Alignment columns", dtype="int", unit="bp"),
            Column("compared", "Comparable columns", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=tuple(per_gene),
        scope=summary_scope,
        key="gene",
    )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"substitutions_summary": summary, "substitutions_by_gene": genes},
        notes=(
            "CDS scope: the genes both genomes share, one copy each (the inverted repeat is "
            "therefore counted once), aligned gene by gene — the scope PROJECT.md §6 fixed.",
            "A gene one genome does not carry is not in this table: this scope compares shared "
            "coding sequence, and a missing gene is a different finding (Phase 3's).",
            "Ts/Tv is blank when no transversion was found: it is undefined, not infinite.",
        ),
        citations=("mafft",),
    )


def _by_name(genome: Genome, aliases) -> dict[str, Any]:
    """The genome's CDS by canonical name, one feature per name (the IR copy is the same gene)."""
    names = normalize_genome(genome, aliases)
    by_name: dict[str, Any] = {}
    for gene in cds_features(genome):
        by_name.setdefault(names.normalized(gene).canonical, gene)
    return by_name


def _counts_in_columns(alignment: Alignment, row: int, columns: tuple[int, ...]) -> dict[str, int]:
    reference = "".join(alignment.sequences[0][column] for column in columns)
    query = "".join(alignment.sequences[row][column] for column in columns)
    return substitution_counts(reference, query)


def _ratio(counts: dict[str, int]) -> dict[str, Any]:
    transversions = counts["transversions"]
    return {
        "compared": counts["compared"],
        "substitutions": counts["substitutions"],
        "transitions": counts["transitions"],
        "transversions": transversions,
        "ts_tv": (counts["transitions"] / transversions) if transversions else None,
    }
