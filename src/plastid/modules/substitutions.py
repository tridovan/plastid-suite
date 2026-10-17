"""Substitutions — pairwise counts against a reference genome (PROJECT.md §13.2, ROADMAP M2.6).

Every count is stated against the genome it was measured from: row 0 of the MAFFT alignment is
the reference, and each other row is compared to it. Columns where either row is a gap or an
ambiguous base are skipped and reported, never scored as a match or a substitution.

Scope matters more here than anywhere: §6's golden Ts/Tv was measured over **CDS only, shared
genes, IR excluded** — 2.6917 over 79 genes since the shared set was re-measured by locus
(2026-10-04) — which is narrower than this module's default whole-alignment comparison. The two
must never be quoted as if they were the same number. Ts/Tv is left blank when there are no
transversions — it is undefined, not infinite.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from plastid.core.align import (
    SUBSTITUTION_TYPES,
    Alignment,
    reference_positions,
    substitution_class,
    substitution_counts,
)
from plastid.core.gene_names import cds_by_canonical_name
from plastid.core.ir import detect_ir
from plastid.core.model import Genome, coding_sequence
from plastid.core.normalize import normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "substitutions"
MODULE_VERSION = (
    "0.3.0"  # 0.3.0 adds substitutions_by_type; 0.2.0 changed behaviour 2026-10-04 (§9.C-35)
)


# The two scopes this module can report, and why the difference is not cosmetic: §6's golden
# Ts/Tv is measured over CDS only, the genes the genomes share, IR excluded. The whole-alignment
# default is a different number over different bases, and the two must never be quoted as one
# another (PROJECT.md §6, §8.25, §9.B-14). The label says *which scope*, not *what number* — a
# number written into a label goes stale the moment the shared set is re-measured.
SCOPES: tuple[str, ...] = ("whole", "cds")
DEFAULT_SCOPE = "whole"


@module(
    MODULE_ID,
    group=GROUP_SHARED,
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
    type_rows: list[dict[str, Any]] = []
    for index in range(1, alignment.n_rows):
        counts = substitution_counts(alignment.sequences[0], alignment.sequences[index])
        rows.append(
            {
                "genome": alignment.names[index],
                "group": labels.get(alignment.names[index], ""),
                **_ratio(counts),
            }
        )
        type_rows.extend(
            _type_rows(alignment.names[index], labels.get(alignment.names[index], ""), counts)
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
            # The report's prose prints this cell rather than asserting a scope of its own:
            # measured 10/10/2026 (P2B-03), the template said "this is the whole aligned
            # plastome" over a run made with --set scope=cds, so the document described
            # numbers it was not showing.
            Column("scope", "Scope (this run)"),
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
                "scope": "the whole aligned plastome against the reference, gap columns skipped",
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
            "substitutions_by_type": _by_type_table(
                type_rows,
                Scope(
                    label=(
                        "each genome against the reference; every one of the twelve directed "
                        "changes is a row, a zero included; gap and ambiguous columns skipped"
                    ),
                    coord_system="normalized",
                ),
            ),
        },
        notes=(
            "Ts/Tv is blank when no transversion was found: it is undefined, not infinite.",
            "§6's golden Ts/Tv has the narrower scope CDS only · 76 shared genes · IR excluded.",
            "substitutions_by_type holds the twelve directed changes the transitions and "
            "transversions totals are summed from; a manuscript's six complementary pairs "
            "(A↔G, C↔T, A↔C, T↔G, C↔G, A↔T) are two rows each.",
        ),
        citations=("mafft",),
    )


def _cds_scope(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    """§6's method: the shared CDS, one copy per gene, IR excluded, aligned gene by gene.

    This is the scope PROJECT.md §6 fixed for the golden ratio — CDS only, shared genes, IR
    excluded — and the scope in which the toolkit's Ts/Tv is comparable to that table at all.
    The numbers it reproduces today (2026-10-04: 323 Ts / 120 Tv = 2.6917 over 79 shared genes)
    live in PROJECT.md, not here.
    """
    reference = genomes[0]
    # No table is passed: the record's own organelle picks the vocabulary (aliases_for) —
    # measured 2026-10-04, a mitogenome whose CDS are named by product sentence needs the
    # mitochondrial one, and passing the plastid table explicitly found no shared gene at all.
    ref_genes = _by_name(reference, None)
    # A record with no annotation cannot share a CDS with anything — it is excluded from the
    # comparison and NAMED, rather than quietly emptying the table for every other genome.
    # Measured 04/10/2026 on the pinned panel: 3 of its 11 records carry no CDS, and the first
    # version of this scope returned an empty table for the whole run.
    queries = [(genome, _by_name(genome, None)) for genome in genomes[1:]]
    annotated = [(genome, genes) for genome, genes in queries if genes]
    unannotated = [genome.accession or genome.name for genome, genes in queries if not genes]
    if not annotated:
        raise ValueError(
            "the CDS scope has nothing to compare: no genome other than the reference carries "
            "an annotated CDS — use the whole-alignment scope, or pass annotated records"
        )
    shared = [
        name for name in sorted(ref_genes) if all(genes.get(name) for _genome, genes in annotated)
    ]
    if not shared:
        raise ValueError(
            "the CDS scope found no gene shared by the reference and every annotated genome — "
            "the records may be too distant, or their gene names unresolved (see the gene_names "
            "module); the whole-alignment scope is the honest comparison for them"
        )

    totals = {"compared": 0, "substitutions": 0, "transitions": 0, "transversions": 0}
    per_gene: list[dict[str, Any]] = []
    per_query: dict[str, dict[str, int]] = {
        genome.accession or genome.name: {
            "genes": 0,
            "columns": 0,
            "compared": 0,
            "substitutions": 0,
            "transitions": 0,
            "transversions": 0,
        }
        for genome, _genes in annotated
    }
    # The per-type counters live beside `per_query`, not inside it: that dict's keys are the
    # `substitutions_by_genome` columns, and a nested dict there would be a column the table
    # does not declare.
    per_query_types: dict[str, dict[str, int]] = {
        name: dict.fromkeys(SUBSTITUTION_TYPES, 0) for name in per_query
    }
    for name in shared:
        rows = [
            (reference.accession or reference.name, coding_sequence(reference, ref_genes[name])),
            *[
                (genome.accession or genome.name, coding_sequence(genome, genes[name]))
                for genome, genes in annotated
            ],
        ]
        alignment = align(rows, ctx)
        counts = _zero_counts()
        for index in range(1, alignment.n_rows):
            # The counters now carry the twelve directed types as well as the two aggregates, and
            # this scope sums both: a gene's per-type counts are what the type table is built from.
            row_counts = substitution_counts(alignment.sequences[0], alignment.sequences[index])
            for key in counts:
                counts[key] += row_counts[key]
            # The per-pair numbers a manuscript table needs: one gene, one query genome.
            entry = per_query[rows[index][0]]
            entry["genes"] += 1
            entry["columns"] += alignment.length
            for key in ("compared", "substitutions", "transitions", "transversions"):
                entry[key] += row_counts[key]
            for type_name in SUBSTITUTION_TYPES:
                per_query_types[rows[index][0]][type_name] += row_counts[type_name]
        for key in totals:
            totals[key] += counts[key]
        per_gene.append({"gene": name, "columns": alignment.length, **_ratio(counts)})

    summary_scope = Scope(
        label=(
            "CDS only; the genes both genomes share, one copy each (IR excluded); each gene "
            "aligned on its own — the scope PROJECT.md §6 fixed for the golden Ts/Tv"
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
            Column("scope", "Scope (this run)"),
        ),
        rows=(
            {
                "accession": reference.accession or "",
                "organism": reference.name,
                "n_genomes": 1 + len(annotated),
                "n_genes": len(shared),
                "alignment_length": sum(row["columns"] for row in per_gene),
                "compared_columns": totals["compared"],
                "substitutions": totals["substitutions"],
                "transitions": totals["transitions"],
                "transversions": totals["transversions"],
                "ts_tv": _ratio(totals)["ts_tv"],
                "tool": "MAFFT",
                "scope": (
                    "CDS only; the genes both genomes share, one copy each (IR excluded), each "
                    "gene aligned on its own"
                ),
            },
        ),
        scope=summary_scope,
    )
    genomes_table = Table(
        id="substitutions_by_genome",
        title="Substitutions by genome (CDS scope)",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("genes", "Shared genes", dtype="int"),
            Column("columns", "Alignment columns", dtype="int", unit="bp"),
            Column("compared", "Comparable columns", dtype="int"),
            Column("substitutions", "Substitutions", dtype="int"),
            Column("transitions", "Transitions", dtype="int"),
            Column("transversions", "Transversions", dtype="int"),
            Column("ts_tv", "Ts/Tv", dtype="float", precision=4),
        ),
        rows=tuple(
            {
                "genome": key,
                "group": (params.get("groups") or {}).get(key, ""),
                **{name: value for name, value in counts.items()},
                "ts_tv": _ratio(counts)["ts_tv"],
            }
            for key, counts in per_query.items()
        ),
        scope=summary_scope,
        key="genome",
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
    notes = [
        "CDS scope: the genes the genomes share, one copy each (the inverted repeat is "
        "therefore counted once), aligned gene by gene — the scope PROJECT.md §6 fixed.",
        "A gene one genome does not carry is not in this table: this scope compares shared "
        "coding sequence, and a missing gene is a different finding (Phase 3's).",
        "Ts/Tv is blank when no transversion was found: it is undefined, not infinite.",
    ]
    if unannotated:
        notes.append(
            f"{len(unannotated)} record(s) carry no annotated CDS and cannot take part in a "
            f"CDS comparison; they are left out and named: {', '.join(unannotated)}."
        )
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "substitutions_summary": summary,
            "substitutions_by_genome": genomes_table,
            "substitutions_by_gene": genes,
            "substitutions_by_type": _by_type_table(
                [
                    row
                    for key, counts in per_query_types.items()
                    for row in _type_rows(key, (params.get("groups") or {}).get(key, ""), counts)
                ],
                summary_scope,
            ),
        },
        notes=tuple(notes),
        citations=("mafft",),
    )


def _zero_counts() -> dict[str, int]:
    """A fresh accumulator: the four aggregates plus every directed type, all at zero."""
    return {
        "compared": 0,
        "substitutions": 0,
        "transitions": 0,
        "transversions": 0,
        **dict.fromkeys(SUBSTITUTION_TYPES, 0),
    }


def _type_rows(genome: str, group: str, counts: dict[str, int]) -> list[dict[str, Any]]:
    """One row per directed substitution type for one genome — zero counts included.

    A zero is a measurement ("this pair never occurred in the comparison"), so every one of the
    twelve is written rather than only the non-zero ones: a reader must not have to tell an
    absent row from a zero.
    """
    return [
        {
            "row_id": f"{genome}|{type_name}",
            "genome": genome,
            "group": group,
            "substitution_type": type_name,
            "class": substitution_class(type_name),
            "count": counts[type_name],
        }
        for type_name in SUBSTITUTION_TYPES
    ]


def _by_type_table(rows: list[dict[str, Any]], scope: Scope) -> Table:
    """The ``substitutions_by_type`` table — one shape in both scopes, one rule for a figure."""
    return Table(
        id="substitutions_by_type",
        title="Substitutions by type",
        columns=(
            Column("row_id", "Row"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("substitution_type", "Substitution type"),
            Column("class", "Class"),
            Column("count", "Count", dtype="int"),
        ),
        rows=tuple(rows),
        scope=scope,
        key="row_id",
    )


def _by_name(genome: Genome, aliases) -> dict[str, Any]:
    """The genome's CDS by canonical name — ``core.gene_names`` owns the rule (both IR copies
    are one gene); this name is kept because the module reads it in several places."""
    return cds_by_canonical_name(genome, aliases)


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
