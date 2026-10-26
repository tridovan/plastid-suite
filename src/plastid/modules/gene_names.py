"""Gene-name check across a set of records (M3.2, P3-02).

One record can be normalized alone; a *panel* can do one thing more: a bare tRNA name
(``trnK``, no anticodon) is resolved by 100 % sequence identity against the named copies in the
other records. Measured 2026-10-03: reading the anticodon at a fixed offset matched only 106 of
293 named tRNAs, so sequence identity is the only honest way to do it — and when it does not
resolve, the name is reported, never guessed.

What the module reports is deliberately narrow: the names that *changed* (syntax), the names
that *are one locus under two nomenclatures* (semantic), and the names that could not be
resolved at all. A name that needed nothing is not repeated back — the summary counts it.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from typing import Any

from plastid.core.gene_names import (
    IDENTICAL_RULE,
    SEMANTIC_RULE,
    SEQUENCE_RULE,
    SYNTAX_RULE,
    GenomeNames,
    load_gene_aliases,
    normalize_genome,
    references_from,
)
from plastid.core.model import GENE_KINDS, Gene, Genome
from plastid.core.organelle import GROUP_SHARED
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "gene_names"
MODULE_VERSION = "0.1.0"

SCOPE_LABEL = "every counted gene of every record in the run"

_RULE_MEANING = {
    IDENTICAL_RULE: "already canonical",
    SYNTAX_RULE: "the same name, written differently",
    SEMANTIC_RULE: "one locus under two nomenclatures (D3-02)",
    SEQUENCE_RULE: "bare tRNA name resolved by 100 % sequence identity",
}


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Gene-name check & normalization",
    requires_ir=False,  # names are not a position-based analysis
    arity=2,  # two records minimum; the pool is what resolves a bare name
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    aliases = load_gene_aliases()
    references = references_from(genomes, aliases)
    normalized = [_genome_names(genome, aliases, references) for genome in genomes]
    keys = [_key(genome, index) for index, genome in enumerate(genomes, start=1)]
    _reject_repeated_keys(keys)

    panel = Table(
        id="panel_summary",
        title="Gene names across the run",
        columns=(
            Column("n_records", "Records", dtype="int"),
            Column("n_genes", "Counted genes", dtype="int"),
            Column("n_distinct_names", "Distinct names as they arrived", dtype="int"),
            Column("n_distinct_canonical", "Distinct canonical names", dtype="int"),
            Column("n_renamed", "Written differently", dtype="int"),
            Column("n_aliased", "One locus, two nomenclatures", dtype="int"),
            Column("n_resolved_by_sequence", "Resolved by sequence", dtype="int"),
            Column("n_unresolved", "Unresolved", dtype="int"),
        ),
        rows=(_panel_row(genomes, normalized),),
        scope=Scope(label=SCOPE_LABEL),
    )

    summary = Table(
        id="name_summary",
        title="Gene names by record",
        columns=(
            Column("genome", "Record"),
            Column("organism", "Organism"),
            Column("n_genes", "Counted genes", dtype="int"),
            Column("n_distinct_names", "Distinct names as they arrived", dtype="int"),
            Column("n_distinct_canonical", "Distinct canonical names", dtype="int"),
            Column("n_renamed", "Written differently", dtype="int"),
            Column("n_aliased", "Two nomenclatures", dtype="int"),
            Column("n_resolved_by_sequence", "Resolved by sequence", dtype="int"),
            Column("n_unresolved", "Unresolved", dtype="int"),
        ),
        rows=tuple(
            _summary_row(key, genome, names)
            for key, genome, names in zip(keys, genomes, normalized, strict=True)
        ),
        scope=Scope(label=SCOPE_LABEL),
        key="genome",
    )

    tables: dict[str, Table] = {panel.id: panel, summary.id: summary}
    mappings = _mapping_rows(keys, genomes, normalized)
    if mappings:
        tables["name_mappings"] = Table(
            id="name_mappings",
            title="Names that were normalized (or carry a synonym)",
            columns=(
                Column("name_id", "Name id"),
                Column("genome", "Record"),
                Column("name", "Name as it arrived"),
                Column("canonical", "Canonical name"),
                Column("printed", "Printed as"),
                Column("rule", "Why"),
                Column("product", "Product in the record"),
            ),
            rows=mappings,
            scope=Scope(label=SCOPE_LABEL),
            key="name_id",
        )
    unresolved = _unresolved_rows(keys, genomes, normalized)
    if unresolved:
        tables["unresolved_names"] = Table(
            id="unresolved_names",
            title="Names that could not be resolved — reported, not guessed",
            columns=(
                Column("name_id", "Name id"),
                Column("genome", "Record"),
                Column("name", "Name as it arrived"),
                Column("product", "Product in the record"),
                Column("reason", "Why it is unresolved"),
            ),
            rows=unresolved,
            scope=Scope(label=SCOPE_LABEL),
            key="name_id",
        )

    notes = [
        "A name is only changed by the normalizer in core/gene_names.py — no module renames one.",
        "A canonical name is the one this toolkit compares and prints; the other name of the "
        "same locus is printed in parentheses beside it (D3-02).",
    ]
    total_unresolved = sum(len(names.unresolved) for names in normalized)
    if total_unresolved:
        notes.append(
            f"{total_unresolved} name(s) could not be resolved and are listed in "
            f"unresolved_names; they keep the name they arrived with."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
    )


def _key(genome: Genome, index: int) -> str:
    return genome.accession or genome.name or f"record {index}"


def _reject_repeated_keys(keys: Sequence[str]) -> None:
    """Every row of this report is addressed by its record's key, so two rows cannot share one.

    Two records sharing an accession means the same file was passed twice — a panel report of
    one genome against itself answers nothing, and the Table would otherwise raise a message
    about columns rather than about the input.
    """
    seen: set[str] = set()
    for key in keys:
        if key in seen:
            raise ValueError(
                f"two records in this run share the key {key!r} — every row of the panel report "
                f"is addressed by accession, so pass each record once"
            )
        seen.add(key)


def _genome_names(genome: Genome, aliases, references) -> GenomeNames:
    return normalize_genome(genome, aliases, references=references)


def _counted(genome: Genome) -> list[Gene]:
    return [gene for gene in genome.genes if gene.kind in GENE_KINDS]


def _panel_row(genomes: Sequence[Genome], normalized: Sequence[GenomeNames]) -> dict[str, Any]:
    """The headline numbers, counted over *distinct names* — one name counted once, not once
    per record that used it, so the row answers "how many names needed work", not "how many
    gene copies are in the run"."""
    by_name: dict[str, str] = {}  # arriving name -> rule
    canonical: set[str] = set()
    unresolved: Counter[str] = Counter()
    for genome, names in zip(genomes, normalized, strict=True):
        for gene in _counted(genome):
            entry = names.normalized(gene)
            by_name.setdefault(gene.name, entry.rule)
            canonical.add(entry.canonical)
    for names in normalized:
        for item in names.unresolved:
            unresolved[item.name] += 1
    rules = Counter(by_name.values())
    return {
        "n_records": len(genomes),
        "n_genes": sum(len(_counted(genome)) for genome in genomes),
        "n_distinct_names": len(by_name),
        "n_distinct_canonical": len(canonical),
        "n_renamed": rules[SYNTAX_RULE],
        "n_aliased": rules[SEMANTIC_RULE],
        "n_resolved_by_sequence": rules[SEQUENCE_RULE],
        "n_unresolved": len(unresolved),
    }


def _summary_row(key: str, genome: Genome, names: GenomeNames) -> dict[str, Any]:
    counted = _counted(genome)
    rules = Counter(names.normalized(gene).rule for gene in counted)
    return {
        "genome": key,
        "organism": genome.name,
        "n_genes": len(counted),
        "n_distinct_names": len({gene.name for gene in counted}),
        "n_distinct_canonical": len({names.canonical(gene) for gene in counted}),
        "n_renamed": rules[SYNTAX_RULE],
        "n_aliased": rules[SEMANTIC_RULE],
        "n_resolved_by_sequence": rules[SEQUENCE_RULE],
        "n_unresolved": len(names.unresolved),
    }


def _mapping_rows(
    keys: Sequence[str], genomes: Sequence[Genome], normalized: Sequence[GenomeNames]
) -> tuple[dict[str, Any], ...]:
    """One row per (record, name) that was rewritten or that has a second name."""
    rows: list[dict[str, Any]] = []
    for key, genome, names in zip(keys, genomes, normalized, strict=True):
        seen: set[str] = set()
        for gene in _counted(genome):
            entry = names.normalized(gene)
            if gene.name in seen:
                continue
            if not entry.changed and not entry.synonyms:
                continue
            seen.add(gene.name)
            rows.append(
                {
                    "name_id": f"{key}:{gene.name}",
                    "genome": key,
                    "name": gene.name,
                    "canonical": entry.canonical,
                    "printed": entry.label,
                    "rule": entry.rule,
                    "product": gene.product or "",
                }
            )
    return tuple(rows)


def _unresolved_rows(
    keys: Sequence[str], genomes: Sequence[Genome], normalized: Sequence[GenomeNames]
) -> tuple[dict[str, Any], ...]:
    rows: list[dict[str, Any]] = []
    for key, names in zip(keys, normalized, strict=True):
        for item in names.unresolved:
            rows.append(
                {
                    "name_id": f"{key}:{item.name}",
                    "genome": key,
                    "name": item.name,
                    "product": item.gene.product or "",
                    "reason": item.reason,
                }
            )
    return tuple(rows)
