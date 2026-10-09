"""Reference-driven gene naming — the step that makes one vocabulary out of many (D5-09, P3B-04).

Every GenBank that enters the pipeline is written back out named after **the project's naming
reference** — one record, chosen by the end user from the audit report (D5-09 §7.1) — plus the
tables that say what changed and why. Two runs of one module, two expectations (D5-09 §2):

* ``mode=normalize`` on the reference DB: renames happen; the DB comes out in vocabulary A;
* ``mode=verify`` on an annotated sample: **nothing should change** — the sample inherited A
  through annotation, so a rename here means something named a gene on its own.

**A and B never mix.** ``genomes[0]`` is the naming reference (vocabulary **A**, what every result
table prints); ``genomes[1:]`` are the records being normalized (vocabulary **B**). B's spellings
appear in exactly one place — the ``pairing`` table, in their own column — because A is the join
key every later table joins on (D5-09 §3).

**The name proposes, the sequence disposes** (D5-09 §4). Four tests decide, in order: the name
(syntax, then the alias table oriented by the reference), the sequence ranking, the margin between
the best and the runner-up, and the neighbours inside one region. Which test decided is the
``rule`` column, and every gene carries a ``status`` — ``measured`` / ``not measured`` /
``outside vocabulary`` — so "there is no counterpart" and "we could not tell" never look alike
(D5-09 §5).

**What the measurement changed** (PROJECT.md §8.41). The sequence tests can decide for a CDS but
**not** for a tRNA: true and wrong tRNA loci overlap in identity (53.76–100 % against up to
73.97 %), and the ranking prefers another serine tRNA in 26 of 1,182 comparisons. A tRNA is
therefore never renamed on a ranking — its name or its neighbours decide — and the nucleotide
comparison is recorded as evidence rather than used as a verdict.

Everything that decides anything lives in :mod:`plastid.core.gene_names`; this module runs the
searches, arranges the tables and writes the artifact. It is not the place a name is chosen.
"""

from __future__ import annotations

import tempfile
from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

from plastid.core.genbank import to_genbank
from plastid.core.gene_names import (
    OUTCOME_CONFIRMED,
    STATUS_MEASURED,
    STATUS_NOT_MEASURED,
    STATUS_OUTSIDE,
    AliasTable,
    GeneDecision,
    Hit,
    NamingConflictError,
    decide_gene,
    duplicate_names,
    load_gene_aliases,
    neighbours_of,
    nucleotide_identity,
    oriented_aliases,
    protein_fasta_rows,
    vocabulary_of,
)
from plastid.core.ir import IRDetectionError, IRGateError, IRResult, detect_ir, enforce_ir_gate
from plastid.core.model import GENE_KINDS, Gene, Genome, coding_sequence
from plastid.core.organelle import GROUP_PLASTID
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "gene_normalize"
MODULE_VERSION = "0.1.0"

MODE_NORMALIZE = "normalize"
MODE_VERIFY = "verify"
MODES = (MODE_NORMALIZE, MODE_VERIFY)

# The conflict tiers' default: stop. `on_conflict=continue` is the second of two invocations —
# never a prompt, so the step stays runnable under Docker and in scripts (D5-09 §5).
ON_CONFLICT_STOP = "stop"
ON_CONFLICT_CONTINUE = "continue"

SCOPE_LABEL = "one record against the project's naming reference (D5-09, vocabulary A)"

# How many nucleotide candidates the evidence column keeps. The ranking never decides for a tRNA
# (§8.41), so the tail past the first handful is not read by anything.
_NUCLEOTIDE_CANDIDATES = 8


@module(
    MODULE_ID,
    group=GROUP_PLASTID,
    version=MODULE_VERSION,
    title="Reference-driven gene naming (A/B pairing)",
    requires_ir=False,  # names are not a position-based analysis; the IR frame only feeds test 4
    arity=2,  # the naming reference plus at least one record to normalize
    reference_param="naming_reference",  # D5-09 §7.2: recorded, not spelled `reference`
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    reference, records = genomes[0], tuple(genomes[1:])
    mode = str(params.get("mode", MODE_NORMALIZE))
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    on_conflict = str(params.get("on_conflict", ON_CONFLICT_STOP))
    if on_conflict not in (ON_CONFLICT_STOP, ON_CONFLICT_CONTINUE):
        raise ValueError(
            f"on_conflict must be {ON_CONFLICT_STOP!r} or {ON_CONFLICT_CONTINUE!r}, "
            f"got {on_conflict!r}"
        )

    # The alias table's direction follows the naming reference (D5-09 §8): the file says which
    # names are one locus, the reference says which of them the project prints.
    aliases = oriented_aliases(load_gene_aliases(), reference)
    vocabulary = vocabulary_of(reference)
    reference_ir = _frame_of(reference)

    notes: list[str] = []
    for record in records:
        already = duplicate_names(record, _frame_of(record))
        if already:
            notes.append(
                f"{_label(record)} already carries more than one gene under {_names(already)} "
                f"— reported, not stopped on: the record arrived that way, and a duplicate the "
                f"normalization did not create is a finding about the record (D5-09 §4)"
            )
    duplicates = duplicate_names(reference, reference_ir)
    if duplicates:
        notes.append(
            f"the naming reference spells more than one gene with {_names(duplicates)} — a "
            f"pairing for a duplicated name points at the first of them, and a rename that would "
            f"copy the duplicate into a record is a conflict, not a merge (D5-09 §4)"
        )

    queries = _protein_queries(records)
    protein_hits = _search(queries, reference, params, ctx) if queries else {}
    if not protein_hits:
        notes.append(
            "no CDS carried a /translation on both sides, so the sequence tests could not run — "
            "every such gene reports 'not measured' rather than being assumed correct"
        )

    decisions: list[tuple[int, GeneDecision]] = []
    for index, record in enumerate(records):
        decisions.extend(
            (index, decision)
            for decision in _decide_record(
                record,
                index=index,
                reference=reference,
                reference_ir=reference_ir,
                aliases=aliases,
                vocabulary=vocabulary,
                protein_hits=protein_hits,
            )
        )

    conflicts = _conflicts(records, decisions)
    if conflicts and on_conflict != ON_CONFLICT_CONTINUE:
        raise NamingConflictError(
            f"gene_normalize: {len(conflicts)} conflict(s) between the records and "
            f"{_label(reference)} — the names and the sequences disagree, so nothing was "
            f"written. Fix the record, or re-run with --set on_conflict=continue to write the "
            f"output and the conflicts table.\n" + "\n".join(conflicts)
        )
    if conflicts:
        notes.append(
            f"{len(conflicts)} conflict(s) were carried through on on_conflict=continue — the "
            f"conflicts table names each one with its coordinates"
        )

    artifacts = _artifacts(mode, records, decisions, notes)

    tables = {
        "normalize_summary": _summary(mode, reference, records, decisions, conflicts, artifacts),
        "pairing": _pairing(records, decisions),
    }
    differences = _differences(reference, records, decisions)
    if differences is not None:
        tables["differences"] = differences
    if conflicts:
        tables["conflicts"] = _conflicts_table(records, decisions, conflicts)

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        artifacts=artifacts,
    )


# ---------------------------------------------------------------- helpers


def _label(genome: Genome) -> str:
    return genome.accession or genome.name


def _names(names: Sequence[str]) -> str:
    return ", ".join(repr(name) for name in names)


def _location(gene: Gene) -> str:
    start_1b, end_1b = gene.span
    return f"({start_1b}..{end_1b}{'+' if gene.strand >= 0 else '-'})"


def _frame_of(genome: Genome) -> IRResult | None:
    """The record's IR, or None when it has none — test 4 then reports itself not measured.

    The gate is applied, never bypassed: a record whose IR copies disagree is exactly the record
    whose gene order cannot be trusted, so its neighbours are not evidence either. Naming is not a
    position-based analysis (``requires_ir=False``); this is a test's input, not a coordinate.
    """
    try:
        return enforce_ir_gate(detect_ir(genome))
    except (IRGateError, IRDetectionError):
        return None


# ---------------------------------------------------------------- the searches


def _protein_queries(records: Sequence[Genome]) -> tuple[tuple[str, str, int, int], ...]:
    """``(query id, protein, record index, gene index)`` for every CDS that carries a protein.

    The id is positional and carries no name: BLAST cuts an id at the first whitespace, and gene
    names contain spaces in some records ("16S ribosomal RNA"), so a name-based id would arrive
    truncated and pair two genes to one hit.
    """
    rows: list[tuple[str, str, int, int]] = []
    for index, record in enumerate(records):
        for gene_number, gene in enumerate(record.genes):
            if gene.kind == "CDS" and gene.translation:
                rows.append((f"r{index}g{gene_number}", gene.translation, index, gene_number))
    return tuple(rows)


def _search(
    queries: Sequence[tuple[str, str, int, int]],
    reference: Genome,
    params: dict[str, Any],
    ctx: RunContext,
) -> dict[tuple[int, int], tuple[Hit, ...]]:
    """One blastp of every record's proteins against the reference's own — the ranking of test 2.

    One database for the run, not one per record: the same reference proteins have to rank the
    same way for every record, and a rebuild per record would be a second, silently different
    database. The reference's own ``/translation`` is the subject — never a self-translation
    (D5-08 §8.21) — written with the three-field subject key the mito path already uses, with the
    vocabulary field saying ``A``.
    """
    from plastid.external.blast import run_blastp

    rows = protein_fasta_rows(reference, vocabulary="A")
    if not rows:
        return {}
    with tempfile.TemporaryDirectory(prefix="gene_normalize_") as work:
        fasta = Path(work) / "naming_reference.faa"
        fasta.write_text("".join(f">{key}\n{protein}\n" for key, protein in rows), encoding="utf-8")
        result = run_blastp(
            [(query_id, protein) for query_id, protein, _, _ in queries],
            fasta,
            ctx,
            directory=params.get("tools"),
        )
    return {
        (index, gene_number): _by_locus(result.for_query(query_id))
        for query_id, _, index, gene_number in queries
    }


def _by_locus(hits: Sequence[Any]) -> tuple[Hit, ...]:
    """The best hit per **locus**, best first — one row per candidate, not one per feature.

    The reference's inverted repeat carries two copies of a gene under one name, and BLAST returns
    both with nearly equal scores; kept as two hits, the runner-up would be the same locus as the
    winner and the margin of test 3 would read zero on every IR gene (measured 08/10/2026: ycf2
    against the golden came out a *conflict* with the name and the sequence agreeing). The margin
    is between the best match and the best **different** match — which is what "runner-up" means.
    """
    best: dict[str, Hit] = {}
    for hit in hits:
        candidate = Hit(
            locus=hit.subject, identity=hit.identity, coverage=hit.coverage, bits=hit.bits
        )
        current = best.get(candidate.locus)
        if current is None or (candidate.bits or 0) > (current.bits or 0):
            best[candidate.locus] = candidate
    return tuple(sorted(best.values(), key=lambda hit: (-(hit.bits or 0), hit.locus)))


def _nucleotide_hits(gene: Gene, record: Genome, reference: Genome) -> tuple[Hit, ...]:
    """The tRNA/rRNA ranking as evidence — identity only, never a verdict (§8.41)."""
    sequence = coding_sequence(record, gene).upper()
    hits = [
        Hit(
            locus=candidate.name,
            identity=nucleotide_identity(sequence, coding_sequence(reference, candidate).upper()),
            coverage=1.0,
        )
        for candidate in reference.genes
        if candidate.kind == gene.kind
    ]
    return tuple(sorted(hits, key=lambda hit: (-hit.identity, hit.locus))[:_NUCLEOTIDE_CANDIDATES])


# ---------------------------------------------------------------- the decisions


def _decide_record(
    record: Genome,
    *,
    index: int,
    reference: Genome,
    reference_ir: IRResult | None,
    aliases: AliasTable,
    vocabulary: dict[str, Gene],
    protein_hits: dict[tuple[int, int], tuple[Hit, ...]],
) -> list[GeneDecision]:
    record_ir = _frame_of(record)
    decisions: list[GeneDecision] = []
    for gene_number, gene in enumerate(record.genes):
        if gene.kind not in GENE_KINDS:
            continue
        if gene.kind == "CDS":
            hits = protein_hits.get((index, gene_number), ())
        else:
            hits = _nucleotide_hits(gene, record, reference)
        coupled = vocabulary.get(_hypothesis(gene, aliases))
        decisions.append(
            decide_gene(
                gene,
                molecules="protein" if gene.kind == "CDS" else "nucleotide",
                vocabulary=vocabulary,
                hits=hits,
                aliases=aliases,
                neighbours=(
                    neighbours_of(record, gene, record_ir, aliases)
                    if record_ir is not None
                    else None
                ),
                reference_neighbours=(
                    neighbours_of(reference, coupled, reference_ir, aliases)
                    if reference_ir is not None and coupled is not None
                    else None
                ),
            )
        )
    return decisions


def _hypothesis(gene: Gene, aliases: AliasTable) -> str:
    from plastid.core.gene_names import normalize_name

    return normalize_name(gene.name, aliases).canonical


# ---------------------------------------------------------------- conflicts and the output


def _conflicts(
    records: Sequence[Genome], decisions: Sequence[tuple[int, GeneDecision]]
) -> list[str]:
    """Every case D5-09 §5's stop tier covers, each with the coordinates it names.

    Two kinds only, and both are about *disagreement* rather than about biology: a name whose
    sequence says another locus, and a rename that would give one record two different genes under
    one name. Missing and extra genes are deliberately not here — a lost gene is a result, not a
    tool failure, and stopping on it would make a real finding indistinguishable from a broken
    input.
    """
    lines: list[str] = []
    for index, decision in decisions:
        if decision.conflict:
            lines.append(
                f"{_label(records[index])} — {decision.gene.name} "
                f"{_location(decision.gene)}: {decision.reason}"
            )
    for index, record in enumerate(records):
        renamed = _renamed_genome(record, index, decisions)
        created = {
            decision.canonical
            for row_index, decision in decisions
            if row_index == index and decision.renamed and decision.canonical
        }
        for name in duplicate_names(renamed, _frame_of(renamed)):
            if name not in created:
                continue  # the record arrived that way: reported as a note, not stopped on
            where = ", ".join(_location(gene) for gene in renamed.genes if gene.name == name)
            lines.append(
                f"{_label(record)} — {name} at {where}: the normalized record would carry one "
                f"name for more than one gene, and a rename that creates a duplicate is a "
                f"conflict, not a merge"
            )
    return lines


def _renamed_genome(
    record: Genome, index: int, decisions: Sequence[tuple[int, GeneDecision]]
) -> Genome:
    """The record with every applied rename in place — the artifact, and the duplicate check."""
    by_gene = {
        id(decision.gene): decision for row_index, decision in decisions if row_index == index
    }
    genes = []
    for gene in record.genes:
        decision = by_gene.get(id(gene))
        if decision is not None and decision.renamed and decision.canonical:
            genes.append(replace(gene, name=decision.canonical))
        else:
            genes.append(gene)
    return replace(record, genes=tuple(genes))


def _artifacts(
    mode: str,
    records: Sequence[Genome],
    decisions: Sequence[tuple[int, GeneDecision]],
    notes: list[str],
) -> dict[str, str]:
    """The normalized records, or — in verify mode — none, and a note saying what it found.

    The artifact keeps the record's own ``/translation``: that protein is what the annotator
    searches (D5-08 §8.21), and a normalization that dropped it would quietly disarm the step
    standing behind it.
    """
    if mode == MODE_VERIFY:
        changed = [decision for _, decision in decisions if decision.renamed]
        notes.append(
            f"verify mode: {len(changed)} name(s) would change — a sample annotated against a "
            f"normalized reference should need none, so each one is a finding"
            if changed
            else "verify mode: nothing would change — the record already carries vocabulary A"
        )
        return {}
    return {
        f"{_label(record)}.normalized.gb": to_genbank(
            _renamed_genome(record, index, decisions), translations=True
        )
        for index, record in enumerate(records)
    }


# ---------------------------------------------------------------- the tables


def _summary(
    mode: str,
    reference: Genome,
    records: Sequence[Genome],
    decisions: Sequence[tuple[int, GeneDecision]],
    conflicts: Sequence[str],
    artifacts: dict[str, str],
) -> Table:
    counts = Counter(decision.status for _, decision in decisions)
    return Table(
        id="normalize_summary",
        title="Reference-driven naming",
        columns=(
            Column("mode", "Mode"),
            Column("reference", "Naming reference (A)"),
            Column("reference_genes", "Genes in the reference", dtype="int"),
            Column("records", "Records normalized", dtype="int"),
            Column("genes", "Genes read", dtype="int"),
            Column("renamed", "Renamed", dtype="int"),
            Column("confirmed", "Confirmed", dtype="int"),
            Column("conflicts", "Conflicts", dtype="int"),
            Column("measured", "Status: measured", dtype="int"),
            Column("not_measured", "Status: not measured", dtype="int"),
            Column("outside", "Status: outside vocabulary", dtype="int"),
            Column("artifacts", "Normalized records written", dtype="int"),
        ),
        rows=(
            {
                "mode": mode,
                "reference": _label(reference),
                "reference_genes": len(reference.genes),
                "records": len(records),
                "genes": len(decisions),
                "renamed": sum(1 for _, decision in decisions if decision.renamed),
                "confirmed": sum(
                    1 for _, decision in decisions if decision.outcome == OUTCOME_CONFIRMED
                ),
                "conflicts": len(conflicts),
                "measured": counts.get(STATUS_MEASURED, 0),
                "not_measured": counts.get(STATUS_NOT_MEASURED, 0),
                "outside": counts.get(STATUS_OUTSIDE, 0),
                "artifacts": len(artifacts),
            },
        ),
        scope=Scope(label=SCOPE_LABEL),
    )


def _pairing(records: Sequence[Genome], decisions: Sequence[tuple[int, GeneDecision]]) -> Table:
    """The one table where vocabulary B appears — in its own column, beside A's (D5-09 §3)."""
    width = len(str(len(decisions)))
    rows = []
    for number, (index, decision) in enumerate(decisions, start=1):
        gene = decision.gene
        start_1b, end_1b = gene.span
        rows.append(
            {
                "pair_id": f"p{number:0{width}d}",
                "record": _label(records[index]),
                "gene": gene.name,
                "kind": gene.kind,
                "canonical": decision.canonical,
                "status": decision.status,
                "outcome": decision.outcome,
                "rule": decision.rule,
                "identity_pct": decision.identity,
                "coverage": decision.coverage,
                "margin_bits": decision.margin_bits,
                "best_locus": decision.best_locus,
                "start_1b": start_1b,
                "end_1b": end_1b,
                "strand": "+" if gene.strand >= 0 else "-",
                "reason": decision.reason,
            }
        )
    return Table(
        id="pairing",
        title="Every gene of every record, paired with the naming reference",
        columns=(
            Column("pair_id", "Pair id"),
            Column("record", "Record"),
            Column("gene", "Name as it arrived (B)"),
            Column("kind", "Kind"),
            Column("canonical", "Name in the reference (A)"),
            Column("status", "Status"),
            Column("outcome", "Outcome"),
            Column("rule", "Decided by"),
            Column("identity_pct", "Identity", dtype="float", unit="%", precision=2),
            Column("coverage", "Coverage", dtype="float", precision=2),
            Column("margin_bits", "Margin", dtype="float", unit="bits", precision=1),
            Column("best_locus", "Best match in the reference"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("end_1b", "End", dtype="int", unit="bp"),
            Column("strand", "Strand"),
            Column("reason", "Why"),
        ),
        rows=tuple(rows),
        scope=Scope(label=SCOPE_LABEL),
        key="pair_id",
    )


def _differences(
    reference: Genome, records: Sequence[Genome], decisions: Sequence[tuple[int, GeneDecision]]
) -> Table | None:
    """Extra, missing, /product, length and part differences — reported, never stopping.

    Only the rows where something differs: a locus that agrees is not a finding, and ninety
    agreeing rows are a table nobody reads. What is compared is the **pairing**, so a locus the
    reference names and a record does not carry shows up as missing rather than absent, and a gene
    the record carries alone shows up as extra — the two cases D5-09 §5 keeps out of the stop tier
    on purpose.
    """
    reference_genes: dict[str, Gene] = {}
    for gene in reference.genes:
        if gene.kind in GENE_KINDS:
            reference_genes.setdefault(gene.name, gene)

    rows: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        paired: dict[str, Gene] = {}
        paired_names: set[str] = set()
        for row_index, decision in decisions:
            if row_index != index or not decision.canonical:
                continue
            paired.setdefault(decision.canonical, decision.gene)
            paired_names.add(decision.gene.name)
        for locus, reference_gene in sorted(reference_genes.items()):
            record_gene = paired.get(locus)
            if record_gene is None:
                rows.append(_difference_row(_label(record), locus, None, reference_gene, "missing"))
                continue
            differs = []
            if record_gene.length != reference_gene.length:
                differs.append("length")
            if len(record_gene.parts) != len(reference_gene.parts):
                differs.append("parts")
            if (
                record_gene.product
                and reference_gene.product
                and record_gene.product != reference_gene.product
            ):
                differs.append("/product")
            if differs:
                rows.append(
                    _difference_row(
                        _label(record), locus, record_gene, reference_gene, ", ".join(differs)
                    )
                )
        for gene in record.genes:
            if gene.kind not in GENE_KINDS or gene.name in paired_names:
                continue
            if gene.name in reference_genes:
                continue  # the reference spells this locus the same way; nothing differs
            rows.append(_difference_row(_label(record), None, gene, None, "extra"))
    if not rows:
        return None
    return Table(
        id="differences",
        title="Where the record and the naming reference differ (reported, never stopping)",
        columns=(
            Column("record", "Record"),
            Column("locus", "Locus (A)"),
            Column("in_reference", "In the reference"),
            Column("in_record", "In the record"),
            Column("length_record", "Record length", dtype="int", unit="bp"),
            Column("length_reference", "Reference length", dtype="int", unit="bp"),
            Column("product_record", "Record /product"),
            Column("product_reference", "Reference /product"),
            Column("parts_record", "Record parts", dtype="int"),
            Column("parts_reference", "Reference parts", dtype="int"),
            Column("differs", "Differs in"),
        ),
        rows=tuple(rows),
        scope=Scope(label=SCOPE_LABEL),
    )


def _difference_row(
    record: str,
    locus: str | None,
    record_gene: Gene | None,
    reference_gene: Gene | None,
    differs: str,
) -> dict[str, Any]:
    return {
        "record": record,
        "locus": locus,
        "in_reference": "yes" if reference_gene is not None else "no",
        "in_record": "yes" if record_gene is not None else "no",
        "length_record": record_gene.length if record_gene is not None else None,
        "length_reference": reference_gene.length if reference_gene is not None else None,
        "product_record": record_gene.product if record_gene is not None else None,
        "product_reference": reference_gene.product if reference_gene is not None else None,
        "parts_record": len(record_gene.parts) if record_gene is not None else None,
        "parts_reference": len(reference_gene.parts) if reference_gene is not None else None,
        "differs": differs,
    }


def _conflicts_table(
    records: Sequence[Genome],
    decisions: Sequence[tuple[int, GeneDecision]],
    conflicts: Sequence[str],
) -> Table:
    """Written only when the run was told to carry on past its conflicts (D5-09 §5)."""
    rows: list[dict[str, Any]] = []
    for number, (index, decision) in enumerate(
        ((row_index, d) for row_index, d in decisions if d.conflict), start=1
    ):
        gene = decision.gene
        start_1b, end_1b = gene.span
        rows.append(
            {
                "conflict_id": f"c{number:03d}",
                "record": _label(records[index]),
                "gene": gene.name,
                "canonical": decision.canonical,
                "kind": gene.kind,
                "start_1b": start_1b,
                "end_1b": end_1b,
                "strand": "+" if gene.strand >= 0 else "-",
                "reason": decision.reason,
            }
        )
    for offset, line in enumerate(conflicts[len(rows) :], start=len(rows) + 1):
        record, _, reason = line.partition(" — ")
        rows.append(
            {
                "conflict_id": f"c{offset:03d}",
                "record": record,
                "gene": "",
                "canonical": "",
                "kind": "",
                "start_1b": None,
                "end_1b": None,
                "strand": "",
                "reason": reason or line,
            }
        )
    return Table(
        id="conflicts",
        title="Conflicts carried through on on_conflict=continue",
        columns=(
            Column("conflict_id", "Conflict id"),
            Column("record", "Record"),
            Column("gene", "Gene (B)"),
            Column("canonical", "Name in the reference (A)"),
            Column("kind", "Kind"),
            Column("start_1b", "Start", dtype="int", unit="bp"),
            Column("end_1b", "End", dtype="int", unit="bp"),
            Column("strand", "Strand"),
            Column("reason", "Why it is a conflict"),
        ),
        rows=tuple(rows),
        scope=Scope(label=SCOPE_LABEL),
        key="conflict_id",
    )


__all__ = ["MODULE_ID", "MODULE_VERSION", "run"]
