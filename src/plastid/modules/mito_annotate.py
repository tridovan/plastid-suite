"""Annotating a plant mitogenome from a curated reference set (P3A-03, prototype).

The module is wiring, not algorithm: it resolves parameters, runs the two BLAST searches through
``external/blast.py``, hands the hits to ``core/mito_caller.py``, and turns the result into the
tables a reader checks and a GenBank file the rest of the pipeline can consume. The calling rules
and every measurement behind them live in the caller's own docstring.

Four things this module refuses to do, each for a measured reason:

* **merge with the input's own annotation.** The artifact carries the calls and nothing else. Run
  this on an already-annotated record and the obvious alternative — appending — would double the
  gene count of the very record the caller is checked against.
* **invent a product string.** Features are named by their symbol (``/gene=nad1``); the product
  wording belongs to a submitter, and ``data/mito_gene_names.yaml``'s synonyms were quoted *from*
  the answer-key record, so writing them back would be the annotation grading its own homework.
* **report success on nothing.** A search that finds no hit at all is a failed run, not an
  annotation of a genome with no genes (the rule the plastid path states the same way).
* **let the artifact go out unread.** The GenBank text is read back with our own reader before it
  is returned; a file that does not re-read is not a result (P3A-01's rule).

The output enters the same pipeline as a downloaded record, with no second code path: the calls
become a ``Genome``, and the ORF profile in the tables is ``core/orf.check_genome``'s — the same
function that scores the answer key.

**Scope (D3A-05), stated because a paper has to:** annotation and submission only. No de-novo gene
discovery, no intron model, no RNA-editing prediction, no repeats. Multi-exon calls report the
exons the alignment supports and never claim a join they cannot show (D3A-04). tRNA is P3A-04's
job (ARAGORN) and is not in this module's output. A reference gene this toolkit has no symbol for
appears in ``unnamed_hits`` — reported, never renamed (H-06).
"""

from __future__ import annotations

import tempfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plastid.core.genbank import read_genbank, to_genbank
from plastid.core.gene_names import load_mito_aliases
from plastid.core.mito_caller import CallResult, ReferenceIndex, call_cds, call_rrna
from plastid.core.mito_reference import (
    RRNA_LENGTH_BANDS,
    build_proteins,
    build_rrna,
    read_reference,
    render_fasta,
)
from plastid.core.model import Gene, Genome, count_genes
from plastid.core.orf import ORF_STATUSES, check_genome
from plastid.core.organelle import GROUP_MITO
from plastid.core.resources import data_file
from plastid.external.aragorn import AragornResult, run_aragorn
from plastid.external.blast import (
    BLAST_VERSION,
    MAX_HSPS,
    run_blastn,
    run_blastx,
)
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults
from plastid.run.tools import ToolError

MODULE_ID = "mito_annotate"
MODULE_VERSION = "0.1.0"

SCOPE_LABEL = (
    "CDS and rRNA called by transfer from the pinned mito reference set; one call per vocabulary "
    "locus; GenBank coordinates"
)
NO_IR_NOTE = (
    "No inverted-repeat gate applies: a plant mitogenome has no canonical IR, so this module "
    "never asks for one — the same rule modules/mito_genes.py follows."
)

# The measured separation on the answer-key record (core/mito_caller.py carries the evidence).
DEFAULT_MIN_IDENTITY = 50.0
DEFAULT_MIN_COVERAGE = 0.5


@dataclass(frozen=True, slots=True)
class _Reference:
    """The two FASTAs a run searches, and where they came from."""

    proteins: Path
    rrna: Path
    label: str
    exclusions: tuple[tuple[str, str, str], ...] = ()


@contextmanager
def _reference_set(params: dict[str, Any], vocabulary: Sequence[str]) -> Iterator[_Reference]:
    """The packaged reference set, or one rebuilt from a directory of GenBank records.

    A caller-supplied directory goes through P3A-02's own builders, and **what they exclude is
    printed** (``ReferenceBuild.excluded``): a reference set that quietly dropped half its genes
    would make every "not called" row a lie.
    """
    given = params.get("reference")
    if given is None:
        yield _Reference(
            proteins=data_file("mito_reference_proteins.faa"),
            rrna=data_file("mito_reference_rrna.fa"),
            label="the packaged reference set (data/mito_reference_proteins.faa + _rrna.fa)",
        )
        return
    directory = Path(given)
    if not directory.is_dir():
        raise ValueError(f"--set reference={given} is not a directory of GenBank records")
    records = sorted(directory.glob("*.gb"))
    if not records:
        raise ValueError(f"--set reference={given} holds no *.gb records to build a set from")
    proteins = build_proteins(records, vocabulary)
    rrna = build_rrna(records)
    with tempfile.TemporaryDirectory(prefix="mito_reference_") as work:
        proteins_path = Path(work) / "proteins.faa"
        rrna_path = Path(work) / "rrna.fa"
        proteins_path.write_text(render_fasta(proteins.entries), encoding="utf-8")
        rrna_path.write_text(render_fasta(rrna.entries), encoding="utf-8")
        yield _Reference(
            proteins=proteins_path,
            rrna=rrna_path,
            label=f"rebuilt from {len(records)} GenBank record(s) in {directory}",
            exclusions=(*proteins.excluded, *rrna.excluded),
        )


@module(
    MODULE_ID,
    group=GROUP_MITO,
    version=MODULE_VERSION,
    title="Mitochondrial annotation (reference transfer, prototype)",
    requires_ir=False,  # deliberate: see NO_IR_NOTE,
    params=(
        ParamSpec("reference", "str", "directory of GenBank records to transfer from"),
        ParamSpec("min_identity", "float", "minimum identity for a transferred call"),
        ParamSpec("min_coverage", "float", "minimum query coverage for a transferred call"),
        ParamSpec("evalue", "float", "maximum BLAST E-value"),
        ParamSpec("transl_table", "int", "translation table stamped on called CDS"),
        ParamSpec("tools", "str", "directory holding BLAST+/ARAGORN when they are not on PATH"),
    ),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    defaults = load_defaults()["mito"]
    min_identity = float(params.get("min_identity", defaults["min_identity"]))
    min_coverage = float(params.get("min_coverage", defaults["min_coverage"]))
    evalue = float(params.get("evalue", defaults["evalue"]))
    table = int(params.get("transl_table", defaults["translation_table"]))
    tools_dir = params.get("tools")
    seq_id = genome.accession or genome.name

    vocabulary = tuple(locus.canonical for locus in load_mito_aliases().loci)
    rrna_vocabulary = tuple(sorted(RRNA_LENGTH_BANDS))

    with _reference_set(params, vocabulary) as reference:
        proteins = read_reference(reference.proteins)
        rrna = read_reference(reference.rrna)
        protein_index = ReferenceIndex.from_entries(proteins, vocabulary=vocabulary)
        rrna_index = ReferenceIndex.from_entries(rrna, vocabulary=rrna_vocabulary)

        protein_run = run_blastx(
            genome,
            reference.proteins,
            ctx,
            seq_id=seq_id,
            directory=tools_dir,
            evalue=evalue,
            max_hsps=MAX_HSPS,
        )
        rrna_run = run_blastn(
            genome,
            reference.rrna,
            ctx,
            seq_id=seq_id,
            directory=tools_dir,
            evalue=evalue,
            max_hsps=MAX_HSPS,
        )

        trna_run = run_aragorn(genome, ctx, seq_id=seq_id, directory=tools_dir)
        if not trna_run.hits:
            findings_trna = (
                "ARAGORN found no tRNA gene at all — for a plant mitogenome that is worth a look "
                "before the annotation is trusted"
            )
        else:
            findings_trna = ""

        if not protein_run.hits and not rrna_run.hits:
            raise ToolError(
                "mito_annotate: the reference searches returned no hit at all — an empty "
                "annotation is a failed annotation, not a genome without genes. Check that the "
                "input is a plant mitochondrion and that the reference set is the pinned one."
            )

        cds_result = call_cds(
            genome,
            protein_index,
            protein_run.hits,
            min_identity=min_identity,
            min_coverage=min_coverage,
            transl_table=table,
            max_hsps_per_subject=MAX_HSPS,
        )
        rrna_result = call_rrna(
            genome,
            rrna_index,
            rrna_run.hits,
            min_identity=min_identity,
            min_coverage=min_coverage,
            max_hsps_per_subject=MAX_HSPS,
        )

        if not cds_result.calls and not rrna_result.calls:
            raise ToolError(
                "mito_annotate: the searches found hits but no locus passed the floors — an "
                "annotation with no genes is a failed annotation. Every rejected locus is named "
                "in `not_called` with the numbers it failed on; check that the reference set "
                "belongs to this plant before loosening a floor."
            )

        annotated = _annotated(genome, cds_result, rrna_result, trna_run.genes, table)
        text = to_genbank(annotated)
        _verify_round_trip(text, annotated, seq_id)

    return _result(
        genome=annotated,
        cds=cds_result,
        rrna=rrna_result,
        trna=trna_run,
        reference=reference,
        floors=(min_identity, min_coverage),
        evalue=evalue,
        table=table,
        text=text,
        raw=(protein_run.raw, rrna_run.raw, trna_run.raw),
        notes=(
            *cds_result.findings,
            *rrna_result.findings,
            *([findings_trna] if findings_trna else []),
        ),
    )


def _annotated(
    genome: Genome,
    cds: CallResult,
    rrna: CallResult,
    trna: Sequence[Gene],
    table: int,
) -> Genome:
    """The calls as a genome: the calls only, never the input's own features (module docstring).

    ``table`` is unused here but named in the signature on purpose — every CDS' table travels on
    the call itself (:meth:`GeneCall.to_gene`), so this function cannot stamp a different one.
    """
    genes = (*[call.to_gene() for call in (*cds.calls, *rrna.calls)], *trna)
    return Genome(
        name=genome.name,
        sequence=genome.sequence,
        accession=genome.accession,
        topology=genome.topology,
        genes=genes,
        repeat_regions=genome.repeat_regions,
        source=dict(genome.source),
    )


def _verify_round_trip(text: str, annotated: Genome, seq_id: str) -> None:
    """Read the artifact back with our own reader — a file that does not re-read is not a result."""
    with tempfile.TemporaryDirectory(prefix="mito_annotate_check_") as work:
        path = Path(work) / f"{seq_id}.gb"
        path.write_text(text, encoding="utf-8")
        back = read_genbank(path)
    before, after = count_genes(annotated), count_genes(back)
    if after != before or back.sequence != annotated.sequence:
        raise ValueError(
            f"mito_annotate: the GenBank it wrote does not read back as the genome it wrote "
            f"({before.total} genes in, {after.total} out) — the export is broken, and an "
            f"artifact nobody can re-read is not a result"
        )


def _result(
    *,
    genome: Genome,
    cds: CallResult,
    rrna: CallResult,
    trna: AragornResult,
    reference: _Reference,
    floors: tuple[float, float],
    evalue: float,
    table: int,
    text: str,
    raw: tuple[str, ...],
    notes: Sequence[str],
) -> Result:
    calls = (*cds.calls, *rrna.calls)
    not_called = (*cds.not_called, *rrna.not_called)
    dropped = (*cds.dropped, *rrna.dropped)
    unnameable = (*cds.unnameable, *rrna.unnameable)
    competing = (*cds.competing, *rrna.competing)
    min_identity, min_coverage = floors
    checks = check_genome(genome, default_table=table)
    # one status per called CDS: `check_genome` walks the same genes in the same order, but a
    # lookup by name cannot silently pair a status with the wrong feature
    statuses = {check.gene.name: check.status for check in checks}

    summary = Table(
        id="mito_annotation_summary",
        title="Mitochondrial annotation — what was called, from what, with which floors",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("n_calls", "Genes called", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_not_called", "Vocabularies not called", dtype="int"),
            Column("n_cds_failing_orf", "Called CDS failing the ORF checks", dtype="int"),
            Column("n_dropped", "Hits not used", dtype="int"),
            Column("n_unnameable", "Hits this vocabulary cannot name", dtype="int"),
            Column("min_identity_pct", "Identity floor", dtype="float", unit="%"),
            Column("min_coverage_pct", "Coverage floor", dtype="float", unit="%"),
            Column("evalue", "E-value ceiling"),
            Column("transl_table", "Translation table", dtype="int"),
            Column("blast", "BLAST"),
            Column("aragorn", "ARAGORN"),
            Column("reference", "Reference set"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "length_bp": len(genome),
                "n_calls": len(calls) + len(trna.hits),
                "n_cds": sum(1 for call in calls if call.kind == "CDS"),
                "n_rrna": sum(1 for call in calls if call.kind == "rRNA"),
                "n_trna": len(trna.hits),
                "n_not_called": len(not_called),
                "n_cds_failing_orf": sum(1 for check in checks if check.is_failure),
                "n_dropped": len(dropped),
                "n_unnameable": len(unnameable),
                "min_identity_pct": min_identity,
                "min_coverage_pct": min_coverage,
                "evalue": f"{evalue:g}",
                "transl_table": table,
                "blast": f"BLAST+ {BLAST_VERSION}",
                "aragorn": f"ARAGORN {trna.version}",
                "reference": reference.label,
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    called = Table(
        id="called_genes",
        title="Called genes — one row per vocabulary locus",
        columns=(
            Column("gene_id", "Gene id"),
            Column("gene", "Gene"),
            Column("kind", "Kind"),
            Column("strand", "Strand"),
            Column("start", "Start", dtype="int"),
            Column("end", "End", dtype="int"),
            Column("n_exons", "Exons", dtype="int"),
            Column("exons", "Exons (transcript order)"),
            Column("length_bp", "Aligned length", dtype="int", unit="bp"),
            Column("identity_pct", "Identity", dtype="float", unit="%", precision=1),
            Column("coverage_pct", "Reference covered", dtype="float", unit="%", precision=0),
            Column("reference_accession", "Reference"),
            Column("reference_gene", "Reference's own name"),
            Column("partial_start", "5' end"),
            Column("partial_end", "3' end"),
            Column("orf_status", "ORF status"),
            Column("notes", "Notes"),
        ),
        rows=tuple(_call_row(call, statuses) for call in calls),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="gene_id",
    )

    trna_table = Table(
        id="trna_calls",
        title="tRNA genes — ARAGORN's calls, named through the shared normalizer",
        columns=(
            Column("gene_id", "Gene id"),
            Column("gene", "Gene"),
            Column("as_reported", "As ARAGORN named it"),
            Column("anticodon", "Anticodon (RNA)"),
            Column("start", "Start", dtype="int"),
            Column("end", "End", dtype="int"),
            Column("strand", "Strand"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("score", "Score", dtype="int"),
            Column("crosses_origin", "Crosses the origin"),
        ),
        rows=tuple(
            {
                "gene_id": f"{hit.name}:{hit.start}-{hit.end}",
                "gene": hit.name,
                "as_reported": hit.source_name,
                "anticodon": hit.anticodon,
                "start": hit.start,
                "end": hit.end,
                "strand": "+" if hit.strand > 0 else "-",
                "length_bp": sum(part.length for part in hit.parts(len(genome))),
                "score": hit.score,
                "crosses_origin": "yes" if hit.wrapped else "no",
            }
            for hit in trna.hits
        ),
        scope=Scope(
            label=(
                "tRNA genes called by ARAGORN on this record, both strands; a gene crossing the "
                "origin is written as the two parts a join() needs"
            ),
            coord_system="genbank",
        ),
        key="gene_id",
    )

    tables: dict[str, Table] = {
        summary.id: summary,
        called.id: called,
        trna_table.id: trna_table,
    }
    if not_called:
        tables["not_called"] = Table(
            id="not_called",
            title="Vocabulary loci with no call — and the best evidence that was rejected",
            columns=(
                Column("gene", "Gene"),
                Column("reason", "Why"),
                Column("best_accession", "Best reference"),
                Column("best_identity_pct", "Its identity", dtype="float", unit="%", precision=1),
                Column("best_coverage_pct", "Its coverage", dtype="float", unit="%", precision=0),
                Column("best_bitscore", "Its bitscore", dtype="float", precision=0),
            ),
            rows=tuple(
                {
                    "gene": row.symbol,
                    "reason": row.reason,
                    "best_accession": row.best_accession or "",
                    "best_identity_pct": row.best_identity,
                    "best_coverage_pct": (
                        None if row.best_coverage is None else row.best_coverage * 100
                    ),
                    "best_bitscore": row.best_bitscore,
                }
                for row in not_called
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="gene",
        )
    if competing:
        tables["competing_references"] = Table(
            id="competing_references",
            title="Every reference considered for a locus, and which one was chosen",
            columns=(
                Column("reference_id", "Reference (accession|gene|length)"),
                Column("gene", "Locus"),
                Column("reference_gene", "Reference's own name"),
                Column("length_aa", "Reference length", dtype="int", unit="aa"),
                Column("identity_pct", "Identity", dtype="float", unit="%", precision=1),
                Column("coverage_pct", "Coverage", dtype="float", unit="%", precision=0),
                Column("bitscore", "Bitscore", dtype="float", precision=0),
                Column("passes_floors", "Passes the floors"),
                Column("chosen", "Chosen"),
            ),
            rows=tuple(
                {
                    # the length belongs in the key: one record can carry two sequences under
                    # one symbol (measured: Arabidopsis atp6 is 385 and 349 aa), and each is a
                    # candidate of its own
                    "reference_id": (
                        f"{row.reference_accession}|{row.reference_gene}|{row.reference_length}"
                    ),
                    "gene": row.symbol,
                    "reference_gene": row.reference_gene,
                    "length_aa": row.reference_length,
                    "identity_pct": row.identity,
                    "coverage_pct": row.coverage * 100,
                    "bitscore": row.bitscore,
                    "passes_floors": "yes" if row.passes_floors else "no",
                    "chosen": "yes" if row.chosen else "no",
                }
                for row in competing
            ),
            scope=Scope(
                label="the reference sequences this run searched; BLAST coordinates",
                coord_system="genbank",
            ),
            key="reference_id",
        )
    if dropped:
        tables["dropped_hits"] = Table(
            id="dropped_hits",
            title="Hits this run did not use — every one with its reason",
            columns=(
                Column("hit_id", "Hit"),
                Column("gene", "Locus"),
                Column("query_start", "Query start", dtype="int"),
                Column("query_end", "Query end", dtype="int"),
                Column("identity_pct", "Identity", dtype="float", unit="%", precision=1),
                Column("bitscore", "Bitscore", dtype="float", precision=0),
                Column("reason", "Why it was not used"),
            ),
            rows=tuple(
                {
                    # the subject belongs in the key: the same query span is a separate hit
                    # against every reference sequence that aligns to it
                    "hit_id": f"{row.hsp.subject_key}:{row.hsp.qlo}-{row.hsp.qhi}",
                    "gene": row.symbol,
                    "query_start": row.hsp.qlo,
                    "query_end": row.hsp.qhi,
                    "identity_pct": row.hsp.pident,
                    "bitscore": row.hsp.bitscore,
                    "reason": row.reason,
                }
                for row in dropped
            ),
            scope=Scope(
                label="BLAST hits to the reference set, before chaining; query coordinates",
                coord_system="genbank",
            ),
            key="hit_id",
        )
    if unnameable:
        tables["unnamed_hits"] = Table(
            id="unnamed_hits",
            title=(
                "Genes found whose symbol this toolkit does not hold (H-06) — reported, not named"
            ),
            columns=(
                Column("gene_id", "Reference name"),
                Column("reference_accession", "Reference"),
                Column("start", "Start", dtype="int"),
                Column("end", "End", dtype="int"),
                Column("strand", "Strand"),
                Column("identity_pct", "Identity", dtype="float", unit="%", precision=1),
                Column("coverage_pct", "Coverage", dtype="float", unit="%", precision=0),
            ),
            rows=tuple(
                {
                    # the accession belongs in the key: two records can give the same region
                    # two different unnameable symbols, and both are worth reporting
                    "gene_id": (
                        f"{call.reference.accession}|{call.symbol}:{call.span[0]}-{call.span[1]}"
                    ),
                    "reference_accession": call.reference.accession,
                    "start": call.span[0],
                    "end": call.span[1],
                    "strand": _strand(call.strand),
                    "identity_pct": call.identity,
                    "coverage_pct": call.coverage * 100,
                }
                for call in unnameable
            ),
            scope=Scope(
                label=(
                    "hits to reference genes with no entry in data/mito_gene_names.yaml — they "
                    "are deliberately outside the annotation until H-06 rules on the vocabulary"
                ),
                coord_system="genbank",
            ),
            key="gene_id",
        )

    status_counts = {status: 0 for status in ORF_STATUSES}
    for check in checks:
        status_counts[check.status] += 1
    tables["orf_status"] = Table(
        id="orf_status",
        title="ORF status of the called CDS — the same checks a downloaded record gets",
        columns=(Column("status", "Status"), Column("n", "CDS", dtype="int")),
        rows=tuple({"status": status, "n": status_counts[status]} for status in ORF_STATUSES),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="status",
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(_notes(genome, cds, rrna, trna, reference, floors, evalue, table, notes)),
        artifacts={"mito_annotation.gb": text, "blast_hits.tsv": _hit_lines(*raw)},
    )


def _notes(
    genome: Genome,
    cds: CallResult,
    rrna: CallResult,
    trna: AragornResult,
    reference: _Reference,
    floors: tuple[float, float],
    evalue: float,
    table: int,
    findings: Sequence[str],
) -> Iterator[str]:
    min_identity, min_coverage = floors
    yield (
        f"Prototype: a reference-transfer annotation, not a submission. The calls come from "
        f"BLAST+ hits against {reference.label}, and this module does no de-novo gene discovery, "
        f"no intron model and no RNA-editing prediction (D3A-05)."
    )
    yield (
        f"Floors stated: identity ≥ {min_identity:.1f} %, coverage ≥ {min_coverage:.0%} of the "
        f"reference, E-value ≤ {evalue:g}. A locus that fails them is listed in `not_called` with "
        f"the numbers it failed on; nothing is rounded up to a call."
    )
    yield (
        f"CDS are checked against translation table {table} and every call carries it — the "
        f"initiator test that decides a partial start and the ORF checks downstream read the same "
        f"number, so they cannot disagree."
    )
    yield (
        "Boundaries come from the alignment to the reference protein, never from ORF-walking "
        "(D3A-02): the stop codon is taken only when the sequence confirms it, one codon past the "
        "alignment, and an end that is not established is reported as partial rather than "
        "extended to the next stop."
    )
    yield (
        "Exons are reported where the alignment breaks and in the transcript order the reference "
        "gives; a join is never claimed (D3A-04). A trans-spliced gene therefore appears as its "
        "exons, not as one joined CDS."
    )
    yield (
        "No `/product` is written: the feature carries its symbol (`/gene=nad1`), because the "
        "product wording is the submitter's and this toolkit's mito vocabulary was derived from "
        "the very record the calls are checked against."
    )
    yield (
        f"tRNA genes come from ARAGORN {trna.version} (`-w -gc1`, circular or linear as the "
        f"record says), and their names go through `core/gene_names.normalize_name` — the same "
        f"normalizer a downloaded record gets, so there is no second naming rule. A tRNA whose "
        f"coordinates wrap the origin is written as the two parts a join() needs, and a gene "
        f"ARAGORN reports that is not a tRNA (a tmRNA, say) is listed in the notes rather than "
        f"parsed as one."
    )
    if trna.others:
        yield (
            f"{len(trna.others)} gene line(s) from ARAGORN were not tRNAs and are reported here "
            f"rather than dropped: {'; '.join(trna.others)}"
        )
    if trna.declared is not None and trna.declared != len(trna.hits):
        yield (
            f"ARAGORN's own summary says {trna.declared} genes and this reader kept "
            f"{len(trna.hits)} tRNAs — the difference is a finding, not a rounding"
        )
    yield NO_IR_NOTE
    yield (
        f"Gene count here: {len(cds.calls)} CDS + {len(trna.hits)} tRNA + {len(rrna.calls)} "
        f"rRNA. The counting rule is the fixed one (CDS + tRNA + rRNA); the numbers are this "
        f"record's, not a constant."
    )
    if reference.exclusions:
        sample = "; ".join(
            f"{accession} {name} ({why})" for accession, name, why in reference.exclusions[:5]
        )
        yield (
            f"{len(reference.exclusions)} reference feature(s) were excluded while building the "
            f"set and are named rather than dropped in silence: {sample}"
            + (" …" if len(reference.exclusions) > 5 else "")
        )
    if cds.unnameable or rrna.unnameable:
        yield (
            f"{len(cds.unnameable) + len(rrna.unnameable)} gene(s) were found whose reference "
            f"symbol this toolkit has no vocabulary entry for — listed in `unnamed_hits`, not "
            f"annotated: renaming them here would change naming for every mito record (H-06)."
        )
    if findings:
        yield f"{len(findings)} finding(s) are reported in the notes below rather than hidden:"
    yield from findings


def _hit_lines(*raw: str) -> str:
    """The raw BLAST+ output this run consumed — the evidence behind every row above.

    Kept verbatim rather than re-rendered from the parsed records: a reader who doubts a call can
    check it against what the search actually printed, and the run's own tool logs hold the same
    text beside the argv that produced it.
    """
    return "".join(text for text in raw)


def _call_row(call: Any, statuses: dict[str, str]) -> dict[str, Any]:
    return {
        "gene_id": f"{call.symbol}:{call.span[0]}-{call.span[1]}",
        "gene": call.symbol,
        "kind": call.kind,
        "strand": _strand(call.strand),
        "start": call.span[0],
        "end": call.span[1],
        "n_exons": len(call.parts),
        "exons": " ".join(
            f"{part.to_1based()[0]}-{part.to_1based()[1]}({'+' if part.strand > 0 else '-'})"
            for part in call.parts
        ),
        "length_bp": sum(part.length for part in call.parts),
        "identity_pct": call.identity,
        "coverage_pct": call.coverage * 100,
        "reference_accession": call.reference.accession,
        "reference_gene": call.reference.source_gene,
        "partial_start": "partial" if call.partial_start else "complete",
        "partial_end": "partial" if call.partial_end else "complete",
        "orf_status": statuses.get(call.symbol, "not checked"),
        "notes": "; ".join(call.notes),
    }


def _strand(strand: int) -> str:
    return {1: "+", -1: "-"}.get(strand, "mixed")
