"""Annotating a plastome from a curated reference set (P3B-08, D5-08's engine).

The module is wiring, not algorithm: it resolves the reference set, runs the searches through
``external/blast.py`` and ``external/aragorn.py``, hands the hits to ``core/plastid_caller.py``,
and turns the result into the tables a reader checks and a GenBank file the rest of the pipeline
consumes. The calling rules and every measurement behind them live in ``core/caller.py``.

**One call per locus, and the losers are listed.** The caller chooses the highest-identity
candidate that passes the floors and reports every other reference sequence it chose over
(``competing_references``, one row each, with the numbers it lost on). PGA pooled every
reference's sequences and annotated each independently — measured: 147 genes where the golden has
132, four loci called twice under two spellings (§8.15). Here a locus is one row.

**Overlap is not error by default.** Measured on the golden: four overlapping pairs are *in the
contract* (psbC×psbD 53 bp · ndhC×ndhK 10 · atpB×atpE 4 · rps12×rps12 114). A rule that dropped
the weaker overlap would delete real genes, so both calls stand and the pair is reported.

**Four things this module refuses to do**, each for a measured reason:

* **merge with the input's own annotation.** An assembled plastome arrives bare; if it carries
  features, they are not the answer and appending them would double the gene count.
* **invent a product string.** Features are named by their symbol; the product wording belongs to
  a submitter, and the vocabulary's synonyms were quoted *from* the records this run is checked
  against.
* **report success on nothing.** A search that finds no hit at all is a failed run, not an
  annotation of a genome with no genes.
* **let the artifact go out unread.** The GenBank text is read back with our own reader before it
  is returned (P3A-01's rule), and the ``source`` carries ``/organelle`` (H-11) — the one thing
  the mito artifact was measured missing.

**The inverted repeat is read behind the gate, never instead of it.** ``ctx.require_ir()`` stops
the run when the two copies are not ≥ 99.9 % identical (§9.A: naive detection got 2 of 9 published
plastomes wrong, one silently). The self-search is a second opinion and a disagreement is a
finding (``core/plastid_caller.repeat_agreement``); the structural frame is what the artifact
writes and what positions are normalized against.

**Scope (D5-08's non-goals), stated because a paper has to:** reference transfer only. No de-novo
gene discovery, no intron model beyond the alignment's own evidence, no RNA-editing prediction.
A reference gene this vocabulary cannot name is reported (``unnamed_hits``), never renamed.
"""

from __future__ import annotations

import tempfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plastid.core.caller import render_fasta
from plastid.core.genbank import read_genbank, to_genbank
from plastid.core.model import Gene, Genome, Interval, count_genes
from plastid.core.orf import check_genome
from plastid.core.plastid_caller import (
    DEFAULT_MIN_COVERAGE,
    DEFAULT_MIN_IDENTITY,
    PLASTID_TRANSL_TABLE,
    TruncatedAnnotation,
    read_inverted_repeat,
    reference_entries,
    reference_index,
    repeat_agreement,
    truncated_annotations,
    vocabulary_of,
)
from plastid.core.plastid_caller import (
    call_cds as call_plastid_cds,
)
from plastid.core.plastid_caller import (
    call_rrna as call_plastid_rrna,
)
from plastid.core.plastid_reference import read_manifest
from plastid.core.resources import data_file
from plastid.core.submission import to_fasta
from plastid.external.aragorn import AragornResult, run_aragorn
from plastid.external.blast import (
    BLAST_VERSION,
    MAX_HSPS,
    BlastResult,
    run_blastn,
    run_blastx,
)
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults
from plastid.run.tools import ToolError

MODULE_ID = "plastid_annotate"
MODULE_VERSION = "0.1.0"

# The organelle this run serves, declared on the artifact's source (H-11) and in the scope label.
ORGANELLE = "plastid"

SCOPE_LABEL = (
    "CDS and rRNA called by transfer from the pinned plastid reference set; tRNA called de novo "
    "by ARAGORN; one call per vocabulary locus; GenBank coordinates"
)
REFERENCE_MANIFEST = "plastid_reference.json"


@dataclass(frozen=True, slots=True)
class _Reference:
    """The records a run searches, where they came from, and which of them were left out."""

    records: tuple[Genome, ...]
    label: str
    excluded: tuple[str, ...] = ()
    truncated: tuple[TruncatedAnnotation, ...] = ()


class ReferenceSetMissing(ToolError):
    """The pinned reference set is not on this machine — named, with the way to supply one."""


def _reference_directory(params: dict[str, Any], target: str) -> tuple[Path, tuple[str, ...]]:
    """Where the reference records are, and which accessions may annotate ``target``.

    The pinned set lives outside the repository (``.panel/`` is gitignored), so a run has to be
    told where it is; the recorded manifest supplies both the default directory and the
    **leave-one-out** rule — a record must never annotate a copy of itself, which is the corpus's
    own measured trap (§8.31: the PGA set carried the golden's sequence and read 85 CDS where the
    contract says 86).
    """
    recorded = data_file(REFERENCE_MANIFEST)
    manifest = read_manifest(recorded) if recorded.is_file() else None
    given = params.get("reference")
    directory = Path(given) if given is not None else Path(manifest.directory if manifest else "")
    if not directory.is_dir():
        raise ReferenceSetMissing(
            f"plastid_annotate: the reference set is not at {directory} — the pinned plastid "
            f"records live outside the repository (.panel/ is gitignored), so a run is pointed at "
            f"them with `--set reference=<directory of GenBank records>`"
        )
    allowed = manifest.references_for(target) if manifest is not None else ()
    return directory, allowed


@contextmanager
def _reference_set(params: dict[str, Any], target: str) -> Iterator[_Reference]:
    """The records to search, read from disk, with the target's own copies left out."""
    directory, allowed = _reference_directory(params, target)
    paths = sorted(directory.glob("*.gb"))
    if not paths:
        raise ReferenceSetMissing(
            f"plastid_annotate: {directory} holds no *.gb records to transfer an annotation from"
        )
    records, excluded = [], []
    for path in paths:
        genome = read_genbank(path)
        accession = genome.accession or genome.name
        if allowed and accession not in allowed:
            excluded.append(accession)
            continue
        records.append(genome)
    if not records:
        raise ReferenceSetMissing(
            f"plastid_annotate: every record in {directory} is excluded for {target} — "
            f"leave-one-out leaves nothing to search"
        )
    truncated = (
        *truncated_annotations(records, kind="CDS"),
        *truncated_annotations(records, kind="rRNA"),
    )
    label = f"{len(records)} GenBank record(s) in {directory}"
    if excluded:
        label += f", {len(excluded)} left out by leave-one-out ({', '.join(sorted(excluded))})"
    if truncated:
        label += f", {len(truncated)} truncated duplicate(s) not searched"
    yield _Reference(
        records=tuple(records),
        label=label,
        excluded=tuple(sorted(excluded)),
        truncated=truncated,
    )


def _search_files(records: Sequence[Genome], work: Path) -> tuple[Path, Path]:
    """The two FASTA files the searches run against — one per kind, never one pooled database."""
    proteins = work / "proteins.faa"
    rrna = work / "rrna.fa"
    proteins.write_text(render_fasta(reference_entries(records, kind="CDS")), encoding="utf-8")
    rrna.write_text(render_fasta(reference_entries(records, kind="rRNA")), encoding="utf-8")
    return proteins, rrna


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="Plastid annotation (reference transfer)",
    requires_ir=True,
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    defaults = load_defaults().get("plastid", {})
    min_identity = float(
        params.get("min_identity", defaults.get("min_identity", DEFAULT_MIN_IDENTITY))
    )
    min_coverage = float(
        params.get("min_coverage", defaults.get("min_coverage", DEFAULT_MIN_COVERAGE))
    )
    evalue = float(params.get("evalue", defaults.get("evalue", 1.0e-5)))
    table = int(params.get("transl_table", defaults.get("translation_table", PLASTID_TRANSL_TABLE)))
    tools_dir = params.get("tools")
    seq_id = genome.accession or genome.name

    with (
        _reference_set(params, seq_id) as reference,
        tempfile.TemporaryDirectory(prefix="plastid_annotate_") as work,
    ):
        workdir = Path(work)
        proteins, rrna = _search_files(reference.records, workdir)
        # The self-search's subject is the record itself, so its header has to be the same
        # three-field shape every other reference header is (measured 08/10/2026: a bare
        # accession comes back as a one-field id and the parser refuses it by name — rightly,
        # since a key that cannot match an entry is worse than an error).
        self_fasta = workdir / "self.fna"
        self_fasta.write_text(to_fasta(genome, f"{seq_id}|self|self"), encoding="utf-8")

        protein_run = run_blastx(
            genome,
            proteins,
            ctx,
            seq_id=seq_id,
            directory=tools_dir,
            evalue=evalue,
            max_hsps=MAX_HSPS,
        )
        rrna_run = run_blastn(
            genome, rrna, ctx, seq_id=seq_id, directory=tools_dir, evalue=evalue, max_hsps=MAX_HSPS
        )
        repeat_run = run_blastn(
            genome,
            self_fasta,
            ctx,
            seq_id=seq_id,
            directory=tools_dir,
            evalue=evalue,
            max_hsps=MAX_HSPS,
        )
        trna_run = run_aragorn(genome, ctx, seq_id=seq_id, directory=tools_dir)

        if not protein_run.hits and not rrna_run.hits:
            raise ToolError(
                "plastid_annotate: the reference searches returned no hit at all — an empty "
                "annotation is a failed annotation, not a genome without genes. Check that the "
                "input is a plastome and that the reference set belongs to its family."
            )

        cds_vocabulary = vocabulary_of(reference.records, kind="CDS")
        rrna_vocabulary = vocabulary_of(reference.records, kind="rRNA")
        cds = call_plastid_cds(
            genome,
            reference_index(reference.records, kind="CDS", vocabulary=cds_vocabulary),
            protein_run.hits,
            min_identity=min_identity,
            min_coverage=min_coverage,
            transl_table=table,
            max_hsps_per_subject=MAX_HSPS,
        )
        rrna_calls = call_plastid_rrna(
            genome,
            reference_index(reference.records, kind="rRNA", vocabulary=rrna_vocabulary),
            rrna_run.hits,
            min_identity=min_identity,
            min_coverage=min_coverage,
            max_hsps_per_subject=MAX_HSPS,
        )
        if not cds.calls and not rrna_calls.calls:
            raise ToolError(
                "plastid_annotate: the searches found hits but no locus passed the floors — an "
                "annotation with no genes is a failed annotation. Every rejected locus is named "
                "in `not_called` with the numbers it failed on; check that the reference set "
                "belongs to this plant before loosening a floor."
            )

        searched_repeat = read_inverted_repeat(repeat_run.hits)
        annotated = _annotated(genome, cds, rrna_calls, trna_run.genes, ir)
        text = to_genbank(annotated)
        _verify_round_trip(text, annotated, seq_id)

    return _result(
        genome=annotated,
        cds=cds,
        rrna=rrna_calls,
        trna=trna_run,
        reference=reference,
        protein_run=protein_run,
        rrna_run=rrna_run,
        repeat_run=repeat_run,
        repeat_note=repeat_agreement(ir, searched_repeat),
        floors=(min_identity, min_coverage),
        evalue=evalue,
        table=table,
        text=text,
    )


def _annotated(
    genome: Genome,
    cds: Any,
    rrna: Any,
    trna: Sequence[Gene],
    ir: Any,
) -> Genome:
    """The calls as a genome — the calls only, never the input's own features.

    ``source`` carries ``/organelle``: this toolkit writes the record it annotated, and a record
    that does not say which organelle it is cannot be checked by the modules behind it (H-11 —
    the mito artifact was measured copying an empty ``source`` through, and failing
    ``submission``'s first pass 24 times). ``repeat_regions`` are the gate's own two IR spans.
    """
    genes = _ir_copies([call.to_gene() for call in (*cds.calls, *rrna.calls)], ir) + tuple(trna)
    source = dict(genome.source)
    source["organelle"] = ORGANELLE
    return Genome(
        name=genome.name,
        sequence=genome.sequence,
        accession=genome.accession,
        topology=genome.topology,
        genes=genes,
        repeat_regions=(Interval(ir.ira.start, ir.ira.end), Interval(ir.irb.start, ir.irb.end)),
        source=source,
    )


def _ir_copies(genes: Sequence[Gene], ir: Any) -> tuple[Gene, ...]:
    """Every gene the inverted repeat carries, written **twice** — one copy per side.

    The repeat is two copies of the same sequence in opposite orientations, so a plastome's
    published annotation lists each IR gene twice: the golden counts **86 CDS** where its distinct
    loci are 79, and all four rRNA loci are inside the repeat (§8.34's 132 = 86 + 38 + 8). The
    caller makes **one call per locus** on purpose; the copies are the *export's* job, and they
    are the same fact written in the record's own coordinates rather than a second measurement.

    The mapping is the inversion itself: a part at ``[s, e]`` inside one copy sits at the mirrored
    offset inside the other, on the opposite strand. A call that straddles an IR boundary is left
    exactly as it was called and **named in the notes** — a trans-spliced gene whose exons span
    the single-copy region and the repeat (the golden's ``rps12``) is not copied by a rule that
    only knows how to mirror a gene that sits inside one side.
    """
    ira, irb = ir.ira, ir.irb
    if ira.length != irb.length:
        raise ValueError(
            f"plastid_annotate: the two inverted-repeat copies are {ira.length} and {irb.length} "
            f"bp — a copy cannot be placed inside a repeat whose sides differ in length"
        )
    inside: list[Gene] = []
    for gene in genes:
        where = _ir_side(gene, ira, irb)
        if where is None:
            inside.append(gene)
            continue
        start, end = ira.start, irb.start
        if where == "IRb":
            start, end = irb.start, ira.start
        mirrored = Gene(
            name=gene.name,
            kind=gene.kind,
            parts=tuple(
                Interval(
                    start=start + (end - part.end),
                    end=start + (end - part.start),
                    strand=-part.strand,
                    partial_start=part.partial_end,
                    partial_end=part.partial_start,
                )
                for part in gene.parts
            ),
            transl_table=gene.transl_table,
            translation=gene.translation,
        )
        inside.append(gene)
        inside.append(mirrored)
    return tuple(sorted(inside, key=lambda gene: gene.parts[0].start))


def _ir_side(gene: Gene, ira: Any, irb: Any) -> str | None:
    """Which copy a gene sits wholly inside, or ``None`` when it straddles or sits outside."""
    for label, region in (("IRa", ira), ("IRb", irb)):
        if all(region.start <= part.start and part.end <= region.end for part in gene.parts):
            return label
    return None


def _verify_round_trip(text: str, annotated: Genome, seq_id: str) -> None:
    """Read the artifact back with our own reader — a file that does not re-read is not a result."""
    with tempfile.TemporaryDirectory(prefix="plastid_annotate_check_") as back_work:
        path = Path(back_work) / f"{seq_id}.gb"
        path.write_text(text, encoding="utf-8")
        back = read_genbank(path)
    before, after = count_genes(annotated), count_genes(back)
    if after != before or back.sequence != annotated.sequence:
        raise ValueError(
            f"plastid_annotate: the GenBank it wrote does not read back as the genome it wrote "
            f"({before.total} genes in, {after.total} out) — the export is broken, and an "
            f"artifact nobody can re-read is not a result"
        )
    if back.source.get("organelle") != ORGANELLE:
        raise ValueError(
            f"plastid_annotate: the artifact re-reads without /organelle={ORGANELLE} "
            f"(source={back.source!r}) — the record this run wrote has to say which organelle it "
            f"is (H-11)"
        )


def _result(
    *,
    genome: Genome,
    cds: Any,
    rrna: Any,
    trna: AragornResult,
    reference: _Reference,
    protein_run: BlastResult,
    rrna_run: BlastResult,
    repeat_run: BlastResult,
    repeat_note: str,
    floors: tuple[float, float],
    evalue: float,
    table: int,
    text: str,
) -> Result:
    calls = (*cds.calls, *rrna.calls)
    not_called = (*cds.not_called, *rrna.not_called)
    dropped = (*cds.dropped, *rrna.dropped)
    unnameable = (*cds.unnameable, *rrna.unnameable)
    competing = (*cds.competing, *rrna.competing)
    min_identity, min_coverage = floors
    checks = check_genome(genome, default_table=table)
    statuses = {check.gene.name: check.status for check in checks}
    overlaps = list(_overlaps(calls))

    summary = Table(
        id="plastid_annotation_summary",
        title="Plastid annotation — what was called, from what, with which floors",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("organelle", "Organelle"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("n_calls", "Genes called", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_not_called", "Vocabulary loci not called", dtype="int"),
            Column("n_cds_failing_orf", "Called CDS failing the ORF checks", dtype="int"),
            Column("n_dropped", "Hits not used", dtype="int"),
            Column("n_unnameable", "Hits this vocabulary cannot name", dtype="int"),
            Column("n_overlapping_pairs", "Overlapping call pairs", dtype="int"),
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
                "organelle": ORGANELLE,
                "length_bp": len(genome),
                "n_calls": len(calls) + len(trna.hits),
                "n_cds": sum(1 for call in calls if call.kind == "CDS"),
                "n_rrna": sum(1 for call in calls if call.kind == "rRNA"),
                "n_trna": len(trna.hits),
                "n_not_called": len(not_called),
                "n_cds_failing_orf": sum(1 for check in checks if check.is_failure),
                "n_dropped": len(dropped),
                "n_unnameable": len(unnameable),
                "n_overlapping_pairs": len(overlaps),
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
        title="Called genes — one row per vocabulary locus, plus the unnamed hits",
        columns=(
            Column("gene", "Gene"),
            Column("kind", "Kind"),
            Column("start", "Start", dtype="int", unit="bp"),
            Column("end", "End", dtype="int", unit="bp"),
            Column("strand", "Strand"),
            Column("exons", "Exons", dtype="int"),
            Column("identity_pct", "Identity", dtype="float", unit="%", precision=2),
            Column("coverage", "Coverage of the reference", dtype="float", precision=3),
            Column("reference", "Reference"),
            Column("partial_start", "Partial 5'"),
            Column("partial_end", "Partial 3'"),
            Column("orf", "ORF check"),
            Column("notes", "Notes"),
        ),
        rows=tuple(
            {
                "gene": call.symbol,
                "kind": call.kind,
                "start": call.span[0],
                "end": call.span[1],
                "strand": _strand(call.strand),
                "exons": len(call.parts),
                "identity_pct": call.identity,
                "coverage": call.coverage,
                "reference": call.reference.accession,
                "partial_start": "yes" if call.partial_start else "no",
                "partial_end": "yes" if call.partial_end else "no",
                "orf": statuses.get(call.symbol, "") if call.kind == "CDS" else "",
                "notes": "; ".join(call.notes),
            }
            for call in (*calls, *unnameable)
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="gene",
    )

    optional: dict[str, Table] = {}
    notes: list[str] = [
        SCOPE_LABEL,
        f"the inverted repeat: {repeat_note}",
        "one call per vocabulary locus; every other reference sequence considered for it is "
        "listed in competing_references with the numbers it lost on",
        *_truncated_notes(reference.truncated),
        *cds.findings,
        *rrna.findings,
    ]

    if not_called:
        optional["not_called"] = Table(
            id="not_called",
            title="Vocabulary loci not called — with the best evidence that was rejected",
            columns=(
                Column("gene", "Gene"),
                Column("reason", "Reason"),
                Column("best_accession", "Best reference"),
                Column("best_identity_pct", "Best identity", dtype="float", unit="%", precision=2),
                Column("best_coverage", "Best coverage", dtype="float", precision=3),
            ),
            rows=tuple(
                {
                    "gene": row.symbol,
                    "reason": row.reason,
                    "best_accession": row.best_accession or "",
                    "best_identity_pct": row.best_identity,
                    "best_coverage": row.best_coverage,
                }
                for row in not_called
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="gene",
        )
    else:
        notes.append(
            "every locus of the reference vocabulary was called — there is no not_called table"
        )

    if competing:
        optional["competing_references"] = Table(
            id="competing_references",
            title="Every reference sequence considered for a locus, and which one was chosen",
            columns=(
                Column("gene", "Gene"),
                Column("reference_accession", "Reference"),
                Column("reference_gene", "Its own spelling"),
                Column("reference_length", "Length", dtype="int", unit="aa"),
                Column("identity_pct", "Identity", dtype="float", unit="%", precision=2),
                Column("coverage", "Coverage", dtype="float", precision=3),
                Column("bitscore", "Bits", dtype="float", precision=1),
                Column("passes_floors", "Passes the floors"),
                Column("chosen", "Chosen"),
            ),
            rows=tuple(
                {
                    "gene": row.symbol,
                    "reference_accession": row.reference_accession,
                    "reference_gene": row.reference_gene,
                    "reference_length": row.reference_length,
                    "identity_pct": row.identity,
                    "coverage": row.coverage,
                    "bitscore": row.bitscore,
                    "passes_floors": "yes" if row.passes_floors else "no",
                    "chosen": "yes" if row.chosen else "no",
                }
                for row in competing
            ),
            scope=Scope(
                label=(
                    "one row per reference sequence that aligned to a vocabulary locus; exactly "
                    "one row per locus is marked chosen — the losing candidates are the evidence "
                    "the choice rests on"
                ),
                coord_system="genbank",
            ),
        )
    else:
        notes.append("no reference sequence aligned to any vocabulary locus — nothing was chosen")

    if dropped:
        optional["rejected_hits"] = Table(
            id="rejected_hits",
            title="Hits not used — one reason each",
            columns=(
                Column("locus", "Locus"),
                Column("span", "Query span", unit="bp"),
                Column("identity_pct", "Identity", dtype="float", unit="%", precision=2),
                Column("aligned", "Aligned", dtype="int"),
                Column("reason", "Reason"),
            ),
            rows=tuple(
                {
                    "locus": hit.symbol or "—",
                    "span": hit.hsp.span_label,
                    "identity_pct": hit.hsp.pident,
                    "aligned": hit.hsp.length,
                    "reason": hit.reason,
                }
                for hit in dropped
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        )
    else:
        notes.append("no hit was rejected: every alignment the searches returned was used")

    if overlaps:
        optional["overlapping_calls"] = Table(
            id="overlapping_calls",
            title="Called genes that overlap",
            columns=(
                Column("first", "Call"),
                Column("second", "Call"),
                Column("overlap_bp", "Overlap", dtype="int", unit="bp"),
                Column("strands", "Strands"),
            ),
            rows=tuple(overlaps),
            scope=Scope(
                label=(
                    "overlap between two calls is reported, never resolved by dropping one — the "
                    "golden itself carries four overlapping pairs in its published annotation"
                ),
                coord_system="genbank",
            ),
        )
    else:
        notes.append("no two calls overlap")

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "plastid_annotation_summary": summary,
            "called_genes": called,
            **optional,
        },
        artifacts={"plastid_annotation.gb": text},
        notes=tuple(notes),
        citations=("blast", "aragorn"),
    )


def _truncated_notes(truncated: Sequence[TruncatedAnnotation]) -> tuple[str, ...]:
    """One line per locus a record annotated twice — the feature the run did **not** search.

    Reported because a reference set that quietly lost a feature would make every "not called"
    row a lie, and because it is the reason a call can land on the complete gene rather than on
    the truncated copy that scores higher (H-19, measured on the golden's ``ycf1``).
    """
    return tuple(
        f"{row.accession} annotates {row.locus} twice: the {row.kept}-residue feature is used "
        f"and its {row.dropped}-residue copy ({row.ratio:.0%} of it) is not searched — structure "
        f"decides before any score (D5-08 §3)"
        for row in truncated
    )


def _overlaps(calls: Sequence[Any]) -> Iterator[dict[str, Any]]:
    """Every pair of calls that share bases — kept and reported, never resolved.

    **Parts, not envelopes.** A trans-spliced call's envelope spans everything between its first
    and last exon — measured on the golden, ``rps12`` is 68,976–98,016, so an envelope test made
    it "overlap" 27 other genes and buried the four pairs that are real (§8.34: psbC×psbD 53 bp ·
    rps12×rps12 114 · ndhC×ndhK 10 · atpB×atpE 4). Two genes overlap when their *exons* do.
    """
    for index, left in enumerate(calls):
        for right in calls[index + 1 :]:
            for left_part in left.parts:
                for right_part in right.parts:
                    lo = max(left_part.start + 1, right_part.start + 1)
                    hi = min(left_part.end, right_part.end)
                    if lo <= hi:
                        yield {
                            "first": f"{left.symbol} {left.span[0]}-{left.span[1]}",
                            "second": f"{right.symbol} {right.span[0]}-{right.span[1]}",
                            "overlap_bp": hi - lo + 1,
                            "strands": (
                                "opposite" if left_part.strand * right_part.strand < 0 else "same"
                            ),
                        }


def _strand(strand: int) -> str:
    return {1: "+", -1: "-"}.get(strand, "mixed")


__all__ = ["MODULE_ID", "MODULE_VERSION", "ORGANELLE", "ReferenceSetMissing", "run"]
