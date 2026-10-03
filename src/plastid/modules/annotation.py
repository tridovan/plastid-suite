"""Annotation of a plastome (M3.1, P3-03; annotator per D3-01).

The module's job is to make an annotator's output *enter the same pipeline as a downloaded
record*: the annotator writes its features, our own reader reads them back, and the gene set, its
ORF status (P3-01) and its name mapping (P3-02) come from the same core functions every other
module uses. There is no second code path — an annotation that skipped the checks would be the one
output nobody verified.

Two annotators can answer, and the report says which one did:

* ``chloe`` (default) — offline, angiosperm plastomes, its own curated references;
* ``pga`` — homology against a directory of GenBank references the caller supplies, which is the
  only option for a genome Chloe does not cover (it annotates flowering plants only).

Both are measured rather than trusted: PGA returned 147 genes where the golden record has 132 and
failed silently three ways (§8.15), and Chloe returns 132 = 86 + 38 + 8 with every CDS translating
— but only when it is told not to transform the sequence (§8.17, ``external/chloe.py``).

The annotators are external and optional: a machine with neither still runs everything else.
``PTS_REQUIRE_TOOLS=1`` turns a missing tool into a failure for a machine that is supposed to
have it.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from plastid.core.gene_groups import assign, load_gene_groups
from plastid.core.gene_names import GenomeNames, load_gene_aliases, normalize_genome
from plastid.core.model import GENE_KINDS, Genome, count_genes
from plastid.core.orf import ORF_STATUSES, STATUS_PSEUDOGENE, check_genome
from plastid.external.chloe import run_chloe
from plastid.external.pga import run_pga
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults
from plastid.run.layout import slugify

MODULE_ID = "annotation"
MODULE_VERSION = "0.2.0"

SCOPE_LABEL = "the annotator's gene set, read back with our reader; GenBank coordinates"

ANNOTATORS: tuple[str, ...] = ("chloe", "pga")
DEFAULT_ANNOTATOR = "chloe"

# Plastids translate with NCBI table 11, and D3-03 measured that all 1,099 CDS in this project's
# dataset declare it. Chloe's GFF declares no table at all, so without this default its CDS would
# be checked against table 1 — which knows fewer initiation codons, and would report a curated
# GTG or ACG start as a failure of the annotator. The value is stated in the run's notes either
# way, and --set transl_table= overrides it (mitochondria use table 1).
PLASTID_TRANSL_TABLE = 11


@dataclass(frozen=True, slots=True)
class _Annotated:
    """What an annotator handed back, in the one shape this module analyses."""

    genome: Genome
    text: str
    source_name: str
    extension: str  # the annotator's own format, for the artifact the run keeps
    annotator: str  # the label the report prints
    reference: str  # what the annotation rests on
    citation: str
    notes: tuple[str, ...] = field(default=())


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="Annotation",
    requires_ir=False,  # annotating is not a position-based analysis
    citations=("chloe", "pga"),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    table = int(params.get("transl_table", PLASTID_TRANSL_TABLE))
    produced = _annotate(genome, params, ctx, table=table)

    annotated = produced.genome
    counts = count_genes(annotated)
    # Every CDS a silent output produced was stamped with `table` by its adapter (Chloe's GFF
    # declares none); a record that declares its own keeps it (PGA's GenBank carries
    # /transl_table). The config default is what is left for a feature that declares nothing at
    # all — the same fallback D3-03 records for a downloaded record.
    checks = check_genome(
        annotated, default_table=int(load_defaults()["orf"]["default_translation_table"])
    )
    names = normalize_genome(annotated, load_gene_aliases())
    groups = load_gene_groups()
    rows = assign(annotated, groups, canonical_of=names.canonical)
    unclassified = sorted({names.label(row.gene) for row in rows if row.group is None})

    status_counts = Counter(check.status for check in checks)
    summary = Table(
        id="annotation_summary",
        title="Annotation summary",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("n_genes", "Genes", dtype="int"),
            Column("n_cds", "CDS", dtype="int"),
            Column("n_trna", "tRNA", dtype="int"),
            Column("n_rrna", "rRNA", dtype="int"),
            Column("n_orf_failed", "CDS failing the ORF checks", dtype="int"),
            Column("n_aliased_names", "Names with a synonym", dtype="int"),
            Column("n_unclassified", "Unclassified", dtype="int"),
            Column("annotator", "Annotator"),
            Column("reference", "References"),
        ),
        rows=(
            {
                "organism": annotated.name,
                "accession": annotated.accession or "",
                "length_bp": len(annotated),
                "n_genes": counts.total,
                "n_cds": counts.cds,
                "n_trna": counts.trna,
                "n_rrna": counts.rrna,
                "n_orf_failed": sum(1 for check in checks if check.is_failure),
                "n_aliased_names": sum(
                    1
                    for gene in annotated.genes
                    if gene.kind in GENE_KINDS and names.normalized(gene).synonyms
                ),
                "n_unclassified": len(unclassified),
                "annotator": produced.annotator,
                "reference": produced.reference,
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    # An annotator can call one locus twice under two names — measured with PGA 2026-10-03: the
    # same clpP/clpP1 span, because some references use the classic names and some the functional
    # ones. That is a finding about the annotation, so it is reported, never merged away and never
    # a crash: the row id carries the name it arrived with, and duplicates get their own table.
    occurrences: Counter[tuple[str, str, int, int]] = Counter()
    gene_rows = []
    for row in rows:
        gene = row.gene
        identity = (gene.kind, gene.name, *gene.span)
        occurrences[identity] += 1
        gene_rows.append(_gene_row(row, names, occurrences[identity]))
    duplicates = _duplicate_rows(rows, names)
    gene_set = Table(
        id="gene_set",
        title="Annotated gene set",
        columns=(
            Column("gene_id", "Gene id"),
            Column("gene", "Gene"),
            Column("kind", "Kind"),
            Column("category", "Category"),
            Column("group", "Group of genes"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("product", "Product"),
        ),
        rows=gene_rows,
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="gene_id",
    )

    status_table = Table(
        id="orf_status",
        title="ORF status of the annotated CDS",
        columns=(
            Column("status", "Status"),
            Column("n", "CDS", dtype="int"),
        ),
        rows=tuple({"status": status, "n": status_counts[status]} for status in ORF_STATUSES),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="status",
    )

    tables: dict[str, Table] = {
        summary.id: summary,
        gene_set.id: gene_set,
        status_table.id: status_table,
    }
    if duplicates:
        tables["duplicate_loci"] = Table(
            id="duplicate_loci",
            title="Loci the annotator called more than once",
            columns=(
                Column("canonical", "Canonical name"),
                Column("names", "Names it arrived under"),
                Column("n", "Calls", dtype="int"),
                Column("groups", "Functional group(s)"),
            ),
            rows=duplicates,
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="canonical",
        )

    if unclassified:
        tables["unclassified"] = Table(
            id="unclassified",
            title="Genes no group claims",
            columns=(Column("gene", "Gene"),),
            rows=tuple({"gene": name} for name in unclassified),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        )

    notes = [
        "The annotator's output was read back with our own reader and then went through the same "
        "ORF checks and name normalization as a downloaded record.",
        f"Raw annotator output kept beside the tables as {produced.source_name}.",
        f"CDS were checked against translation table(s) "
        f"{', '.join(str(number) for number in sorted({check.transl_table for check in checks}))}"
        f": a record that declares its own (PGA writes /transl_table) keeps it, and an output "
        f"that is silent (Chloe's GFF says nothing) is stamped with {table} — pass "
        f"--set transl_table= to change that.",
        *produced.notes,
    ]
    if duplicates:
        notes.append(
            f"{len(duplicates)} locus/loci were annotated more than once "
            f"({', '.join(row['canonical'] for row in duplicates)}) — reported in "
            f"duplicate_loci, not merged: a second call under a second name is a statement about "
            f"the annotator's references, not about the genome."
        )
    if status_counts[STATUS_PSEUDOGENE]:
        notes.append(
            f"{status_counts[STATUS_PSEUDOGENE]} CDS were annotated as pseudogenes and are kept "
            f"out of the failure count."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        citations=(produced.citation,),
        # the annotator's own file, byte for byte, under the format it writes: the report names
        # it, so any number above can be traced back to the line that produced it
        artifacts={f"annotation.{produced.extension}": produced.text},
    )


def _annotate(genome: Genome, params: dict[str, Any], ctx: RunContext, *, table: int) -> _Annotated:
    """Run the requested annotator and hand back one shape, whichever one answered.

    Each annotator's own failures stay its own (they are raised by its adapter, with the fix in
    the message); what this function decides is only *which* one runs and how its provenance is
    worded for the report.
    """
    annotator = str(params.get("annotator", DEFAULT_ANNOTATOR)).lower()
    if annotator not in ANNOTATORS:
        raise ValueError(
            f"annotation: unknown annotator {annotator!r} — available: {', '.join(ANNOTATORS)}"
        )
    seq_id = genome.accession or slugify(genome.name)

    if annotator == "chloe":
        produced = run_chloe(
            genome,
            ctx,
            seq_id=seq_id,
            script=params.get("chloe"),
            julia=params.get("julia"),
            transl_table=table,
        )
        return _Annotated(
            genome=produced.genome,
            text=produced.gff_text,
            source_name=produced.source_name,
            extension="gff",
            annotator=f"Chloe {produced.version}",
            reference=f"its own curated set (Chloe_refs {produced.reference_tree[:12]})",
            citation="chloe",
            notes=(
                "The annotator was run with --no-transform, so its coordinates are this "
                "sequence's own: by default Chloe flips the single-copy region into its standard "
                "configuration and writes a GFF in *that* frame, which no other module of this "
                "run could line up with (measured 2026-10-03 — the gene set is identical either "
                "way, so a gene count alone cannot see the difference).",
            ),
        )

    reference = params.get("reference")
    if not reference:
        raise ValueError(
            "annotation: PGA needs reference plastomes — pass --set "
            "reference=<directory of GenBank files from the same or a close family> (it "
            "annotates by homology), or use the default annotator, which brings its own"
        )
    pga = run_pga(
        genome,
        ctx,
        reference_dir=str(reference),
        seq_id=seq_id,
        script=params.get("pga"),
    )
    notes = (
        (
            "The annotator produced nothing on its first attempt and the run was retried once: "
            "PGA's temporary file names can contain a character the shell it calls blastn "
            "through expands away (measured, about one run in five). The retry succeeded.",
        )
        if pga.attempts > 1
        else ()
    )
    return _Annotated(
        genome=pga.genome,
        text=pga.genbank_text,
        source_name=pga.source_name,
        extension="gb",
        annotator=f"PGA {pga.source_name}",
        reference=str(reference),
        citation="pga",
        notes=notes,
    )


def _gene_row(row, names: GenomeNames, occurrence: int) -> dict[str, Any]:
    """One annotated gene. The id carries the name it arrived with, so two calls of one locus
    under two names stay two addressable rows instead of colliding on a canonical key; the two
    copies of an inverted-repeat gene share a name and a span, so the second gets a suffix."""
    gene = row.gene
    start, end = gene.span
    suffix = f"#{occurrence}" if occurrence > 1 else ""
    return {
        "gene_id": f"{gene.kind}:{gene.name}:{start}-{end}{suffix}",
        "gene": names.label(gene),
        "kind": gene.kind,
        "group": row.group.name if row.group else "",
        "category": row.group.category if row.group else "",
        "length_bp": gene.length,
        "product": gene.product or "",
    }


def _duplicate_rows(rows, names: GenomeNames) -> tuple[dict[str, Any], ...]:
    """Loci the annotator called under more than one name — the honest report of over-calling.

    Measured 2026-10-03 with PGA and seven references: it emitted the same clpP locus twice, once
    as ``clpP`` (references using the classic name) and once as ``clpP1`` (references using the
    functional one). Merging those rows would hide a real property of the annotation; dropping one
    would hide the other name. So they are reported — and only those: the same name twice is the
    two copies of an inverted-repeat gene, which every plastome has.
    """
    by_canonical: dict[str, list] = {}
    for row in rows:
        by_canonical.setdefault(names.canonical(row.gene), []).append(row)
    return tuple(
        {
            "canonical": canonical,
            "names": ", ".join(sorted({row.gene.name for row in calls})),
            "n": len(calls),
            "groups": ", ".join(sorted({row.group.name for row in calls if row.group})),
        }
        for canonical, calls in sorted(by_canonical.items())
        # only *different names* for one locus: a plastome carries two copies of every IR gene,
        # and reporting those as duplicates would bury the finding that matters (measured
        # 2026-10-03: 25 "duplicates", of which 21 were the two inverted-repeat copies)
        if len({row.gene.name for row in calls}) > 1
    )
