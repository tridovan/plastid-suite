"""RNA editing — predicted C-to-U sites, PREP-Cp-compatible (M4.5).

The published PREP method (Mower 2005) run by this toolkit against reference B′, its rebuild of
the PREP-Cp training reference (Mower 2009): each supported CDS is translated with the standard
code, added to the gene's fixed alignment of mature reference proteins with MAFFT, and every codon
that contains a C is asked whether a C-to-U edit would agree with more reference residues in its
column. The rule, the per-locus loop and the mapping live in ``core.rna_editing``; this file wires
the aligner and writes the tables. The specification and the acceptance set are
``docs/tools/prep_cp.md``.

What a reader can rely on:

* **Every site is a prediction**, labelled ``predicted · PREP-Cp-compatible`` with B′'s digest and
  the cutoff in every table's scope. Nothing here is original PREP-Cp output, a measured editing
  frequency, or a statement about the input record: the genome, its /translation and its
  annotation are read, never changed.
* **Every CDS gets one status.** ``evaluated`` (its site count may be 0), ``not supported`` (not one
  of the 35 PREP-Cp genes), ``absent from the record`` (a supported gene with no CDS), or ``not
  evaluated`` with the reason. None of the last three is written as 0 sites.
* **Both IR copies are evaluated and mapped**, each as its own locus (``ndhB#1``, ``ndhB#2``):
  prediction is gene-based and needs no IR gate, but a site's genome position is a fact about one
  copy.
* **RNA-seq is not an input.** It is the independent check in the acceptance set, never a step of
  the prediction.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from fractions import Fraction
from typing import Any

from plastid.core.gene_names import aliases_for
from plastid.core.model import Genome
from plastid.core.organelle import GROUP_PLASTID
from plastid.core.rna_editing import (
    DEFAULT_CUTOFF,
    METHOD_LABEL,
    SCOPE_LABEL,
    STATUS_ABSENT,
    STATUS_EVALUATED,
    STATUS_NOT_EVALUATED,
    STATUS_NOT_SUPPORTED,
    TIE_RULE,
    EditingReference,
    LocusOutcome,
    evaluate_genome,
    load_reference,
    location_text,
    locus_segments,
    parse_cutoff,
    reference_set_entry,
)
from plastid.external import mafft
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.params import ParamValueError

MODULE_ID = "rna_editing"
MODULE_VERSION = "0.1.0"

ALIGNER = f"MAFFT {mafft.MAFFT_VERSION} {' '.join(mafft.PROTEIN_ARGV)} --add"


@module(
    MODULE_ID,
    group=GROUP_PLASTID,
    version=MODULE_VERSION,
    title="RNA editing — predicted C-to-U sites (PREP-Cp-compatible)",
    # Gene-based: no position in the inverted repeat enters a score, and each copy is mapped on
    # its own coordinates.
    requires_ir=False,
    citations=("prep_mt", "prep_suite", "mafft", "biopython"),
    reference_set=reference_set_entry,
    params=(
        ParamSpec(
            "cutoff",
            "float",
            "PREP cutoff C in [0, 1]: a site is reported when its edited state scores >= C "
            "(default 0.8, Mower 2009's optimum for seed-plant chloroplasts)",
        ),
    ),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    try:
        cutoff = parse_cutoff(params.get("cutoff", DEFAULT_CUTOFF))
    except ValueError as error:
        raise ParamValueError(f"{MODULE_ID}: {error}") from None
    reference = load_reference()
    cutoff_text = _decimal(cutoff)
    scope = Scope(
        label=(
            f"{SCOPE_LABEL} · reference B′ {reference.digest} · C ≥ {cutoff_text} · {ALIGNER} · "
            f"35 PREP-Cp genes, each CDS copy on its own coordinates"
        ),
        include=("CDS",),
    )
    sample = genome.accession or genome.name
    sequence_sha256 = hashlib.sha256(genome.sequence.encode("ascii")).hexdigest()

    outcomes = evaluate_genome(
        genome,
        reference,
        cutoff,
        aliases_for(genome),
        lambda rows, query: mafft.add_protein(rows, query, ctx),
    )
    coverage = [_coverage_row(outcome) for outcome in outcomes]
    query_rows = [
        _query_row(sample, sequence_sha256, outcome)
        for outcome in outcomes
        if outcome.prepared is not None
    ]
    site_rows = [
        row
        for outcome in outcomes
        if outcome.prediction is not None
        for row in _site_rows(sample, outcome, cutoff_text, reference.digest)
    ]
    artifacts = {
        f"alignments/{_file_stem(outcome.locus)}.fasta": "".join(
            f">{name}\n{row}\n" for name, row in outcome.alignment
        )
        for outcome in outcomes
        if outcome.alignment
    }

    statuses = Counter(outcome.status for outcome in outcomes)
    tables = {
        "rna_editing_summary": _summary_table(
            genome, sample, reference, cutoff_text, statuses, site_rows, scope
        ),
        "rna_editing_coverage": _coverage_table(coverage, scope),
    }
    if site_rows:
        tables["rna_editing_sites"] = _sites_table(site_rows, scope)
    if query_rows:
        tables["rna_editing_query"] = _query_table(query_rows, scope)

    notes = [
        f"Every site is predicted and {METHOD_LABEL}: the published PREP method (Mower 2005) run "
        f"by this toolkit against reference B′, never output of the original PREP-Cp, whose "
        f"service ended in 2022.",
        f"Reference B′ {reference.digest}: {reference.label}.",
        f"Cutoff C = {cutoff_text}; a site is reported when its edited state scores >= C. Scores "
        f"are reference residues matching the state's amino acid over reference residues in the "
        f"column, gaps excluded.",
        f"Ties: {TIE_RULE}.",
        f"Alignment: the query protein is added to B′'s fixed alignment with {ALIGNER}; the "
        f"reference rows are checked unchanged on every call.",
        "RNA-seq is not used. The genome, its /translation and its annotation are not changed.",
        *_edit_list_notes(reference),
    ]
    if not site_rows:
        notes.append("No site was predicted at this cutoff; rna_editing_coverage gives each gene.")
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        citations=("prep_mt", "prep_suite", "mafft", "biopython"),
        artifacts=artifacts,
    )


def _coverage_row(outcome: LocusOutcome) -> dict[str, Any]:
    gene = outcome.gene
    row: dict[str, Any] = {
        "locus": outcome.locus,
        "name_in_record": outcome.name_in_record,
        "gene_a": outcome.gene_a,
        "prep_gene": outcome.prep_gene or "",
        "status": outcome.status,
        "reason": outcome.reason,
        "partial": "" if gene is None else ("yes" if gene.partial else "no"),
        "location": "" if gene is None else location_text(gene),
    }
    prediction = outcome.prediction
    if prediction is not None:
        # Only an evaluated locus has counts; every other status leaves them blank, never 0.
        row.update(
            {
                "n_sites": len(prediction.sites),
                "n_c_codons": prediction.c_codons,
                "n_scored": prediction.scored,
                "unscored_codons": _numbers(prediction.unscored),
                "ambiguous_codons": _numbers(prediction.ambiguous),
            }
        )
    return row


def _site_rows(
    sample: str, outcome: LocusOutcome, cutoff_text: str, digest: str
) -> list[dict[str, Any]]:
    prepared, prediction = outcome.prepared, outcome.prediction
    if prepared is None or prediction is None:
        return []
    rows = []
    for site in prediction.sites:
        position, strand = prepared.positions[site.cds_position - 1]
        call = site.call
        rows.append(
            {
                "site": f"{outcome.locus}:{site.cds_position}",
                "sample": sample,
                "locus": outcome.locus,
                "gene_a": outcome.gene_a,
                "prep_gene": outcome.prep_gene,
                "cds_position": site.cds_position,
                "aa_position": site.codon_number,
                "codon_position": site.codon_position,
                "alignment_column": site.column,
                "genome_position": position + 1,
                "strand": "+" if strand > 0 else "-",
                # A minus-strand C-to-U is a G-to-A on the forward strand; both are stated.
                "genome_base": "C" if strand > 0 else "G",
                "forward_change": "C>T" if strand > 0 else "G>A",
                "codon_before": call.unedited.codon,
                "codon_after": call.state.codon,
                "aa_before": call.unedited.amino_acid,
                "aa_after": call.state.amino_acid,
                "effect": (
                    f"{call.unedited.codon} ({call.unedited.amino_acid}) => "
                    f"{call.state.codon} ({call.state.amino_acid})"
                ),
                "score": float(call.score),
                "matching": call.matching,
                "compared": call.compared,
                "cutoff": cutoff_text,
                "reference_digest": digest,
                "status": SCOPE_LABEL,
            }
        )
    return rows


def _query_row(sample: str, sequence_sha256: str, outcome: LocusOutcome) -> dict[str, Any]:
    gene, prepared = outcome.gene, outcome.prepared
    if gene is None or prepared is None:
        raise ValueError(f"{outcome.locus}: a query row needs an extracted CDS")
    return {
        "locus": outcome.locus,
        "sample": sample,
        "sequence_sha256": sequence_sha256,
        "name_in_record": outcome.name_in_record,
        "gene_a": outcome.gene_a,
        "prep_gene": outcome.prep_gene,
        "location": location_text(gene),
        "strand": {1: "+", -1: "-"}.get(gene.strand, "mixed"),
        "codon_start": gene.codon_start,
        "transl_table": "" if gene.transl_table is None else str(gene.transl_table),
        "partial": "yes" if gene.partial else "no",
        "cds_length": len(prepared.cds),
        "cds_sha256": hashlib.sha256(prepared.cds.encode("ascii")).hexdigest(),
        "protein_length": len(prepared.protein),
        "mapping": locus_segments(prepared.positions),
    }


def _summary_table(
    genome: Genome,
    sample: str,
    reference: EditingReference,
    cutoff_text: str,
    statuses: Counter[str],
    site_rows: list[dict[str, Any]],
    scope: Scope,
) -> Table:
    return Table(
        id="rna_editing_summary",
        title="RNA editing — predicted C-to-U sites (PREP-Cp-compatible)",
        columns=(
            Column("sample", "Sample"),
            Column("organism", "Organism"),
            Column("method", "Method"),
            Column("reference", "Reference"),
            Column("reference_digest", "Reference digest"),
            Column("cutoff", "Cutoff C"),
            Column("aligner", "Aligner"),
            Column("tie_rule", "Tie rule"),
            Column("n_evaluated", "Loci evaluated", dtype="int"),
            Column("n_not_evaluated", "Loci not evaluated", dtype="int"),
            Column("n_absent", "Supported genes absent", dtype="int"),
            Column("n_not_supported", "CDS not supported", dtype="int"),
            Column("n_sites", "Predicted sites", dtype="int"),
            Column("n_loci_with_sites", "Loci with a predicted site", dtype="int"),
        ),
        rows=(
            {
                "sample": sample,
                "organism": genome.name,
                "method": SCOPE_LABEL,
                "reference": reference.label,
                "reference_digest": reference.digest,
                "cutoff": cutoff_text,
                "aligner": ALIGNER,
                "tie_rule": TIE_RULE,
                "n_evaluated": statuses[STATUS_EVALUATED],
                "n_not_evaluated": statuses[STATUS_NOT_EVALUATED],
                "n_absent": statuses[STATUS_ABSENT],
                "n_not_supported": statuses[STATUS_NOT_SUPPORTED],
                "n_sites": len(site_rows),
                "n_loci_with_sites": len({row["locus"] for row in site_rows}),
            },
        ),
        scope=scope,
    )


def _coverage_table(rows: list[dict[str, Any]], scope: Scope) -> Table:
    return Table(
        id="rna_editing_coverage",
        title="RNA editing — every CDS and every supported gene, with its status",
        columns=(
            Column("locus", "Locus"),
            Column("name_in_record", "Name in the record"),
            Column("gene_a", "Name (vocabulary A)"),
            Column("prep_gene", "PREP-Cp gene"),
            Column("status", "Status"),
            Column("reason", "Reason"),
            Column("n_sites", "Predicted sites", dtype="int"),
            Column("n_c_codons", "Codons with a C", dtype="int"),
            Column("n_scored", "Codons scored", dtype="int"),
            Column("unscored_codons", "C codons with no reference residue"),
            Column("ambiguous_codons", "Codons with an ambiguity code (not evaluated)"),
            Column("partial", "Partial"),
            Column("location", "Location (1-based, transcript order)"),
        ),
        rows=tuple(rows),
        scope=scope,
        key="locus",
    )


def _sites_table(rows: list[dict[str, Any]], scope: Scope) -> Table:
    return Table(
        id="rna_editing_sites",
        title="RNA editing — predicted C-to-U sites",
        columns=(
            Column("site", "Site"),
            Column("sample", "Sample"),
            Column("locus", "Locus"),
            Column("gene_a", "Name (vocabulary A)"),
            Column("prep_gene", "PREP-Cp gene"),
            Column("cds_position", "CDS position (nt)", dtype="int"),
            Column("aa_position", "Amino-acid position", dtype="int"),
            Column("codon_position", "Codon position", dtype="int"),
            Column("alignment_column", "Alignment column", dtype="int"),
            Column("genome_position", "Genome position (1-based)", dtype="int"),
            Column("strand", "Strand"),
            Column("genome_base", "Base on the forward strand"),
            Column("forward_change", "Change on the forward strand"),
            Column("codon_before", "Codon"),
            Column("codon_after", "Edited codon"),
            Column("aa_before", "Amino acid"),
            Column("aa_after", "Edited amino acid"),
            Column("effect", "Effect"),
            Column("score", "Score", dtype="float", precision=2),
            Column("matching", "Reference residues matching", dtype="int"),
            Column("compared", "Reference residues in the column", dtype="int"),
            Column("cutoff", "Cutoff C"),
            Column("reference_digest", "Reference digest"),
            Column("status", "Status"),
        ),
        rows=tuple(rows),
        scope=scope,
        key="site",
    )


def _query_table(rows: list[dict[str, Any]], scope: Scope) -> Table:
    return Table(
        id="rna_editing_query",
        title="RNA editing — the query manifest: what was extracted from each supported CDS",
        columns=(
            Column("locus", "Locus"),
            Column("sample", "Sample"),
            Column("sequence_sha256", "Sequence SHA-256"),
            Column("name_in_record", "Name in the record"),
            Column("gene_a", "Name (vocabulary A)"),
            Column("prep_gene", "PREP-Cp gene"),
            Column("location", "Parts (1-based, transcript order)"),
            Column("strand", "Strand"),
            Column("codon_start", "codon_start", dtype="int"),
            Column("transl_table", "transl_table"),
            Column("partial", "Partial"),
            Column("cds_length", "CDS length (nt)", dtype="int"),
            Column("cds_sha256", "CDS SHA-256"),
            Column("protein_length", "Protein length (aa)", dtype="int"),
            Column("mapping", "CDS → genome (1-based)"),
        ),
        rows=tuple(rows),
        scope=scope,
        key="locus",
    )


def _edit_list_notes(reference: EditingReference) -> list[str]:
    """One note per B′ taxon: how many sourced edit sites it carries against 2009's count."""
    counts = Counter(site["accession"] for site in reference.payload["edit_sites"])
    notes = []
    for taxon in reference.payload["taxa"]:
        accession = taxon["accession"]
        notes.append(
            f"B′ edit list, {taxon['organism']} ({accession}): {counts[accession]} sourced sites "
            f"against {taxon['known_2009']} in PREP-Cp's own reference (Mower 2009, Table 2) — "
            f"{taxon['edit_list']}"
        )
    return notes


def _file_stem(locus_id: str) -> str:
    return locus_id.replace("#", "_")


def _numbers(values: tuple[int, ...]) -> str:
    return ",".join(str(value) for value in values)


def _decimal(value: Fraction) -> str:
    """A cutoff as the decimal it was given as (4/5 → 0.8, 1 → 1)."""
    text = f"{value.numerator / value.denominator:.6f}".rstrip("0").rstrip(".")
    return text or "0"
