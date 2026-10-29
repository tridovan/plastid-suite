"""ORF checks — report, never repair (M3.3, P3-01).

Every ``CDS`` feature is checked against its own declared translation table and reported with
one headline status; ``gene`` features no CDS covers (the truncated-ycf1
shape) get their own table rather than a failure. Nothing here rewrites a CDS to make it pass:
a sequence that disagrees with its annotation is the finding.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from typing import Any

from plastid.core.model import Gene, Genome, count_genes
from plastid.core.orf import (
    ORF_STATUSES,
    STATUS_ALTERNATIVE_START,
    STATUS_INTERNAL_STOP,
    STATUS_NO_START,
    STATUS_NO_STOP,
    STATUS_NOT_MULTIPLE_OF_THREE,
    STATUS_OK,
    STATUS_PSEUDOGENE,
    OrfCheck,
    annotated_not_coding,
    check_genome,
)
from plastid.core.organelle import GROUP_SHARED
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults

MODULE_ID = "orf"
MODULE_VERSION = "0.2.0"

SCOPE_LABEL = "every CDS feature, both IR copies, pseudogenes included; GenBank coordinates"

_STATUS_MEANING = {
    STATUS_NOT_MULTIPLE_OF_THREE: "length is not a multiple of three — no frame to translate",
    STATUS_NO_START: "first codon is not an initiation codon of the declared table",
    STATUS_NO_STOP: "last codon is not a stop codon",
    STATUS_INTERNAL_STOP: "in-frame stop codon before the end",
    STATUS_ALTERNATIVE_START: "non-ATG initiation codon — flagged, not a failure",
    STATUS_PSEUDOGENE: "annotated pseudogene — kept out of the failure counts",
    STATUS_OK: "ATG start, stop codon at the end, no internal stop",
}


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="ORF checks",
    requires_ir=False,  # reading frames come from the sequence, not from IR positions
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    default_table = int(load_defaults()["orf"]["default_translation_table"])
    checks = check_genome(genome, default_table=default_table)
    orphans = annotated_not_coding(genome)
    cds_ids = _cds_ids(checks)
    counts = count_genes(genome)

    status_counts = Counter(check.status for check in checks)
    failed = sum(1 for check in checks if check.is_failure)
    tables_used = sorted({check.transl_table for check in checks})

    summary = Table(
        id="orf_summary",
        title="ORF check summary",
        columns=(
            Column("organism", "Organism"),
            Column("accession", "Accession"),
            Column("n_cds", "CDS checked", dtype="int"),
            Column("n_failed", "Failed", dtype="int"),
            Column("n_pseudogene", "Pseudogenes (excluded from failures)", dtype="int"),
            Column("n_annotated_not_coding", "Annotated, not coding", dtype="int"),
            Column("n_duplicated_cds", "CDS features sharing a name and span", dtype="int"),
            Column("translation_tables", "Translation table(s) used"),
            Column("coord_system", "Coordinates"),
        ),
        rows=(
            {
                "organism": genome.name,
                "accession": genome.accession or "",
                "n_cds": len(checks),
                "n_failed": failed,
                "n_pseudogene": status_counts[STATUS_PSEUDOGENE],
                "n_annotated_not_coding": len(orphans),
                "n_duplicated_cds": sum(1 for cds_id in cds_ids if "#" in cds_id),
                "translation_tables": ", ".join(str(table) for table in tables_used),
                "coord_system": "genbank",
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )

    status_table = Table(
        id="orf_status",
        title="CDS by ORF status",
        columns=(
            Column("status", "Status"),
            Column("meaning", "Meaning"),
            Column("n", "CDS", dtype="int"),
        ),
        rows=tuple(
            {
                "status": status,
                "meaning": _STATUS_MEANING[status],
                "n": status_counts[status],
            }
            for status in ORF_STATUSES
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="status",
    )

    tables: dict[str, Table] = {
        summary.id: summary,
        status_table.id: status_table,
    }
    if checks:
        tables["orf_checks"] = Table(
            id="orf_checks",
            title="ORF check per CDS",
            columns=(
                Column("cds_id", "CDS"),
                Column("gene", "Gene"),
                Column("start_1b", "Start", dtype="int"),
                Column("end_1b", "End", dtype="int"),
                Column("strand", "Strand"),
                Column("length_bp", "Length", dtype="int", unit="bp"),
                Column("transl_table", "Table", dtype="int"),
                Column("codon_start", "Codon start", dtype="int"),
                Column("start_codon", "Start codon"),
                Column("start_class", "Start class"),
                Column("stop_codon", "Stop codon"),
                Column("internal_stops", "Internal stops", dtype="int"),
                Column("internal_stop_codons", "Internal stop positions (codon)"),
                Column("status", "Status"),
                Column("pseudo", "Pseudogene"),
                Column("partial", "Partial CDS"),
            ),
            rows=tuple(
                _check_row(check, cds_id) for check, cds_id in zip(checks, cds_ids, strict=True)
            ),
            scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
            key="cds_id",
        )
    # a record with no CDS is a legitimate input — an unannotated draft, or the target of an
    # annotation run. It gets a summary of zeros and the note below, and no per-CDS table: an
    # empty one cannot be built (Table refuses zero rows) and could not be read if it could.
    # Measured 04/10/2026: constructing the table unconditionally made the module raise *before*
    # this branch, so the intended handling of a bare sequence had never actually run (P4a-04).
    if orphans:
        tables["annotated_not_coding"] = _orphan_table(orphans)

    tables["orf_reconciliation"] = Table(
        id="orf_reconciliation",
        title="ORF checks reconciled against the genome-wide CDS count",
        columns=(
            Column("kind", "Kind"),
            Column("in_table", "Checked here", dtype="int"),
            Column("counted", "count_genes()", dtype="int"),
            Column("agreement", "Agreement"),
        ),
        rows=(
            {
                "kind": "CDS",
                "in_table": len(checks),
                "counted": counts.cds,
                "agreement": "agrees" if len(checks) == counts.cds else "DISAGREES",
            },
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
        key="kind",
    )

    notes = ["Pseudogenes are reported as their own class and excluded from the failure count."]
    if counts.cds == 0:
        notes.append("This record carries no CDS feature: sequence-level analyses only.")
    duplicated = _duplicated_ids(cds_ids)
    if duplicated:
        notes.append(
            "This record annotates the same gene more than once at one span, and both are "
            "checked: " + ", ".join(duplicated) + " (ids suffixed with #2, #3, …)"
        )
    non_canonical = sorted(
        {check.gene.name for check in checks if check.start_class == "non_canonical"}
    )
    if non_canonical:
        notes.append(
            f"Non-canonical start codon (not in the declared table, usually edited to ATG in "
            f"the transcript): {', '.join(non_canonical)}"
        )
    if orphans:
        notes.append(
            "Annotated without a CDS (reported, not counted as a failure): "
            + ", ".join(sorted({gene.name for gene in orphans}))
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
    )


def _cds_ids(checks: Sequence[OrfCheck]) -> list[str]:
    """One id per CDS **feature**, unique even when a record annotates one gene twice.

    A record can carry two CDS features that collapse to the same span — *Libertia pulchella*
    NC_072575.1 annotates ``rps16`` twice, at ``5741-5780,4684-4895`` and ``5737-5780,4684-4891``
    (measured 04/10/2026, P4a-03). Keying on name and span alone made those two rows share an id,
    and ``Table`` refuses a repeated key, so the whole run died on a defect that belongs to the
    record. Suffixing keeps both rows — a duplicated annotation is a finding to report, not to
    drop — and the notes name it.
    """
    ids: list[str] = []
    seen: Counter[str] = Counter()
    for check in checks:
        start, end = check.gene.span
        base = f"{check.gene.name}:{start}-{end}"
        seen[base] += 1
        ids.append(base if seen[base] == 1 else f"{base}#{seen[base]}")
    return ids


def _duplicated_ids(cds_ids: Sequence[str]) -> list[str]:
    """The bases that needed a suffix — i.e. the features the record states more than once."""
    return sorted({cds_id.split("#", 1)[0] for cds_id in cds_ids if "#" in cds_id})


def _check_row(check: OrfCheck, cds_id: str) -> dict[str, Any]:
    start, end = check.gene.span
    return {
        "cds_id": cds_id,
        "gene": check.gene.name,
        "start_1b": start,
        "end_1b": end,
        "strand": "+" if check.gene.strand >= 0 else "-",
        "length_bp": check.length_bp,
        "transl_table": check.transl_table,
        "codon_start": check.codon_start,
        "start_codon": check.start_codon,
        "start_class": check.start_class,
        "stop_codon": check.stop_codon,
        "internal_stops": len(check.internal_stop_positions),
        "internal_stop_codons": ", ".join(
            str(position) for position in check.internal_stop_positions
        ),
        "status": check.status,
        "pseudo": "yes" if check.gene.pseudo else "no",
        "partial": "yes" if check.gene.partial else "no",
    }


def _orphan_table(orphans: tuple[Gene, ...]) -> Table:
    """One row per ``gene`` feature no CDS covers — the truncated-ycf1 shape."""
    return Table(
        id="annotated_not_coding",
        title="Annotated as a gene, not coding",
        columns=(
            Column("name", "Name"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("strand", "Strand"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("product", "Product"),
            Column("pseudo", "Pseudogene"),
        ),
        rows=tuple(
            {
                "name": gene.name,
                "start_1b": gene.span[0],
                "end_1b": gene.span[1],
                "strand": "+" if gene.strand >= 0 else "-",
                "length_bp": gene.length,
                "product": gene.product or "",
                "pseudo": "yes" if gene.pseudo else "no",
            }
            for gene in orphans
        ),
        scope=Scope(label=SCOPE_LABEL, coord_system="genbank"),
    )
