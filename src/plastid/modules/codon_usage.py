"""Codon usage — RSCU · ENC · CAI · CBI · Fop (PROJECT.md §13.2, ROADMAP M2.4).

The CDS set is the counted one: every ``CDS`` feature with ``pseudo`` excluded, both IR copies
included — the same rule as ``count_genes`` (golden genome: 86). Sequences come from
``coding_sequence``, which reads each part on its own strand in annotation order (P2-03 measured
that this is the only reading of the trans-spliced rps12 that yields a canonical ORF).

Translation table: NCBI table 1 (decision D2-03). RSCU and ENC are defined by the sequences
alone; **CAI · CBI · Fop need a reference**, so the reference used is named in the output — here
the analysed CDS set itself. A bias index without its reference is not a number anyone can use.
"""

from __future__ import annotations

from typing import Any

from plastid.core.codon import (
    FAMILIES,
    SENSE_FAMILIES,
    cai,
    cbi,
    codon_counts,
    enc,
    fop,
    rscu,
    translate,
)
from plastid.core.model import Genome, cds_features, coding_sequence
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "codon_usage"
MODULE_VERSION = "0.1.0"

TRANSLATION_TABLE = "NCBI table 1 (standard)"  # D2-03; defaults.yaml records the same choice
REFERENCE_LABEL = "the analysed CDS set itself (self-reference)"
RSCU_TOLERANCE = 1e-9


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="Codon usage",
    requires_ir=True,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ctx.require_ir()
    cdss = cds_features(genome)

    counts = dict.fromkeys((codon for family in SENSE_FAMILIES.values() for codon in family), 0)
    internal_stops = 0
    for gene in cdss:
        sequence = coding_sequence(genome, gene)
        for codon, value in codon_counts(sequence).items():
            counts[codon] += value
        internal_stops += translate(sequence)[:-1].count("*")

    total_codons = sum(counts.values())
    values = rscu(counts)
    family_rows = _family_rows(counts, values)
    unobserved = [row["amino_acid"] for row in family_rows if row["status"] != "ok"]

    summary = Table(
        id="codon_usage",
        title="Codon usage summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS in the set", dtype="int"),
            Column("n_codons", "Codons counted", dtype="int"),
            Column("enc", "ENC", dtype="float", precision=2),
            Column("cai", "CAI", dtype="float", precision=4),
            Column("cbi", "CBI", dtype="float", precision=4),
            Column("fop", "Fop", dtype="float", precision=4),
            Column("internal_stops", "Internal stop codons", dtype="int"),
            Column("translation_table", "Translation table"),
            Column("reference", "Reference for CAI/CBI/Fop"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "n_cds": len(cdss),
                "n_codons": total_codons,
                # raw values: the export layer formats by Column.precision, JSON keeps them exact
                "enc": enc(counts),
                "cai": cai(counts, counts),
                "cbi": cbi(counts, counts),
                "fop": fop(counts, counts),
                "internal_stops": internal_stops,
                "translation_table": TRANSLATION_TABLE,
                "reference": REFERENCE_LABEL,
            },
        ),
        scope=Scope(
            label=(
                "CDS only, pseudogenes excluded, both IR copies included; "
                "stops and partial codons not counted"
            )
        ),
    )

    rscu_table = Table(
        id="rscu",
        title="Relative synonymous codon usage",
        columns=(
            Column("codon", "Codon"),
            Column("amino_acid", "Amino acid"),
            Column("count", "Count", dtype="int"),
            Column("rscu", "RSCU", dtype="float", precision=4),
        ),
        rows=tuple(
            {
                "codon": codon,
                "amino_acid": amino_acid,
                "count": counts[codon],
                "rscu": values[codon],
            }
            for amino_acid, members in SENSE_FAMILIES.items()
            for codon in members
        ),
        scope=Scope(label="CDS only, pseudogenes excluded; stops excluded"),
        key="codon",
    )

    families = Table(
        id="rscu_families",
        title="RSCU family sums (the check: sum of RSCU == family size)",
        columns=(
            Column("amino_acid", "Amino acid"),
            Column("family_size", "Family size", dtype="int"),
            Column("observations", "Codons observed", dtype="int"),
            Column("rscu_sum", "Sum of RSCU", dtype="float", precision=6),
            Column("status", "Check"),
        ),
        rows=tuple(family_rows),
        scope=Scope(label="CDS only, pseudogenes excluded; one row per amino acid"),
        key="amino_acid",
    )

    notes = [
        f"CAI · CBI · Fop reference: {REFERENCE_LABEL}.",
        f"Translation table: {TRANSLATION_TABLE} (decision D2-03).",
    ]
    if unobserved:
        notes.append(f"No observations for: {', '.join(unobserved)} — RSCU undefined there.")
    if internal_stops:
        notes.append(f"{internal_stops} internal stop codon(s) inside CDS — check the annotation.")

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"codon_usage": summary, "rscu": rscu_table, "rscu_families": families},
        notes=tuple(notes),
        citations=("biopython",),
    )


def _family_rows(counts: dict[str, int], values: dict[str, float]) -> list[dict[str, Any]]:
    """One row per amino acid with the family-sum check, which is the point of the table."""
    rows: list[dict[str, Any]] = []
    for amino_acid, members in FAMILIES.items():
        if amino_acid == "*":
            continue
        observations = sum(counts[codon] for codon in members)
        total = sum(values[codon] for codon in members)
        if observations == 0:
            status = "no observations"
        else:
            status = "ok" if abs(total - len(members)) < RSCU_TOLERANCE else "FAIL"
        rows.append(
            {
                "amino_acid": amino_acid,
                "family_size": len(members),
                "observations": observations,
                "rscu_sum": total,
                "status": status,
            }
        )
    return rows
