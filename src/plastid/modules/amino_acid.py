"""Amino-acid composition, GRAVY and aromaticity (ROADMAP M2.5).

The CDS set comes from ``cds_features`` — the same source P2-03 uses, so the two modules can
never disagree about which sequences they describe. Every translation is built with
``coding_sequence`` and NCBI table 1, compiled here (D2-03) — and the output names the table the
record's own CDS declare (plastid records declare 11; the two share every amino acid and differ
only in initiation codons, and the initiator is read as its own amino acid).

**No inverted repeat is needed.** Composition is a property of the translated CDS set, not a
position; the module no longer waits behind the IR gate (M11.5). A frame or declared table this
compiled code cannot read is refused by name, never translated silently.

Cross-check that must hold (and is asserted in the golden tests): the residue total equals the
sense-codon total of `codon_usage`, because every sense codon contributes exactly one residue.
"""

from __future__ import annotations

from typing import Any

from plastid.core.codon import require_translatable, translate, translation_table_label
from plastid.core.model import Genome, cds_features, coding_sequence
from plastid.core.organelle import GROUP_SHARED
from plastid.core.protein import (
    KYTE_DOOLITTLE,
    RESIDUES,
    aromaticity,
    composition,
    gravy,
    ignored_symbols,
)
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "amino_acid"
MODULE_VERSION = "0.2.0"  # 0.2.0: no IR requirement; unsupported frame/table refused (M11.5)


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Amino-acid composition",
    # CDS sequences only: no position in the inverted repeat enters any number here (M11.5).
    requires_ir=False,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    cdss = cds_features(genome)
    require_translatable(cdss, module=MODULE_ID)
    table_label = translation_table_label(cdss)

    translation = "".join(translate(coding_sequence(genome, gene)) for gene in cdss)
    counts = composition(translation)
    residues = sum(counts.values())
    ignored = ignored_symbols(translation)

    summary = Table(
        id="amino_acid",
        title="Amino-acid summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("n_cds", "CDS in the set", dtype="int"),
            Column("n_residues", "Residues counted", dtype="int"),
            Column("ignored_symbols", "Stops and ambiguous symbols", dtype="int"),
            Column("gravy", "GRAVY (Kyte–Doolittle)", dtype="float", precision=4),
            Column("aromaticity", "Aromaticity (F+W+Y)", dtype="float", precision=4),
            Column("translation_table", "Translation table"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "n_cds": len(cdss),
                "n_residues": residues,
                "ignored_symbols": ignored,
                # raw values: the export layer formats by Column.precision, JSON keeps them exact
                "gravy": gravy(counts),
                "aromaticity": aromaticity(counts),
                "translation_table": table_label,
            },
        ),
        scope=Scope(
            label=(
                "CDS only, pseudogenes excluded, both IR copies included; "
                "stop codons are not residues"
            )
        ),
    )

    table = Table(
        id="amino_acid_composition",
        title="Amino-acid composition",
        columns=(
            Column("residue", "Residue"),
            Column("count", "Count", dtype="int"),
            Column("percent", "Percent", dtype="float", unit="%", precision=2),
            Column("hydropathy", "Kyte–Doolittle hydropathy", dtype="float", precision=1),
        ),
        rows=tuple(
            {
                "residue": residue,
                "count": counts[residue],
                "percent": counts[residue] / residues * 100 if residues else 0.0,
                "hydropathy": KYTE_DOOLITTLE[residue],
            }
            for residue in RESIDUES
        ),
        scope=Scope(label="CDS only, pseudogenes excluded; percents over the residues counted"),
        key="residue",
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"amino_acid": summary, "amino_acid_composition": table},
        notes=(
            f"Translation table: {table_label} (decision D2-03). Table 1 and the plastid table 11 "
            f"share every amino-acid assignment (they differ only in initiation codons), and the "
            f"initiator codon is read as its own amino acid — a CDS starting GTG/ACG/ATC "
            f"contributes V/T/I where its published protein begins with M.",
            f"GRAVY uses the Kyte–Doolittle scale; {ignored} stop/ambiguous symbol(s) excluded.",
        ),
        citations=("biopython",),
    )
