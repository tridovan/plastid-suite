"""Nucleotide diversity π — windows, regions, genes, introns and intergenic spacers.

π is measured on a MAFFT alignment of the input genomes; row 0 is the reference and its
annotation supplies the regions, genes and spans. Gapped columns are **dropped, never scored**,
and the number of usable columns is published next to every π — a value from a mostly-empty
window must not pass for a well-supported one.

The intron and spacer spans come from :mod:`plastid.core.spans`, the same derivation the
``introns`` module publishes, so the two tables join row for row and cannot disagree about where
a boundary is (D-g). The spacer set is the single-copy regions only: the IR's two copies are one
sequence annotated twice, and — measured on the golden 10/10/2026 — their *spacers* do not always
mirror exactly, so a spacer cannot be deduplicated the way a gene can.

Window and step (decision D2-02, the mVISTA-style convention) come from params and are printed
with the result: a π curve without its window means nothing.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from typing import Any

from plastid.core.align import (
    Alignment,
    columns_for_span,
    nucleotide_diversity,
    reference_positions,
    window_diversity,
)
from plastid.core.model import GENE_KINDS, Genome
from plastid.core.normalize import NormalizedFrame, normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.core.spans import IntronSpan, SpacerSet, introns_of, spacers_of
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "diversity"
MODULE_VERSION = "0.2.0"  # 0.2.0: π per intron and per intergenic spacer (P2B-01)

DEFAULT_WINDOW_BP = 600  # D2-02
DEFAULT_STEP_BP = 200


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Nucleotide diversity (π)",
    requires_ir=True,
    arity=2,
    citations=("mafft",),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    window = int(params.get("window", DEFAULT_WINDOW_BP))
    step = int(params.get("step", DEFAULT_STEP_BP))
    reference = genomes[0]
    alignment = align(
        [(genome.accession or genome.name, genome.sequence) for genome in genomes], ctx
    )

    whole_pi, usable = nucleotide_diversity(alignment, range(alignment.length))
    windows = window_diversity(alignment, window=window, step=step)

    ir = ctx.require_ir()  # seam 4: the gate runs first, and it is the run's IR, not a re-detect
    frame = normalize(reference, ir)
    introns = introns_of(reference, ir)
    spacers = spacers_of(reference, ir)
    positions = reference_positions(alignment)
    window_rows = []
    for entry in windows:
        touched = [
            position
            for position in positions[entry.start_1b - 1 : entry.end_1b]
            if position is not None
        ]
        window_rows.append(
            {
                "window_start_1b": entry.start_1b,
                "window_end_1b": entry.end_1b,
                "usable_columns": entry.usable_columns,
                "pi": entry.pi,
                "region": _region_name(frame, touched[len(touched) // 2]) if touched else "",
            }
        )

    summary = Table(
        id="diversity_summary",
        title="Nucleotide diversity summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes aligned", dtype="int"),
            Column("alignment_length", "Alignment length", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π (whole alignment)", dtype="float", precision=6),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("n_introns", "Introns (reference)", dtype="int"),
            Column("n_spacers", "Intergenic spacers (reference)", dtype="int"),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference.accession or "",
                "organism": reference.name,
                "n_genomes": alignment.n_rows,
                "alignment_length": alignment.length,
                "usable_columns": usable,
                "pi": whole_pi,
                "window_bp": window,
                "step_bp": step,
                "n_introns": len(introns.introns),
                "n_spacers": len(spacers.igs),
                "tool": "MAFFT",
            },
        ),
        scope=Scope(
            label=(
                "whole aligned plastome, gaps dropped column-wise; positions are alignment "
                "columns, the region column follows the reference genome"
            ),
            coord_system="normalized",
        ),
    )

    windows_table = Table(
        id="diversity_windows",
        title="Nucleotide diversity windows",
        columns=(
            Column("window_start_1b", "Start (alignment column)", dtype="int"),
            Column("window_end_1b", "End (alignment column)", dtype="int"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
            Column("region", "Region (reference)"),
        ),
        rows=tuple(window_rows),
        scope=Scope(
            label=f"{window} bp window / {step} bp step over alignment columns; gaps dropped",
            coord_system="normalized",
        ),
    )

    region_rows = _region_rows(alignment, frame, positions)
    regions = Table(
        id="diversity_regions",
        title="Nucleotide diversity by region",
        columns=(
            Column("region", "Region"),
            Column("columns", "Alignment columns", dtype="int"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=tuple(region_rows),
        scope=Scope(
            label=(
                "reference genome's normalized regions mapped onto alignment columns; gaps dropped"
            ),
            coord_system="normalized",
        ),
        key="region",
    )

    genes = Table(
        id="diversity_genes",
        title="Nucleotide diversity by gene (reference annotation)",
        columns=(
            Column("gene", "Gene"),
            Column("kind", "Kind"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=tuple(_gene_rows(reference, alignment)),
        scope=Scope(
            label="reference genome's features mapped onto alignment columns; gaps dropped",
            coord_system="normalized",
        ),
        key="gene",
    )

    intron_rows = [
        {
            "intron_id": span.intron_id,
            "gene": span.gene,
            "region": span.region,
            "length_bp": span.length,
            **_span_diversity(alignment, span.start, span.end),
        }
        for span in introns.introns
    ]
    intron_table = Table(
        id="diversity_introns",
        title="Nucleotide diversity by intron (reference annotation)",
        columns=(
            Column("intron_id", "Intron (gene.copy.n)"),
            Column("gene", "Gene"),
            Column("region", "Region"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=tuple(intron_rows),
        scope=Scope(
            label=(
                "the introns the `introns` module publishes — both IR copies, whose pair "
                "measures one locus — mapped onto alignment columns; gaps dropped"
            ),
        ),
        key="intron_id",
    )

    spacer_rows = [
        {
            "igs_id": span.igs_id,
            "left": span.left,
            "right": span.right,
            "region": span.region,
            "length_bp": span.length,
            **_span_diversity(alignment, span.start, span.end),
        }
        for span in spacers.igs
    ]
    spacer_table = Table(
        id="diversity_igs",
        title="Nucleotide diversity by intergenic spacer (reference annotation)",
        columns=(
            Column("igs_id", "Spacer (left–right)"),
            Column("left", "Left gene"),
            Column("right", "Right gene"),
            Column("region", "Region"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
        ),
        rows=tuple(spacer_rows),
        scope=Scope(
            label=(
                "gaps between the exons of two neighbouring called genes inside LSC or SSC; the "
                "IR is excluded (its two copies are one sequence annotated twice) — see the "
                "notes for the gaps that are not spacers; gaps dropped"
            ),
        ),
        key="igs_id",
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "diversity_summary": summary,
            "diversity_windows": windows_table,
            "diversity_regions": regions,
            "diversity_genes": genes,
            "diversity_introns": intron_table,
            "diversity_igs": spacer_table,
        },
        notes=(
            f"π over {alignment.n_rows} aligned genomes; {usable} of {alignment.length} "
            f"alignment columns were usable.",
            "π = mean over sites of (n/(n−1))·(1 − Σp²); a gapped column is dropped, not scored.",
            *_span_notes(introns.introns, spacers),
        ),
        citations=("mafft",),
    )


def _span_diversity(alignment: Alignment, start: int, end: int) -> dict[str, Any]:
    """π over one reference span, as the two columns every span table publishes."""
    pi, usable = nucleotide_diversity(alignment, columns_for_span(alignment, 0, (start, end)))
    return {"usable_columns": usable, "pi": pi}


def _span_notes(introns: Sequence[IntronSpan], spacers: SpacerSet) -> tuple[str, ...]:
    """What the span tables measured, and every gap that was looked at and is not one.

    Named rather than counted away: a reader has to be able to see which gap was left out and
    why, or two runs that disagree about the scope cannot be told apart from two that disagree
    about the genome (the pitfall D-g exists for).
    """
    notes = [
        f"{len(introns)} introns and {len(spacers.igs)} intergenic spacers on the reference; "
        f"every π is computed over the reference's own span, mapped onto alignment columns."
    ]
    if spacers.in_ir:
        notes.append(
            f"{len(spacers.in_ir)} spacer(s) lie inside the IR and are not reported: the two "
            f"copies are one sequence annotated twice."
        )
    if spacers.crossing:
        notes.append(
            "gaps crossing a plastid-region boundary are not spacers: "
            + ", ".join(f"{gap.left}–{gap.right} ({gap.length} bp)" for gap in spacers.crossing)
        )
    if spacers.no_spacer:
        notes.append(
            "gene pairs leaving no gap to measure: "
            + ", ".join(f"{gap.left}·{gap.right} ({gap.length} bp)" for gap in spacers.no_spacer)
        )
    return tuple(notes)


def _region_name(frame: NormalizedFrame, position: int) -> str:
    """The reference's region at one of its **own** positions — the frames are not the same.

    Measured 10/10/2026: this module compared GenBank-frame positions against normalized-frame
    bounds until P2B-01. The golden's rotation offset is 0, so the two agreed by accident; the
    first record whose origin is not the LSC's first base would have been labelled with another
    genome's regions. ``region_at`` does the conversion.
    """
    region = frame.region_at(position)
    return region.name if region else ""


def _region_rows(
    alignment: Alignment, frame: NormalizedFrame, positions: tuple[int | None, ...]
) -> list[dict[str, Any]]:
    """π per normalized region, with each region's columns taken from the reference positions."""
    rows: list[dict[str, Any]] = []
    for region in frame.regions:
        columns = tuple(
            index
            for index, position in enumerate(positions)
            if position is not None and _region_name(frame, position) == region.name
        )
        pi, usable = nucleotide_diversity(alignment, columns)
        rows.append(
            {"region": region.name, "columns": len(columns), "usable_columns": usable, "pi": pi}
        )
    return rows


def _gene_rows(reference: Genome, alignment: Alignment) -> list[dict[str, Any]]:
    """π per annotated feature, built from its **exons** — an intron is not part of its gene."""
    copies: Counter[str] = Counter()
    rows: list[dict[str, Any]] = []
    for gene in reference.genes:
        if gene.kind not in GENE_KINDS:
            continue
        copies[gene.name] += 1
        columns: list[int] = []
        for part in gene.parts:
            columns.extend(columns_for_span(alignment, 0, (part.start, part.end)))
        pi, usable = nucleotide_diversity(alignment, columns)
        rows.append(
            {
                "gene": f"{gene.name}.{copies[gene.name]}",
                "kind": gene.kind,
                "usable_columns": usable,
                "pi": pi,
            }
        )
    return rows
