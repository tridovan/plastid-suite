"""Nucleotide diversity π — windows, regions and genes (PROJECT.md §13.2, ROADMAP M2.9).

π is measured on a MAFFT alignment of the input genomes; row 0 is the reference and its
annotation supplies the regions and genes. Gapped columns are **dropped, never scored**, and the
number of usable columns is published next to every π — a value from a mostly-empty window must
not pass for a well-supported one.

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
from plastid.core.ir import detect_ir
from plastid.core.model import GENE_KINDS, Genome
from plastid.core.normalize import NormalizedFrame, normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "diversity"
MODULE_VERSION = "0.1.0"

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

    frame = normalize(reference, detect_ir(reference))
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
                "region": _region_of(touched[len(touched) // 2], frame) if touched else "",
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

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "diversity_summary": summary,
            "diversity_windows": windows_table,
            "diversity_regions": regions,
            "diversity_genes": genes,
        },
        notes=(
            f"π over {alignment.n_rows} aligned genomes; {usable} of {alignment.length} "
            f"alignment columns were usable.",
            "π = mean over sites of (n/(n−1))·(1 − Σp²); a gapped column is dropped, not scored.",
        ),
        citations=("mafft",),
    )


def _region_of(position: int, frame: NormalizedFrame) -> str:
    for region in frame.regions:
        if region.start <= position < region.end:
            return region.name
    return ""


def _region_rows(
    alignment: Alignment, frame: NormalizedFrame, positions: tuple[int | None, ...]
) -> list[dict[str, Any]]:
    """π per normalized region, with each region's columns taken from the reference positions."""
    rows: list[dict[str, Any]] = []
    for region in frame.regions:
        columns = tuple(
            index
            for index, position in enumerate(positions)
            if position is not None and region.start <= position < region.end
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
