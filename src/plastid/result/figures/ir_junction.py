"""The two IR junctions: which genes sit on each boundary.

One panel per junction — LSC/IRa and SSC/IRb — showing **3 kb of the single-copy region before the
boundary and the whole IR copy after it**, with every gene whose span falls in it, drawn on its own
strand. The whole IR is shown because the question the figure answers is not only where the
junction is: it is which genes straddle it and which loci are duplicated across the copy, and a
3 kb slice would cut the duplicated loci that sit further inside the IR. The panel width therefore
varies with the IR's length (15–47 kb over the 80-record panel), so the figure's physical width
grows with the window rather than squeezing a 46 kb panel into the same 7.4 inches.

Everything is drawn in the **normalized** frame (LSC → IRa → SSC → IRb, the frame ``structure``
exports), so the genes — which ``genome_stats`` exports in the GenBank frame — are converted with
the rotation offset the summary reports. A feature's parts come from the genes table's
``segments`` cell: a gene annotated across the record origin is drawn as its parts, not as the
whole-genome arc its outer bounds would give, and a part that crosses the normalized origin is
split again. The figure reads no sequence: no module stores one in the run directory, so a junction
diagram drawn from a run directory alone can only be made of spans.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.figures.cells import gene_segments
from plastid.result.figures.frames import Frame
from plastid.result.result import Result

_TYPED_KINDS = ("CDS", "tRNA", "rRNA")
# How much of the single-copy region the panel shows before the boundary. The other side is the
# whole IR copy, so the window is `_SINGLE_COPY_CONTEXT_BP + IR length`, not `2 * _SINGLE_COPY_...`.
_SINGLE_COPY_CONTEXT_BP = 3000
_BASE_WIDTH_IN = 7.4  # the figure width while the window fits in _BASE_WINDOW_BP
_BASE_WINDOW_BP = 30_000
_WIDTH_PER_10KB_IN = 1.0  # one extra inch per 10 kb of window above the base
_MAX_WIDTH_IN = 10.0
_PANEL_HEIGHT_IN = 4.6
_LANES = {"+": 1.0, "-": 0.0}
_GENE_COLOURS = {"CDS": "#2f5d8a", "tRNA": "#b4472e", "rRNA": "#3f7d52"}


@dataclass(frozen=True, slots=True)
class PlacedGene:
    name: str
    kind: str
    start_1b: int  # normalized frame, clipped to the shown window
    end_1b: int
    strand: str


@dataclass(frozen=True, slots=True)
class Junction:
    name: str  # "LSC / IRa"
    left_region: str
    right_region: str
    boundary_1b: int  # normalized: the first base of the right-hand region
    window_start_1b: int
    window_end_1b: int
    genes: tuple[PlacedGene, ...]


@dataclass(frozen=True, slots=True)
class JunctionSeries:
    accession: str
    organism: str
    junctions: tuple[Junction, ...]

    @property
    def n_genes(self) -> int:
        return sum(len(junction.genes) for junction in self.junctions)


def series(results: Mapping[str, Result]) -> JunctionSeries:
    """Read the region table and the gene table into the two junction panels."""
    stats = results.get("genome_stats")
    structure = results.get("structure")
    if stats is None or structure is None:
        raise ValueError(
            f"ir_junction needs the 'genome_stats' and 'structure' results; have {sorted(results)}"
        )
    for name in ("genome_stats", "genes"):
        if name not in stats.tables:
            raise ValueError(f"ir_junction needs table genome_stats.{name}")
    for name in ("regions", "structure_summary"):
        if name not in structure.tables:
            raise ValueError(f"ir_junction needs table structure.{name}")

    summary = stats.tables["genome_stats"]
    frame = Frame(
        length_bp=summary.value("length_bp"),
        rotation_offset=structure.tables["structure_summary"].value("rotation_offset"),
    )
    genes = [
        PlacedGene(
            name=row["name"],
            kind=row["kind"],
            start_1b=start,
            end_1b=end,
            strand=row["strand"] if row["strand"] in _LANES else "+",
        )
        for row in stats.tables["genes"].rows
        if row["kind"] in _TYPED_KINDS
        # A feature annotated across the record origin has no single linear interval; its outer
        # bounds are the whole genome, and drawing those would fill the panel. Each part is
        # placed where the annotation says it is, then split again if it crosses the
        # *normalized* origin (M11.7).
        for segment_start, segment_end in gene_segments(row)
        for start, end in frame.span_to_normalized(segment_start, segment_end)
    ]
    regions = {row["region"]: row for row in structure.tables["regions"].rows}
    junctions: list[Junction] = []
    for left, right in (("LSC", "IRa"), ("SSC", "IRb")):
        if left not in regions or right not in regions:
            raise ValueError(f"ir_junction: the region table has no {left!r}/{right!r} pair")
        boundary = regions[right]["start_1b"]
        # 3 kb of the single-copy region, never reaching back past that region's own start; the
        # right side is the whole IR copy, which is what makes the duplicated loci visible.
        start = max(regions[left]["start_1b"], boundary - _SINGLE_COPY_CONTEXT_BP)
        end = min(frame.length_bp, regions[right]["end_1b"])
        shown = []
        for gene in genes:
            if gene.end_1b < start or gene.start_1b > end:
                continue
            shown.append(
                PlacedGene(
                    name=gene.name,
                    kind=gene.kind,
                    start_1b=max(gene.start_1b, start),
                    end_1b=min(gene.end_1b, end),
                    strand=gene.strand,
                )
            )
        junctions.append(
            Junction(
                name=f"{left} / {right}",
                left_region=left,
                right_region=right,
                boundary_1b=boundary,
                window_start_1b=start,
                window_end_1b=end,
                genes=tuple(shown),
            )
        )
    return JunctionSeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        junctions=tuple(junctions),
    )


def caption(data: JunctionSeries) -> str:
    """State the rule the panels actually use, and every panel's own width.

    Measured 28/10/2026: the old caption said "±3000 bp" while the code showed the whole IR copy —
    29,214 bp on the golden genome. A caption is a claim; it names the window each panel got.
    """
    widths = " / ".join(
        f"{junction.window_end_1b - junction.window_start_1b + 1:,}" for junction in data.junctions
    )
    return (
        f"IR junctions — {data.organism} ({data.accession}) · "
        f"{_SINGLE_COPY_CONTEXT_BP:,} bp of the single-copy region + the whole IR copy · "
        f"windows {widths} bp · {data.n_genes} gene spans in the two windows"
    )


def _figure_width(window_bp: int) -> float:
    """The physical width for the widest panel: one inch per extra 10 kb above the base (M11.7).

    A 3 kb + whole-IR window is 15–47 kb over the 80-record panel. At one fixed 7.4 in the largest
    record draws 6.3 kb per inch and its gene labels collide; widening the figure keeps the gene
    density roughly constant. The cap keeps a pathological IR from producing an unprintable page.
    """
    extra = max(0, window_bp - _BASE_WINDOW_BP) / 10_000
    return min(_MAX_WIDTH_IN, _BASE_WIDTH_IN + extra * _WIDTH_PER_10KB_IN)


def _build(data: JunctionSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        panels = figure.subplots(len(data.junctions), 1)
        figure.suptitle(title, x=0.02, ha="left", fontsize=THEME.font_size, color=THEME.line)
        for axes, junction in zip(panels, data.junctions, strict=True):
            for colour, from_left in (
                (THEME.region_fill(junction.left_region), True),
                (THEME.region_fill(junction.right_region), False),
            ):
                if from_left:
                    axes.axvspan(
                        junction.window_start_1b,
                        junction.boundary_1b,
                        color=colour,
                        linewidth=0,
                        zorder=0,
                    )
                else:
                    axes.axvspan(
                        junction.boundary_1b,
                        junction.window_end_1b,
                        color=colour,
                        linewidth=0,
                        zorder=0,
                    )
            axes.axvline(junction.boundary_1b, color=THEME.accent, linewidth=1.0, zorder=4)
            axes.text(
                junction.boundary_1b,
                1.36,
                f"{junction.left_region} | {junction.right_region} at {junction.boundary_1b}",
                ha="center",
                va="bottom",
                fontsize=THEME.font_size - 1,
                color=THEME.accent,
            )
            for gene in junction.genes:
                lane = _LANES[gene.strand]
                length = gene.end_1b - gene.start_1b + 1
                axes.arrow(
                    gene.start_1b,
                    lane,
                    length,
                    0,
                    head_width=0.18,
                    head_length=min(length * 0.4, _SINGLE_COPY_CONTEXT_BP * 0.12),
                    length_includes_head=True,
                    linewidth=0,
                    color=_GENE_COLOURS[gene.kind],
                    zorder=3,
                )
                axes.text(
                    (gene.start_1b + gene.end_1b) / 2,
                    lane + 0.26,
                    gene.name,
                    ha="center",
                    va="bottom",
                    rotation=45,
                    fontsize=THEME.font_size - 2,
                    color=THEME.line,
                )
            axes.set_xlim(junction.window_start_1b, junction.window_end_1b)
            axes.set_ylim(-0.4, 1.5)
            axes.set_yticks([_LANES["-"], _LANES["+"]])
            axes.set_yticklabels(["− strand", "+ strand"], fontsize=THEME.font_size - 1)
            axes.set_xlabel("normalized position (bp)")
            axes.grid(axis="x", zorder=1)
            axes.set_axisbelow(True)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    data = series(results)
    text = caption(data)
    widest = max(
        junction.window_end_1b - junction.window_start_1b + 1 for junction in data.junctions
    )
    return export_figure(
        figure_id="ir_junction",
        title=text,
        source_tables=(
            "genome_stats.genome_stats",
            "genome_stats.genes",
            "structure.regions",
            "structure.structure_summary",
        ),
        build=_build(data, text),
        figsize=(_figure_width(widest), _PANEL_HEIGHT_IN),
    )


SPEC = FigureSpec(
    id="ir_junction",
    title="IR junction diagram",
    requires=(
        ("genome_stats", "genome_stats"),
        ("genome_stats", "genes"),
        ("structure", "regions"),
        ("structure", "structure_summary"),
    ),
    build=build,
)
