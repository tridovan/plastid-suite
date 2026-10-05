"""The two IR junctions: which genes sit on each boundary.

One panel per junction — LSC/IRa and SSC/IRb — showing a window of ±3 kb around the boundary with
every gene whose span falls in it, drawn on its own strand. This is what the IR boundary question
is read off: a gene that straddles the junction, the two copies of a duplicated locus, and the
exact base where LSC becomes IRa.

Everything is drawn in the **normalized** frame (LSC → IRa → SSC → IRb, the frame ``structure``
exports), so the genes — which ``genome_stats`` exports in the GenBank frame — are converted with
the rotation offset the summary reports. The figure reads no sequence: no module stores one in the
run directory, so a junction diagram drawn from a run directory alone can only be made of spans.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.figures.frames import Frame
from plastid.result.result import Result

_TYPED_KINDS = ("CDS", "tRNA", "rRNA")
_WINDOW_BP = 3000  # how much of each side of a boundary the panel shows
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
            start_1b=frame.to_normalized(row["start_1b"]),
            end_1b=frame.to_normalized(row["end_1b"]),
            strand=row["strand"] if row["strand"] in _LANES else "+",
        )
        for row in stats.tables["genes"].rows
        if row["kind"] in _TYPED_KINDS
    ]
    regions = {row["region"]: row for row in structure.tables["regions"].rows}
    junctions: list[Junction] = []
    for left, right in (("LSC", "IRa"), ("SSC", "IRb")):
        if left not in regions or right not in regions:
            raise ValueError(f"ir_junction: the region table has no {left!r}/{right!r} pair")
        boundary = regions[right]["start_1b"]
        start = max(1, boundary - _WINDOW_BP)
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
    return (
        f"IR junctions — {data.organism} ({data.accession}) · "
        f"±{_WINDOW_BP} bp shown · {data.n_genes} gene spans in the two windows"
    )


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
                    head_length=min(length * 0.4, _WINDOW_BP * 0.12),
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
        figsize=(7.4, 4.6),
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
