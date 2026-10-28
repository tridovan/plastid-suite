"""The circular plastome map: regions, genes by strand, and the GC-skew ring.

Three rings, read from three tables — ``structure.regions`` (the region spans and their GC),
``genome_stats.genes`` (every gene with its strand) and ``structure.gc_skew_windows`` (the skew
ring). Regions are exported in the normalized frame and genes in the GenBank frame, so the regions
are converted with the rotation offset the summary reports (:mod:`plastid.result.figures.frames`).

Only the typed gene rows are drawn (``CDS``/``tRNA``/``rRNA``): ``genome_stats.genes`` also carries
the ``gene`` feature rows, and drawing both would show a plastome of 265 genes where the contract
says 132 = 86 + 38 + 8. How many rows were dropped is stated in the caption, never silently.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.figures.cells import gene_segments
from plastid.result.figures.frames import Frame
from plastid.result.result import Result

_TYPED_KINDS = ("CDS", "tRNA", "rRNA")

# Radii of the three rings, in axes fractions. One place, so the rings cannot drift apart, and
# each label sits in the gap between two rings rather than on top of one.
_REGION_RING = (1.0, 1.045)
_REGION_LABEL = 1.11
_GENE_RING = {"+": (0.90, 0.96), "-": (0.82, 0.88)}
_GENE_LABEL = {"+": 0.98, "-": 0.79}
_SKEW_RING = (0.70, 0.78)

# A gene label is drawn only when its arc is wide enough to hold it, and only when the previous
# label on the same strand is far enough around the circle — 132 labels on one circle is a smear,
# not a map. Both are layout rules; neither changes a number.
_LABEL_MIN_DEGREES = 4.0
_LABEL_MIN_GAP_DEGREES = 8.0

_GENE_COLOURS = {"CDS": "#2f5d8a", "tRNA": "#b4472e", "rRNA": "#3f7d52"}


@dataclass(frozen=True, slots=True)
class Gene:
    name: str
    kind: str
    start_1b: int  # GenBank frame, 1-based inclusive, outer span of the feature
    end_1b: int
    strand: str


@dataclass(frozen=True, slots=True)
class RegionArc:
    region: str
    start_1b: int  # GenBank frame after conversion
    end_1b: int
    length_bp: int


@dataclass(frozen=True, slots=True)
class SkewWindow:
    region: str
    start_1b: int  # normalized frame
    end_1b: int
    gc_skew: float | None


@dataclass(frozen=True, slots=True)
class CircularMapSeries:
    """Everything the map draws — every field a table cell or read off one."""

    accession: str
    organism: str
    length_bp: int
    gc_pct: str  # the summary cell's own rendering
    n_genes: str
    frame: Frame
    regions: tuple[RegionArc, ...]
    genes: tuple[Gene, ...]
    skew: tuple[SkewWindow, ...]
    n_gene_features_skipped: int

    def theta(self, position_1b: int) -> float:
        """A 1-based position as an angle in radians, 0 at the origin of the frame."""
        return 2 * math.pi * (position_1b - 1) / self.length_bp


def series(results: Mapping[str, Result]) -> CircularMapSeries:
    """Read the three tables into the plain data the map draws."""
    stats = results.get("genome_stats")
    structure = results.get("structure")
    if stats is None or structure is None:
        raise ValueError(
            f"circular_map needs the 'genome_stats' and 'structure' results; have {sorted(results)}"
        )
    for name in ("genome_stats", "genes"):
        if name not in stats.tables:
            raise ValueError(f"circular_map needs table genome_stats.{name}")
    for name in ("regions", "structure_summary", "gc_skew_windows"):
        if name not in structure.tables:
            raise ValueError(f"circular_map needs table structure.{name}")

    summary = stats.tables["genome_stats"]
    frame = Frame(
        length_bp=summary.value("length_bp"),
        rotation_offset=structure.tables["structure_summary"].value("rotation_offset"),
    )
    typed_rows = [row for row in stats.tables["genes"].rows if row["kind"] in _TYPED_KINDS]
    typed = [
        Gene(
            name=row["name"],
            kind=row["kind"],
            start_1b=start,
            end_1b=end,
            strand=row["strand"],
        )
        for row in typed_rows
        # One arc per part: a feature annotated across the record origin has no single linear
        # interval, and its outer bounds would draw one arc the long way round (M11.7).
        for start, end in gene_segments(row)
    ]
    regions: list[RegionArc] = []
    for row in structure.tables["regions"].rows:
        for start, end in frame.span_to_genbank(row["start_1b"], row["end_1b"]):
            regions.append(
                RegionArc(
                    region=row["region"], start_1b=start, end_1b=end, length_bp=row["length_bp"]
                )
            )
    return CircularMapSeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        length_bp=summary.value("length_bp"),
        gc_pct=summary.format("gc_pct"),
        n_genes=summary.format("n_genes_total"),
        frame=frame,
        regions=tuple(regions),
        genes=tuple(typed),
        skew=tuple(
            SkewWindow(
                region=row["region"],
                start_1b=row["window_start_1b"],
                end_1b=row["window_end_1b"],
                gc_skew=row["gc_skew"],
            )
            for row in structure.tables["gc_skew_windows"].rows
        ),
        n_gene_features_skipped=len(stats.tables["genes"].rows) - len(typed_rows),
    )


def caption(data: CircularMapSeries) -> str:
    """The line the figure carries — no number of its own, every one a cell's rendering."""
    text = (
        f"{data.organism} ({data.accession}) — {data.length_bp} bp · GC {data.gc_pct} % · "
        f"{data.n_genes} genes"
    )
    if data.n_gene_features_skipped:
        text += (
            f" · {data.n_gene_features_skipped} 'gene' feature row(s) not drawn (see genes table)"
        )
    return text


def _ring(
    axes: Any,
    data: CircularMapSeries,
    start_1b: int,
    end_1b: int,
    radii: tuple[float, float],
    colour: str,
    *,
    zorder: int = 2,
) -> None:
    """One arc as a polar bar — split in two when it crosses the origin of the frame."""
    theta1 = data.theta(start_1b)
    theta2 = data.theta(end_1b + 1)
    if theta2 <= theta1:
        theta2 += 2 * math.pi
    height = radii[1] - radii[0]
    axes.bar(
        (theta1 + theta2) / 2,
        height,
        width=theta2 - theta1,
        bottom=radii[0],
        color=colour,
        linewidth=0,
        zorder=zorder,
    )


def _build(data: CircularMapSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        axes = figure.add_subplot(projection="polar")
        axes.set_theta_offset(math.pi / 2)  # position 1 at the top
        axes.set_theta_direction(-1)  # and the genome running clockwise, as plastome maps do
        axes.set_axis_off()
        axes.set_ylim(0, 1.28)
        axes.set_title(title, loc="left", fontsize=THEME.font_size, color=THEME.line, pad=10)

        for region in data.regions:
            _ring(
                axes,
                data,
                region.start_1b,
                region.end_1b,
                _REGION_RING,
                THEME.region_fill(region.region),
                zorder=1,
            )
            middle = data.theta((region.start_1b + region.end_1b) / 2)
            axes.text(
                middle,
                _REGION_LABEL,
                region.region,
                ha="center",
                va="center",
                fontsize=THEME.font_size - 1,
                color=THEME.muted,
            )

        for gene in data.genes:
            radii = _GENE_RING.get(gene.strand, _GENE_RING["+"])
            _ring(axes, data, gene.start_1b, gene.end_1b, radii, _GENE_COLOURS[gene.kind], zorder=3)

        labelled: dict[str, float] = {}
        gap = math.radians(_LABEL_MIN_GAP_DEGREES)
        for gene in sorted(data.genes, key=lambda item: item.start_1b):
            span = 360.0 * (gene.end_1b + 1 - gene.start_1b) / data.length_bp
            if span < _LABEL_MIN_DEGREES:
                continue
            middle = data.theta((gene.start_1b + gene.end_1b) / 2)
            previous = labelled.get(gene.strand)
            if previous is not None and abs(middle - previous) < gap:
                continue
            labelled[gene.strand] = middle
            axes.text(
                middle,
                _GENE_LABEL.get(gene.strand, _GENE_LABEL["+"]),
                gene.name,
                ha="center",
                va="center",
                fontsize=THEME.font_size - 2,
                color=THEME.line,
            )

        # The ring is a shape, not a scale: bars are sized against the largest skew this genome
        # actually shows, so a genome whose skew is uniformly small still reads as a ring.
        peak = max((abs(w.gc_skew) for w in data.skew if w.gc_skew is not None), default=0.0)
        for window in data.skew:
            if window.gc_skew is None or peak == 0.0:
                continue  # an undefined skew is a hole in the ring, never a zero
            middle = data.theta((window.start_1b + window.end_1b) / 2)
            # A bar's width is the window's whole arc: half of it would leave the ring dotted.
            width = 2 * math.pi * (window.end_1b - window.start_1b + 1) / data.length_bp
            base = (_SKEW_RING[0] + _SKEW_RING[1]) / 2
            reach = (_SKEW_RING[1] - _SKEW_RING[0]) / 2
            axes.bar(
                middle,
                (window.gc_skew / peak) * reach,
                width=width,
                bottom=base,
                color=THEME.accent if window.gc_skew >= 0 else THEME.line,
                linewidth=0,
                zorder=2,
            )
        axes.plot([], [], color=THEME.accent, linewidth=3, label="GC skew > 0")
        axes.plot([], [], color=THEME.line, linewidth=3, label="GC skew < 0")
        for kind, colour in _GENE_COLOURS.items():
            axes.plot([], [], color=colour, linewidth=3, label=kind)
        axes.legend(
            loc="lower center",
            bbox_to_anchor=(0.5, -0.06),
            ncol=len(_GENE_COLOURS) + 2,
            frameon=False,
            fontsize=THEME.font_size - 1,
        )

    return build


def build(results: Mapping[str, Result]) -> Figure:
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="circular_map",
        title=text,
        source_tables=(
            "genome_stats.genome_stats",
            "genome_stats.genes",
            "structure.regions",
            "structure.structure_summary",
            "structure.gc_skew_windows",
        ),
        build=_build(data, text),
        figsize=(5.6, 6.0),
    )


SPEC = FigureSpec(
    id="circular_map",
    title="Circular plastome map",
    requires=(
        ("genome_stats", "genome_stats"),
        ("genome_stats", "genes"),
        ("structure", "regions"),
        ("structure", "structure_summary"),
        ("structure", "gc_skew_windows"),
    ),
    build=build,
)
