"""One page per genome — the sheet a supplementary file carries, PDF only (P5-04).

A single portrait page: the accession and its summary numbers, the four regions drawn to scale,
and the GC-skew profile. Every number on it is ``Column.format`` of the cell it comes from, so the
page and the tables cannot disagree — which is the whole reason the supplement is generated rather
than laid out.

It is PDF only on purpose: a page is a document, and an SVG or an EPS of one is a file nobody
asked for. A run that asks for the page gets the PDF, and the figure pass says so.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result

_REGIONS = "regions"
_SUMMARY = "structure_summary"
_SKEW = "gc_skew_windows"
_STATS = "genome_stats"
_PAGE_FORMATS = ("pdf",)


@dataclass(frozen=True, slots=True)
class RegionSpan:
    region: str
    start_1b: int  # normal_ized frame, 1-based inclusive
    end_1b: int
    length: str  # the cell's own rendering
    gc: str


@dataclass(frozen=True, slots=True)
class SkewPoint:
    position: float  # window midpoint, normalized frame
    gc_skew: float | None


@dataclass(frozen=True, slots=True)
class GenomePageSeries:
    """Everything the page prints, read from the two tables — the D5-03 audit trail."""

    accession: str
    organism: str
    length_bp: int
    length: str
    gc: str
    genes: str
    cds: str
    trna: str
    rrna: str
    ir_length: str
    ir_identity: str
    ir_source: str
    regions: tuple[RegionSpan, ...]
    skew: tuple[SkewPoint, ...]

    def stat_lines(self) -> tuple[str, ...]:
        """The summary block, one line per quantity — every value a rendered cell."""
        return (
            f"Length {self.length} bp · GC {self.gc} %",
            f"Genes {self.genes} = {self.cds} CDS + {self.trna} tRNA + {self.rrna} rRNA",
            f"Inverted repeat {self.ir_length} bp · copies {self.ir_identity} % identical "
            f"({self.ir_source})",
        )


def series(results: Mapping[str, Result]) -> GenomePageSeries:
    """Read the genome-statistics and structure tables into the page's plain data."""
    stats = results.get("genome_stats")
    structure = results.get("structure")
    if stats is None or structure is None:
        raise ValueError(
            f"genome_page needs the 'genome_stats' and 'structure' results; have {sorted(results)}"
        )
    for name in (_STATS,):
        if name not in stats.tables:
            raise ValueError(f"genome_page needs table genome_stats.{name}")
    for name in (_REGIONS, _SUMMARY, _SKEW):
        if name not in structure.tables:
            raise ValueError(f"genome_page needs table structure.{name}")

    summary = stats.tables[_STATS]
    region_table = structure.tables[_REGIONS]
    skew_table = structure.tables[_SKEW]

    def cell(table: Any, row_key: str, column: str) -> str:
        return table.format(f"{row_key}.{column}")

    return GenomePageSeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        length_bp=summary.value("length_bp"),
        length=summary.format("length_bp"),
        gc=summary.format("gc_pct"),
        genes=summary.format("n_genes_total"),
        cds=summary.format("n_cds"),
        trna=summary.format("n_trna"),
        rrna=summary.format("n_rrna"),
        ir_length=summary.format("ir_length_bp"),
        ir_identity=summary.format("ir_identity_pct"),
        ir_source=summary.format("ir_source"),
        regions=tuple(
            RegionSpan(
                region=row["region"],
                start_1b=row["start_1b"],
                end_1b=row["end_1b"],
                length=cell(region_table, row["region"], "length_bp"),
                gc=cell(region_table, row["region"], "gc_pct"),
            )
            for row in region_table.rows
        ),
        skew=tuple(
            SkewPoint(
                position=(row["window_start_1b"] + row["window_end_1b"]) / 2,
                gc_skew=row["gc_skew"],
            )
            for row in skew_table.rows
        ),
    )


def caption(data: GenomePageSeries) -> str:
    return f"{data.organism} ({data.accession}) — {data.length} bp, GC {data.gc} %"


def _build(data: GenomePageSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        figure.suptitle(
            title, x=0.08, y=0.95, ha="left", fontsize=THEME.font_size + 3, color=THEME.line
        )
        bar = figure.add_axes((0.10, 0.66, 0.82, 0.09))
        for span in data.regions:
            bar.barh(
                0,
                span.end_1b - span.start_1b + 1,
                left=span.start_1b,
                height=0.6,
                color=THEME.region_fill(span.region),
                edgecolor=THEME.muted,
                linewidth=0.4,
            )
            bar.text(
                (span.start_1b + span.end_1b) / 2,
                0,
                f"{span.region}\n{span.length} bp\nGC {span.gc} %",
                ha="center",
                va="center",
                fontsize=THEME.font_size - 1,
                color=THEME.line,
            )
        bar.set_xlim(1, data.length_bp)
        bar.set_ylim(-0.6, 0.6)
        bar.set_yticks([])
        for spine in ("top", "right", "left"):
            bar.spines[spine].set_visible(False)
        bar.tick_params(axis="x", labelsize=THEME.font_size - 1)

        skew = figure.add_axes((0.10, 0.14, 0.82, 0.44))
        xs = [point.position for point in data.skew]
        ys = [float("nan") if point.gc_skew is None else point.gc_skew for point in data.skew]
        skew.plot(xs, ys, color=THEME.line, linewidth=THEME.line_width)
        skew.axhline(0.0, color=THEME.grid, linewidth=0.8)
        skew.set_xlim(1, data.length_bp)
        skew.set_xlabel("normalized position (bp)")
        skew.set_ylabel("GC skew (G−C)/(G+C)")
        skew.grid(axis="y", zorder=1)
        skew.set_axisbelow(True)

        # Above the drawings, not beside them: a page whose summary sits on the plot is a page
        # nobody can read — the first render had exactly that (§8.39).
        for index, line in enumerate(data.stat_lines()):
            figure.text(
                0.10, 0.89 - index * 0.032, line, fontsize=THEME.font_size, color=THEME.line
            )

    return build


def build(results: Mapping[str, Result]) -> Figure:
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="genome_page",
        title=text,
        source_tables=(
            "genome_stats.genome_stats",
            "structure.regions",
            "structure.gc_skew_windows",
        ),
        build=_build(data, text),
        figsize=(8.27, 11.69),  # A4 portrait
        formats=_PAGE_FORMATS,
    )


SPEC = FigureSpec(
    id="genome_page",
    title="One-page genome sheet (PDF)",
    requires=(
        ("genome_stats", _STATS),
        ("structure", _REGIONS),
        ("structure", _SUMMARY),
        ("structure", _SKEW),
    ),
    build=build,
    formats=_PAGE_FORMATS,
)
