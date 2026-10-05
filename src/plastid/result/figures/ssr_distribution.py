"""SSR distribution — counts per plastid region, and every SSR's position per repeat unit.

Sources: ``ssr.ssr_by_region`` (stacked counts), ``ssr.ssrs`` (positions) and ``ssr.ssr_summary``
(the caption's total and thresholds). The right panel's x axis is the **GenBank frame** (1-based,
inclusive) because ``ssrs`` carries GenBank coordinates; the left panel is counts, not
coordinates, so the panels never share an axis and the label says which frame the strip is in. An
SSR spans far less than a pixel at genome scale, so it is drawn as a tick at the centre of its
span, not as a drawn interval: the tick marks the position honestly, and start..end stay in the
series (1-based, inclusive). A blank count stays a gap — no width, no number.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result

_BY_REGION = "ssr_by_region"
_SSRS = "ssrs"
_SUMMARY = "ssr_summary"
_SEGMENT_NAMES = ("mono", "di", "tri", "other")
_REGION_COLUMNS = ("region", "total", *_SEGMENT_NAMES)
_SSR_COLUMNS = ("unit", "start_1b", "end_1b")
_SUMMARY_COLUMNS = ("accession", "total", "thresholds")

# Repeat-unit strips in unit-size order. The labels are MISA's own, but this figure must not
# import the external-tool layer to learn their order, so the order is restated here; a label the
# tuple does not know keeps its table order, after the known ones.
_UNIT_ORDER = ("mono", "di", "tri", "tetra", "penta", "hexa")

# Stack colours: mono/di are the theme's line/accent; the two extra hexes keep every pair
# separable under colour-vision deficiency and at least 3:1 against the white surface (measured
# with a palette checker: worst all-pairs dE 11.6 deutan, 15.9 normal vision).
_COLOUR_TRI = "#b8860b"  # darkgoldenrod
_COLOUR_OTHER = "#2f7da8"  # steel blue
_SEGMENT_STYLES: dict[str, tuple[str, str]] = {  # name -> (fill, the ink that reads on it)
    "mono": (THEME.line, THEME.background),
    "di": (THEME.accent, THEME.background),
    "tri": (_COLOUR_TRI, THEME.line),
    "other": (_COLOUR_OTHER, THEME.background),
}
_FALLBACK_STYLE = (THEME.muted, THEME.line)
_POSITION_PAD = 0.02  # axis headroom, a layout constant — not a value read off the data


@dataclass(frozen=True, slots=True)
class RegionCounts:
    """One region row: its total and its unit segments, in the table's column order."""

    region: str
    total: int | None
    segments: tuple[tuple[str, int | None, str], ...]  # (unit column, count, the cell's rendering)


@dataclass(frozen=True, slots=True)
class SsrSpan:
    """One SSR — GenBank frame, 1-based inclusive, exactly the cells' values."""

    unit: str
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class SsrDistributionSeries:
    """Everything the figure draws, read from the tables — the D5-03 audit trail."""

    accession: str
    total: str
    thresholds: str
    regions: tuple[RegionCounts, ...]
    spans: tuple[SsrSpan, ...]

    @property
    def units(self) -> tuple[str, ...]:
        """The strip labels, smallest repeat unit first; unknown labels keep table order, last."""
        seen = list(dict.fromkeys(span.unit for span in self.spans))
        known = [unit for unit in _UNIT_ORDER if unit in seen]
        return tuple(known + [unit for unit in seen if unit not in _UNIT_ORDER])


def series(results: Mapping[str, Result]) -> SsrDistributionSeries:
    """Read the SSR tables into the plain data the figure draws."""
    result = results.get("ssr")
    if result is None:
        raise ValueError(
            f"ssr_distribution needs the 'ssr' module's result; have {sorted(results)}"
        )
    missing = [name for name in (_BY_REGION, _SSRS, _SUMMARY) if name not in result.tables]
    if missing:
        raise ValueError(f"ssr_distribution needs table(s) {missing} from 'ssr'")
    region_table = result.tables[_BY_REGION]
    features = result.tables[_SSRS]
    summary = result.tables[_SUMMARY]
    checks = (
        (region_table, _REGION_COLUMNS),
        (features, _SSR_COLUMNS),
        (summary, _SUMMARY_COLUMNS),
    )
    for table, columns in checks:
        absent = [name for name in columns if name not in {c.name for c in table.columns}]
        if absent:
            raise ValueError(f"ssr_distribution: {table.id} has no column(s) {absent}")

    return SsrDistributionSeries(
        accession=summary.format("accession"),
        total=summary.format("total"),
        thresholds=summary.format("thresholds"),
        regions=tuple(
            RegionCounts(
                region=row["region"],
                total=row["total"],
                segments=tuple(
                    (name, row[name], region_table.format(f"{row['region']}.{name}"))
                    for name in _SEGMENT_NAMES
                ),
            )
            for row in region_table.rows
        ),
        spans=tuple(SsrSpan(row["unit"], row["start_1b"], row["end_1b"]) for row in features.rows),
    )


def caption(data: SsrDistributionSeries) -> str:
    """The line the figure carries — the summary cells' own renderings, verbatim."""
    return f"SSR distribution — {data.accession} · {data.total} SSRs · thresholds {data.thresholds}"


def _style_panel(axes: Any, title: str, xlabel: str) -> None:
    """The axis furniture the two panels share, so they read as one figure."""
    axes.set_xlabel(xlabel)
    axes.set_title(title, loc="left", fontsize=THEME.font_size - 1, color=THEME.muted, pad=4)
    axes.grid(axis="x", color=THEME.grid)
    axes.set_axisbelow(True)


def _regions_panel(axes: Any, data: SsrDistributionSeries) -> None:
    """Stacked counts per region; every segment carries its own cell's rendering."""
    from matplotlib.patches import Patch

    for index, region in enumerate(data.regions):
        offset = 0.0
        for name, count, text in region.segments:
            if count is None:
                continue  # a blank is a gap: no width, no number — never drawn as a zero
            fill, ink = _SEGMENT_STYLES.get(name, _FALLBACK_STYLE)
            axes.barh(
                index,
                count,
                left=offset,
                height=0.62,
                color=fill,
                edgecolor=THEME.background,  # a surface gap between neighbouring fills
                linewidth=0.8,
            )
            if count:  # a zero-width segment has nowhere to print its number
                axes.text(
                    offset + count / 2,
                    index,
                    text,
                    ha="center",
                    va="center",
                    fontsize=THEME.font_size - 3,
                    color=ink,
                )
            offset += count
    handles = [
        Patch(facecolor=_SEGMENT_STYLES.get(name, _FALLBACK_STYLE)[0], label=name)
        for name in _SEGMENT_NAMES
    ]
    axes.legend(handles=handles, loc="center right", frameon=False, fontsize=THEME.font_size - 1)
    axes.set_yticks(range(len(data.regions)), labels=[region.region for region in data.regions])
    axes.set_ylim(len(data.regions) - 0.5, -0.5)  # first row at the top, one axis unit per row
    _style_panel(axes, "counts by region", "SSRs (count)")


def _positions_panel(axes: Any, data: SsrDistributionSeries) -> None:
    """One strip per repeat unit; each SSR is a tick at the centre of its start..end span."""
    units = data.units
    for index, unit in enumerate(units):
        spans = [span for span in data.spans if span.unit == unit]
        fill = _SEGMENT_STYLES.get(unit, _SEGMENT_STYLES["other"])[0]  # the stack's own colour
        # A vertical stroke, not hlines: the interval is sub-pixel and matplotlib draws a
        # near-zero-length stroke as nothing at all (measured) — the tick's width is the
        # linewidth, its height is in data units, and the span's midpoint is where it sits.
        axes.vlines(
            [(span.start + span.end) / 2 for span in spans],
            index - 0.3,
            index + 0.3,
            color=fill,
            linewidth=1.8,
        )
    starts = [span.start for span in data.spans]
    ends = [span.end for span in data.spans]
    pad = (max(ends) - min(starts)) * _POSITION_PAD or 1.0  # a lone SSR has a zero-width range
    axes.set_xlim(min(starts) - pad, max(ends) + pad)
    axes.set_yticks(range(len(units)), labels=list(units))
    axes.set_ylim(len(units) - 0.5, -0.5)  # first unit at the top, one axis unit per strip
    _style_panel(axes, "positions by unit", "GenBank position (1-based, inclusive)")


def _build(data: SsrDistributionSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        left, right = figure.subplots(1, 2, width_ratios=(2, 3))
        figure.suptitle(title, x=0.02, ha="left", fontsize=THEME.font_size, color=THEME.line)
        _regions_panel(left, data)
        _positions_panel(right, data)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the tables, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="ssr_distribution",
        title=text,
        source_tables=(f"ssr.{_BY_REGION}", f"ssr.{_SSRS}", f"ssr.{_SUMMARY}"),
        build=_build(data, text),
        figsize=(7.2, 3.6),
    )


SPEC = FigureSpec(
    id="ssr_distribution",
    title="SSR distribution by region and unit",
    requires=(("ssr", _BY_REGION), ("ssr", _SSRS), ("ssr", _SUMMARY)),
    build=build,
)
