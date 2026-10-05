"""π along the alignment with the plastid regions shaded — M5.1's sliding-window figure.

The source is ``diversity.diversity_windows`` and nothing else (plus the module's summary for the
two numbers the axis labels print). Every window already carries **the reference region it falls
in**, so the shading needs no mapping from alignment columns back to reference coordinates — which
is exactly the mapping no module exports. A window with no usable column stays a **gap**: a
missing π is a hole in the curve, never a zero.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result

# Headroom above the tallest π so the curve does not touch the frame — a layout constant, not a
# value read off the data.
_HEADROOM = 1.15
_WINDOWS = "diversity_windows"
_SUMMARY = "diversity_summary"
_COLUMNS = ("window_start_1b", "window_end_1b", "pi", "region")


@dataclass(frozen=True, slots=True)
class Window:
    """One alignment window, as the table states it: 1-based inclusive alignment columns."""

    start: int
    end: int
    pi: float | None
    region: str

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2


def _boundary(left: Window, right: Window) -> float:
    """Where one region's shading ends and the next begins — midway between the two windows."""
    return (left.midpoint + right.midpoint) / 2


@dataclass(frozen=True, slots=True)
class PiWindowSeries:
    """Everything the figure draws, read from the tables — the D5-03 audit trail."""

    accession: str
    organism: str
    window_bp: int
    step_bp: int
    alignment_length: int
    pi_whole: str  # the summary cell's own rendering — the figure prints no number of its own
    windows: tuple[Window, ...]

    @property
    def n_windows_without_value(self) -> int:
        """Windows whose π is blank — drawn as gaps; the count is worth stating in the caption."""
        return sum(1 for window in self.windows if window.pi is None)

    @property
    def bands(self) -> tuple[tuple[str, int, int], ...]:
        """The regions as shaded spans: ``(region, first column, last column)``.

        Windows overlap (600 bp window, 200 bp step), so taking each band from its windows' outer
        edges would paint neighbouring regions on top of each other. A boundary is placed **midway
        between the last window of one region and the first of the next** instead — the resolution
        the windowing actually has, which the figure does not pretend to beat.
        """
        runs: list[list[Window]] = []
        for window in self.windows:
            if not window.region:
                continue
            if runs and runs[-1][0].region == window.region:
                runs[-1].append(window)
            else:
                runs.append([window])
        bands: list[tuple[str, int, int]] = []
        for index, run in enumerate(runs):
            start = run[0].start if index == 0 else _boundary(runs[index - 1][-1], run[0])
            end = run[-1].end if index == len(runs) - 1 else _boundary(run[-1], runs[index + 1][0])
            bands.append((run[0].region, round(start), round(end)))
        return tuple(bands)


def series(results: Mapping[str, Result]) -> PiWindowSeries:
    """Read the diversity tables into the plain data the figure draws."""
    result = results.get("diversity")
    if result is None:
        raise ValueError(f"pi_window needs the 'diversity' module's result; have {sorted(results)}")
    missing = [name for name in (_WINDOWS, _SUMMARY) if name not in result.tables]
    if missing:
        raise ValueError(f"pi_window needs table(s) {missing} from 'diversity'")
    windows_table = result.tables[_WINDOWS]
    summary = result.tables[_SUMMARY]
    absent = [name for name in _COLUMNS if name not in {c.name for c in windows_table.columns}]
    if absent:
        raise ValueError(f"pi_window: {_WINDOWS} has no column(s) {absent}")

    return PiWindowSeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        window_bp=summary.value("window_bp"),
        step_bp=summary.value("step_bp"),
        alignment_length=summary.value("alignment_length"),
        pi_whole=summary.format("pi"),
        windows=tuple(
            Window(
                start=row["window_start_1b"],
                end=row["window_end_1b"],
                pi=row["pi"],
                region=row["region"] or "",
            )
            for row in windows_table.rows
        ),
    )


def caption(data: PiWindowSeries) -> str:
    """The line the figure carries, built from the series — and from no number of its own."""
    text = (
        f"π along the alignment — {data.organism} ({data.accession}) · "
        f"whole-alignment π {data.pi_whole}"
    )
    if data.n_windows_without_value:
        text += (
            f" · {data.n_windows_without_value} window(s) without a usable column (shown as gaps)"
        )
    return text


def _build(data: PiWindowSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        axes = figure.subplots()
        axes.set_title(title, loc="left", fontsize=THEME.font_size, color=THEME.line, pad=6)
        for region, start, end in data.bands:
            axes.axvspan(start, end, color=THEME.region_fill(region), linewidth=0, zorder=0)
        axes.plot(
            [window.midpoint for window in data.windows],
            [float("nan") if window.pi is None else window.pi for window in data.windows],
            color=THEME.line,
            linewidth=THEME.line_width,
            solid_joinstyle="round",
            zorder=3,
        )
        axes.set_xlabel(f"alignment column — {data.window_bp} bp window, {data.step_bp} bp step")
        axes.set_ylabel("π (per site)")
        axes.set_xlim(1, data.alignment_length)
        values = [window.pi for window in data.windows if window.pi is not None]
        if values:
            axes.set_ylim(0, max(values) * _HEADROOM)
        # Region names, drawn inside their band at the top — the shading is the label's legend.
        for region, start, end in data.bands:
            axes.text(
                (start + end) / 2,
                0.96,
                region,
                transform=axes.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=THEME.font_size - 1,
                color=THEME.muted,
            )
        axes.grid(axis="y", zorder=1)
        axes.set_axisbelow(True)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the tables, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="pi_window",
        title=text,
        source_tables=(f"diversity.{_WINDOWS}", f"diversity.{_SUMMARY}"),
        build=_build(data, text),
        figsize=(7.2, 2.8),
    )


SPEC = FigureSpec(
    id="pi_window",
    title="π sliding window with LSC/SSC/IR shading",
    requires=(("diversity", _WINDOWS), ("diversity", _SUMMARY)),
    build=build,
)
