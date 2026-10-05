"""π along the reference with the candidate barcode loci marked — M4.3's picture.

The source is ``barcodes``: the windowed π, the candidate loci, the run's own numbers. The figure
recomputes nothing — the threshold line prints the summary cell's own rendering (D5-03).

Choices worth stating:

* x is the **reference** coordinate: a window sits at the midpoint of its own reference span, and a
  window the run could not place is a **gap** (``nan``), like a blank π;
* a candidate whose span is blank gets **no band** — a band at zero would be a wrong claim;
* a locus the reference does not carry (``found == "no"``) still gets its band, in the muted
  colour: absence is a result, not a reason to drop the row;
* a locus that crosses position 1 (``reference_start`` > ``reference_end``) is drawn as **two**
  bands out to the range the cells place numbers in — one band would shade the locus's complement
  (measured on the golden genome's trnH-psbA spacer).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result, Table

_WINDOWS = "barcodes_windows"
_CANDIDATES = "barcodes_candidates"
_SUMMARY = "barcodes_summary"
_WINDOW_COLUMNS = ("reference_start", "reference_end", "pi")
_CANDIDATE_COLUMNS = ("locus", "found", "reference_start", "reference_end", "note")
_SUMMARY_COLUMNS = (
    "accession",
    "organism",
    "n_genomes",
    "threshold_pi",
    "threshold_rule",
    "n_candidates",
    "n_found",
)

# (reference_start, reference_end, pi): one window's three cells, blanks kept as None — a blank is
# a gap, never a zero.
BarcodeWindow = tuple[int | None, int | None, float | None]
# (locus, found, reference_start, reference_end, note): ``found`` is the table's own
# "yes"/"no" cell, and ``note`` is what says whether a span that runs backwards wraps the
# origin — see `_bands`.
BarcodeCandidate = tuple[str, str, int | None, int | None, str]

_HEADROOM = 1.15
_BAND_TINT = 0.18  # how far a band's colour is blended into the background


@dataclass(frozen=True, slots=True)
class BarcodeAssessmentSeries:
    """Everything the figure draws, read from the tables — the D5-03 audit trail."""

    accession: str
    organism: str
    n_genomes: int
    n_candidates: int
    n_found: int
    threshold_pi: float | None
    threshold_label: str  # the cell's own rendering; "" when the threshold cell is blank
    threshold_rule: str
    reference_limits: tuple[int, int] | None  # from the cells, never an invented span
    windows: tuple[BarcodeWindow, ...]
    candidates: tuple[BarcodeCandidate, ...]


def _require_columns(table: Table, names: tuple[str, ...]) -> None:
    absent = [name for name in names if name not in {column.name for column in table.columns}]
    if absent:
        raise ValueError(f"barcode_assessment: {table.id} has no column(s) {absent}")


def _pale(colour: str) -> str:
    """Blend a theme colour into the background — a solid tint, because EPS drops ``alpha``.

    matplotlib's PostScript backend warns and renders partially transparent artists opaque; a band
    in the π line's own colour would then hide the curve in the EPS copy.
    """
    base = bytes.fromhex(THEME.background.lstrip("#"))
    tint = bytes.fromhex(colour.lstrip("#"))
    blend = (round(b + (t - b) * _BAND_TINT) for b, t in zip(base, tint, strict=True))
    return "#" + bytes(blend).hex()


_ORIGIN_NOTE = "spans the origin"


def _bands(start: int, end: int, note: str, limits: tuple[int, int]) -> tuple[tuple[int, int], ...]:
    """Where to shade for a locus whose two coordinates may run either way.

    ``reference_start > reference_end`` is **not** by itself a wrap: measured over both Phase-4
    deliverable runs, the one such candidate is ``rps16`` — a two-part gene (its intron) whose
    parts descend — and its ``note`` is empty. Shading the complement for it painted almost the
    whole axis. Only the module's own note says a span crosses position 1, and the module writes
    that note precisely because the two cases are indistinguishable from the numbers alone.
    """
    if start <= end:
        return ((start, end),)
    if note.startswith(_ORIGIN_NOTE):
        return ((start, limits[1]), (limits[0], end))
    return ((end, start),)


def series(results: Mapping[str, Result]) -> BarcodeAssessmentSeries:
    """Read the barcode tables into the plain data the figure draws."""
    result = results.get("barcodes")
    if result is None:
        raise ValueError(
            f"barcode_assessment needs the 'barcodes' module's result; have {sorted(results)}"
        )
    missing = [name for name in (_WINDOWS, _CANDIDATES, _SUMMARY) if name not in result.tables]
    if missing:
        raise ValueError(f"barcode_assessment needs table(s) {missing} from 'barcodes'")
    windows = result.tables[_WINDOWS]
    candidates = result.tables[_CANDIDATES]
    summary = result.tables[_SUMMARY]
    _require_columns(windows, _WINDOW_COLUMNS)
    _require_columns(candidates, _CANDIDATE_COLUMNS)
    _require_columns(summary, _SUMMARY_COLUMNS)

    # The x limits are the coordinates the cells actually carry — no axis span is invented.
    coordinates = [
        row[name]
        for row in (*windows.rows, *candidates.rows)
        for name in ("reference_start", "reference_end")
        if row[name] is not None
    ]
    return BarcodeAssessmentSeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        n_genomes=summary.value("n_genomes"),
        n_candidates=summary.value("n_candidates"),
        n_found=summary.value("n_found"),
        threshold_pi=summary.value("threshold_pi"),
        threshold_label=summary.format("threshold_pi"),
        threshold_rule=summary.value("threshold_rule"),
        reference_limits=(min(coordinates), max(coordinates)) if coordinates else None,
        windows=tuple(
            (row["reference_start"], row["reference_end"], row["pi"]) for row in windows.rows
        ),
        candidates=tuple(
            (
                row["locus"],
                row["found"],
                row["reference_start"],
                row["reference_end"],
                row["note"],
            )
            for row in candidates.rows
        ),
    )


def caption(data: BarcodeAssessmentSeries) -> str:
    """The line the figure carries, built from the series — and from no number of its own."""
    return (
        f"Barcode candidate assessment — {data.organism} ({data.accession}) · "
        f"{data.n_candidates} candidate(s), {data.n_found} found · "
        f"{data.n_genomes} genome(s) aligned"
    )


def _build(data: BarcodeAssessmentSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        axes = figure.subplots()
        axes.set_title(title, loc="left", fontsize=THEME.font_size, color=THEME.line, pad=6)

        # The candidate bands first, so the π curve draws over them; the patch colours are the
        # legend, and the locus name is written inside the band.
        handles: dict[str, Any] = {}
        limits = data.reference_limits
        for locus, found_cell, start, end, note in data.candidates:
            if start is None or end is None or limits is None:
                continue
            found = found_cell == "yes"
            colour = THEME.line if found else THEME.muted
            bands = _bands(start, end, note, limits)
            for index, (low, high) in enumerate(bands):
                patch = axes.axvspan(low, high, color=_pale(colour), zorder=0)
                if index == 0:
                    handles.setdefault("found" if found else "not found", patch)
            low, high = max(bands, key=lambda band: band[1] - band[0])
            axes.text(
                (low + high) / 2,
                0.97,
                locus,
                transform=axes.get_xaxis_transform(),
                ha="center",
                va="top",
                rotation=90,
                fontsize=THEME.font_size - 2,
                color=colour,
            )

        # x is the midpoint of the window's own reference span — an unplaced window is a gap.
        xs = [
            float("nan") if start is None or end is None else (start + end) / 2
            for start, end, _ in data.windows
        ]
        ys = [float("nan") if pi is None else pi for _, _, pi in data.windows]
        (line,) = axes.plot(xs, ys, color=THEME.line, linewidth=THEME.line_width, zorder=3)
        if data.threshold_pi is not None:
            axes.axhline(data.threshold_pi, color=THEME.accent, linestyle="--", zorder=2)
            axes.text(
                0.99,
                data.threshold_pi,
                f"π = {data.threshold_label} · {data.threshold_rule}",
                transform=axes.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=THEME.font_size - 1,
                color=THEME.accent,
            )
        axes.set_xlabel("reference position (bp, normalized frame)")
        axes.set_ylabel("π (per site)")
        if limits is not None:
            axes.set_xlim(*limits)
        heights = [pi for _, _, pi in data.windows if pi is not None]
        if data.threshold_pi is not None:
            heights.append(data.threshold_pi)
        if heights:
            axes.set_ylim(0, max(heights) * _HEADROOM)
        axes.legend(
            [line, *handles.values()],
            ["π", *handles],
            loc="upper right",
            frameon=False,
            fontsize=THEME.font_size - 1,
        )
        axes.grid(axis="y", zorder=1)
        axes.set_axisbelow(True)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the tables, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="barcode_assessment",
        title=text,
        source_tables=(f"barcodes.{_WINDOWS}", f"barcodes.{_CANDIDATES}", f"barcodes.{_SUMMARY}"),
        build=_build(data, text),
        figsize=(7.2, 3.2),
    )


SPEC = FigureSpec(
    id="barcode_assessment",
    title="Barcode candidate assessment",
    requires=(("barcodes", _WINDOWS), ("barcodes", _CANDIDATES), ("barcodes", _SUMMARY)),
    build=build,
)
