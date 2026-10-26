"""mVISTA-style % identity along the reference, with the synteny blocks beneath it.

The source is ``synteny`` and nothing else: the windowed % identity series, the alignment blocks,
the run's own parameters. Every number drawn is a table cell — the figure recomputes nothing
(D5-03). Window coordinates are the **reference's normalized frame**, the same frame the block
table's ``target_*`` columns use, so the two panels share one x axis; and a window with no aligned
pair has no identity, so it is drawn as a **gap** (``nan``), never as zero.

The per-genome series is carried as plain tuples in table order — one
``(genome, group, ((window_start, window_end, identity), ...))`` per genome — so a test can
compare the figure's own data with the table cell by cell (the P5-01 pattern, and the shape the
N-scaling measurement reads).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result, Table

_IDENTITY = "synteny_identity"
_BLOCKS = "synteny_blocks"
_SUMMARY = "synteny_summary"
_IDENTITY_COLUMNS = ("genome", "group", "window_start", "window_end", "identity")
_BLOCK_COLUMNS = ("genome", "target_start", "target_end", "orientation")
_SUMMARY_COLUMNS = (
    "accession",
    "organism",
    "reference_length",
    "n_genomes",
    "window_bp",
    "step_bp",
    "tool",
)

# (window_start, window_end, identity): one window's three cells, a blank identity kept as None.
IdentityWindow = tuple[int, int, float | None]
# (genome, group, windows): one per genome, in the table's own order.
GenomeIdentity = tuple[str, str, tuple[IdentityWindow, ...]]
# (genome, target_start, target_end, orientation): one alignment block, the reference's frame.
SyntenyBlock = tuple[str, int, int, str]

# One colour per series, cycled past the third; the theme owns the palette, the figure picks.
_SERIES_COLOURS = (THEME.line, THEME.accent, THEME.muted)
_BLOCK_COLOURS = {"forward": THEME.line, "reverse": THEME.accent}
_BAR_HEIGHT = 0.6

# Measured, not guessed (P5-02): the same figure at N = 2, 8 and 34. At 8 the seven
# series and their legend read; at 34 the legend covered the panel and 33 lines merged into one
# grey mass — an unreadable figure is a failed figure, not a dense one (D5-04). Above the cap the
# panel switches to the span the genomes cover, whose two edges are still cell values.
_SERIES_CAP = 8
_ENVELOPE_FILL = "#c9d2dc"  # THEME.line blended into the background — no alpha (EPS drops it)


@dataclass(frozen=True, slots=True)
class IdentitySyntenySeries:
    """Everything the figure draws, read from the tables — the D5-03 audit trail."""

    accession: str
    organism: str
    reference_length: int
    n_genomes: int
    window_bp: int
    step_bp: int
    tool: str
    per_genome: tuple[GenomeIdentity, ...]
    blocks: tuple[SyntenyBlock, ...]


def _require_columns(table: Table, names: tuple[str, ...]) -> None:
    absent = [name for name in names if name not in {column.name for column in table.columns}]
    if absent:
        raise ValueError(f"identity_synteny: {table.id} has no column(s) {absent}")


def series(results: Mapping[str, Result]) -> IdentitySyntenySeries:
    """Read the synteny tables into the plain data the figure draws."""
    result = results.get("synteny")
    if result is None:
        raise ValueError(
            f"identity_synteny needs the 'synteny' module's result; have {sorted(results)}"
        )
    missing = [name for name in (_IDENTITY, _BLOCKS, _SUMMARY) if name not in result.tables]
    if missing:
        raise ValueError(f"identity_synteny needs table(s) {missing} from 'synteny'")
    identity = result.tables[_IDENTITY]
    blocks = result.tables[_BLOCKS]
    summary = result.tables[_SUMMARY]
    _require_columns(identity, _IDENTITY_COLUMNS)
    _require_columns(blocks, _BLOCK_COLUMNS)
    _require_columns(summary, _SUMMARY_COLUMNS)

    # Group rows by genome, first encounter first: the table is ordered by genome then by window.
    grouped: dict[str, list[IdentityWindow]] = {}
    groups: dict[str, str] = {}
    for row in identity.rows:
        genome = row["genome"]
        groups.setdefault(genome, row["group"] or "")
        grouped.setdefault(genome, []).append(
            (row["window_start"], row["window_end"], row["identity"])
        )

    return IdentitySyntenySeries(
        accession=summary.value("accession"),
        organism=summary.value("organism"),
        reference_length=summary.value("reference_length"),
        n_genomes=summary.value("n_genomes"),
        window_bp=summary.value("window_bp"),
        step_bp=summary.value("step_bp"),
        tool=summary.value("tool"),
        per_genome=tuple((name, groups[name], tuple(windows)) for name, windows in grouped.items()),
        blocks=tuple(
            (row["genome"], row["target_start"], row["target_end"], row["orientation"])
            for row in blocks.rows
        ),
    )


def caption(data: IdentitySyntenySeries) -> str:
    """The line the figure carries, built from the series — and from no number of its own."""
    return (
        f"mVISTA-style % identity along the reference — {data.organism} ({data.accession}) · "
        f"{data.n_genomes} genome(s) · {data.window_bp} bp window / {data.step_bp} bp step · "
        f"{data.tool}"
    )


def _draw_series(top: Any, data: IdentitySyntenySeries) -> None:
    """One line per genome, up to the cap — see ``_SERIES_CAP``."""
    for index, (genome, group, windows) in enumerate(data.per_genome):
        top.plot(
            [(start + end) / 2 for start, end, _ in windows],
            [identity if identity is not None else float("nan") for _, _, identity in windows],
            color=_SERIES_COLOURS[index % len(_SERIES_COLOURS)],
            linewidth=THEME.line_width,
            label=f"{genome} ({group})" if group else genome,
        )
    top.legend(
        loc="lower left",
        ncols=min(4, len(data.per_genome)),
        frameon=False,
        fontsize=THEME.font_size - 1,
    )


def _draw_envelope(top: Any, data: IdentitySyntenySeries) -> None:
    """The span the genomes cover, window by window — above the cap, one line per genome is ash.

    Both edges are **cell values**: the lowest and the highest identity any genome reports at
    that window. No number is invented, and the table still carries every series — the figure
    has merely stopped pretending 33 overlapping lines can be told apart.
    """
    columns: dict[int, list[float]] = {}
    mids: dict[int, float] = {}
    for _, _, windows in data.per_genome:
        for start, end, identity in windows:
            mids[start] = (start + end) / 2
            if identity is not None:
                columns.setdefault(start, []).append(identity)
    xs = sorted(columns)
    if not xs:
        return
    middle = [mids[start] for start in xs]
    lows = [min(columns[start]) for start in xs]
    highs = [max(columns[start]) for start in xs]
    # A solid tint rather than alpha: the PostScript backend warns and renders translucent
    # artists opaque, which would make the EPS copy differ from the SVG and the PDF.
    top.fill_between(
        middle,
        lows,
        highs,
        color=_ENVELOPE_FILL,
        linewidth=0,
        label=f"span of {len(data.per_genome)} genomes (lowest to highest)",
    )


def _build(data: IdentitySyntenySeries, title: str) -> Any:
    def build(figure: Any) -> None:
        rows = max(len({name for name, _, _ in data.per_genome} | {b[0] for b in data.blocks}), 1)
        top, bottom = figure.subplots(
            2,
            1,
            sharex=True,
            # The block panel is one row per genome, so it takes the space its rows need rather
            # than a fixed third of the canvas (measured at N = 34).
            height_ratios=(2, max(1.0, rows / 8)),
            gridspec_kw={"hspace": 0.08},
        )
        top.set_title(title, loc="left", fontsize=THEME.font_size, color=THEME.line, pad=6)
        if len(data.per_genome) <= _SERIES_CAP:
            _draw_series(top, data)
        else:
            _draw_envelope(top, data)
        top.set_ylabel("% identity")
        top.set_xlim(1, data.reference_length)
        values = [value for _, _, windows in data.per_genome for _, _, value in windows]
        values = [value for value in values if value is not None]
        if values:
            padding = max((max(values) - min(values)) * 0.1, 1.0)
            top.set_ylim(max(0.0, min(values) - padding), min(100.0, max(values) + padding))
        top.grid(axis="y", zorder=1)
        top.set_axisbelow(True)

        # One bar row per genome — the series' own genome order, extended when a genome appears
        # in the block table but has no window there.
        order = [name for name, _, _ in data.per_genome]
        for block in data.blocks:
            if block[0] not in order:
                order.append(block[0])
        rows = {name: index for index, name in enumerate(order)}
        for orientation in dict.fromkeys(block[3] for block in data.blocks):
            selected = [block for block in data.blocks if block[3] == orientation]
            bottom.barh(
                [rows[block[0]] for block in selected],
                [block[2] - block[1] + 1 for block in selected],
                left=[block[1] for block in selected],
                height=_BAR_HEIGHT,
                color=_BLOCK_COLOURS.get(orientation, THEME.muted),
                label=orientation,
            )
        bottom.set_yticks(range(len(order)))
        bottom.set_yticklabels(order)
        bottom.set_ylim(-0.6, len(order) - 0.4)
        bottom.set_ylabel("genome")
        bottom.set_xlabel("reference position (bp, normalized frame)")
        bottom.legend(loc="lower right", frameon=False, fontsize=THEME.font_size - 1)
        bottom.grid(axis="x", zorder=1)
        bottom.set_axisbelow(True)

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the tables, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="identity_synteny",
        title=text,
        source_tables=(f"synteny.{_IDENTITY}", f"synteny.{_BLOCKS}", f"synteny.{_SUMMARY}"),
        build=_build(data, text),
        # The block panel is one row per genome; a fixed height would smear 33 labels into a
        # line, so the canvas grows with them (measured at N = 34).
        figsize=(7.2, 3.0 + 0.12 * max(len(data.per_genome), 1)),
    )


SPEC = FigureSpec(
    id="identity_synteny",
    title="mVISTA-style % identity with the synteny blocks",
    requires=(("synteny", _IDENTITY), ("synteny", _BLOCKS), ("synteny", _SUMMARY)),
    build=build,
)
