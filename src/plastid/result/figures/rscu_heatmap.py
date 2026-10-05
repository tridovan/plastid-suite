"""RSCU heatmap — one row per amino-acid family, one column per codon of that family.

Sources: ``codon_usage.rscu`` (the 61 sense codons, in family order) and the module's one-row
summary for the caption's values. Every number drawn is a table cell — a cell's RSCU text is the
column's own ``format`` rendering, so the figure cannot round a value twice — and a blank RSCU
stays blank, never a zero.

RSCU 1.0 is parity (a family's codons used equally often), so the scale is **diverging**, centred
there: cool below, warm above, near-white at parity (two hues, neutral midpoint, written out from
the theme's line/accent rather than a named map). Each cell also carries its codon — codons are
packed left within a family, so a bare grid of numbers would not say which codon it draws.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from plastid.result.figure import THEME, Figure, FigureSpec, export_figure
from plastid.result.result import Result

_RSCU = "rscu"
_SUMMARY = "codon_usage"
_RSCU_COLUMNS = ("codon", "amino_acid", "count", "rscu")
_SUMMARY_COLUMNS = ("accession", "organism", "n_cds", "enc")

# The diverging ramp: THEME.line (cool) -> tints -> near-white parity -> THEME.accent (warm).
_RSCU_STOPS = (THEME.line, "#8ba4bc", "#f2f4f6", "#dfa591", THEME.accent)


@dataclass(frozen=True, slots=True)
class RscuEntry:
    """One sense codon, exactly as the table's row states it."""

    amino_acid: str
    codon: str
    count: int | None
    rscu: float | None
    text: str  # the rscu column's own rendering — the figure prints no number of its own


@dataclass(frozen=True, slots=True)
class RscuHeatmapSeries:
    """Everything the figure draws, read from the tables — the D5-03 audit trail."""

    accession: str
    organism: str
    n_cds: str
    enc: str
    entries: tuple[RscuEntry, ...]

    @property
    def codons(self) -> tuple[tuple[str, str, int | None, float | None], ...]:
        """The rows as plain tuples, so a test can compare them with the table they came from."""
        return tuple(
            (entry.amino_acid, entry.codon, entry.count, entry.rscu) for entry in self.entries
        )

    @property
    def families(self) -> tuple[tuple[str, tuple[RscuEntry, ...]], ...]:
        """Entries grouped per amino acid, families in first-appearance order and codons in table
        order — the table lists each family's codons consecutively."""
        grouped: dict[str, list[RscuEntry]] = {}
        for entry in self.entries:
            grouped.setdefault(entry.amino_acid, []).append(entry)
        return tuple((name, tuple(members)) for name, members in grouped.items())


def series(results: Mapping[str, Result]) -> RscuHeatmapSeries:
    """Read the codon-usage tables into the plain data the figure draws."""
    result = results.get("codon_usage")
    if result is None:
        raise ValueError(
            f"rscu_heatmap needs the 'codon_usage' module's result; have {sorted(results)}"
        )
    missing = [name for name in (_RSCU, _SUMMARY) if name not in result.tables]
    if missing:
        raise ValueError(f"rscu_heatmap needs table(s) {missing} from 'codon_usage'")
    rscu_table = result.tables[_RSCU]
    summary = result.tables[_SUMMARY]
    for table, columns in ((rscu_table, _RSCU_COLUMNS), (summary, _SUMMARY_COLUMNS)):
        absent = [name for name in columns if name not in {c.name for c in table.columns}]
        if absent:
            raise ValueError(f"rscu_heatmap: {table.id} has no column(s) {absent}")

    return RscuHeatmapSeries(
        accession=summary.format("accession"),
        organism=summary.format("organism"),
        n_cds=summary.format("n_cds"),
        enc=summary.format("enc"),
        entries=tuple(
            RscuEntry(
                amino_acid=row["amino_acid"],
                codon=row["codon"],
                count=row["count"],
                rscu=row["rscu"],
                text=rscu_table.format(f"{row['codon']}.rscu"),
            )
            for row in rscu_table.rows
        ),
    )


def caption(data: RscuHeatmapSeries) -> str:
    """The line the figure carries, built from the summary cells — and from no number of its own."""
    return (
        f"RSCU by amino-acid family — {data.organism} ({data.accession}) · "
        f"{data.n_cds} CDS · ENC {data.enc}"
    )


def _luminance(colour: str) -> float:
    """WCAG relative luminance of an ``#rrggbb`` colour."""
    red, green, blue = (
        value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
        for value in (int(colour[index : index + 2], 16) / 255 for index in (1, 3, 5))
    )
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def _ink(fill: str) -> str:
    """The theme ink with the higher WCAG contrast ratio on ``fill`` — small text on a dark cell
    must stay legible, so the ramp's dark ends get background ink and its middle gets line ink."""

    def contrast(other: str) -> float:
        other_luminance = _luminance(other)
        fill_luminance = _luminance(fill)
        return (max(other_luminance, fill_luminance) + 0.05) / (
            min(other_luminance, fill_luminance) + 0.05
        )

    return max((THEME.background, THEME.line), key=contrast)


def _build(data: RscuHeatmapSeries, title: str) -> Any:
    def build(figure: Any) -> None:
        from matplotlib.colors import LinearSegmentedColormap, to_hex

        axes = figure.subplots()
        families = data.families
        values = [entry.rscu for entry in data.entries if entry.rscu is not None]
        # Symmetric limits put RSCU 1.0 exactly at the ramp's midpoint, whatever the data range.
        span = max((abs(value - 1.0) for value in values), default=0.0) or 1e-9
        low, high = 1.0 - span, 1.0 + span
        ramp = LinearSegmentedColormap.from_list("rscu", _RSCU_STOPS)
        ramp.set_bad(THEME.background)  # a blank RSCU is a gap, and so is an unused trailing slot
        width = max(len(members) for _, members in families)
        grid = [
            [float("nan") if entry.rscu is None else entry.rscu for entry in members]
            + [float("nan")] * (width - len(members))
            for _, members in families
        ]
        image = axes.imshow(
            grid, cmap=ramp, vmin=low, vmax=high, aspect="auto", interpolation="nearest"
        )
        for row, (_, members) in enumerate(families):
            for column, entry in enumerate(members):
                fill = (
                    to_hex(ramp((entry.rscu - low) / (high - low)))
                    if entry.rscu is not None
                    else THEME.background
                )
                axes.text(
                    column,
                    row,
                    f"{entry.codon} {entry.text}".strip(),
                    ha="center",
                    va="center",
                    fontsize=THEME.font_size - 3,
                    color=_ink(fill),
                )
        axes.set_title(title, loc="left", fontsize=THEME.font_size, color=THEME.line, pad=6)
        axes.set_yticks(range(len(families)), labels=[name for name, _ in families])
        axes.set_xticks([])
        axes.set_xticks([column - 0.5 for column in range(width + 1)], minor=True)
        axes.set_yticks([row - 0.5 for row in range(len(families) + 1)], minor=True)
        axes.tick_params(which="both", length=0)
        axes.grid(which="minor", color=THEME.background, linewidth=1.0)
        axes.set_axisbelow(True)
        axes.set_xlabel("codon within amino-acid family")
        axes.set_ylabel("amino acid")
        axes.set_xlim(-0.5, width - 0.5)
        figure.colorbar(image, ax=axes, fraction=0.035, pad=0.02).set_label("RSCU")

    return build


def build(results: Mapping[str, Result]) -> Figure:
    """The figure, end to end — series from the tables, then the shared vector export."""
    data = series(results)
    text = caption(data)
    return export_figure(
        figure_id="rscu_heatmap",
        title=text,
        source_tables=(f"codon_usage.{_RSCU}", f"codon_usage.{_SUMMARY}"),
        build=_build(data, text),
        figsize=(7.2, 5.8),
    )


SPEC = FigureSpec(
    id="rscu_heatmap",
    title="RSCU heatmap by amino-acid family",
    requires=(("codon_usage", _RSCU), ("codon_usage", _SUMMARY)),
    build=build,
)
