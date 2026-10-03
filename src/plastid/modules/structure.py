"""Genome structure — region boundaries, %GC and GC-skew per region (PROJECT.md §13.2).

Golden numbers (§6): LSC 83,127 · IRa = IRb = 26,214 · SSC 18,169 — the four regions tile the
genome exactly (153,724 bp). Boundaries come from the IR pair and are reported in the
**normalized** frame (LSC → IRa → SSC → IRb); ``rotation_offset`` maps back to GenBank.

Both skew column headers state their convention — ``(G − C)/(G + C)``, and a cumulative column
that is the running sum *within* its region — because a silent sign flip is this phase's
classic quiet error (TASKS.md P2-02).
"""

from __future__ import annotations

from typing import Any

from plastid.core.model import Genome
from plastid.core.normalize import Region, normalize
from plastid.core.sequence import gc_percent, gc_skew
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "structure"
MODULE_VERSION = "0.1.0"

DEFAULT_WINDOW_BP = 1000
DEFAULT_STEP_BP = 1000


@module(
    MODULE_ID,
    version=MODULE_VERSION,
    title="Genome structure",
    requires_ir=True,
    citations=("biopython",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    window = int(params.get("window", DEFAULT_WINDOW_BP))
    step = int(params.get("step", DEFAULT_STEP_BP))
    if window <= 0 or step <= 0:
        raise ValueError(f"window and step must be positive, got window={window}, step={step}")
    if step > window:
        raise ValueError(
            f"step {step} exceeds window {window} — sequence between windows is skipped"
        )

    frame = normalize(genome, ir)
    region_rows: list[dict[str, Any]] = []
    window_rows: list[dict[str, Any]] = []
    for region in frame.regions:
        sequence = frame.genome.sequence[region.start : region.end]
        region_rows.append(
            {
                "region": region.name,
                "start_1b": region.start + 1,
                "end_1b": region.end,
                "length_bp": region.length,
                "gc_pct": round(gc_percent(sequence), 4),
                "gc_skew": round(gc_skew(sequence), 6),
            }
        )
        window_rows.extend(_window_rows(region, sequence, window, step))

    summary = Table(
        id="structure_summary",
        title="Structure summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("n_windows", "Windows", dtype="int"),
            Column("rotation_offset", "Rotation offset", dtype="int", unit="bp"),
            Column("coord_system", "Coordinates"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "window_bp": window,
                "step_bp": step,
                "n_windows": len(window_rows),
                "rotation_offset": frame.rotation_offset,
                "coord_system": "normalized",
            },
        ),
        scope=Scope(
            label="one row; the parameters the windowed table was built with",
            coord_system="normalized",
        ),
    )

    regions = Table(
        id="regions",
        title="Regions",
        columns=(
            Column("region", "Region"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("gc_pct", "GC", dtype="float", unit="%", precision=2),
            Column("gc_skew", "GC-skew (G−C)/(G+C)", dtype="float", precision=4),
        ),
        rows=tuple(region_rows),
        scope=Scope(
            label="normalized frame; boundaries from the IR pair; whole-region values",
            coord_system="normalized",
        ),
        key="region",
    )

    windows = Table(
        id="gc_skew_windows",
        title="GC-skew windows",
        columns=(
            Column("region", "Region"),
            Column("window_start_1b", "Start", dtype="int"),
            Column("window_end_1b", "End", dtype="int"),
            Column("window_length_bp", "Window length", dtype="int", unit="bp"),
            Column("gc_skew", "Window GC-skew (G−C)/(G+C)", dtype="float", precision=4),
            Column(
                "cumulative_gc_skew",
                "Cumulative GC-skew (running sum within region)",
                dtype="float",
                precision=4,
            ),
        ),
        rows=tuple(window_rows),
        scope=Scope(
            label=(
                f"normalized frame; {window} bp window / {step} bp step; "
                "the last window of a region may be shorter; cumulative resets per region"
            ),
            coord_system="normalized",
        ),
    )

    total = sum(row["length_bp"] for row in region_rows)
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={"structure_summary": summary, "regions": regions, "gc_skew_windows": windows},
        notes=(
            f"Regions tile the genome exactly: {total} bp over {len(region_rows)} regions.",
            "GC-skew is (G − C)/(G + C) per window; the cumulative column is the running sum.",
        ),
        citations=("biopython",),
    )


def _window_rows(region: Region, sequence: str, window: int, step: int) -> list[dict[str, Any]]:
    """One row per window inside ``region``; the cumulative skew restarts at each region."""
    rows: list[dict[str, Any]] = []
    cumulative = 0.0
    offset = 0
    while offset < region.length:
        end = min(offset + window, region.length)
        value = gc_skew(sequence[offset:end])
        cumulative += value
        rows.append(
            {
                "region": region.name,
                "window_start_1b": region.start + offset + 1,
                "window_end_1b": region.start + end,
                "window_length_bp": end - offset,
                "gc_skew": round(value, 6),
                "cumulative_gc_skew": round(cumulative, 6),
            }
        )
        if end == region.length:
            break
        offset += step
    return rows
