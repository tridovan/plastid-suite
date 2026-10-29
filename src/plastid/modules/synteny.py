"""Synteny and % identity along the reference — the mVISTA-style comparison (P4-01, M4.2).

Two products, one alignment: a **block table** (one row per alignment, with the tool's own
``%IDY``) and an **% identity series** along the reference, windowed. The series is computed from
the delta file's own gap list rather than re-aligned here, so the table and the curve cannot
describe different alignments — and ``external/mummer.py`` refuses the run if they do.

Scope, stated because two numbers in this module could be mistaken for one another:

* the block table's ``%IDY`` is MUMmer's: ``100 · (1 − errors / alignment length)``, where an
  indel counts as an error;
* the series' identity compares **base pairs only** — a base the other genome does not carry is
  a gap, not a mismatch — so a block rich in indels reads slightly higher in the series than in
  the table. Both are labelled where they are exported.

The reference is aligned in its **normalized** frame (LSC → IRa → SSC → IRb, Phase 1), so a
window's coordinates mean the same thing in every pairwise plot against it; the IR gate is
required, as it is for every position-based analysis.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from plastid.core.delta import aligned_totals, window_identity
from plastid.core.model import Genome
from plastid.core.normalize import normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.external.mummer import run_nucmer
from plastid.modules.registry import ParamSpec, module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults

MODULE_ID = "synteny"
MODULE_VERSION = "0.1.0"


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Synteny and identity",
    requires_ir=True,
    arity=2,
    citations=("mummer",),
    params=(
        ParamSpec("window", "int", "window length in bp for the %identity series"),
        ParamSpec("step", "int", "window step in bp"),
        ParamSpec("min_match", "int", "nucmer -l: minimum exact match length"),
        ParamSpec("min_cluster", "int", "nucmer -c: minimum cluster length"),
        ParamSpec("filter", "str", "delta-filter mode (the toolkit default is -1)"),
        ParamSpec("mummer", "str", "directory holding nucmer/delta-filter/show-coords"),
        ParamSpec("groups", "str", "sample-group label for each input (supplied by --groups)"),
    ),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    reference = normalize(genomes[0], ir).genome
    reference_name = reference.accession or reference.name

    defaults = load_defaults()["synteny"]
    window = int(params.get("window", defaults["window_bp"]))
    step = int(params.get("step", defaults["step_bp"]))
    min_match = int(params.get("min_match", defaults["min_match"]))
    min_cluster = int(params.get("min_cluster", defaults["min_cluster"]))
    filter_mode = str(params.get("filter", defaults["filter"]))
    directory = params.get("mummer")
    labels: dict[str, str] = params.get("groups", {})

    block_rows: list[dict[str, Any]] = []
    genome_rows: list[dict[str, Any]] = []
    series_rows: list[dict[str, Any]] = []
    artifacts: dict[str, str] = {}

    for genome in genomes[1:]:
        name = genome.accession or genome.name
        result = run_nucmer(
            (reference_name, reference.sequence),
            (name, genome.sequence),
            ctx,
            directory=directory,
            min_match=min_match,
            min_cluster=min_cluster,
            filter_mode=filter_mode,
        )
        compared, mismatches = aligned_totals(
            result.delta_blocks, target=reference.sequence, query=genome.sequence
        )
        for index, row in enumerate(result.blocks, start=1):
            block_rows.append(
                {
                    "block": f"{name}#{index}",
                    "genome": name,
                    "group": labels.get(name, ""),
                    "target_start": row.target_start,
                    "target_end": row.target_end,
                    "query_start": row.query_start,
                    "query_end": row.query_end,
                    "length_target": row.length_target,
                    "length_query": row.length_query,
                    "identity": row.identity,
                    "orientation": "reverse" if row.reverse else "forward",
                }
            )
        identities = [row.identity for row in result.blocks]
        genome_rows.append(
            {
                "genome": name,
                "group": labels.get(name, ""),
                "n_blocks": len(result.blocks),
                "aligned_bp": compared,
                "coverage_pct": 100.0 * compared / len(reference),
                "min_identity": min(identities) if identities else None,
                "max_identity": max(identities) if identities else None,
                "mean_identity": (sum(identities) / len(identities)) if identities else None,
                "aligned_identity": (100.0 * (1 - mismatches / compared)) if compared else None,
            }
        )
        for entry in window_identity(
            result.delta_blocks,
            target=reference.sequence,
            query=genome.sequence,
            window=window,
            step=step,
        ):
            series_rows.append(
                {
                    "genome": name,
                    "group": labels.get(name, ""),
                    "window_start": entry.start_1b,
                    "window_end": entry.end_1b,
                    "compared_bp": entry.compared_bp,
                    "identity": entry.identity_pct,
                }
            )
        tag = _safe(name)
        artifacts[f"mummer_{tag}.coords"] = result.coords
        artifacts[f"mummer_{tag}.delta"] = result.delta

    if not block_rows:
        raise ValueError(
            f"no alignment block between {reference_name} and "
            f"{', '.join(genome.accession or genome.name for genome in genomes[1:])} — "
            f"an empty comparison is reported as zero blocks, and there are none to report"
        )

    block_scope = Scope(
        label=(
            "nucmer --maxmatch alignments filtered to one-to-one; target coordinates are the "
            "reference genome's normalized frame, query coordinates the query's own; %IDY is "
            "the tool's own (indels count as errors)"
        ),
        coord_system="normalized",
    )
    summary_scope = Scope(
        label=(
            "the reference and the parameters of this comparison; per-genome numbers are in "
            "synteny_genomes"
        ),
        coord_system="normalized",
    )
    genome_scope = Scope(
        label=(
            "each genome aligned against the reference with nucmer --maxmatch; identity is "
            "MUMmer's %IDY per block, mean is the unweighted mean over blocks; coverage is "
            "aligned base pairs over the reference's length"
        ),
        coord_system="normalized",
    )
    series_scope = Scope(
        label=(
            f"{window} bp window / {step} bp step along the reference's normalized coordinates; "
            "only aligned base pairs count (gaps and ambiguous bases are dropped, not scored), "
            "so the series reads slightly higher than %IDY in indel-rich blocks; a window with "
            "no aligned pair has no value"
        ),
        coord_system="normalized",
    )

    summary = Table(
        id="synteny_summary",
        title="Synteny summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("reference_length", "Reference length", dtype="int", unit="bp"),
            Column("n_genomes", "Genomes aligned", dtype="int"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference_name,
                "organism": reference.name,
                "reference_length": len(reference),
                "n_genomes": len(genomes) - 1,
                "window_bp": window,
                "step_bp": step,
                "tool": "MUMmer nucmer",
            },
        ),
        scope=summary_scope,
    )

    genomes_table = Table(
        id="synteny_genomes",
        title="Synteny by genome",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("n_blocks", "Alignment blocks", dtype="int"),
            Column("aligned_bp", "Aligned base pairs", dtype="int", unit="bp"),
            Column("coverage_pct", "Coverage", dtype="float", unit="%", precision=2),
            Column("min_identity", "Min %IDY", dtype="float", unit="%", precision=2),
            Column("max_identity", "Max %IDY", dtype="float", unit="%", precision=2),
            Column("mean_identity", "Mean %IDY", dtype="float", unit="%", precision=2),
            Column(
                "aligned_identity",
                "Identity over base pairs",
                dtype="float",
                unit="%",
                precision=2,
            ),
        ),
        rows=tuple(genome_rows),
        scope=genome_scope,
        key="genome",
    )

    blocks = Table(
        id="synteny_blocks",
        title="Alignment blocks against the reference",
        columns=(
            Column("block", "Block"),
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("target_start", "Target start", dtype="int", unit="bp"),
            Column("target_end", "Target end", dtype="int", unit="bp"),
            Column("query_start", "Query start", dtype="int", unit="bp"),
            Column("query_end", "Query end", dtype="int", unit="bp"),
            Column("length_target", "Target length", dtype="int", unit="bp"),
            Column("length_query", "Query length", dtype="int", unit="bp"),
            Column("identity", "%IDY", dtype="float", unit="%", precision=2),
            Column("orientation", "Orientation"),
        ),
        rows=tuple(block_rows),
        scope=block_scope,
        key="block",
    )

    series = Table(
        id="synteny_identity",
        title="% identity along the reference",
        columns=(
            Column("genome", "Genome"),
            Column("group", "Group"),
            Column("window_start", "Window start", dtype="int", unit="bp"),
            Column("window_end", "Window end", dtype="int", unit="bp"),
            Column("compared_bp", "Aligned base pairs", dtype="int", unit="bp"),
            Column("identity", "% identity", dtype="float", unit="%", precision=2),
        ),
        rows=tuple(series_rows),
        scope=series_scope,
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "synteny_summary": summary,
            "synteny_genomes": genomes_table,
            "synteny_blocks": blocks,
            "synteny_identity": series,
        },
        notes=(
            "One row per alignment block as nucmer reports it; the series is read from the same "
            "delta file, so the table and the curve describe one alignment.",
            f"A window with no aligned base pair has no identity — {window} bp window, "
            f"{step} bp step.",
            "The pipeline is nucmer --maxmatch -> delta-filter -> "
            "show-coords; the parameters of this run are in the manifest's argv.",
        ),
        citations=("mummer",),
        artifacts=artifacts,
    )


def _safe(name: str) -> str:
    """An artifact-name fragment: letters, digits, dot, dash and underscore survive."""
    return "".join(
        character if character.isalnum() or character in ".-_" else "_" for character in name
    )
