"""Barcode candidates and hypervariable regions (P4-04, M4.3).

Three outputs, one π: the window series (the same quantity :mod:`plastid.modules.diversity`
reports — ``core.align.window_diversity``, same window and step), the **hypervariable regions**
called from it, and the **candidate barcode loci** the reference actually carries, each with the
π of its own span.

Reuse, not re-computation, is the rule here: π comes from ``core.align.nucleotide_diversity`` and
the series from ``core.align.window_diversity`` — the *same* functions the diversity module calls
— so a window's π in this table and in ``diversity_windows`` are the same number by construction,
and ``tests/golden/test_tier_b.py`` checks exactly that on a real run rather than trusting it.

What the module refuses to invent:

* the region rule is stated and exported (``mean + k·sd`` over this run's own windows, k a
  parameter; a region is at least ``min_windows`` consecutive windows), because "hypervariable
  region" without its rule is a shape somebody eyeballed;
* the candidate list is **data** (``data/barcode_loci.yaml``), each locus with the source of its
  claim, and a locus the reference does not carry is reported as **not found** with the reason —
  never dropped from the table;
* coordinates are the reference genome's normalized frame, and the scope label says so.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from plastid.core.align import (
    Alignment,
    columns_for_span,
    nucleotide_diversity,
    reference_positions,
    window_diversity,
)
from plastid.core.barcodes import (
    DEFAULT_MIN_WINDOWS,
    DEFAULT_THRESHOLD_SD,
    hypervariable_regions,
    load_barcode_loci,
    threshold_for,
)
from plastid.core.gene_names import aliases_for, normalize_genome
from plastid.core.model import GENE_KINDS, Gene, Genome
from plastid.core.normalize import normalize
from plastid.core.organelle import GROUP_PLASTID
from plastid.external.mafft import align
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext
from plastid.run.defaults import load_defaults

MODULE_ID = "barcodes"
MODULE_VERSION = "0.2.0"  # behaviour changed 2026-10-04 after the first runs (see §9.C-35)


@module(
    MODULE_ID,
    group=GROUP_PLASTID,
    version=MODULE_VERSION,
    title="Barcode candidates and hypervariable regions",
    # positions along the reference: the same gate every position-based module passes
    requires_ir=True,
    arity=2,
    citations=("mafft", "cbol"),
)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    reference = normalize(genomes[0], ir).genome
    reference_name = reference.accession or reference.name

    # The same defaults diversity uses: the two tables are read side by side, and a different
    # window would invite comparing a peak in one with a trough in the other.
    defaults = load_defaults()["diversity"]
    window = int(params.get("window", defaults["window_bp"]))
    step = int(params.get("step", defaults["step_bp"]))
    threshold_k = float(params.get("threshold_sd", DEFAULT_THRESHOLD_SD))
    min_windows = int(params.get("min_windows", DEFAULT_MIN_WINDOWS))

    alignment = align(
        [(genome.accession or genome.name, genome.sequence) for genome in genomes], ctx
    )
    positions = reference_positions(alignment)
    series = window_diversity(alignment, window=window, step=step)

    # Where each window sits on the reference (its first and last non-gap reference position).
    spans: list[tuple[int | None, int | None]] = []
    for entry in series:
        touched = [
            position
            for position in positions[entry.start_1b - 1 : entry.end_1b]
            if position is not None
        ]
        spans.append((touched[0] + 1, touched[-1] + 1) if touched else (None, None))

    threshold = threshold_for(series, k=threshold_k)
    regions = (
        hypervariable_regions(series, threshold=threshold, min_windows=min_windows)
        if threshold is not None
        else ()
    )
    region_of: dict[int, str] = {}
    for index, region in enumerate(regions, start=1):
        for entry in series:
            if region.start_1b <= entry.start_1b and entry.end_1b <= region.end_1b:
                region_of.setdefault(entry.start_1b, f"R{index}")

    ranked = sorted(
        (entry for entry in series if entry.pi is not None),
        key=lambda entry: (-entry.pi, entry.start_1b),
    )
    ranks = {entry.start_1b: rank for rank, entry in enumerate(ranked, start=1)}

    window_rows = [
        {
            "window": index,
            "alignment_start": entry.start_1b,
            "alignment_end": entry.end_1b,
            "reference_start": spans[index - 1][0],
            "reference_end": spans[index - 1][1],
            "usable_columns": entry.usable_columns,
            "pi": entry.pi,
            "rank": ranks.get(entry.start_1b),
            "region": region_of.get(entry.start_1b, ""),
        }
        for index, entry in enumerate(series, start=1)
    ]

    region_rows = []
    for index, region in enumerate(sorted(regions, key=lambda item: -item.mean_pi), start=1):
        first = next(
            span
            for entry, span in zip(series, spans, strict=True)
            if entry.start_1b == region.start_1b
        )
        last = next(
            span for entry, span in zip(series, spans, strict=True) if entry.end_1b == region.end_1b
        )
        region_rows.append(
            {
                "region": f"R{index}",
                "rank": index,
                "reference_start": first[0],
                "reference_end": last[1],
                "span_bp": (
                    last[1] - first[0] + 1 if first[0] is not None and last[1] is not None else None
                ),
                "n_windows": region.n_windows,
                "mean_pi": region.mean_pi,
                "max_pi": region.max_pi,
            }
        )

    candidate_rows = _candidates(reference, alignment)

    values = [entry.pi for entry in series if entry.pi is not None]
    summary = Table(
        id="barcodes_summary",
        title="Barcode and variability summary",
        columns=(
            Column("accession", "Reference accession"),
            Column("organism", "Reference organism"),
            Column("n_genomes", "Genomes aligned", dtype="int"),
            Column("window_bp", "Window", dtype="int", unit="bp"),
            Column("step_bp", "Step", dtype="int", unit="bp"),
            Column("n_windows", "Windows with a value", dtype="int"),
            Column("mean_pi", "Mean π", dtype="float", precision=6),
            Column("threshold_pi", "Region threshold", dtype="float", precision=6),
            Column("threshold_rule", "Threshold rule"),
            Column("n_regions", "Hypervariable regions", dtype="int"),
            Column("n_candidates", "Candidate loci", dtype="int"),
            Column("n_found", "Candidates found", dtype="int"),
            Column("tool", "Aligner"),
        ),
        rows=(
            {
                "accession": reference_name,
                "organism": reference.name,
                "n_genomes": alignment.n_rows,
                "window_bp": window,
                "step_bp": step,
                "n_windows": len(values),
                "mean_pi": (sum(values) / len(values)) if values else None,
                "threshold_pi": threshold,
                "threshold_rule": f"mean + {threshold_k:g} sd of this run's windows",
                "n_regions": len(regions),
                "n_candidates": len(candidate_rows),
                "n_found": sum(1 for row in candidate_rows if row["found"] == "yes"),
                "tool": "MAFFT",
            },
        ),
        scope=Scope(
            label=(
                "the reference genome's normalized frame; π is the same quantity the diversity "
                "module reports, computed by the same function over the same alignment"
            ),
            coord_system="normalized",
        ),
    )

    windows_table = Table(
        id="barcodes_windows",
        title="π by window, ranked",
        columns=(
            Column("window", "Window", dtype="int"),
            Column("alignment_start", "Alignment start", dtype="int"),
            Column("alignment_end", "Alignment end", dtype="int"),
            Column("reference_start", "Reference start", dtype="int", unit="bp"),
            Column("reference_end", "Reference end", dtype="int", unit="bp"),
            Column("usable_columns", "Usable columns", dtype="int"),
            Column("pi", "π", dtype="float", precision=6),
            Column("rank", "Rank by π", dtype="int"),
            Column("region", "Region"),
        ),
        rows=tuple(window_rows),
        scope=Scope(
            label=(
                f"{window} bp window / {step} bp step over alignment columns; rank 1 is the most "
                f"variable window; reference coordinates follow the reference genome; a window "
                f"with no usable column has no π and no rank"
            ),
            coord_system="normalized",
        ),
    )

    regions_table = Table(
        id="barcodes_regions",
        title="Hypervariable regions",
        columns=(
            Column("region", "Region"),
            Column("rank", "Rank by mean π", dtype="int"),
            Column("reference_start", "Reference start", dtype="int", unit="bp"),
            Column("reference_end", "Reference end", dtype="int", unit="bp"),
            Column("span_bp", "Span", dtype="int", unit="bp"),
            Column("n_windows", "Windows", dtype="int"),
            Column("mean_pi", "Mean π", dtype="float", precision=6),
            Column("max_pi", "Max π", dtype="float", precision=6),
        ),
        rows=tuple(region_rows),
        scope=Scope(
            label=(
                f"maximal runs of at least {min_windows} consecutive windows above "
                f"{threshold_k:g} sd of this run's own mean π, found at {window}/{step} bp; a "
                f"window with no usable column breaks a run rather than joining it"
            ),
            coord_system="normalized",
        ),
        key="region",
    )

    candidates_table = Table(
        id="barcodes_candidates",
        title="Candidate barcode loci",
        columns=(
            Column("locus", "Locus"),
            Column("kind", "Kind"),
            Column("members", "Members"),
            Column("found", "Found"),
            Column("reference_start", "Reference start", dtype="int", unit="bp"),
            Column("reference_end", "Reference end", dtype="int", unit="bp"),
            Column("span_bp", "Span", dtype="int", unit="bp"),
            Column("pi", "π", dtype="float", precision=6),
            Column("rank", "Rank by π", dtype="int"),
            Column("note", "Note"),
        ),
        rows=tuple(candidate_rows),
        scope=Scope(
            label=(
                "the loci in data/barcode_loci.yaml, each checked against the reference "
                "genome's own annotation; π is over the locus's aligned columns"
            ),
            coord_system="normalized",
        ),
        key="locus",
    )

    not_found = [row["locus"] for row in candidate_rows if row["found"] == "no"]
    notes = [
        "π is computed by the core function the diversity module uses, over the same alignment "
        "and the same window and step — a window's value here is that window's value there.",
        "A barcode locus the reference does not carry keeps its row and says so"
        + (f": {', '.join(not_found)}." if not_found else "."),
        "Each candidate locus carries the source of its claim in data/barcode_loci.yaml; a locus "
        "not in that file is not a barcode this toolkit knows about.",
    ]
    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            "barcodes_summary": summary,
            "barcodes_windows": windows_table,
            "barcodes_regions": regions_table,
            "barcodes_candidates": candidates_table,
        },
        notes=tuple(notes),
        citations=("mafft", "cbol"),
    )


def _candidates(reference: Genome, alignment: Alignment) -> list[dict[str, Any]]:
    """Every locus of the shipped list, found or not, with the π of the columns it covers."""
    names = normalize_genome(reference, aliases_for(reference))
    copies: dict[str, list[Gene]] = {}
    features: list[Gene] = []
    for gene in reference.genes:
        if gene.kind not in GENE_KINDS:
            continue
        features.append(gene)
        copies.setdefault(names.normalized(gene).canonical, []).append(gene)

    rows: list[dict[str, Any]] = []
    for locus in load_barcode_loci():
        row: dict[str, Any] = {
            "locus": locus.name,
            "kind": locus.kind,
            "members": ", ".join(locus.members),
            "found": "no",
            "reference_start": None,
            "reference_end": None,
            "span_bp": None,
            "pi": None,
            "rank": None,
            "note": "",
        }
        span = _locate(locus, copies, features, reference)
        if span is None:
            missing = [member for member in locus.members if not copies.get(member)]
            if locus.kind == "nuclear":
                row["note"] = "not a plastid locus — a plastome does not carry it"
            elif missing:
                row["note"] = f"not carried by this record: {', '.join(missing)}"
            else:
                # Both members are here but nothing joins them: say that, rather than "absent".
                # Measured on the golden genome: trnH-GUG and psbA are split by rps19 at the
                # IRb/LSC junction, so the classic spacer does not exist as a spacer there.
                row["note"] = (
                    "present but not adjacent in this record's annotation — another feature "
                    "lies between the two members"
                )
            rows.append(row)
            continue
        row["found"] = "yes"
        row["reference_start"] = span.start
        row["reference_end"] = span.end
        row["span_bp"] = span.span
        if span.wrapped:
            (first_start, first_end), (second_start, second_end) = span.parts
            row["note"] = (
                f"spans the origin: {first_start}–{first_end} and {second_start}–{second_end}"
            )
        columns: list[int] = []
        for part_start, part_end in span.parts:
            columns.extend(columns_for_span(alignment, 0, (part_start - 1, part_end)))
        pi, usable = nucleotide_diversity(alignment, columns)
        row["pi"] = pi
        if usable == 0:
            row["note"] = (row["note"] + "; " if row["note"] else "") + (
                "no usable column falls inside this locus's alignment"
            )
        rows.append(row)

    ranked = sorted(
        (row for row in rows if row["pi"] is not None),
        key=lambda row: (-row["pi"], row["locus"]),
    )
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank
    return rows


@dataclass(frozen=True, slots=True)
class _LocusSpan:
    """Where a locus sits: one 1-based inclusive part, or two when it crosses the origin.

    A plastome is circular and its record starts wherever the submitter started it, so a locus
    may run through position 1 — measured on the golden genome: its ``trnH-psbA`` spacer is
    ``153273–153724`` plus ``1–71``, because the record begins inside that spacer. Reporting it
    as "not carried" (which the first version of this code did) is the exact failure this
    project exists to avoid.
    """

    parts: tuple[tuple[int, int], ...]

    @property
    def start(self) -> int:
        return self.parts[0][0]

    @property
    def end(self) -> int:
        return self.parts[-1][1]

    wrapped: bool = False  # a second part is an intron or a spacer across position 1 — only the
    #                        caller knows which, so it says so rather than the count guessing

    @property
    def span(self) -> int:
        return sum(part_end - part_start + 1 for part_start, part_end in self.parts)


def _locate(
    locus: Any, copies: dict[str, list[Gene]], features: Sequence[Gene], reference: Genome
) -> _LocusSpan | None:
    """The span of a locus in the reference's normalized frame, or ``None`` when it is absent.

    A gene locus is the first copy of that name (ycf1's truncated second copy is why "first" is
    said out loud); its span is its **parts**, so a feature the normalization split at the origin
    keeps both, and π is taken over the union of their columns — the same reading the diversity
    module gives its gene rows.

    An intergenic locus is the spacer between its two members: it counts when *some* pair of
    copies is adjacent in this record's own annotation with nothing between them — along the
    **circle** when the record is circular, which is how the golden genome's ``trnH-psbA`` spacer
    turns out to exist (it crosses position 1). An overlap is not a spacer.
    """
    if locus.kind == "nuclear":
        return None
    length = len(reference)
    if locus.kind == "gene":
        found = copies.get(locus.members[0])
        if not found:
            return None
        parts = tuple((part.start + 1, part.end) for part in found[0].parts)
        wrapped = any(part.start == 0 for part in found[0].parts) and any(
            part.end == length for part in found[0].parts
        )
        return _LocusSpan(parts, wrapped=wrapped)

    circular = reference.topology == "circular"
    upstreams, downstreams = (copies.get(member, []) for member in locus.members)
    for upstream in upstreams:
        for downstream in downstreams:
            direct = _spacer(upstream, downstream, features)
            if direct is not None:
                return direct
            if not circular:
                continue
            for first, second in ((upstream, downstream), (downstream, upstream)):
                wrapped = _wrapped_spacer(first, second, features, length)
                if wrapped is not None:
                    return wrapped
    return None


def _spacer(upstream: Gene, downstream: Gene, features: Sequence[Gene]) -> _LocusSpan | None:
    """The gap between two features when it really is a gap — start < end and nothing inside."""
    start, end = upstream.span[1] + 1, downstream.span[0] - 1
    if start > end:
        return None
    if any(
        gene is not upstream
        and gene is not downstream
        and start <= gene.span[0]
        and gene.span[1] <= end
        for gene in features
    ):
        return None
    return _LocusSpan(((start, end),))


def _wrapped_spacer(
    upstream: Gene, downstream: Gene, features: Sequence[Gene], length: int
) -> _LocusSpan | None:
    """The gap that crosses position 1: ``upstream`` late in the record, ``downstream`` early."""
    if upstream.span[1] <= downstream.span[0]:
        return None  # that ordering is the direct case, not the wrapped one
    tail_start = upstream.span[1] + 1
    head_end = downstream.span[0] - 1
    # nothing may sit after the upstream copy's end, nor before the downstream copy's start
    if any(
        gene is not upstream and gene is not downstream and gene.span[0] >= tail_start
        for gene in features
    ):
        return None
    if any(
        gene is not upstream and gene is not downstream and gene.span[1] <= head_end
        for gene in features
    ):
        return None
    parts: list[tuple[int, int]] = []
    if tail_start <= length:
        parts.append((tail_start, length))
    if head_end >= 1:
        parts.append((1, head_end))
    return _LocusSpan(tuple(parts), wrapped=True) if parts else None
