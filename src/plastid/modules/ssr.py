"""SSR (microsatellite) analysis — tier B: needs the bundled MISA + perl.

Golden numbers (PROJECT.md §6): 37 = mono 31 · di 4 · tri 2, compound 3.
Counts are MISA's; the module adds the plastid region (LSC/IRa/SSC/IRb) and the gene
context (Gene/Intron/IGS) for every SSR, both scope-labelled.
"""

from __future__ import annotations

from typing import Any

from plastid.core.ir import IRResult
from plastid.core.model import Genome, Interval
from plastid.core.organelle import GROUP_SHARED
from plastid.external.misa import UNIT_LABELS, MisaResult, run_misa
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "ssr"
# 0.3.0 (10/10/2026): the compound flag now marks every component of a cluster. MISA writes a
# cluster two ways and the parser read only one, so a component on its own line said "no" —
# 0 of PX413282's 16 components while MISA counted 8 clusters (P8B-01, found by the Crataegus
# answer-key corpus). The run cache keys on this version, not on the source (H-08).
# 0.2.0 (10/10/2026): the Context column changed — a trans-spliced gene's outer span no longer
# makes everything inside it an "Intron" (§9.M-1, found by the CPStools cross-check §8.61).
MODULE_VERSION = "0.3.0"

_THRESHOLDS = "(1/10) (2/6) (3/5) (4/5) (5/5) (6/5)"  # MISA defaults, PROJECT.md §10.2
_GENE_KINDS = ("CDS", "tRNA", "rRNA")


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="SSR (microsatellites)",
    requires_ir=True,
    citations=("misa",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    misa = run_misa(genome, ctx)

    rows = []
    for feature in misa.features:
        span = Interval(feature.start - 1, feature.end)
        rows.append(
            {
                "motif": feature.motif,
                "unit": feature.unit_label,
                "start_1b": feature.start,
                "end_1b": feature.end,
                "length_bp": span.length,
                "region": _plastid_region(span, ir, len(genome)),
                "context": _gene_context(span, genome, ir),
                "compound": "yes" if feature.compound else "no",
            }
        )

    summary = Table(
        id="ssr_summary",
        title="SSR summary",
        columns=(
            Column("accession", "Accession"),
            Column("total", "Total SSRs", dtype="int"),
            Column("mono", "Monomeric", dtype="int"),
            Column("di", "Dinucleotide", dtype="int"),
            Column("tri", "Trinucleotide", dtype="int"),
            Column("compound", "In compound formation", dtype="int"),
            Column("tool", "Tool"),
            Column("thresholds", "Unit/repeat thresholds"),
            Column("reconciliation", "GFF vs statistics"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "total": misa.total,
                "mono": misa.by_unit.get(1, 0),
                "di": misa.by_unit.get(2, 0),
                "tri": misa.by_unit.get(3, 0),
                "compound": misa.compound,
                "tool": f"MISA {_misa_version(ctx)}",
                "thresholds": _THRESHOLDS,
                "reconciliation": _reconciliation(misa, rows),
            },
        ),
        scope=Scope(label="whole plastome; MISA defaults; both IR copies included"),
    )

    by_unit = Table(
        id="ssr_by_unit",
        title="SSRs by repeat unit",
        columns=(
            Column("unit", "Unit"),
            Column("unit_size", "Unit size", dtype="int", unit="bp"),
            Column("n", "Count", dtype="int"),
        ),
        rows=tuple(
            {"unit": feature_label(size), "unit_size": size, "n": count}
            for size, count in sorted(misa.by_unit.items())
        ),
        scope=Scope(label="whole plastome; counts from MISA .statistics"),
        key="unit",
    )

    features_table = Table(
        id="ssrs",
        title="SSRs (positions from the MISA GFF)",
        columns=(
            Column("motif", "Motif"),
            Column("unit", "Unit"),
            Column("start_1b", "Start", dtype="int"),
            Column("end_1b", "End", dtype="int"),
            Column("length_bp", "Length", dtype="int", unit="bp"),
            Column("region", "Region"),
            Column("context", "Context"),
            Column("compound", "In compound formation"),
        ),
        rows=tuple(rows),
        scope=Scope(
            label=(
                "GenBank coordinates; region split from the IR pair; a compound cluster's "
                "components share the cluster's span"
            ),
            coord_system="genbank",
        ),
    )

    region_table = Table(
        id="ssr_by_region",
        title="SSRs by plastid region",
        columns=(
            Column("region", "Region"),
            Column("total", "Total", dtype="int"),
            Column("mono", "Monomeric", dtype="int"),
            Column("di", "Dinucleotide", dtype="int"),
            Column("tri", "Trinucleotide", dtype="int"),
            Column("other", "Other unit sizes", dtype="int"),
        ),
        rows=tuple(_region_row(region, rows) for region in ("LSC", "IRa", "SSC", "IRb")),
        scope=Scope(
            label="whole plastome; counts from the per-SSR table, so the two always reconcile"
        ),
        key="region",
    )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables={
            summary.id: summary,
            by_unit.id: by_unit,
            region_table.id: region_table,
            features_table.id: features_table,
        },
        notes=[
            f"MISA {_misa_version(ctx)} with its default thresholds {_THRESHOLDS}; "
            f"{misa.compound} SSRs in compound formation — MISA's own count, **one per cluster**.",
            "The per-SSR table flags **every component** of a cluster, so its `yes` rows number "
            "more than MISA's count (a two-component cluster contributes 2 rows to 1): the two "
            "numbers measure different things and neither is derived from the other.",
            "A compound cluster's components share the cluster's coordinates — MISA does not "
            "separate them, and a component on its own line carries the `compound_repeat` note.",
        ],
        citations=("misa",),
    )


def _region_row(region: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """One region's counts, taken from the per-SSR rows so the summary cannot drift from them."""
    in_region = [row for row in rows if row["region"] == region]
    return {
        "region": region,
        "total": len(in_region),
        "mono": sum(1 for row in in_region if row["unit"] == "mono"),
        "di": sum(1 for row in in_region if row["unit"] == "di"),
        "tri": sum(1 for row in in_region if row["unit"] == "tri"),
        "other": sum(1 for row in in_region if row["unit"] not in ("mono", "di", "tri")),
    }


def _reconciliation(misa: MisaResult, rows: list[dict[str, Any]]) -> str:
    """MISA's own statistics against the SSRs actually parsed out of its GFF.

    They disagreed once (2026-10-03): the parser read only ``monomeric_repeat`` lines and
    silently lost every longer unit, leaving 30 of 37 SSRs. This column exists so that can
    never pass unnoticed again.
    """
    problems: list[str] = []
    if len(rows) != misa.total:
        problems.append(f"GFF has {len(rows)} SSRs, statistics say {misa.total}")
    for size, count in sorted(misa.by_unit.items()):
        label = UNIT_LABELS.get(size, f"{size}-unit")
        parsed = sum(1 for row in rows if row["unit"] == label)
        if parsed != count:
            problems.append(f"{label}: GFF {parsed} vs statistics {count}")
    if problems:
        return "DISAGREES — " + "; ".join(problems)
    return f"agrees — {len(rows)} SSRs, every unit class matches"


def _misa_version(ctx: RunContext) -> str:
    for record in ctx.tools:
        if record.name == "misa" and record.version:
            return record.version
    return "unknown"


def feature_label(unit_size: int) -> str:
    return UNIT_LABELS.get(unit_size, f"{unit_size}-unit")


def _plastid_region(span: Interval, ir: IRResult, genome_length: int) -> str:
    for name, copy in (("IRa", ir.ira), ("IRb", ir.irb)):
        if span.overlaps(copy):
            return name
    first, second = sorted((ir.ira, ir.irb), key=lambda interval: interval.start)
    between = Interval(first.end, second.start)
    between_length = between.length
    wrapping_length = (genome_length - second.end) + first.start
    if between_length >= wrapping_length:
        return "LSC" if span.overlaps(between) else "SSC"
    return "SSC" if span.overlaps(between) else "LSC"


def _gene_context(span: Interval, genome: Genome, ir: IRResult) -> str:
    """Gene · Intron · IGS for one span, by the rule ``core/spans.py`` already states.

    A gene occupies its **exons**, never the space between them. Read by outer span instead, the
    golden's trans-spliced ``rps12`` — exon 1 in the LSC, exons 2 and 3 in the IR, all three
    listed on each CDS feature — is one gene covering 68,975–139,634, so every SSR in between
    came out "Intron" whichever gene it actually sits in. Measured 10/10/2026 in the CPStools
    cross-check (§8.61): 11 of the golden's 37 rows, and cpstools' table named the genes they sit
    in (``ycf1`` · ``ndhA`` · the ``rpl32`` spacer).

    So a gap between two parts of one gene is an intron only when **both ends sit in one region**
    — the same rule ``introns_of`` splits introns from trans-splicing junctions by. Exon
    containment is checked across every gene before a gap is believed, because one gene's intron
    can carry another gene (``matK`` is annotated inside ``trnK``'s), and that inner gene's exon
    is a "Gene", not an intron.
    """
    length = len(genome)
    in_intron = False
    for gene in genome.genes:
        if gene.kind not in _GENE_KINDS:
            continue
        if any(part.start <= span.start and span.end <= part.end for part in gene.parts):
            return "Gene"
        parts = sorted(gene.parts, key=lambda part: part.start)
        for first, second in zip(parts, parts[1:], strict=False):
            if not (first.end <= span.start and span.end <= second.start):
                continue
            head = _plastid_region(Interval(first.end, first.end + 1), ir, length)
            tail = _plastid_region(Interval(second.start - 1, second.start), ir, length)
            if head == tail:
                in_intron = True
    return "Intron" if in_intron else "IGS"


__all__ = ["MODULE_ID", "MODULE_VERSION", "MisaResult", "run"]
