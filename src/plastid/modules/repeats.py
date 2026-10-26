"""Repeats — maximal exact matches at or above a length floor (M2.12, D-f, H-16).

The engine is MUMmer's ``repeat-match``, and the two class names are the tool's own: **dispersed**
(a direct repeat) and **palindromic** (an inverted one). Vmatch 2.3.1 is the independent second
engine the count is checked against — measured 10/10/2026, both give the golden's **9** matches at
30 bp or more, the same coordinates and the same two classes — the way MISA 2.1 and cpstools 3.0
both give its 37 SSRs (D-f).

Every match carries both copies' coordinates **and the plastid region of each copy**, the
inverted-repeat pair among them: on the golden the longest match *is* the IR pair, 26,214 bp
between IRa and IRb, and that is a test rather than something a reader checks by eye.

The count is measured per run and never inherited: the "9 matches ≥30 bp" the phase's board
carried came from the 26/9 tool survey, measured on other data. It happens to be the same
number here — recorded as a coincidence, not as a reason.

A **tandem** array is a different class and is not this module: ``repeat-match`` reports its
copies as overlapping maximal matches, not as one tandem row, and the phrase "tandem repeat" in
this toolkit means ``mreps``/``trf`` ground (M2.12's other half, still open).
"""

from __future__ import annotations

from typing import Any

from plastid.core.model import Genome
from plastid.core.normalize import NormalizedFrame, normalize
from plastid.core.organelle import GROUP_SHARED
from plastid.external.repeatmatch import MIN_LENGTH_BP, run_repeat_match
from plastid.modules.registry import module
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.context import RunContext

MODULE_ID = "repeats"
MODULE_VERSION = "0.1.0"

DEFAULT_MIN_LENGTH_BP = MIN_LENGTH_BP


@module(
    MODULE_ID,
    group=GROUP_SHARED,
    version=MODULE_VERSION,
    title="Repeats (exact)",
    requires_ir=True,
    citations=("mummer",),
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result:
    ir = ctx.require_ir()
    frame = normalize(genome, ir)
    min_length = int(params.get("min_length_bp", DEFAULT_MIN_LENGTH_BP))
    found = run_repeat_match(genome, ctx, min_length=min_length, directory=params.get("repeats"))

    # Longest first: a manuscript table reads that way, and the order is fixed here rather than
    # inherited from the tool, so two runs of two builds can be diffed.
    ranked = sorted(found.matches, key=lambda match: (-match.length, match.start1, match.start2))
    rows = []
    crossing = 0
    for rank, match in enumerate(ranked, start=1):
        start1, end1, start2, end2 = match.to_1based()
        row = {
            "repeat_id": f"R{rank:02d}",
            "class": "palindromic" if match.palindromic else "dispersed",
            "length_bp": match.length,
            "start1_1b": start1,
            "end1_1b": end1,
            "region1": _region_name(frame, match.start1),
            "start2_1b": start2,
            "end2_1b": end2,
            "region2": _region_name(frame, match.start2),
        }
        if _crosses(frame, match.start1, match.end1) or _crosses(frame, match.start2, match.end2):
            crossing += 1
        rows.append(row)

    summary = Table(
        id="repeat_summary",
        title="Exact repeat summary",
        columns=(
            Column("accession", "Accession"),
            Column("organism", "Organism"),
            Column("min_length_bp", "Minimum match length", dtype="int", unit="bp"),
            Column("n_matches", "Matches", dtype="int"),
            Column("n_dispersed", "Dispersed", dtype="int"),
            Column("n_palindromic", "Palindromic", dtype="int"),
            Column("longest_bp", "Longest match", dtype="int", unit="bp"),
            Column("longest_regions", "Longest match regions"),
            Column("tool", "Tool"),
        ),
        rows=(
            {
                "accession": genome.accession or "",
                "organism": genome.name,
                "min_length_bp": min_length,
                "n_matches": len(rows),
                "n_dispersed": sum(1 for row in rows if row["class"] == "dispersed"),
                "n_palindromic": sum(1 for row in rows if row["class"] == "palindromic"),
                "longest_bp": rows[0]["length_bp"] if rows else 0,
                "longest_regions": (f"{rows[0]['region1']}–{rows[0]['region2']}" if rows else ""),
                "tool": f"repeat-match (MUMmer {found.version})",
            },
        ),
        scope=Scope(
            label=(
                f"maximal **exact** matches (dispersed and palindromic — no mismatches, no "
                f"tandem class) at least {min_length} bp, both strands, both IR copies included"
            )
        ),
    )

    # A genome with no match gets no match table — Table refuses zero rows, and an empty table
    # could not be read if it existed (the rule reference_audit states). The summary still
    # reports the zero, so "found nothing" cannot be told apart from "did not run".
    tables = {"repeat_summary": summary}
    if rows:
        tables["repeats"] = Table(
            id="repeats",
            title="Exact repeats",
            columns=(
                Column("repeat_id", "Repeat"),
                Column("class", "Class"),
                Column("length_bp", "Length", dtype="int", unit="bp"),
                Column("start1_1b", "Copy 1 start", dtype="int"),
                Column("end1_1b", "Copy 1 end", dtype="int"),
                Column("region1", "Copy 1 region"),
                Column("start2_1b", "Copy 2 start", dtype="int"),
                Column("end2_1b", "Copy 2 end", dtype="int"),
                Column("region2", "Copy 2 region"),
            ),
            rows=tuple(rows),
            scope=Scope(
                label=(
                    "GenBank coordinates, 1-based inclusive; region = the plastid region of the "
                    "copy's first base; copy 1 is the upstream copy"
                )
            ),
            key="repeat_id",
        )

    notes = [
        f"{len(rows)} exact match(es) of at least {min_length} bp: "
        f"{summary.rows[0]['n_dispersed']} dispersed, "
        f"{summary.rows[0]['n_palindromic']} palindromic.",
        "Maximal exact matches only: a tandem array arrives as overlapping matches rather than "
        "one row, so the count is of matches, not of repeat loci.",
    ]
    if not rows:
        notes.append(
            f"no exact match of {min_length} bp or more was found — the table of matches is "
            f"absent rather than empty, and this note is what says the search ran."
        )
    if (
        rows
        and rows[0]["length_bp"] == ir.length
        and {rows[0]["region1"], rows[0]["region2"]} == {"IRa", "IRb"}
    ):
        notes.append(
            f"the longest match is the inverted-repeat pair itself: {ir.length} bp between "
            f"IRa and IRb."
        )
    if crossing:
        notes.append(
            f"{crossing} match(es) span a plastid-region boundary — the region columns give each "
            f"copy's first base, not the whole span."
        )

    return Result(
        module=MODULE_ID,
        module_version=MODULE_VERSION,
        tables=tables,
        notes=tuple(notes),
        citations=("mummer",),
    )


def _region_name(frame: NormalizedFrame, position: int) -> str:
    region = frame.region_at(position)
    return region.name if region else ""


def _crosses(frame: NormalizedFrame, start: int, end: int) -> bool:
    """True when a span's two ends sit in different plastid regions."""
    return _region_name(frame, start) != _region_name(frame, end - 1)
