"""The pipeline diagram has one canonical copy, and ROADMAP's coverage tables cite real milestones.

`PROJECT.md` §13.0 is the canonical diagram; `ROADMAP.md` ("The pipeline") and
`docs/D5-09_reference_naming.md` §1 keep copies so each document reads on its own. Both copies
state the rule — *edit the canonical one first, then copy it back* — and this test is what makes
the rule checkable: a copy that drifts silently is a reader being handed a different pipeline.

The other half of ROADMAP's "How to use" rule 5 — every box, and every tool inside a box, names
the milestone that owns it — cannot be tested, because ownership is a human judgement. What *is*
mechanical is the table's own integrity: a milestone id written in a table must exist as a
milestone in the file. Measured 07/10/2026: the first pass ran at box level and missed three
unowned tools (`OATK`/`raven`, `BWA`/`minimap2`, `Sibelia`) and one team artifact (the `.sbt`
template) — so the rule now runs at tool level, and this guard keeps the tables honest about the
ids they cite.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).parents[2]
CANONICAL = ROOT / "PROJECT.md"
COPIES = (ROOT / "ROADMAP.md", ROOT / "docs" / "D5-09_reference_naming.md")

MARKER = "NGUỒN THAM CHIẾU"  # the diagram's first line — one fenced block per file
DEFINED = re.compile(r"^- \[[ x~]\] \*\*(M[0-9A-Za-z]+\.\d+)", re.MULTILINE)
CITED = re.compile(r"M[0-9A-Za-z]+\.\d+")


def _diagram(text: str, path: Path) -> str:
    """The fenced block holding the pipeline diagram; exactly one per file."""
    hits = [part for part in text.split("```") if MARKER in part]
    assert len(hits) == 1, f"{path}: {len(hits)} blocks contain {MARKER!r}, expected 1"
    return hits[0]


def test_the_three_diagram_copies_are_identical() -> None:
    canonical = _diagram(CANONICAL.read_text(), CANONICAL)
    for copy in COPIES:
        assert _diagram(copy.read_text(), copy) == canonical, (
            f"{copy} has drifted from the canonical diagram (PROJECT.md §13.0) — "
            "edit the canonical copy first, then copy it back"
        )


def test_every_milestone_a_table_cites_exists() -> None:
    roadmap = (ROOT / "ROADMAP.md").read_text()
    defined = set(DEFINED.findall(roadmap))
    assert defined, "no milestone definitions found in ROADMAP.md"
    dangling = sorted(
        {
            cited
            for line in roadmap.splitlines()
            if line.startswith("|")
            for cited in CITED.findall(line)
        }
        - defined
    )
    assert not dangling, (
        f"ROADMAP.md tables cite milestones that do not exist: {', '.join(dangling)} — "
        "a coverage row pointing at nothing is the silent gap rule 5 exists to catch"
    )
