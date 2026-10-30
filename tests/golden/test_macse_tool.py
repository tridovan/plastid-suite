"""The MACSE adapter against the real jar (M11.21) — golden, tier B.

One gene, not a supermatrix: 78 MACSE calls are ~23 minutes (measured 28/10/2026), and the suite
needs a wiring check, not a second scale measurement. The module's end-to-end MACSE path is
exercised by the stubbed test in `tests/unit/test_phylogeny_module.py`; this file is what proves
the jar itself runs and returns the alphabet the adapter expects.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from plastid.external.macse import MacseBinary, align_cds, resolve_macse
from plastid.run.tools import ToolError

pytestmark = [pytest.mark.tier_b, pytest.mark.panel]

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"
RECORDS = ("NC_008641.1", "NC_016708.2", "NC_028027.1", "NC_016736.1")


def _real_binary() -> MacseBinary:
    """The jar named by ``PTS_MACSE_JAR``, else the installed package's ``macse`` on ``PATH``.

    No path of one developer's machine is baked in: CI installs the pinned package and this
    resolves it through ``PATH``, while a local tree that keeps only the jar sets the variable.
    A probe of an unreadable candidate used to raise ``PermissionError`` from ``Path.is_file``
    instead of skipping, which is why the resolution is the adapter's, not the test's.
    """
    explicit = os.environ.get("PTS_MACSE_JAR")
    try:
        return resolve_macse(explicit) if explicit else resolve_macse()
    except ToolError as error:
        reason = f"MACSE is not installed: {error}"
        if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
            pytest.fail(reason)
        pytest.skip(reason)


def test_the_real_jar_aligns_one_gene_across_four_records(tmp_path: Path) -> None:
    """One gene, real jar: equal-length rows, ACGT/-, and the input order preserved."""
    binary = _real_binary()
    if not PANEL.is_dir():
        pytest.skip("the IR panel is not in this tree (.panel/ is gitignored)")
    from plastid.core.genbank import read_genbank
    from plastid.core.gene_names import cds_by_canonical_name
    from plastid.core.model import coding_sequence
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    rows = []
    for accession in RECORDS:
        genome = read_genbank(PANEL / f"{accession}.gb")
        rows.append((accession, coding_sequence(genome, cds_by_canonical_name(genome)["rbcL"])))

    layout = RunLayout.create(tmp_path, "macse", fingerprint="test")
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="phylogeny", module_order=1)
    alignment = align_cds(rows, ctx, binary=binary)

    assert alignment.names == tuple(name for name, _ in rows)
    assert len({len(sequence) for sequence in alignment.sequences}) == 1
    assert set("".join(alignment.sequences)) <= set("ACGTN-")
