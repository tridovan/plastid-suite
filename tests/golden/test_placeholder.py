"""Placeholder for the golden-number tests — the harness arrives in P1-06.

Kept as a loud skip instead of an empty directory: `pytest` exits 5 (failure)
when it collects nothing, and a silently empty golden suite is exactly the
kind of quiet hole this project forbids.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.skip(reason="golden harness arrives in P1-06 (ROADMAP.md M1.6)")


def test_golden_numbers() -> None:
    # When the harness lands, this asserts the §6 numbers:
    # 153,724 bp · GC 37.85 % · 132 = 86 CDS + 38 tRNA + 8 rRNA
    # IR 26,214 bp @ 100.00 % · LSC 83,127 · SSC 18,169 · 37 SSR (31/4/2).
    raise AssertionError("golden harness not built yet")
