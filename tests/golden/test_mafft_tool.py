"""The MAFFT adapter against the pinned build — golden, tier B.

Every strategy name this toolkit documents must be an argv the pinned MAFFT accepts. The names
lived in a dict and a unit test checked their shape, so `--fftns` and `--nwns` stayed in the table
while MAFFT 7 had removed both: `--set mafft_strategy=fftns` reached the binary and exited 1 with
`Unknown option: --fftns` (measured 28/10/2026). Only a call to the real tool can answer this, so
the check lives here, on a four-sequence input that costs milliseconds.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import require_tools
from plastid.external.mafft import MAFFT_STRATEGIES, MAFFT_VERSION, align
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

pytestmark = pytest.mark.tier_b

ROWS = (
    ("a", "ATGGC TTCAA ACCCG GGTTT TAAGC TTGAC CGGAT CCGTA AGCCT TGGCA".replace(" ", "")),
    ("b", "ATGGCTTCAAACCCGGGTTTTAAGCTTGACCGGATCCGTAAGCCTTGGCA"),
    ("c", "ATGGCTTCAAGCCCGGGCTTTAAGCTTGACCGGATCCGTAAGCTTTGGCG"),
    ("d", "ATGGCCTCAAGCCCAGGTTTTAAGCTTGACCGGTTCCGTAAGCCTTGGCG"),
)


def _ctx(root: Path) -> RunContext:
    layout = RunLayout.create(root, "mafft", fingerprint="strategy")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="mafft_check", module_order=1)


def test_the_pinned_build_is_the_one_the_table_was_written_for(
    tools: dict[str, str | None],
) -> None:
    """The banner MAFFT prints is `v7.526 (2024/Apr/26)`; the pin is the `7.526` in it."""
    require_tools(tools, "mafft")
    from plastid.external.mafft import tool_version

    banner, _source = tool_version("mafft", "binary", tools["mafft"])
    assert banner is not None and MAFFT_VERSION in banner, banner


@pytest.mark.parametrize("strategy", sorted(MAFFT_STRATEGIES))
def test_every_documented_strategy_runs(
    strategy: str, tools: dict[str, str | None], tmp_path: Path
) -> None:
    """A name in the table is a promise that the pinned build accepts its argv."""
    require_tools(tools, "mafft")
    alignment = align(ROWS, _ctx(tmp_path / strategy), strategy=strategy)
    assert alignment.names == tuple(name for name, _ in ROWS)
    assert len({len(sequence) for sequence in alignment.sequences}) == 1
    assert set("".join(alignment.sequences)) <= set("ACGTN-")
