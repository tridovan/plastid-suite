"""H-22 — ARAGORN's intron range is declared, and a locus beyond it is never called intron-free.

Measured 26/10/2026 on ARAGORN 1.2.41. The golden's ``trnK-UUU`` intron is lengthened in steps
inside its own locus: 3,000 bp is reported, 3,001 bp is not, and ``-i0,4000``/``-i4000``/
``-i0,6000`` do not change that — the documented maximum is an internal cap, not a default that
the switch lifts. A caller that asks for a wider search is refused by name, and every result
carries the range it actually searched, so "no intron" and "not assessed" cannot be the same row.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core.model import Genome
from plastid.external.aragorn import ARAGORN_INTRON_LIMITS, run_aragorn
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

pytestmark = pytest.mark.tier_b

# 1-based trnK-UUU: exons 1367–1402 and 4020–4057, intron 1403–4019 (2,618 bp). The locus is
# searched without the rest of the genome; the tRNA structure is what ARAGORN needs.
LOCUS_START_1B = 1_200
LOCUS_END_1B = 4_260
INTRON_MID = ((1_403 - LOCUS_START_1B) + (4_019 - LOCUS_START_1B)) // 2
BASE_INTRON = 2_618


def _context(tmp_path: Path, label: str) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", label, fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="test", module_order=1)


def _locus_with_intron(golden: Genome, extra: int) -> Genome:
    """The trnK-UUU locus with ``extra`` random bases inserted at the intron's midpoint."""
    locus = golden.sequence[LOCUS_START_1B:LOCUS_END_1B]
    random.seed(7)  # deterministic: the probe is about length, not about a lucky insertion
    insertion = "".join(random.choice("ACGT") for _ in range(extra))
    sequence = locus[:INTRON_MID] + insertion + locus[INTRON_MID:]
    return Genome(
        name=f"trnK_{BASE_INTRON + extra}",
        accession=f"trnK_{BASE_INTRON + extra}",
        sequence=sequence,
        topology="linear",
    )


@pytest.fixture(scope="module")
def aragorn(blast_tools: dict[str, str | None]) -> str:
    require_tools(blast_tools, "aragorn")
    assert blast_tools["aragorn"] is not None
    return blast_tools["aragorn"]


def test_the_declared_range_is_the_bound_the_caller_actually_applies() -> None:
    assert ARAGORN_INTRON_LIMITS == (0, 3_000)


def test_an_intron_at_the_maximum_is_found(
    aragorn: str, golden_genome: Genome, tmp_path: Path
) -> None:
    genome = _locus_with_intron(golden_genome, 3_000 - BASE_INTRON)
    ctx = _context(tmp_path, "at_limit")

    result = run_aragorn(genome, ctx, seq_id=genome.accession or "", aragorn=aragorn, introns=True)

    assert result.intron_search == (0, 3_000)
    assert "-i" in ctx.tools[-1].argv
    introns = [hit.intron for hit in result.hits if hit.intron is not None]
    assert [length for _start, length in introns] == [3_000]


def test_one_base_past_the_maximum_is_not_reported(
    aragorn: str, golden_genome: Genome, tmp_path: Path
) -> None:
    """The failure mode H-22 names: the gene is simply absent, so the run must say what it searched.

    The result carries ``(0, 3_000)``; the run's notes (``plastid_annotate``) state that a tRNA
    outside the range is not compared as absent. The test cannot assert a note from here — the
    annotator owns it — but it can assert the evidence the note rests on.
    """
    genome = _locus_with_intron(golden_genome, 3_001 - BASE_INTRON)
    ctx = _context(tmp_path, "over_limit")

    result = run_aragorn(genome, ctx, seq_id=genome.accession or "", aragorn=aragorn, introns=True)

    assert result.intron_search == (0, 3_000)
    assert not [hit for hit in result.hits if hit.intron is not None]


def test_a_wider_request_is_refused_rather_than_silently_ignored(
    aragorn: str, golden_genome: Genome, tmp_path: Path
) -> None:
    """ARAGORN 1.2.41 accepts ``-i0,4000`` and still stops at 3,000; passing it would be a lie."""
    genome = _locus_with_intron(golden_genome, 3_001 - BASE_INTRON)
    ctx = _context(tmp_path, "refused_limit")

    with pytest.raises(ValueError, match=r"does not search a tRNA intron longer than 3,000 bp"):
        run_aragorn(
            genome,
            ctx,
            seq_id=genome.accession or "",
            aragorn=aragorn,
            introns=True,
            intron_limits=(0, 4_000),
        )
