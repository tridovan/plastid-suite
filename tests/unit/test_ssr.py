"""The SSR context column — Gene · Intron · IGS — and the trans-spliced gene that broke it.

The golden's own record annotates ``rps12`` the way the biology is: exon 1 in the LSC, exons 2
and 3 in the IR, and **each** CDS feature lists all three. Read by its outer span, that is one
gene covering 68,975–139,634 — most of the genome — so every SSR between those two numbers came
out "Intron", whatever gene it actually sits in.

Measured 10/10/2026 in the CPStools cross-check (P8-02): **11 of the golden's 37 rows were
wrong**, and cpstools' table named the genes they sit in (``ycf1`` · ``ndhA`` · the ``rpl32``
spacer). The rule that keeps it right is the one ``core/spans.py`` already states — a gap between
two parts of one gene is an intron only when both ends sit in the **same region**; when it
crosses a region boundary it is a trans-splicing junction, and the sequence between the two exons
is not an intron of anything.
"""

from __future__ import annotations

from plastid.core.ir import detect_ir
from plastid.core.model import Genome, Interval
from plastid.modules.ssr import _gene_context


def _context(genome: Genome, start_1b: int, end_1b: int) -> str:
    """The module's own call, from the 1-based inclusive coordinates a record is written in."""
    return _gene_context(Interval(start_1b - 1, end_1b), genome, detect_ir(genome))


def test_an_ssr_inside_a_gene_is_gene_not_intron(golden_genome: Genome) -> None:
    """123,879 sits inside ``ycf1`` — the trans-spliced ``rps12`` span used to shadow it."""
    assert _context(golden_genome, 123879, 123889) == "Gene"


def test_an_ssr_inside_a_genes_intron_is_intron(golden_genome: Genome) -> None:
    """5,083 sits between ``rps16``'s two exons — a real intron, and it must stay one."""
    assert _context(golden_genome, 5083, 5096) == "Intron"


def test_an_ssr_between_the_trans_spliced_exons_is_igs(golden_genome: Genome) -> None:
    """82,711 is in the ``rps12`` trans-splicing junction: LSC sequence, not an intron."""
    assert _context(golden_genome, 82711, 82722) == "IGS"


def test_an_intron_whose_ends_share_one_ir_copy_still_counts(golden_genome: Genome) -> None:
    """97,300 is inside ``rps12``'s second intron — both ends in one IR copy, so it is an intron.

    The region rule must not over-correct: it rejects the junction that crosses LSC→IR, not every
    gap of a gene whose parts are spread across the genome.
    """
    assert _context(golden_genome, 97300, 97310) == "Intron"


def test_an_ssr_in_an_intergenic_spacer_is_igs(golden_genome: Genome) -> None:
    assert _context(golden_genome, 4562, 4572) == "IGS"
