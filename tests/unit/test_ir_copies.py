"""The inverted repeat as the export writes it — the copy rule, on a synthetic plastome.

`_ir_copies` is where the toolkit writes the same gene twice, and the rule has two halves that
were measured separately: a gene **inside** one side is mirrored into the other (a measured
failure caught the arithmetic being one repeat-length out, and no count could see it), and a
gene that **straddles** a boundary — trans-spliced, one exon in the single-copy region and one
inside the repeat — is copied too, with the shared exon written into both features (the
golden's `rps12`, which the published record writes twice).

The counter-case is `ycf1`: it runs from the SSC into the repeat as **one** call, and one call is
not a second copy — the published golden carries one `ycf1` and two `rps12`. A rule that copied
any gene touching the repeat would invent a second `ycf1` on every genome.
"""

from __future__ import annotations

from plastid.core.model import Gene, Interval
from plastid.modules.plastid_annotate import _ir_copies

# A repeat of 1,000 bp on each side, with 2,000 bp of single-copy region either side of it:
# LSC 0..1999 · IRa 2000..2999 · SSC 3000..4999 · IRb 5000..5999 (0-based, half-open).
IRA = Interval(2000, 3000)
IRB = Interval(5000, 6000)


class _Repeat:
    def __init__(self) -> None:
        self.ira = IRA
        self.irb = IRB


def gene(name: str, span: tuple[int, int], *more: tuple[int, int], strand: int = 1) -> Gene:
    """A gene from 1-based inclusive exons, the way a caller reports one."""
    return Gene(
        name=name,
        kind="CDS",
        parts=tuple(Interval(start - 1, end, strand) for start, end in (span, *more)),
    )


def copies(genes: tuple[Gene, ...], name: str) -> list[list[tuple[int, int, int]]]:
    return [
        [(part.start + 1, part.end, part.strand) for part in found.parts]
        for found in genes
        if found.name == name
    ]


def test_a_gene_inside_one_side_is_written_twice_at_the_mirrored_offset() -> None:
    """The inversion as arithmetic: 300 bp inside IRa from its start is 300 bp inside IRb."""
    written = _ir_copies([gene("ndhB", (2300, 2600), strand=1)], _Repeat())
    found = copies(written, "ndhB")

    # 299 bp from IRa's start is 299 bp from IRb's *end* — the inversion, not a translation
    assert found == [
        [(2300, 2600, 1)],
        [(5401, 5701, -1)],
    ]


def test_a_gene_outside_the_repeat_is_left_alone() -> None:
    written = _ir_copies([gene("psbA", (100, 400))], _Repeat())

    assert copies(written, "psbA") == [[(100, 400, 1)]]


def test_a_trans_spliced_gene_across_the_boundary_keeps_its_shared_exon_in_both_copies() -> None:
    """The golden's `rps12`: an exon in the single-copy region and one inside IRa. The record
    writes it as two features sharing that first exon, so the copy keeps it exactly as called."""
    # one call, two exons — the caller makes one call per locus
    written = _ir_copies([gene("rps12", (1500, 1600), (2400, 2700), strand=-1)], _Repeat())
    found = copies(written, "rps12")

    assert found == [
        [(1500, 1600, -1), (2400, 2700, -1)],
        [(1500, 1600, -1), (5301, 5601, 1)],  # the shared exon travels unchanged
    ]


def test_a_call_that_only_reaches_into_the_repeat_gets_no_second_copy() -> None:
    """`ycf1` crosses the SSC/IRb boundary as a single part. The published golden carries one
    `ycf1`, so a rule that copied every gene touching the repeat would invent one."""
    written = _ir_copies([gene("ycf1", (4900, 5300), strand=-1)], _Repeat())

    assert copies(written, "ycf1") == [[(4900, 5300, -1)]]


def test_two_sides_of_different_length_refuse_rather_than_guess() -> None:
    """A copy cannot be placed inside a repeat whose sides differ — the mirror would be off by
    the difference, which is the mirrored-from-the-wrong-end failure all over again."""

    class _Unequal(_Repeat):
        def __init__(self) -> None:
            self.ira = Interval(2000, 3000)
            self.irb = Interval(5000, 6100)

    import pytest

    with pytest.raises(ValueError, match="differ in length"):
        _ir_copies([gene("ndhB", (2300, 2600))], _Unequal())
