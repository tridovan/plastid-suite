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

import random

from plastid.core.ir import measure_ir
from plastid.core.model import Gene, Genome, Interval, coding_sequence
from plastid.core.sequence import reverse_complement
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

    from plastid.core.ir import IRGateError

    # A named refusal (exit 2 at the CLI), not a bare ValueError: without an alignment there is
    # no evidence for where the other copy's bases sit.
    with pytest.raises(IRGateError, match="different length"):
        _ir_copies([gene("ndhB", (2300, 2600))], _Unequal())


# --- unequal copies: the mapping is the gate's own alignment (M11.11) -------------------------
#
# A 1,000 bp repeat pair where copy B carries three inserted bases in the middle. The insertion
# is inside the alignment, not at a terminal end, so the old mirror constant would pair every
# base after it with the wrong partner. The expected coordinates below are found by locating the
# gene's own sequence in the reverse complement of copy B — independent of the mapping code.

_IR = "".join(random.Random(20261026).choices("ACGT", k=1000))
_INSERT_AT = 500
_INSERTED = "GGG"
_COPY_B = reverse_complement(_IR)
_COPY_B = _COPY_B[:_INSERT_AT] + _INSERTED + _COPY_B[_INSERT_AT:]
_SEQUENCE = "ACGT" * 50 + _IR + "TTAA" * 50 + _COPY_B + "GGCC" * 50
_IRA_UNEQUAL = Interval(200, 1200)
_IRB_UNEQUAL = Interval(1400, 1400 + len(_COPY_B))
_IR_UNEQUAL = measure_ir(_SEQUENCE, _IRA_UNEQUAL, _IRB_UNEQUAL, source="annotation")


def _publish_copies() -> list[list[tuple[int, int, int]]]:
    written = _ir_copies(
        [gene("ndhB", (400, 500), strand=1)],
        _IR_UNEQUAL,
    )
    return [[(part.start + 1, part.end, part.strand) for part in found.parts] for found in written]


def test_an_unequal_pair_carries_an_alignment_mapping() -> None:
    """The gate measured this pair end to end, so the result carries the path, not a constant."""
    assert _IR_UNEQUAL.ira.length != _IR_UNEQUAL.irb.length
    assert _IR_UNEQUAL.mapping is not None
    assert _IR_UNEQUAL.mapping.method == "anchored"


def test_a_gene_in_the_shorter_copy_is_placed_by_the_alignment_not_the_constant() -> None:
    """The independent answer: the gene's own sequence, found in copy B's reverse complement.

    The insertion sits after the gene, so the mirror constant would put the copy three bases off
    (and then every later base on the wrong partner). Locating the sequence itself has no such
    assumption.
    """
    written = _publish_copies()
    found = written[1]

    piece = _SEQUENCE[399:500]  # the gene, as copy A carries it (1-based 400..500)
    opponent = reverse_complement(_SEQUENCE[_IRB_UNEQUAL.start : _IRB_UNEQUAL.end])
    position = opponent.index(piece)
    expected = (
        _IRB_UNEQUAL.end - (position + len(piece)) + 1,
        _IRB_UNEQUAL.end - position,
        -1,
    )
    assert found == [expected]


def test_a_gene_spanning_an_insertion_loses_exactly_the_inserted_bases() -> None:
    """A boundary gene is the case a length difference can actually land in.

    The insertion is in copy B; a gene there that spans it has no counterpart for those three
    bases, so its copy in A is three bases shorter and the rest of the gene still matches.
    """
    start = _IRB_UNEQUAL.start + _INSERT_AT - 50
    end = _IRB_UNEQUAL.start + _INSERT_AT + 50 + len(_INSERTED)
    source = Gene(name="rps12", kind="CDS", parts=(Interval(start, end, -1),))
    written = _ir_copies([source], _IR_UNEQUAL)
    assert len(written) == 2
    mapped = next(found for found in written if found.strand == 1)
    assert [part.length for part in mapped.parts] == [end - start - len(_INSERTED)]
    genome = Genome(name="x", sequence=_SEQUENCE)
    source_coding = coding_sequence(genome, source)
    mapped_coding = coding_sequence(genome, mapped)
    # The gene begins 50 bp before the insertion on the minus strand, so the insertion is coding
    # positions 51..53. The three bases are copy B's, in this gene's orientation.
    assert source_coding[50:53] == reverse_complement(_INSERTED)
    assert mapped_coding == source_coding[:50] + source_coding[53:]


def test_a_trans_spliced_gene_on_an_unequal_pair_keeps_its_shared_exon() -> None:
    """The golden's ``rps12`` shape, on a pair whose copies differ: only the IR exon moves."""
    # exon 1 in the single-copy region before IRa, exon 2 inside IRa
    source = Gene(
        name="rps12",
        kind="CDS",
        parts=(Interval(120, 199, -1), Interval(300, 400, -1)),
    )
    written = _ir_copies([source], _IR_UNEQUAL)
    assert len(written) == 2
    assert written[0].parts[0] == Interval(120, 199, -1)
    assert written[1].parts[0] == Interval(120, 199, -1), "the shared exon must travel unchanged"
    mirrored = written[1].parts[1]
    assert _IRB_UNEQUAL.start <= mirrored.start and mirrored.end <= _IRB_UNEQUAL.end
    assert mirrored.strand == 1
