"""The ARAGORN adapter (P3A-04) — parsing its batch output, its names, and its version.

The sample lines are the real ones ARAGORN 1.2.41 printed for the answer-key record on
03/10/2026, including the two shapes that would silently corrupt an annotation: a tRNA whose
coordinates **wrap the origin** (``[258274,73]``) and the anticodon printed as DNA rather than in
the RNA form the published names use.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.core.model import Genome
from plastid.external.aragorn import (
    ARAGORN_VERSION,
    AragornResult,
    aragorn_version,
    parse_batch,
    parse_calls,
    run_aragorn,
    trna_name,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

# Verbatim from the measured run, including its header lines.
BATCH = """>OX459128.1
17 genes found
1   tRNA-Phe                   [258274,73]\t35  \t(gaa)
2   tRNA-Pro                     [329,403]\t36  \t(tgg)
3   tRNA-Met                   [3691,3766]\t36  \t(cat)
8   tRNA-Met              c[107050,107126]\t35  \t(cat)
13  tRNA-His              c[228034,228108]\t35  \t(gtg)
17  tRNA-Lys                 [256837,256909]\t34  \t(ttt)
"""

GENOME = Genome(
    name="OX459128.1",
    accession="OX459128.1",
    # the answer-key record's own length: its coordinates mean nothing against another one
    sequence="ATGC" * 64568 + "AT",
    topology="circular",
)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "aragorn-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="mito_annotate", module_order=1)


def _stub(tmp_path: Path, output: str = BATCH, version: str = ARAGORN_VERSION) -> Path:
    script = tmp_path / "aragorn"
    script.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "-h" ]; then\n'
        f"  echo '----------------------------'\n"
        f"  echo 'ARAGORN v{version} Dean Laslett'\n"
        f"  echo '----------------------------'\n"
        "  exit 0\n"
        "fi\n"
        f"cat <<'HITS'\n{output}HITS\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


# --------------------------------------------------------------------------- the parser


def test_the_batch_lines_are_read_with_their_names_and_anticodons() -> None:
    hits, others = parse_batch(BATCH)

    assert [hit.name for hit in hits] == [
        "trnF-GAA",
        "trnP-UGG",
        "trnM-CAU",
        "trnM-CAU",
        "trnH-GUG",
        "trnK-UUU",
    ]
    assert others == ()
    third = hits[2]
    assert (third.start, third.end, third.strand) == (3691, 3766, 1)
    assert third.source_name == "tRNA-Met"


def test_a_minus_strand_hit_is_marked_by_the_c_before_its_coordinates() -> None:
    hits, _ = parse_batch(BATCH)
    minus = hits[3]

    assert (minus.start, minus.end, minus.strand) == (107050, 107126, -1)
    assert minus.anticodon == "CAU"


def test_a_hit_that_crosses_the_origin_becomes_two_parts_in_transcript_order() -> None:
    """Measured: the record's first tRNA is printed ``[258274,73]`` — a start greater than the
    end. That is not a GenBank coordinate, and writing it as one would put a 74-base gene at
    position 258,274 of a 258,274-base record."""
    hit = parse_batch(BATCH)[0][0]

    assert hit.wrapped is True
    parts = hit.parts(len(GENOME))
    assert [part.to_1based() for part in parts] == [(258274, 258274), (1, 73)]
    assert hit.gene(len(GENOME)).length == 74


def test_a_gene_that_is_not_a_trna_is_reported_not_parsed_as_one() -> None:
    text = BATCH + "18  tmRNA                     [1000,1300]\t80  \t(gga)\n"
    hits, others = parse_batch(text)

    assert len(hits) == 6
    assert len(others) == 1 and "tmRNA" in others[0]


def test_a_trna_whose_anticodon_is_a_stop_codon_is_reported_not_a_crash() -> None:
    """Measured 09/10/2026 on the pinned corpus: ``NC_072561.1`` and ``NC_072580.1`` each bring
    one ``tRNA-Stop`` line at the ``trnK-UUU`` locus — ARAGORN reading the anticodon one base off,
    so that it spells a stop codon. The three records carry the **same** exon-1 sequence.

    The gene is real and the name is what this adapter cannot supply, so the line is *reported*
    (the module's notes carry it, next to ARAGORN's own "N genes found") — neither dropped, which
    would make the gene look absent, nor fatal, which would cost the other thirty-seven tRNAs and
    every CDS of that record with it.
    """
    text = BATCH + "18  tRNA-Stop             c[1399,4136]\t35  \t(tta)i(40,2662)\n"
    hits, others = parse_batch(text)

    assert len(hits) == 6, "a line this adapter cannot name took the whole run with it"
    assert len(others) == 1 and "tRNA-Stop" in others[0]
    # the refusal itself is still a refusal: asking for the name directly raises, and nothing
    # here guesses a residue letter for it
    with pytest.raises(ValueError, match="not one of the standard names"):
        trna_name("Stop", "tta")


def test_an_intron_splits_the_hit_into_the_two_exons_a_join_needs() -> None:
    """ARAGORN was measured to skip every intron-containing tRNA unless ``-i`` is passed (its help
    says so; the golden lost 8 of 38 that way). With ``-i`` its span **covers the intron**, and a
    feature written over that span would claim the intron is coding sequence.

    The intron is reported in the gene's own coordinates — 1-based from its 5' end — so the
    conversion depends on the strand: on the minus strand the 5' end is the *higher* coordinate,
    and walking towards the 3' end means walking down the genome.

    **The second field is the intron's length, not its end** — corrected 10/10/2026, and
    this test asserted the wrong reading until then. Read as an end, ``trnL-UAA``'s ``i(36,406)``
    spliced into 120 nt and ``trnK-UUU``'s ``i(39,2618)`` into 112 — lengths no tRNA has — while
    read as a length they splice into 35+50 = 85 and 38+36 = 74, and the published record's own
    joined exons for ``trnL-UAA`` are **exactly** 35+50. ARAGORN's help documents the intron's
    maximum *length* (``-i<min>,<max>``) and never mentions an end coordinate.
    """
    plus = parse_batch("1   tRNA-Ser                 [8889,9634]\t688 \t(cga)i(32,655)\n")[0][0]
    assert plus.intron == (32, 655)
    # 31 + 60 = 91 nt spliced; the intron is positions 32..686 and is not coding sequence
    assert [part.to_1based() for part in plus.parts(100_000)] == [(8889, 8919), (9575, 9634)]

    minus = parse_batch("1   tRNA-Lys               c[1367,4058]\t34  \t(ttt)i(39,2618)\n")[0][0]
    parts = minus.parts(153_724)
    assert [part.to_1based() for part in parts] == [(4021, 4058), (1367, 1402)]
    assert {part.strand for part in parts} == {-1}


def test_a_spliced_trna_comes_out_the_length_a_trna_can_be() -> None:
    """The measurement that caught the reading: every intron tRNA on the golden, spliced.

    A tRNA is 72–92 nt. The published record writes those eight as two-part joins, and its own
    ``trnL-UAA`` is 35 + 50 — the number this parser has to reproduce, on the strand where the
    5' exon is the high coordinate.
    """
    lines = [
        "1   tRNA-Lys               c[1367,4058]\t34  \t(ttt)i(39,2618)\n",
        "16  tRNA-Leu                 [47453,47943]\t35  \t(taa)i(36,406)\n",
        "27  tRNA-Ala                 [103034,103915]\t34 \t(tgc)i(38,809)\n",
    ]
    spliced = {}
    for line in lines:
        hit = parse_batch(line)[0][0]
        spliced[hit.name] = sum(part.end - part.start for part in hit.parts(153_724))

    assert spliced == {"trnK-UUU": 74, "trnL-UAA": 85, "trnA-UGC": 73}
    assert all(60 <= length <= 95 for length in spliced.values()), spliced


def test_the_declared_count_is_read_from_aragorns_own_line() -> None:
    from plastid.external.aragorn import _DECLARED

    match = _DECLARED.search(BATCH)
    assert match is not None and match.group("count") == "17"


# --------------------------------------------------------------------------- the names


def test_the_name_is_the_published_form_and_comes_from_the_shared_normalizer() -> None:
    """ARAGORN prints an amino acid in words and the anticodon as DNA; the record names the same
    gene ``trnM-CAU``. The conversion ends at ``normalize_name`` — the same function a downloaded
    record goes through, so there is no second naming rule in this toolkit."""
    assert trna_name("Met", "cat") == "trnM-CAU"
    assert trna_name("Phe", "gaa") == "trnF-GAA"
    assert trna_name("His", "gtg") == "trnH-GUG"


def test_an_amino_acid_this_adapter_does_not_know_is_refused() -> None:
    with pytest.raises(ValueError, match="not one of the standard names"):
        trna_name("Xaa", "cat")


# --------------------------------------------------------------------------- the binary


def test_the_version_comes_from_the_help_banner(tmp_path: Path) -> None:
    """ARAGORN has no version flag: ``-V`` exits **0** and prints an error sentence (measured),
    which is why the probe reads ``-h``."""
    assert aragorn_version(str(_stub(tmp_path))) == ARAGORN_VERSION


def test_another_version_is_refused_by_name(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="1.2.38"):
        aragorn_version(str(_stub(tmp_path, version="1.2.38")))


def test_a_missing_binary_names_the_document_that_declares_it(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="docs/tools/aragorn.md"):
        run_aragorn(GENOME, _ctx(tmp_path), seq_id="OX459128.1", directory=tmp_path / "empty")


# --------------------------------------------------------------------------- a run


def test_a_run_returns_the_hits_their_genes_and_what_was_printed(tmp_path: Path) -> None:
    script = _stub(tmp_path)
    ctx = _ctx(tmp_path)

    result = run_aragorn(GENOME, ctx, seq_id="OX459128.1", aragorn=str(script))

    assert isinstance(result, AragornResult)
    assert (result.declared, len(result.hits), len(result.genes)) == (17, 6, 6)
    assert result.version == ARAGORN_VERSION
    assert result.raw == BATCH
    assert [gene.kind for gene in result.genes] == ["tRNA"] * 6
    argv = ctx.tools[-1].argv
    assert "-w" in argv and "-gc1" in argv and "-c" in argv
    assert ctx.tools[-1].version == ARAGORN_VERSION


def test_a_linear_record_is_searched_as_a_linear_one(tmp_path: Path) -> None:
    """The circular flag makes ARAGORN wrap its search around the ends — right for a mitogenome,
    wrong for a fragment, and it changes what is found at the boundaries."""
    linear = Genome(
        name="fragment", accession="FRAG.1", sequence="ATGC" * 64568 + "AT", topology="linear"
    )
    ctx = _ctx(tmp_path)

    run_aragorn(linear, ctx, seq_id="FRAG.1", aragorn=str(_stub(tmp_path)))

    assert "-l" in ctx.tools[-1].argv


def test_an_unnameable_line_keeps_its_coordinates_for_locus_evidence() -> None:
    """H-24: the line this adapter cannot name is still a gene at a place.

    ``parse_calls`` keeps every tRNA line, named or not, so ``plastid_annotate`` can name an
    intron tRNA from the reference set's own annotation at that locus. ``parse_batch``'s public
    contract is unchanged: the line is still reported in ``others``, never silently called, and
    asking for its gene without a name raises rather than inventing one.
    """
    text = BATCH + "18  tRNA-Stop             c[1399,4136]\t35  \t(tta)i(40,2662)\n"

    hits, others = parse_batch(text)
    calls = parse_calls(text)

    assert len(hits) == 6 and len(others) == 1
    assert len(calls) == 7, "the unnameable line vanished from the locus evidence"
    unnamed = calls[-1]
    assert unnamed.name is None
    assert (unnamed.start, unnamed.end, unnamed.strand) == (1399, 4136, -1)
    assert unnamed.intron == (40, 2662)
    assert "tRNA-Stop" in unnamed.raw_line
    with pytest.raises(ValueError, match="has no name"):
        unnamed.gene(len(GENOME))
    assert unnamed.with_name("trnV-UAC").gene(len(GENOME)).name == "trnV-UAC"


def test_the_intron_range_is_forwarded_and_carried_on_the_result(tmp_path: Path) -> None:
    """H-22: a consumer has to know what was searched, not only that a gene was absent.

    ``introns=False`` with an explicit range is not a contradiction: the range *is* the request to
    search. The switch is ``-i<min>,<max>`` and the result records the same pair.
    """
    ctx = _ctx(tmp_path)

    result = run_aragorn(
        GENOME,
        ctx,
        seq_id="OX459128.1",
        aragorn=str(_stub(tmp_path)),
        intron_limits=(500, 2_500),
    )

    assert "-i500,2500" in ctx.tools[-1].argv
    assert result.intron_search == (500, 2_500)


def test_a_search_without_introns_records_that_it_did_not_look(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)

    result = run_aragorn(GENOME, ctx, seq_id="OX459128.1", aragorn=str(_stub(tmp_path)))

    assert "-i" not in ctx.tools[-1].argv
    assert result.intron_search is None


def test_a_range_beyond_the_caller_maximum_is_refused_by_name(tmp_path: Path) -> None:
    """Measured 26/10/2026: `-i0,4000` is accepted by ARAGORN and still stops at 3,000."""
    ctx = _ctx(tmp_path)

    with pytest.raises(ValueError, match=r"does not search a tRNA intron longer than 3,000 bp"):
        run_aragorn(
            GENOME,
            ctx,
            seq_id="OX459128.1",
            aragorn=str(_stub(tmp_path)),
            introns=True,
            intron_limits=(0, 4_000),
        )
