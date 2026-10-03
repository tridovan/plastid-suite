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
