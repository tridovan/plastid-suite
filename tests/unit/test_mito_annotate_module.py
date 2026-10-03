"""The mito annotation module (P3A-03) — the wiring, the tables and the artifact.

BLAST itself is stubbed here: the searches are ``external/blast.py``'s job and its own tests hold
the command line still. What is tested is everything around them, and the four refusals the module
docstring states — no merging with the input's annotation, no invented product text, no success on
an empty search, and no artifact that has not been read back with our own reader.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from plastid.core.model import Gene, Genome, Interval
from plastid.core.sequence import reverse_complement
from plastid.external.blast import BLAST_VERSION
from plastid.modules import mito_annotate
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

SEQ_ID = "SYN_MITO.1"
_PROTEIN = "M" * 9
# The same synthetic mitogenome the caller's tests use: a plus-strand CDS at 1..30, a
# minus-strand one at 61..90, then filler.
GENOME = Genome(
    name="synthetic mito",
    accession=SEQ_ID,
    sequence="ATG"
    + "GCT" * 8
    + "TAA"
    + "GCT" * 9
    + "TTA"
    + reverse_complement("ATG" + "GCT" * 8 + "TAA")
    + "ATGC" * 40,
    topology="circular",
)

# What the stubbed searches print. The ccmB hit covers the reference protein (9 aa) and stops on
# the stop codon the sequence confirms; the rrn5 hit covers its whole reference.
BLASTX_HITS = (
    f"{SEQ_ID}\tNC_999999.1|ccmB|ccmB|Testus testus\t99.0\t9\t0\t0\t1\t27\t1\t9\t1e-50\t200\t1\n"
)
# One tRNA in ARAGORN's measured batch shape (the real lines live in test_aragorn.py).
ARAGORN_HITS = (
    f">{SEQ_ID}\n1 genes found\n1   tRNA-Phe                     [181,250]\t35  \t(gaa)\n"
)
BLASTN_HITS = (
    f"{SEQ_ID}\tNC_999999.1|rrn5|rrn5|Testus testus\t96.0\t116\t0\t0\t61\t176\t1\t116\t1e-60\t300\n"
)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "mito-annotate-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="mito_annotate", module_order=1)


def _stub(directory: Path, name: str, output: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / name
    if name == "aragorn":
        script.write_text(
            "#!/bin/sh\n"
            'if [ "$1" = "-h" ]; then echo "ARAGORN v1.2.41 Dean Laslett"; exit 0; fi\n'
            f"cat <<'HITS'\n{output}HITS\n",
            encoding="utf-8",
        )
    else:
        script.write_text(
            "#!/bin/sh\n"
            'case "$1" in\n'
            f"  --version|-version) echo '{name}: {BLAST_VERSION}+'; exit 0;;\n"
            "esac\n"
            f"cat <<'HITS'\n{output}HITS\n",
            encoding="utf-8",
        )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)


def _binaries(
    tmp_path: Path,
    *,
    blastx: str = BLASTX_HITS,
    blastn: str = BLASTN_HITS,
    aragorn: str = ARAGORN_HITS,
) -> Path:
    directory = tmp_path / "bin"
    _stub(directory, "makeblastdb", "")
    _stub(directory, "blastx", blastx)
    _stub(directory, "blastn", blastn)
    _stub(directory, "aragorn", aragorn)
    return directory


def _reference(tmp_path: Path) -> Path:
    """A two-gene reference record: ccmB (a real translation) and rrn5 (a real rRNA length)."""
    directory = tmp_path / "reference"
    directory.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(Seq("ATGC" * 100), id="NC_999999.1", name="NC_999999.1")
    record.annotations = {"molecule_type": "DNA", "topology": "circular", "organism": "Testus"}
    record.features = [
        SeqFeature(
            SimpleLocation(0, 10, strand=1), type="source", qualifiers={"organism": ["Testus"]}
        ),
        SeqFeature(
            SimpleLocation(0, 30, strand=1),
            type="CDS",
            qualifiers={"gene": ["ccmB"], "translation": [_PROTEIN], "codon_start": ["1"]},
        ),
        SeqFeature(SimpleLocation(100, 216, strand=1), type="rRNA", qualifiers={"gene": ["rrn5"]}),
    ]
    SeqIO.write(record, str(directory / "NC_999999.1.gb"), "genbank")
    return directory


def _params(tmp_path: Path, **extra: object) -> dict[str, object]:
    params: dict[str, object] = {
        "reference": str(_reference(tmp_path)),
        "tools": str(_binaries(tmp_path)),
    }
    params.update(extra)
    return params


def _run(tmp_path: Path, genome: Genome = GENOME, **extra: object):
    return mito_annotate.run(genome, _params(tmp_path, **extra), _ctx(tmp_path))


# --------------------------------------------------------------------------- the tables


def test_the_module_calls_what_the_hits_support_and_lists_what_it_did_not(tmp_path: Path) -> None:
    result = _run(tmp_path)

    called = {row["gene"] for row in result.tables["called_genes"].rows}
    assert called == {"ccmB", "rrn5"}
    # the other 30 vocabulary loci and the other two rRNA have no reference in this tiny set
    assert {row["gene"] for row in result.tables["not_called"].rows} >= {"nad1", "rrn18"}
    summary = result.tables["mito_annotation_summary"].rows[0]
    assert (summary["n_calls"], summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
        3,
        1,
        1,
        1,
    )
    assert [row["gene"] for row in result.tables["trna_calls"].rows] == ["trnF-GAA"]


def test_every_table_carries_a_scope_label(tmp_path: Path) -> None:
    result = _run(tmp_path)
    for table in result.tables.values():
        assert table.scope.label, f"{table.id} has no scope label"
        assert table.scope.coord_system == "genbank"


def test_the_summary_states_the_floors_the_reference_and_the_tool(tmp_path: Path) -> None:
    row = _run(tmp_path).tables["mito_annotation_summary"].rows[0]

    assert (row["min_identity_pct"], row["min_coverage_pct"]) == (50.0, 0.5)
    assert (row["blast"], row["aragorn"]) == (f"BLAST+ {BLAST_VERSION}", "ARAGORN 1.2.41")
    assert "rebuilt from 1 GenBank record" in row["reference"]
    assert row["transl_table"] == 1


def test_a_called_cds_carries_an_orf_status_from_the_shared_checks(tmp_path: Path) -> None:
    """No second code path: the status column is core/orf's, the same function that scores a
    downloaded record."""
    row = next(row for row in _run(tmp_path).tables["called_genes"].rows if row["gene"] == "ccmB")

    assert row["orf_status"] == "ok"
    assert (row["start"], row["end"], row["strand"]) == (1, 30, "+")


def test_a_called_gene_carries_no_product_text(tmp_path: Path) -> None:
    """The vocabulary's synonyms were quoted from the answer-key record; writing them back would
    be this module grading its own homework."""
    result = _run(tmp_path)

    for gene in read_back(result).genes:
        assert gene.product is None
        assert gene.name in {"ccmB", "rrn5", "trnF-GAA"}


def test_a_reference_gene_without_a_vocabulary_entry_is_reported_not_named(tmp_path: Path) -> None:
    """A hit against a reference this toolkit cannot name shows up in its own table and never in
    the annotation — H-06 owns the naming decision."""
    directory = _reference(tmp_path)
    record = next(SeqIO.parse(str(directory / "NC_999999.1.gb"), "genbank"))
    record.features.append(
        SeqFeature(
            SimpleLocation(300, 330, strand=1),
            type="CDS",
            qualifiers={"gene": ["rps12"], "translation": ["M" * 10], "codon_start": ["1"]},
        )
    )
    SeqIO.write(record, str(directory / "NC_999999.1.gb"), "genbank")
    hits = BLASTX_HITS + (
        f"{SEQ_ID}\tNC_999999.1|rps12|-\t95.0\t10\t0\t0\t181\t210\t1\t10\t1e-40\t150\t1\n"
    )
    binaries = _binaries(tmp_path, blastx=hits)

    result = mito_annotate.run(
        GENOME,
        {"reference": str(directory), "tools": str(binaries)},
        _ctx(tmp_path),
    )

    assert [row["gene_id"] for row in result.tables["unnamed_hits"].rows] == [
        "NC_999999.1|rps12:181-210"
    ]
    assert "rps12" not in {gene.name for gene in read_back(result).genes}


# --------------------------------------------------------------------------- the artifact


def read_back(result):
    """The artifact, through our own reader — the same check the module itself performs."""
    import tempfile

    from plastid.core.genbank import read_genbank

    with tempfile.TemporaryDirectory() as work:
        path = Path(work) / "artifact.gb"
        path.write_text(result.artifacts["mito_annotation.gb"], encoding="utf-8")
        return read_genbank(path)


def test_the_artifact_re_reads_and_carries_only_the_calls(tmp_path: Path) -> None:
    """Measured trap: run this on a record that is already annotated and appending would double
    its gene count — which is exactly the record the caller is checked against."""
    annotated = Genome(
        name=GENOME.name,
        sequence=GENOME.sequence,
        accession=GENOME.accession,
        topology=GENOME.topology,
        genes=(Gene("old_annotation", "CDS", (Interval(0, 30),)),),
    )
    result = _run(tmp_path, genome=annotated)
    back = read_back(result)

    assert sorted(gene.name for gene in back.genes) == ["ccmB", "rrn5", "trnF-GAA"]
    assert back.sequence == annotated.sequence
    assert "old_annotation" not in result.artifacts["mito_annotation.gb"]


def test_the_raw_search_output_travels_with_the_run(tmp_path: Path) -> None:
    result = _run(tmp_path)
    assert result.artifacts["blast_hits.tsv"] == BLASTX_HITS + BLASTN_HITS + ARAGORN_HITS


# --------------------------------------------------------------------------- refusals


def test_a_search_that_finds_nothing_is_a_failed_run_not_an_empty_annotation(
    tmp_path: Path,
) -> None:
    directory = _binaries(tmp_path, blastx="", blastn="")

    with pytest.raises(ToolError, match="empty annotation is a failed annotation"):
        mito_annotate.run(
            GENOME,
            {"reference": str(_reference(tmp_path)), "tools": str(directory)},
            _ctx(tmp_path),
        )


def test_the_floors_can_be_overridden_and_the_summary_says_so(tmp_path: Path) -> None:
    """ccmB aligns at 99 % and rrn5 at 96 %: raising the floor to 98 % must drop one and keep the
    other, and the run has to say which floor it used."""
    result = _run(tmp_path, min_identity=98.0)

    assert result.tables["mito_annotation_summary"].rows[0]["min_identity_pct"] == 98.0
    assert {row["gene"] for row in result.tables["called_genes"].rows} == {"ccmB"}
    row = next(row for row in result.tables["not_called"].rows if row["gene"] == "rrn5")
    assert "identity floor" in row["reason"]


def test_a_run_where_nothing_passes_the_floors_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="no locus passed the floors"):
        _run(tmp_path, min_identity=99.9)


def test_a_reference_directory_that_holds_no_record_is_refused(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError, match="no \\*.gb records"):
        mito_annotate.run(
            GENOME, {"reference": str(empty), "tools": str(_binaries(tmp_path))}, _ctx(tmp_path)
        )


def test_a_missing_blast_is_named_the_way_the_cli_reports_it(tmp_path: Path) -> None:
    """Exit 3 with a message that names the fix — the CLI turns ToolError into exactly that."""
    missing = tmp_path / "no-tools"
    with pytest.raises(ToolError, match="docs/tools/blast.md"):
        mito_annotate.run(
            GENOME, {"reference": str(_reference(tmp_path)), "tools": str(missing)}, _ctx(tmp_path)
        )


def test_the_module_does_not_ask_for_an_ir_result(tmp_path: Path) -> None:
    """A plant mitogenome has no canonical IR; both stand-in mitogenomes fail detection, and a
    module that gated here would refuse to describe a genome it can read perfectly well. The
    context below carries no IR result at all."""
    ctx = _ctx(tmp_path)
    assert ctx.ir is None
    result = mito_annotate.run(GENOME, _params(tmp_path), ctx)

    summary = result.tables["mito_annotation_summary"].rows[0]
    assert (summary["n_calls"], summary["n_cds"], summary["n_trna"], summary["n_rrna"]) == (
        3,
        1,
        1,
        1,
    )
    assert [row["gene"] for row in result.tables["trna_calls"].rows] == ["trnF-GAA"]
    assert any("never asks for one" in note for note in result.notes)


def test_the_notes_state_the_scope_and_the_tool_that_is_not_here(tmp_path: Path) -> None:
    notes = " ".join(_run(tmp_path).notes)

    assert "no de-novo gene discovery" in notes
    assert "ORF-walking" in notes
    assert "ARAGORN 1.2.41" in notes and "no second naming rule" in notes
