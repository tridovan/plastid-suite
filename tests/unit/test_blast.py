"""The BLAST+ adapter (P3A-03) — the two halves that must not drift apart: the command line and
the parser for what it prints.

BLAST+ writes tab-separated output whose *column order is whatever ``-outfmt`` asked for*, so the
field list and the parser are one contract. They live in the same module for that reason, and the
argv is asserted here without running BLAST — the same trick ``test_chloe.py`` uses.

Two measured facts are pinned as tests, because both are silent when they go wrong:

* BLAST truncates ``sseqid`` at the first whitespace, so a subject key is the first three ``|``
  fields (the reference headers end with an organism name);
* ``makeblastdb -parse_seqids`` **fails** on the committed reference set (``Error: Duplicate
  seq_ids are found: LCL|NC_006581.1|SDH3|-|NICOTIANA``) and upper-cases the ids it does parse —
  so the flag the obvious naming suggests is exactly the one that breaks the pipeline.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.core.model import Genome
from plastid.external.blast import (
    BLAST_VERSION,
    SEGMENT_OFF,
    blast_version,
    makeblastdb_argv,
    parse_hits,
    resolve_binary,
    run_blastn,
    run_blastx,
    search_argv,
    subject_count,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

TARGET = Genome(name="synthetic mito", accession="SYN_MITO.1", sequence="ATGC" * 100)

# What a real blastx run printed on 03/10/2026 (the answer-key record against the committed
# reference set), kept verbatim in shape: the organism field is already truncated by BLAST.
BLASTX_OUT = (
    "OX459128.1\tNC_057654.1|matR|matR|Scyphiphora\t95.525\t648\t28\t1\t196060\t194117\t1\t647"
    "\t0.0\t1196\t-1\n"
    "OX459128.1\tNC_006581.1|mat-R|matR|Nicotiana\t92.427\t647\t48\t1\t196060\t194120\t1\t646"
    "\t0.0\t1147\t-1\n"
    "OX459128.1\tNC_006581.1|rpl16|-\t42.0\t138\t80\t2\t78383\t78769\t31\t168"
    "\t3.3e-27\t210\t-1\n"
)
BLASTN_OUT = (
    "OX459128.1\tNC_057654.1|rrn18|rrn18|Scyphiphora\t97.783\t1849\t40\t1\t197230\t199073"
    "\t1848\t1\t0.0\t3182\n"
)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "blast-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="mito_annotate", module_order=1)


def _stub(directory: Path, name: str, output: str) -> Path:
    """A stand-in for one BLAST binary: it answers the version probe and prints ``output``."""
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / name
    script.write_text(
        "#!/bin/sh\n"
        'case "$1" in\n'
        "  --version|-version) echo '" + name + ": " + BLAST_VERSION + "+'; exit 0;;\n"
        "esac\n"
        f"cat <<'HITS'\n{output}HITS\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script


# --------------------------------------------------------------------------- the parser


def test_the_parser_reads_the_columns_it_asked_for() -> None:
    hits = parse_hits(BLASTX_OUT, seq_id="OX459128.1", frame=True)

    assert len(hits) == 3
    first = hits[0]
    assert first.subject_key == "NC_057654.1|matR|matR"
    assert (first.pident, first.length, first.bitscore) == (95.525, 648, 1196.0)
    assert (first.qstart, first.qend) == (196060, 194117)
    assert (first.sstart, first.send) == (1, 647)
    assert first.frame == -1


def test_a_subject_id_with_a_space_is_read_as_its_first_three_fields() -> None:
    """Measured: ``NC_006581.1|mat-R|matR|Nicotiana tabacum`` arrives as ``…|Nicotiana`` — the
    organism field never survives, so a key built from four fields would match nothing."""
    hits = parse_hits(BLASTX_OUT, seq_id="OX459128.1", frame=True)

    assert hits[1].subject_key == "NC_006581.1|mat-R|matR"
    assert hits[2].subject_key == "NC_006581.1|rpl16|-", "an unnameable reference keeps its dash"


def test_a_nucleotide_hit_has_no_frame() -> None:
    hit = parse_hits(BLASTN_OUT, seq_id="OX459128.1", frame=False)[0]

    assert hit.frame is None
    assert (hit.qstart, hit.qend) == (197230, 199073)
    assert (hit.sstart, hit.send) == (1848, 1)
    assert hit.strand == -1, "the strand comes from the subject, since a blastn field has none"


def test_a_subject_id_with_fewer_than_three_fields_is_refused() -> None:
    line = "OX459128.1\tNC_057654.1|matR\t95.5\t648\t28\t1\t1\t648\t1\t647\t0.0\t1196\t-1\n"
    with pytest.raises(ValueError, match="three"):
        parse_hits(line, seq_id="OX459128.1", frame=True)


def test_a_hit_for_another_query_is_refused() -> None:
    """One query goes in, so a line for something else means the file is not this run's output."""
    with pytest.raises(ValueError, match="query"):
        parse_hits(BLASTX_OUT, seq_id="AP045469.1", frame=True)


def test_a_malformed_line_names_itself() -> None:
    with pytest.raises(ValueError, match="expected 13"):
        parse_hits("OX459128.1\tNC_057654.1|matR|matR\t95.5\n", seq_id="OX459128.1", frame=True)


def test_an_empty_result_is_an_empty_tuple_not_an_error() -> None:
    """A run that finds nothing is a fact; whether that is acceptable is the caller's decision."""
    assert parse_hits("", seq_id="OX459128.1", frame=True) == ()


# --------------------------------------------------------------------------- the command line


def test_the_search_argv_pins_what_the_measurements_decided() -> None:
    argv = search_argv(
        "blastx",
        binary="/usr/bin/blastx",
        query=Path("/tmp/target.fna"),
        database=Path("/tmp/refdb"),
        evalue=1e-5,
        max_hsps=10,
        max_target_seqs=287,
    )
    joined = " ".join(argv)

    assert argv[0] == "/usr/bin/blastx"
    assert "-query /tmp/target.fna" in joined
    assert "-db /tmp/refdb" in joined
    assert "-evalue 1e-05" in joined
    assert "-max_hsps 10" in joined
    assert "-max_target_seqs 287" in joined
    assert "qframe" in joined, "a protein hit's strand comes from its frame"
    assert " pident " in joined and " sstart " in joined and " bitscore " in joined


def test_the_target_cap_is_the_size_of_the_reference_set(tmp_path: Path) -> None:
    """Measured 08/10/2026 (P3B-07): a **constant** cap truncates in silence.

    The 500 that sufficed for the mito path's 287 sequences met the plastid set's 2,838 and BLAST
    returned 929 HSPs over 467 subjects where the same search with the cap lifted returned 4,382
    over 2,603 — and nothing in the output said a subject had been dropped. So the value is
    counted from the file the database is built from, and this test holds that it is counted
    rather than chosen.
    """
    reference = tmp_path / "reference.faa"
    reference.write_text(
        "".join(f">NC_00000{i}.1|gene{i}|gene{i}\nMAAA\n" for i in range(7)), encoding="utf-8"
    )
    assert subject_count(reference) == 7
    argv = search_argv(
        "blastx",
        binary="blastx",
        query=Path("q"),
        database=Path("d"),
        evalue=1e-5,
        max_hsps=10,
        max_target_seqs=subject_count(reference),
    )
    assert "-max_target_seqs 7" in " ".join(argv)


def test_seg_masking_is_off_and_the_test_says_why() -> None:
    """Measured: SEG masking halves atp9's best alignment (137 → 67 bits) and shortens six other
    real genes — it is built to hide low-complexity *query* sequence, and a plant mitogenome is
    full of hydrophobic membrane proteins. The floor, not the mask, is what keeps junk out."""
    argv = search_argv(
        "blastx",
        binary="blastx",
        query=Path("q"),
        database=Path("d"),
        evalue=1e-5,
        max_hsps=10,
        max_target_seqs=1,
    )
    assert (SEGMENT_OFF, "no") == ("-seg", "no")
    assert "-seg no" in " ".join(argv)


def test_a_nucleotide_search_asks_for_the_sensitive_task() -> None:
    """Megablast and blastn found the same three rRNAs here, but the target may be a divergent
    sample and the database is 25 sequences — sensitivity is nearly free."""
    argv = search_argv(
        "blastn",
        binary="blastn",
        query=Path("q"),
        database=Path("d"),
        evalue=1e-5,
        max_hsps=10,
        max_target_seqs=1,
    )
    assert "-task blastn" in " ".join(argv)
    assert "qframe" not in " ".join(argv)


def test_makeblastdb_is_not_asked_to_parse_the_sequence_ids() -> None:
    """Measured 03/10/2026: ``-parse_seqids`` fails outright on the committed reference set —
    ``Error: Duplicate seq_ids are found: LCL|NC_006581.1|SDH3|-|NICOTIANA`` — and it upper-cases
    the ids it does parse, which would break every lowercase gene symbol in the key."""
    argv = makeblastdb_argv(
        binary="/usr/bin/makeblastdb", reference=Path("/tmp/ref.faa"), database=Path("/tmp/refdb")
    )

    assert "-parse_seqids" not in argv
    assert "-dbtype prot" in " ".join(argv)


# --------------------------------------------------------------------------- the binaries


def test_a_missing_binary_names_the_document_that_declares_it(tmp_path: Path) -> None:
    """The message has to say where the tool is declared and how to point a run at one — Điềm's
    machine keeps BLAST+ outside PATH, so "not found" alone would be a dead end."""
    with pytest.raises(ToolError, match="docs/tools/blast.md"):
        resolve_binary("blastx", explicit=tmp_path / "nowhere" / "blastx")
    with pytest.raises(ToolError, match="docs/tools/blast.md"):
        resolve_binary("blastx", directory=tmp_path / "empty")


def test_a_binary_reporting_another_version_is_refused(tmp_path: Path) -> None:
    script = tmp_path / "blastx"
    script.write_text("#!/bin/sh\necho 'blastx: 2.11.0+'\n", encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)

    with pytest.raises(ToolError, match="2.11.0"):
        blast_version(str(script))


def test_the_version_probe_reads_blasts_own_header(tmp_path: Path) -> None:
    script = _stub(tmp_path, "blastx", "")
    assert blast_version(str(script)) == BLAST_VERSION


# --------------------------------------------------------------------------- an end-to-end stub


def test_a_run_records_both_tools_and_returns_the_parsed_hits(tmp_path: Path) -> None:
    binaries = tmp_path / "bin"
    blastx = _stub(binaries, "blastx", BLASTX_OUT)
    makeblastdb = _stub(binaries, "makeblastdb", "")
    ctx = _ctx(tmp_path)
    reference = tmp_path / "reference.faa"
    reference.write_text(">NC_057654.1|matR|matR|Scyphiphora\nMKK\n", encoding="utf-8")

    result = run_blastx(
        TARGET,
        reference,
        ctx,
        seq_id="OX459128.1",
        blastx=str(blastx),
        makeblastdb=str(makeblastdb),
    )

    assert len(result.hits) == 3
    assert result.raw == BLASTX_OUT
    assert [record.name for record in ctx.tools] == ["makeblastdb", "blastx"]
    assert all("-seg" in record.argv for record in ctx.tools if record.name == "blastx")
    assert result.version == BLAST_VERSION


def test_a_nucleotide_run_parses_a_nucleotide_file(tmp_path: Path) -> None:
    binaries = tmp_path / "bin"
    blastn = _stub(binaries, "blastn", BLASTN_OUT)
    makeblastdb = _stub(binaries, "makeblastdb", "")
    ctx = _ctx(tmp_path)
    reference = tmp_path / "rrna.fa"
    reference.write_text(">NC_057654.1|rrn18|rrn18|Scyphiphora\nACGT\n", encoding="utf-8")

    result = run_blastn(
        TARGET,
        reference,
        ctx,
        seq_id="OX459128.1",
        blastn=str(blastn),
        makeblastdb=str(makeblastdb),
    )

    assert [hit.subject_key for hit in result.hits] == ["NC_057654.1|rrn18|rrn18"]
    assert result.hits[0].frame is None
