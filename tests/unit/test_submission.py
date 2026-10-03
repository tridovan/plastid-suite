"""NCBI submission files (P3-05): the five-column rules, partials, and the round trip.

The format rules here are not invented — they are NCBI's, and each one fails at *upload* time
if it is wrong, after the work is done. The two examples below are NCBI's own (BankIt feature
table help, retrieved 2026-10-03), re-tabulated: the help page renders them with the tab
columns collapsed, and the documented column definition says a qualifier name belongs in
column 4 and its value in column 5.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Gene, Genome, Interval
from plastid.core.submission import (
    parse_tbl,
    submission_problems,
    to_fasta,
    to_tbl,
)
from plastid.run.context import RunContext

# NCBI's documented example, in the documented five-column form (tabs, not spaces).
NCBI_JOIN_EXAMPLE = (
    ">Feature Seq2\n2626\t2590\ttRNA\t\t\n2570\t2535\t\t\t\n\t\t\tproduct\ttRNA-Phe\n"
)

NCBI_PARTIAL_EXAMPLE = (
    ">Feature Seq1\n"
    "<1\t>1050\tgene\t\t\n"
    "\t\t\tgene\tATH1\n"
    "<1\t1009\tCDS\t\t\n"
    "\t\t\tcodon_start\t2\n"
    "<1\t>1050\tmRNA\t\t\n"
)

_PAD = "ACGT" * 20


def _genome(*genes: Gene, length: int | None = None) -> Genome:
    sequence = _PAD if length is None else "A" * length
    return Genome(
        name="synthetic", accession="SYN001.1", sequence=sequence, genes=genes, topology="linear"
    )


def test_a_multi_interval_feature_is_stacked_not_joined() -> None:
    """There is no ``join()`` keyword in this format: intervals stack, key on the first line."""
    gene = Gene("trnK-UUU", "tRNA", (Interval(100, 200), Interval(300, 350)))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[0] == ">Features SYN001.1"
    assert lines[1] == "101\t200\ttRNA\t\t"
    assert lines[2] == "301\t350\t\t\t"
    assert lines[3] == "\t\t\tgene\ttrnK-UUU"


def test_a_minus_strand_feature_writes_its_interval_reversed() -> None:
    gene = Gene("rbcL", "CDS", (Interval(100, 200, -1),), product="rubisco")
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "200\t101\tCDS\t\t"
    assert "\t\t\tproduct\trubisco" in lines


def test_a_mixed_strand_trans_spliced_gene_keeps_each_interval_its_own_way() -> None:
    """rps12 in a plastome is exactly this: one exon on the plus strand, one on the minus."""
    gene = Gene(
        "rps12",
        "CDS",
        (Interval(100, 120, -1), Interval(300, 320, 1)),
        product="ribosomal protein S12",
        trans_splicing=True,
    )
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "120\t101\tCDS\t\t"
    assert lines[2] == "301\t320\t\t\t"
    assert "\t\t\ttrans_splicing\t" in lines


def test_a_partial_coordinate_keeps_its_symbol_in_the_right_column() -> None:
    """``<`` always in column 1, ``>`` always in column 2 — NCBI's rule, whatever the strand."""
    gene = Gene("ycf1", "CDS", (Interval(0, 100, 1, partial_start=True, partial_end=True),))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "<1\t>100\tCDS\t\t"


def test_a_partial_on_a_minus_strand_interval_stays_with_its_own_coordinate() -> None:
    """The coordinate is what is partial, not the column.

    Coordinate 1 of a minus-strand CDS is the feature's **3'** end, and it is the one the record
    marked incomplete. After the interval is reversed (higher coordinate first, as the format
    requires), that partial end sits in column 2 — and NCBI's rule puts ``>`` in column 2 and
    ``<`` in column 1 whatever the strand, so the marker travels with its coordinate and takes
    the column's symbol.
    """
    gene = Gene(
        "ndhF",
        "CDS",
        (Interval(0, 100, -1, partial_start=True, partial_end=False),),
    )
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "100\t>1\tCDS\t\t"


def test_a_partial_five_prime_end_of_a_minus_strand_feature_takes_the_other_symbol() -> None:
    """The mirror case: `complement(1..>100)` — coordinate 100 is the feature's 5' end."""
    gene = Gene("ndhF", "CDS", (Interval(0, 100, -1, partial_start=False, partial_end=True),))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "<100\t1\tCDS\t\t"


def test_a_partial_on_an_internal_interval_is_refused_not_dropped() -> None:
    """NCBI rejects internal partials; writing the feature without the marker would be a
    silently different record, so this refuses instead."""
    gene = Gene(
        "rps12",
        "CDS",
        (
            Interval(100, 120),
            Interval(200, 220, 1, partial_start=True),
            Interval(300, 320),
        ),
    )
    with pytest.raises(ValueError, match="internal interval"):
        to_tbl(_genome(gene), "SYN001.1")


def test_the_fasta_header_is_the_sequence_id_and_the_lines_are_wrapped() -> None:
    text = to_fasta(_genome(length=130), "SYN001.1", width=60)
    lines = text.splitlines()
    assert lines[0] == ">SYN001.1"
    assert [len(line) for line in lines[1:]] == [60, 60, 10]


def test_a_sequence_id_with_a_space_is_refused() -> None:
    """The ID must match between the FASTA and the ``>Features`` line, and a FASTA ID may not
    contain a space — a name like "Iris domestica" would split into two."""
    with pytest.raises(ValueError, match="free of spaces"):
        to_fasta(_genome(length=10), "Iris domestica")


def test_the_parsers_reads_ncbi_documented_examples() -> None:
    joined = parse_tbl(NCBI_JOIN_EXAMPLE)
    assert len(joined) == 1
    assert joined[0].kind == "tRNA"
    assert joined[0].parts == ((2626, 2590, False, False), (2570, 2535, False, False))
    assert joined[0].qualifiers == {"product": "tRNA-Phe"}

    partial = parse_tbl(NCBI_PARTIAL_EXAMPLE)
    assert [feature.kind for feature in partial] == ["gene", "CDS", "mRNA"]
    assert partial[0].parts == ((1, 1050, True, True),)
    assert partial[0].qualifiers == {"gene": "ATH1"}
    assert partial[1].qualifiers == {"codon_start": "2"}


def test_the_parser_reads_the_compact_qualifier_layout_too() -> None:
    """The BankIt help page renders qualifiers with the leading tabs collapsed; a table pasted
    from it still has to read."""
    text = ">Feature Seq2\n2626\t2590\ttRNA\t\t\n2570\t2535\t\t\t\nproduct\ttRNA-Phe\n"
    parsed = parse_tbl(text)
    assert parsed[0].qualifiers == {"product": "tRNA-Phe"}


def _canonical_parts(gene: Gene) -> tuple[tuple[tuple[int, bool], ...], ...]:
    """Each interval as coordinate-sorted (coordinate, partial) pairs — the form in which a
    source feature and a parsed one are comparable without knowing the strand."""
    return tuple(
        tuple(sorted(((part.start + 1, part.partial_start), (part.end, part.partial_end))))
        for part in gene.parts
    )


def test_the_golden_genome_round_trips_through_the_table(data_dir: Path) -> None:
    """P3-05's Done-when: every feature the source has, the table re-reads — same kind, same
    name, same intervals (partials included), same product and translation table. Compared as a
    multiset: both IR copies are features in their own right, not one feature written twice."""
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    parsed = parse_tbl(to_tbl(golden, "NC_050833.1"))

    def from_gene(gene: Gene) -> tuple:
        return (
            gene.kind,
            gene.name,
            gene.product,
            str(gene.transl_table) if gene.transl_table is not None else None,
            _canonical_parts(gene),
        )

    def from_tbl(feature) -> tuple:
        return (
            feature.kind,
            feature.qualifiers.get("gene"),
            feature.qualifiers.get("product"),
            feature.qualifiers.get("transl_table"),
            tuple(tuple(sorted(((start, ps), (end, pe)))) for start, end, ps, pe in feature.parts),
        )

    source = sorted(
        from_gene(gene) for gene in golden.genes if gene.kind in ("gene", "CDS", "tRNA", "rRNA")
    )
    assert len(source) == 265  # 133 gene + 86 CDS + 38 tRNA + 8 rRNA
    assert sorted(from_tbl(feature) for feature in parsed) == source


def test_the_golden_genome_passes_the_first_pass_checks(data_dir: Path) -> None:
    golden = read_genbank(data_dir / "NC_050833.1.gb")
    assert submission_problems(golden) == ()


def test_each_first_pass_check_fires_on_its_own_case() -> None:
    inside = Gene("rbcL", "CDS", (Interval(0, 30),), product="rubisco", transl_table=11)
    outside = Gene("psbA", "CDS", (Interval(0, 500),), product="D1", transl_table=11)
    no_table = Gene("psbB", "CDS", (Interval(50, 80),), product="CP47")
    no_product = Gene("psbC", "CDS", (Interval(80, 110),), transl_table=11)
    genome = _genome(inside, inside, outside, no_table, no_product, length=200)

    problems = submission_problems(genome)
    assert any("outside the sequence" in problem for problem in problems)
    assert any("no /transl_table" in problem and "psbB" in problem for problem in problems)
    assert any("no /product" in problem and "psbC" in problem for problem in problems)
    assert any("identical" in problem and "locus tags" in problem for problem in problems)


# --- the module -------------------------------------------------------------------------------


def _ctx() -> RunContext:
    return RunContext(layout=None, params={}, module_id="submission", module_order=1)  # type: ignore[arg-type]


def test_the_module_writes_both_files_from_the_normalized_set(data_dir: Path) -> None:
    from plastid.modules import submission

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    result = submission.run(golden, {}, _ctx())

    assert set(result.artifacts) == {"submission.fsa", "submission.tbl"}
    assert result.artifacts["submission.fsa"].startswith(">NC_050833.1\n")
    assert result.artifacts["submission.tbl"].startswith(">Features NC_050833.1\n")
    summary = result.tables["submission_summary"].rows[0]
    assert summary["n_problems"] == 0
    assert summary["sqn"] == "not requested"
    assert summary["n_features"] == 265  # 133 gene + 86 CDS + 38 tRNA + 8 rRNA


def test_the_sqn_step_is_a_named_error_when_table2asn_is_absent(data_dir: Path) -> None:
    """D3-04: the .sqn is opt-in, and asking for it without the tool is exit 3, not a silent
    skip — the .fsa/.tbl half stays complete and is what NCBI's web submission takes."""
    from plastid.modules import submission
    from plastid.run.tools import ToolError, which

    if which("table2asn") is not None:
        pytest.skip("table2asn is installed here — this test is about the machine that has none")

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    with pytest.raises(ToolError, match="table2asn not found"):
        submission.run(golden, {"sqn": "true"}, _ctx())


def test_with_the_tool_present_the_missing_template_is_what_stops_the_sqn(
    data_dir: Path, tmp_path: Path
) -> None:
    """This machine has table2asn (installed 2026-10-03, declared in docs/tools/table2asn.md),
    so the path past "tool not found" is testable — and what stops it is the template, which is
    exactly the file only the team can produce."""
    from plastid.modules import submission
    from plastid.run.tools import ToolError, which

    binary = which("table2asn") or (
        "/root/toolbox/bin/table2asn" if Path("/root/toolbox/bin/table2asn").is_file() else None
    )
    if binary is None:
        pytest.skip("table2asn is not installed here — nothing to test past the tool check")

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    with pytest.raises(ToolError, match="needs a submitter template"):
        submission.run(golden, {"sqn": "true", "table2asn": binary}, _ctx())
    with pytest.raises(ToolError, match="does not exist"):
        submission.run(
            golden,
            {"sqn": "true", "table2asn": binary, "template": str(tmp_path / "none.sbt")},
            _ctx(),
        )


def test_the_cli_run_writes_the_files_and_the_manifest_lists_them(
    tmp_path: Path, data_dir: Path
) -> None:
    import json

    from plastid.cli import main

    golden = data_dir / "NC_050833.1.gb"
    assert (
        main(
            [
                "run",
                "-m",
                "submission",
                "-i",
                str(golden),
                "-o",
                str(tmp_path),
                "--dataset",
                "iris",
            ]
        )
        == 0
    )
    run_dir = next((tmp_path / "iris").iterdir())
    module_dir = run_dir / "01_submission"
    assert (module_dir / "submission.fsa").is_file()
    assert (module_dir / "submission.tbl").is_file()

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    listed = {entry["path"] for entry in manifest["files"]}
    assert "01_submission/submission.fsa" in listed
    assert "01_submission/submission.tbl" in listed

    # a second identical run is a cache hit, and the artifacts come back byte-identical
    before = (module_dir / "submission.tbl").read_bytes()
    assert (
        main(
            [
                "run",
                "-m",
                "submission",
                "-i",
                str(golden),
                "-o",
                str(tmp_path),
                "--dataset",
                "iris",
            ]
        )
        == 0
    )
    second = sorted((tmp_path / "iris").iterdir())[-1]
    assert (second / "01_submission" / "submission.tbl").read_bytes() == before
