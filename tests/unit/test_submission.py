"""NCBI submission files (P3-05): the five-column rules, partials, and the round trip.

The format rules here are not invented — they are NCBI's, and each one fails at *upload* time
if it is wrong, after the work is done. The two examples below are NCBI's own (BankIt feature
table help, retrieved 2026-10-03), re-tabulated: the help page renders them with the tab
columns collapsed, and the documented column definition says a qualifier name belongs in
column 4 and its value in column 5.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Gene, Genome, Interval
from plastid.core.submission import (
    parse_tbl,
    partial_start_genes,
    source_qualifiers_left_out,
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
    "<1\t>1050\tgene\n"
    "\t\t\tgene\tATH1\n"
    "<1\t1009\tCDS\n"
    "\t\t\tcodon_start\t2\n"
    "<1\t>1050\tmRNA\n"
)


def _genome(*genes: Gene, length: int = 400) -> Genome:
    """A synthetic genome whose CDS all start with ATG.

    The submission writer marks a CDS whose first codon is not an initiator as partial — NCBI's
    own rule, measured against its validator — and these tests are about the *format*, so
    the sequence is built to have a real start codon. The tests that are about that rule build
    their own genome instead.
    """
    sequence = bytearray(b"A" * length)
    for gene in genes:
        if gene.kind != "CDS" or gene.parts[0].end > length:
            continue
        first = gene.parts[0]
        if first.strand < 0:
            sequence[first.end - 3 : first.end] = b"CAT"  # read back on the minus strand: ATG
        else:
            sequence[first.start : first.start + 3] = b"ATG"
    return Genome(
        name="synthetic",
        accession="SYN001.1",
        sequence=bytes(sequence).decode(),
        genes=genes,
        topology="linear",
    )


def test_a_multi_interval_feature_is_stacked_not_joined() -> None:
    """There is no ``join()`` keyword in this format: intervals stack, key on the first line."""
    gene = Gene("trnK-UUU", "tRNA", (Interval(100, 200), Interval(300, 350)))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[0] == ">Features SYN001.1"
    assert lines[1] == "101\t200\ttRNA"
    assert lines[2] == "301\t350\t"
    assert lines[3] == "\t\t\tgene\ttrnK-UUU"


def test_a_minus_strand_feature_writes_its_interval_reversed() -> None:
    gene = Gene("rbcL", "CDS", (Interval(100, 200, -1),), product="rubisco")
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "200\t101\tCDS"
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
    assert lines[1] == "120\t101\tCDS"
    assert lines[2] == "301\t320\t"
    assert "\t\t\ttrans_splicing\t" in lines


def test_a_partial_coordinate_keeps_its_symbol_in_the_right_column() -> None:
    """``<`` always in column 1, ``>`` always in column 2 — NCBI's rule, whatever the strand."""
    gene = Gene("ycf1", "CDS", (Interval(0, 100, 1, partial_start=True, partial_end=True),))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "<1\t>100\tCDS"


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
    assert lines[1] == "100\t>1\tCDS"


def test_a_partial_five_prime_end_of_a_minus_strand_feature_takes_the_other_symbol() -> None:
    """The mirror case: `complement(1..>100)` — coordinate 100 is the feature's 5' end."""
    gene = Gene("ndhF", "CDS", (Interval(0, 100, -1, partial_start=False, partial_end=True),))
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert lines[1] == "<100\t1\tCDS"


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
    joined = _features(parse_tbl(NCBI_JOIN_EXAMPLE))
    assert len(joined) == 1
    assert joined[0].kind == "tRNA"
    assert joined[0].parts == ((2626, 2590, False, False), (2570, 2535, False, False))
    assert joined[0].qualifiers == {"product": "tRNA-Phe"}

    partial = _features(parse_tbl(NCBI_PARTIAL_EXAMPLE))
    assert [feature.kind for feature in partial] == ["gene", "CDS", "mRNA"]
    assert partial[0].parts == ((1, 1050, True, True),)
    assert partial[0].qualifiers == {"gene": "ATH1"}
    assert partial[1].qualifiers == {"codon_start": "2"}


def test_the_parser_reads_the_compact_qualifier_layout_too() -> None:
    """The BankIt help page renders qualifiers with the leading tabs collapsed; a table pasted
    from it still has to read."""
    text = ">Feature Seq2\n2626\t2590\ttRNA\t\t\n2570\t2535\t\t\t\nproduct\ttRNA-Phe\n"
    parsed = _features(parse_tbl(text))
    assert parsed[0].qualifiers == {"product": "tRNA-Phe"}


def _features(records) -> tuple:
    """The one record a table is expected to hold here — the multi-record case has its own test."""
    assert len(records) == 1, f"expected one record, got {len(records)}"
    return records[0].features


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
    parsed = _features(parse_tbl(to_tbl(golden, "NC_050833.1")))

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

    # The writer adds two things the gene model does not hold, both measured against NCBI's own
    # validator: the record's source feature, and a partial start on the CDS whose first
    # codon is not an initiator (rpl2, both copies — ACG, edited to AUG in the transcript).
    # Everything else must come back identical, which is what the two assertions below say.
    marked = {(gene.kind, gene.name, *gene.span) for gene in partial_start_genes(golden)}
    assert marked  # the golden really has one, so the excepted path is exercised, not assumed

    def identity(feature) -> tuple:
        coordinates = [value for part in feature.parts for value in (part[0], part[1])]
        return (
            feature.kind,
            feature.qualifiers.get("gene"),
            min(coordinates),
            max(coordinates),
        )

    untouched = [
        feature
        for feature in parsed
        if feature.kind != "source" and identity(feature) not in marked
    ]
    untouched_source = sorted(
        from_gene(gene)
        for gene in golden.genes
        if gene.kind in ("gene", "CDS", "tRNA", "rRNA")
        and (gene.kind, gene.name, *gene.span) not in marked
    )
    assert len(untouched_source) == 263  # the 265 features minus the two marked rpl2 copies
    assert sorted(from_tbl(feature) for feature in untouched) == untouched_source

    assert len(parsed) == 266  # those 265 features, plus the record's source feature
    assert sorted(from_tbl(f) for f in parsed if f.kind == "source") == [
        ("source", None, None, None, (((1, False), (153724, False)),))
    ]

    # the marked CDS are exactly the ones carrying the marker, and nothing else moved
    assert {identity(f) for f in parsed if f.kind == "CDS" and _is_partial(f)} == marked


def _is_partial(feature) -> bool:
    """True when any interval of the feature carries a partial marker."""
    return any(partial_start or partial_end for _, _, partial_start, partial_end in feature.parts)


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
    assert summary["validator_status"] == "not requested"
    # an unrun validator is not zero errors: the count is absent, and the status says why
    assert summary["n_validator_errors"] is None
    assert summary["n_validator_warnings"] is None
    assert summary["n_features"] == 265  # 133 gene + 86 CDS + 38 tRNA + 8 rRNA


def test_a_record_without_features_is_refused_not_written_as_an_empty_table() -> None:
    """A submission is built from a checked GenBank, not from a bare sequence (M11.18).

    Writing an empty feature table from a FASTA would turn "this sequence was never annotated"
    into "this sequence has no genes" — a different claim, and one NCBI would receive.
    """
    from plastid.core.submission import SubmissionInputError
    from plastid.modules import submission

    with pytest.raises(SubmissionInputError, match="no gene, CDS, tRNA or rRNA"):
        submission.run(_genome(), {}, _ctx())


def test_the_fsa_and_tbl_carry_the_same_sequence_id(data_dir: Path) -> None:
    """table2asn pairs the two files by id; a mismatch would attach one sample's table to
    another sample's sequence (M11.18)."""
    from plastid.modules import submission

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    result = submission.run(golden, {}, _ctx())
    header = result.artifacts["submission.fsa"].splitlines()[0]
    assert header.startswith(">")
    (record,) = parse_tbl(result.artifacts["submission.tbl"])
    assert record.seq_id == header[1:]
    assert result.tables["submission_summary"].rows[0]["sequence_id"] == record.seq_id


def test_the_sqn_step_is_a_named_error_when_table2asn_is_absent(
    data_dir: Path, submission_tools: dict[str, str | None]
) -> None:
    """D3-04: the .sqn is opt-in, and asking for it without the tool is exit 3, not a silent
    skip — the .fsa/.tbl half stays complete and is what NCBI's web submission takes."""
    from plastid.modules import submission
    from plastid.run.tools import ToolError

    if submission_tools["table2asn"] is not None:
        pytest.skip("table2asn is installed here — this test is about the machine that has none")

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    with pytest.raises(ToolError, match="table2asn not found"):
        submission.run(golden, {"sqn": "true"}, _ctx())


def test_with_the_tool_present_the_missing_template_is_what_stops_the_sqn(
    data_dir: Path, tmp_path: Path, submission_tools: dict[str, str | None]
) -> None:
    """The machine that has table2asn (docs/tools/table2asn.md) can be tested past "tool not
    found" — and what stops it is the template, the file only the team can produce.

    The tool is found the way every other tool-dependent test finds one: through
    ``submission_tools``, which honours ``PTS_TABLE2ASN``. It used to fall back to a hardcoded
    ``/root/toolbox/bin/table2asn`` — a banned hardcoded path, and on a machine where that path
    is present but unreadable (the CI runner, whose user cannot traverse ``/root``) the check
    raised ``PermissionError`` instead of skipping the test. Measured 10/10/2026, P8-05.
    """
    from plastid.modules import submission
    from plastid.run.tools import ToolError

    binary = submission_tools["table2asn"]
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


REAL_NCBI_TABLES = (
    "real_NC_039897.vadr.tbl",
    "real_va-dengue-rpn.vadr.tbl",
    "real_va-noro.vadr.tbl",
)


@pytest.mark.parametrize("name", REAL_NCBI_TABLES)
def test_a_real_ncbi_table_round_trips_through_our_writer(name: str, data_dir: Path) -> None:
    """P3-05's second clause, as far as it can be closed here: the table we write has the same
    shape as tables NCBI's own pipeline produces, field for field.

    These files are not examples from a manual — they are VADR's output for real records
    (provenance in tests/data/PROVENANCE.md). Reading one back and writing it out again must
    reproduce it exactly, record by record: every convention that fails at *upload* time is in
    these bytes, and a writer that disagrees with any of them disagrees with a real table.
    """
    from plastid.core.submission import tbl_text

    text = (data_dir / name).read_text(encoding="utf-8")
    records = parse_tbl(text)
    assert records, "the real table parsed to nothing"

    written = "".join(tbl_text(record.features, record.seq_id) for record in records)
    original_lines = text.rstrip("\n").splitlines()
    written_lines = written.rstrip("\n").splitlines()

    # The header spelling differs between sources: NCBI's documentation says ">Features" (which
    # this writer emits) and VADR's own output says ">Feature". The parser accepts both, because
    # both are in the wild; everything after each header must match byte for byte.
    original_body = [line for line in original_lines if not line.startswith(">")]
    written_body = [line for line in written_lines if not line.startswith(">")]
    assert written_body == original_body
    assert all(line.startswith(">Features ") for line in written_lines if line.startswith(">"))


def test_a_multi_record_table_keeps_its_records_apart(data_dir: Path) -> None:
    """Measured while vendoring these files: NCBI's own output can hold several records in one
    table (one VADR file carried eight). Merging them would attach every record's features to one
    sequence — silently, and to the wrong sequence. So a record boundary is part of the model."""
    noro = (data_dir / "real_va-noro.vadr.tbl").read_text(encoding="utf-8")
    records = parse_tbl(noro)
    assert len(records) == 8
    assert [record.seq_id for record in records][:2] == ["KY887602.1", "KT818729.1"]
    assert all(record.features for record in records)


def test_the_conventions_the_real_tables_pin_down(data_dir: Path) -> None:
    """Spelled out, so the reason these files are in the repository survives a refactor."""
    partial = (data_dir / "real_va-dengue-rpn.vadr.tbl").read_text(encoding="utf-8")
    # a partial marker is a prefix, in column 1 for the start and column 2 for the stop
    assert "<1\t>515\tCDS" in partial
    # a qualifier name sits in column 4, its value in column 5 — three empty fields before it
    assert "\t\t\tcodon_start\t2" in partial
    assert "\t\t\tproduct\tpolyprotein" in partial

    plain = (data_dir / "real_NC_039897.vadr.tbl").read_text(encoding="utf-8")
    assert plain.startswith(">Feature ref|NC_039897.1|\n")
    assert "\t\t\tgene\tORF1" in plain
    # features are written in position order, and only the first interval of a feature carries
    # its key — the shape our writer emits for a multi-interval gene
    starts = [int(line.split("\t")[0]) for line in plain.splitlines()[1:] if line[:1].isdigit()]
    assert starts == sorted(starts)


# --- what a submission needs that a gene model does not hold (measured 2026-10-03) ------
#
# All three below came out of running table2asn on the golden genome's table with a real
# submission template: 13 ERROR-level messages, of which these rules removed four. The remaining
# ones are inherent to the format and are.


def test_the_source_feature_comes_first_with_the_records_own_qualifiers() -> None:
    """Without it NCBI's validator reports ``NoSourceDescriptor`` and ``BioSourceMissing``: a
    submission has to say what the sequence is, and a toolkit that invented an organism would be
    worse than one that refuses."""
    gene = Gene("psbA", "CDS", (Interval(100, 200, 1),), product="D1")
    genome = replace(
        _genome(gene),
        source={
            "organism": "Iris domestica",
            "organelle": "plastid:chloroplast",
            "mol_type": "genomic DNA",
            "specimen_voucher": "HN-2024-001",
        },
    )
    lines = to_tbl(genome, "SYN001.1").splitlines()
    assert lines[1] == "1\t400\tsource"
    assert "\t\t\torganism\tIris domestica" in lines
    assert "\t\t\torganelle\tplastid:chloroplast" in lines
    # the qualifiers a real sample carries travel too — they are the submitter's statement
    assert "\t\t\tspecimen_voucher\tHN-2024-001" in lines
    # and the source is one feature, not one per gene
    assert lines.count("1\t400\tsource") == 1


def test_a_record_that_states_no_organism_is_a_named_problem() -> None:
    gene = Gene("psbA", "CDS", (Interval(100, 200, 1),), product="D1", transl_table=11)
    problems = submission_problems(_genome(gene))
    assert any("no source organism" in problem for problem in problems)


def test_an_exception_travels_with_its_feature() -> None:
    """``/exception`` is the record's own word about a gene whose translation is unusual. Dropping
    it is not neutral: measured, NCBI's validator rejects an RNA-edited start without it."""
    gene = Gene(
        "rpl2",
        "CDS",
        (Interval(100, 200, 1),),
        product="ribosomal protein L2",
        transl_table=11,
        exception="RNA editing",
    )
    lines = to_tbl(_genome(gene), "SYN001.1").splitlines()
    assert "\t\t\texception\tRNA editing" in lines


def _gene_with_start(sequence: str, *, strand: int = 1) -> Genome:
    """A genome whose one CDS sits at 101-200 and whose first codon the test chooses."""
    codon = sequence if strand > 0 else _reverse_complement(sequence)
    body = "A" * 100 + codon + "A" * (300 - 103)
    gene = Gene(
        "psbA",
        "CDS",
        (Interval(100, 200, strand),),
        product="D1",
        transl_table=11,
    )
    return Genome(
        name="synthetic", accession="SYN001.1", sequence=body, genes=(gene,), topology="linear"
    )


def _reverse_complement(codon: str) -> str:
    return codon.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.parametrize("strand", [1, -1])
def test_a_start_codon_the_table_does_not_initiate_is_written_partial(strand: int) -> None:
    """Measured 2026-10-03: ACG — the RNA-edited start of rpl2 — makes table2asn report both
    ``StartCodon`` and ``BadProteinStart``. With a partial start the same CDS validates, which is
    also what the published *Dioscorea alata* plastome submission does for its rpl2 and ndhD."""
    genome = _gene_with_start("ACG", strand=strand)
    assert [gene.name for gene in partial_start_genes(genome)] == ["psbA"]
    lines = to_tbl(genome, "SYN001.1").splitlines()
    if strand > 0:
        assert lines[1] == "<101\t200\tCDS"
    else:
        # a minus-strand feature is written reversed, so its 5' end is in column 1
        assert lines[1] == "<200\t101\tCDS"


def test_a_canonical_start_is_left_exactly_as_the_record_has_it() -> None:
    genome = _gene_with_start("ATG")
    assert partial_start_genes(genome) == ()
    assert to_tbl(genome, "SYN001.1").splitlines()[1] == "101\t200\tCDS"


def test_a_cds_the_toolkit_cannot_translate_is_not_guessed_at() -> None:
    """A table we do not know, or a shifted frame, is not evidence of an edited start — the
    feature is left alone rather than marked on a guess."""
    gene = Gene("psbA", "CDS", (Interval(100, 200, 1),), product="D1", transl_table=4)
    genome = Genome(
        name="synthetic",
        accession="SYN001.1",
        sequence="A" * 400,
        genes=(gene,),
        topology="linear",
    )
    assert partial_start_genes(genome) == ()
    assert to_tbl(genome, "SYN001.1").splitlines()[1] == "101\t200\tCDS"


def test_the_two_source_qualifiers_this_format_has_no_place_for_are_left_out() -> None:
    """Measured 2026-10-03: table2asn calls ``mol_type`` and ``db_xref`` "Unrecognized qualifier
    name" and exits non-zero while still writing a correct ``.sqn`` — the real published plastome
    table carries ``mol_type`` and is noisier for it. Neither is the submitter's to state: NCBI
    derives the molecule type from the submission and the taxon from the organism. They stay on
    the genome, and the list of what was left behind is available rather than silent."""
    gene = Gene("psbA", "CDS", (Interval(100, 200, 1),), product="D1", transl_table=11)
    genome = replace(
        _genome(gene),
        source={
            "organism": "Iris domestica",
            "organelle": "plastid:chloroplast",
            "mol_type": "genomic DNA",
            "db_xref": "taxon:58944",
            "cultivar": "Hong",
        },
    )
    lines = to_tbl(genome, "SYN001.1").splitlines()
    assert "\t\t\torganism\tIris domestica" in lines
    assert "\t\t\tcultivar\tHong" in lines  # a sample qualifier, and the submitter's own to state
    assert not [line for line in lines if "mol_type" in line or "db_xref" in line]
    assert source_qualifiers_left_out(genome) == ("mol_type", "db_xref")
    assert genome.source["mol_type"] == "genomic DNA"  # left out of the table, not off the record


# --- what NCBI's own validator said about the table (measured 2026-10-03) ---------------

_VALIDATOR_STATS = """Total messages:\t\t104

=================================================================
95 WARNING-level messages exist

SEQ_FEAT.GenCodeMismatch:\t86
SEQ_FEAT.CDSgeneRange:\t8
SEQ_FEAT.UnnecessaryException:\t1

=================================================================
9 ERROR-level messages exist

SEQ_INST.BadProteinStart:\t2
SEQ_FEAT.GeneXrefStrandProblem:\t7
"""


def test_the_validator_report_is_read_not_just_kept() -> None:
    """table2asn's ``.stats`` is a file, and a file nobody reads is not a check — the counts go
    into the table, and the per-code breakdown says which class each finding belongs to."""
    from plastid.modules.submission import validator_counts

    counts = validator_counts(_VALIDATOR_STATS)
    assert (counts.errors, counts.warnings) == (9, 95)
    # errors first, and the biggest class first within a level — the order a reader wants
    assert counts.by_class == (
        ("SEQ_FEAT.GeneXrefStrandProblem", "error", 7),
        ("SEQ_INST.BadProteinStart", "error", 2),
        ("SEQ_FEAT.GenCodeMismatch", "warning", 86),
        ("SEQ_FEAT.CDSgeneRange", "warning", 8),
        ("SEQ_FEAT.UnnecessaryException", "warning", 1),
    )


def test_a_stats_file_with_no_findings_is_not_an_error() -> None:
    from plastid.modules.submission import validator_counts

    counts = validator_counts("Total messages:\t\t0\n")
    assert (counts.errors, counts.warnings, counts.by_class) == (0, 0, ())
