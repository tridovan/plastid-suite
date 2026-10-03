"""Tier B golden checks — need external tools; a missing tool is named, never silent.

Run on a tool-equipped machine with ``PTS_REQUIRE_TOOLS=1`` to turn skips into failures.
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.genbank import read_genbank
from plastid.run.manifest import validate
from plastid.run.tools import which

pytestmark = pytest.mark.tier_b


def test_tool_probe_reports_what_is_available(tools: dict[str, str | None]) -> None:
    assert set(tools) == {"perl", "misa.pl", "mafft"}


def test_ssr_reproduces_the_golden_numbers(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    require_tools(tools, "perl", "misa.pl")
    assert (
        main(
            [
                "run",
                "--module",
                "ssr",
                "--input",
                str(data_dir / "NC_050833.1.gb"),
                "--out",
                str(tmp_path),
                "--dataset",
                "golden",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "golden").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)

    tool_records = manifest["modules"][0]["tools"]
    assert [record["name"] for record in tool_records] == ["misa"]
    record = tool_records[0]
    assert record["version"] == "2.1"  # pinned: the vendored MISA
    assert record["interpreter"].startswith("perl ")
    assert record["rc"] == 0
    assert Path(record["stdout_file"]).is_file()

    summary_csv = (run_dir / "01_ssr" / "tables" / "ssr_summary.csv").read_text(encoding="utf-8")
    rows = csv.reader(
        io.StringIO(
            "\n".join(line for line in summary_csv.splitlines() if not line.startswith("#"))
        )
    )
    header = next(rows)
    values = dict(zip(header, next(rows), strict=True))
    assert values["Total SSRs"] == "37"  # PROJECT.md §6 golden number
    assert values["Monomeric"] == "31"
    assert values["Dinucleotide"] == "4"
    assert values["Trinucleotide"] == "2"
    assert values["In compound formation"] == "3"
    # P2-06: the per-SSR table must reconcile with MISA's own statistics. They disagreed once
    # (30 vs 37) until the parser learned MISA's second vocabulary — see tests/unit/test_misa.py.
    assert values["GFF vs statistics"].startswith("agrees")

    per_ssr = (run_dir / "01_ssr" / "tables" / "ssrs.csv").read_text(encoding="utf-8")
    ssr_rows = [line for line in per_ssr.splitlines() if line.strip() and not line.startswith("#")][
        1:
    ]
    assert len(ssr_rows) == 37  # every SSR is listed, compounds included component by component

    by_region = (run_dir / "01_ssr" / "tables" / "ssr_by_region.csv").read_text(encoding="utf-8")
    region_rows = dict(
        (fields[0], fields[1:])
        for fields in (
            line.split(",")
            for line in by_region.splitlines()
            if line.strip() and not line.startswith("#")
        )
    )
    assert set(region_rows) == {"Region", "LSC", "IRa", "SSC", "IRb"}
    region_totals = sum(int(region_rows[name][0]) for name in ("LSC", "IRa", "SSC", "IRb"))
    assert region_totals == 37  # the regions account for every SSR, none lost

    report = (run_dir / "01_ssr" / "report.md").read_text(encoding="utf-8")
    assert "**37 SSRs**" in report
    assert "Draft — needs human review." in report


def test_the_cds_scope_reproduces_the_golden_ts_tv(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """PROJECT.md §6, measured by the toolkit's own path: **2.686 = 317 Ts / 118 Tv** over the
    **76** shared CDS and **66,066** alignment columns, IR excluded.

    This check is why the module's CDS scope exists at all — §6's number is *not* the
    whole-alignment one (measured 03/10/2026: 1.8416 over 152,338 comparable columns), and the two
    must never be quoted as one another.

    It is also the check whose absence hid a real bug: MAFFT 7.526 prints its alignment in
    lowercase, the counters read ``ACGT``, and the first real-data run reported **zero comparable
    columns and a blank ratio** — the unit tests had only ever built uppercase alignments by hand.
    """
    require_tools(tools, "mafft")
    assert (
        main(
            [
                "run",
                "-m",
                "substitutions",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
                "--set",
                "scope=cds",
            ]
        )
        == 0
    )
    summary = _table(tmp_path / "runs", "substitutions_summary")
    row = next(summary)
    assert row["Shared genes"] == "76"
    assert row["Alignment columns (bp)"] == "66066"
    assert (row["Transitions"], row["Transversions"]) == ("317", "118")
    assert row["Ts/Tv"] == "2.6864"


def test_the_whole_alignment_scope_is_a_different_number_and_says_so(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """The default scope is not wrong — it is *wider*, and its label has to keep it apart from
    §6's. Measured 03/10/2026: 1.8416 over the whole aligned plastome."""
    require_tools(tools, "mafft")
    assert (
        main(
            [
                "run",
                "-m",
                "substitutions",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
            ]
        )
        == 0
    )
    summary = _table(tmp_path / "runs", "substitutions_summary")
    row = next(summary)
    assert row["Ts/Tv"] == "1.8416"
    assert "not the CDS-only scope" in _scope_label(tmp_path / "runs", "substitutions_summary")


def test_nucleotide_diversity_runs_on_real_data_and_keeps_its_shape(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P2-07's real-data run, measured 03/10/2026 with MAFFT 7.526: π 0.007181 over the whole
    alignment, and the region split the plastome predicts — the **inverted repeat is the most
    conserved part and its two copies diverge by the same amount** (IRa 0.001682 · IRb 0.001681),
    while the single-copy regions are several times more diverse (LSC 0.009185 · SSC 0.014066).

    There is no §6 golden value for π, so this check holds the *relationships* — which is what
    makes it worth having: a regression that quietly stopped counting columns (the lowercase bug
    of `external/mafft.py`) would make every number zero and every relationship vanish.
    """
    require_tools(tools, "mafft")
    assert (
        main(
            [
                "run",
                "-m",
                "diversity",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
            ]
        )
        == 0
    )
    summary = next(_table(tmp_path / "runs", "diversity_summary"))
    assert float(summary["π (whole alignment)"]) == pytest.approx(0.007181, rel=1e-3)
    assert summary["Usable columns"] != "0"

    regions = {
        row["Region"]: float(row["π"]) for row in _table(tmp_path / "runs", "diversity_regions")
    }
    assert set(regions) == {"LSC", "IRa", "SSC", "IRb"}
    assert regions["IRa"] == pytest.approx(regions["IRb"], rel=0.05), "the IR copies must mirror"
    assert regions["IRa"] < regions["LSC"] < regions["SSC"]


def _table(root: Path, table_id: str) -> csv.DictReader:
    """A module table as rows: the exports carry their provenance as ``#`` comments."""
    path = next(root.rglob(f"{table_id}.csv"))
    body = [
        line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")
    ]
    return csv.DictReader(io.StringIO("\n".join(body)))


def _scope_label(root: Path, table_id: str) -> str:
    path = next(root.rglob(f"{table_id}.csv"))
    return next(line for line in path.read_text(encoding="utf-8").splitlines() if "scope:" in line)


def test_pga_reproduces_the_golden_gene_counts(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P3-03's Done-when, on a machine that has PGA: annotate a de-annotated golden FASTA and
    check the gene counts by the fixed rule. Nothing is claimed here — the test skips with a
    named reason when PGA is absent, and ``PTS_REQUIRE_TOOLS=1`` turns that into a failure on a
    machine that is supposed to have it."""
    # skips with the tool list printed, or fails under PTS_REQUIRE_TOOLS on a machine that is
    # supposed to have PGA
    require_tools({**tools, "PGA.pl": tools.get("PGA.pl") or which("PGA.pl")}, "perl", "PGA.pl")

    from plastid.core.model import count_genes

    golden_path = data_dir / "NC_050833.1.gb"
    golden = read_genbank(golden_path)
    # a de-annotated target: the sequence only, so the annotator has to find the genes itself
    target = tmp_path / "target.gb"
    target.write_text(_strip_features(golden_path.read_text(encoding="utf-8")), encoding="utf-8")

    reference = tmp_path / "refs"
    reference.mkdir()
    (reference / golden_path.name).write_bytes(golden_path.read_bytes())

    assert (
        main(
            [
                "run",
                "--module",
                "annotation",
                "--input",
                str(target),
                "--out",
                str(tmp_path / "runs"),
                "--dataset",
                "pga",
                "--set",
                f"reference={reference}",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "runs" / "pga").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)
    assert manifest["modules"][0]["tools"][0]["name"] == "pga"

    counts = count_genes(read_genbank(run_dir / "01_annotation" / "annotation.gb"))
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)  # the fixed rule, golden value
    assert counts.total == 132
    assert len(golden.genes) > 0  # the source really was annotated


def test_chloe_reproduces_the_golden_gene_counts(
    annotator_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P3-03's Done-when, on a machine that has Chloe: annotate the golden genome's own sequence
    and check the gene counts by the fixed rule. Skips with a named reason when the tool is
    absent; ``PTS_REQUIRE_TOOLS=1`` turns that into a failure on a machine supposed to have it.

    The last assertion is the one that is not a gene count: the run must have passed
    ``--no-transform``. Without it Chloe annotates a *flipped* sequence and returns the same 132
    genes in another frame — measured 2026-10-03, which is why the flag is asserted here rather
    than trusted to stay in the adapter.
    """
    require_tools(annotator_tools, "chloe", "julia")
    golden_path = data_dir / "NC_050833.1.gb"
    target = tmp_path / "target.gb"
    target.write_text(_strip_features(golden_path.read_text(encoding="utf-8")), encoding="utf-8")

    assert (
        main(
            [
                "run",
                "--module",
                "annotation",
                "--input",
                str(target),
                "--out",
                str(tmp_path / "runs"),
                "--dataset",
                "chloe",
                "--set",
                f"chloe={annotator_tools['chloe']}",
                "--set",
                f"julia={annotator_tools['julia']}",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "runs" / "chloe").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)

    record = manifest["modules"][0]["tools"][0]
    assert record["name"] == "chloe"
    assert record["version"] == "0.1.18"  # pinned to the release the numbers were measured with
    assert record["interpreter"].startswith("julia ")
    assert "--no-transform" in record["argv"]
    assert record["rc"] == 0

    from plastid.core.gff import genes_from_gff
    from plastid.core.model import Genome, count_genes

    golden = read_genbank(golden_path)
    annotated = Genome(
        name=golden.name,
        sequence=golden.sequence,
        genes=genes_from_gff((run_dir / "01_annotation" / "annotation.gff").read_text("utf-8")),
    )
    counts = count_genes(annotated)
    assert (counts.cds, counts.trna, counts.rrna) == (86, 38, 8)  # the fixed rule, golden value
    assert counts.total == 132

    from plastid.core.orf import check_genome

    checks = check_genome(annotated, default_table=11)
    assert not [check for check in checks if check.is_failure]


def test_chloe_transforms_the_sequence_unless_it_is_told_not_to(
    annotator_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """The measurement behind ``--no-transform``, against the real tool.

    Run without the flag, Chloe writes a FASTA beside its GFF — the sequence it actually
    annotated, with the single-copy region flipped — and the gene set is *identical*: same 132
    genes, same names. A gene count therefore cannot see it, which is the whole reason this is a
    test rather than a comment in the adapter. What changes is the frame: read against the
    sequence the caller handed over, the transformed GFF leaves a third of the CDS without a
    valid ORF (measured 2026-10-03: 12 of 86 failing, 8 of them `no_start`).
    """
    require_tools(annotator_tools, "chloe", "julia")

    from plastid.core.gff import genes_from_gff
    from plastid.core.model import Genome, count_genes
    from plastid.core.orf import check_genome
    from plastid.core.submission import to_fasta

    golden = read_genbank(data_dir / "NC_050833.1.gb")
    chloe = _script(annotator_tools["chloe"])
    workdir = tmp_path / "transformed"
    workdir.mkdir()
    fasta = workdir / "NC_050833.1.fa"
    fasta.write_text(to_fasta(golden, "NC_050833.1"), encoding="utf-8")

    ctx = _ctx(tmp_path)
    ctx.tool(
        "chloe",
        [
            str(annotator_tools["julia"]),
            f"--project={chloe.parent}",
            str(chloe),
            "annotate",
            str(fasta),
        ],
        kind="binary",
        tool_path=str(chloe),
        pinned_version="0.1.18",
        cwd=workdir,
        timeout=3600,
    )

    written = (workdir / "NC_050833.1.chloe.fa").read_text(encoding="utf-8")
    sequence = "".join(
        line.strip() for line in written.splitlines() if not line.startswith(">")
    ).upper()
    assert sequence != golden.sequence.upper()  # the sequence the GFF's coordinates belong to

    genes = genes_from_gff((workdir / "NC_050833.1.chloe.gff").read_text(encoding="utf-8"))
    annotated = Genome(name=golden.name, sequence=golden.sequence, genes=genes)
    assert count_genes(annotated).total == 132  # the count survives the transform …
    failures = [check for check in check_genome(annotated, default_table=11) if check.is_failure]
    assert failures, "… but the frame does not, so this is what a gene count cannot see"


def _script(value: str) -> Path:
    """``--set chloe=`` takes the clone or the script in it; the tests take both too."""
    path = Path(value)
    return path / "chloe.jl" if path.is_dir() else path


def _ctx(tmp_path: Path):
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    layout = RunLayout.create(tmp_path / "ctx", "chloe-tier-b", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="annotation", module_order=1)


def _strip_features(genbank_text: str) -> str:
    """Blank every feature line, keeping the sequence and the record's own header."""
    lines = genbank_text.splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("FEATURES"))
    end = next(index for index, line in enumerate(lines) if line.startswith("ORIGIN"))
    kept = [*lines[: start + 1], f"     source          1..{_length(lines)}", "", *lines[end:]]
    return "\n".join(kept) + "\n"


def _length(lines: list[str]) -> int:
    """The record's length, from its LOCUS line — the only number left after stripping.

    Read as the token before ``bp``: the field after ``LOCUS`` is the record's *name*, so a
    published accession (``LOCUS  NC_050833  153724 bp``) and a hand-made one
    (``LOCUS  SYN_ANNOT  600 bp``) do not sit in the same column. Measured 2026-10-03 — the
    first version of this helper read position 1, worked on the hand-made record it was written
    against, and raised on the golden one the moment a second annotator's test used it.
    """
    locus = next(line for line in lines if line.startswith("LOCUS")).split()
    return int(locus[locus.index("bp") - 1])


def test_the_sqn_step_runs_and_ncbi_reads_the_table(
    submission_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P3-05's ``.sqn`` half, on a machine that has table2asn: the module builds the pair, the
    tool converts it, and the validator's own verdict is in the run's tables rather than in a file
    nobody opens.

    The template is the repository's **stand-in** — a real submission uses the team's own, which
    NCBI's form generates from cô Trang's details (``tests/data/PROVENANCE.md``). What is checked
    here is the pipeline: the tool exits 0, the ``.sqn`` exists, and every error NCBI reports
    belongs to the one class a duplicated inverted-repeat plastome produces (measured 2026-10-03
    against a published plastome submission, PROJECT.md §8.18) — a new class would fail this test,
    which is the point of asserting it.
    """
    require_tools(submission_tools, "table2asn")

    assert (
        main(
            [
                "run",
                "--module",
                "submission",
                "--input",
                str(data_dir / "NC_050833.1.gb"),
                "--out",
                str(tmp_path / "runs"),
                "--dataset",
                "sqn",
                "--set",
                "sqn=true",
                "--set",
                f"template={data_dir / 'standin_template.sbt'}",
                "--set",
                f"table2asn={submission_tools['table2asn']}",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "runs" / "sqn").iterdir())[0]
    module_dir = run_dir / "01_submission"
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)
    record = manifest["modules"][0]["tools"][0]
    assert record["name"] == "table2asn"
    assert record["rc"] == 0  # the tool accepted the table pair without complaint

    sqn = (module_dir / "submission.sqn").read_text(encoding="utf-8")
    assert sqn.startswith("Seq-submit ::=")
    assert "NC_050833.1" in sqn

    # the report's counts are the validator's own, and the only error class left is the one this
    # format produces for every duplicated IR gene
    summary = (module_dir / "tables" / "submission_summary.csv").read_text(encoding="utf-8")
    rows = csv.reader(
        io.StringIO("\n".join(line for line in summary.splitlines() if not line.startswith("#")))
    )
    values = dict(zip(next(rows), next(rows), strict=True))
    assert int(values["NCBI errors"]) > 0
    val = (module_dir / "submission.val").read_text(encoding="utf-8")
    errors = [line for line in val.splitlines() if line.startswith("Error:")]
    assert len(errors) == int(values["NCBI errors"])
    assert all("GeneXrefStrandProblem" in line for line in errors)

    # and the table itself carries what the record states, plus the one marker NCBI requires
    table = (module_dir / "submission.tbl").read_text(encoding="utf-8")
    assert "\t\t\torganism\tIris domestica" in table
    assert int(values["CDS written with a partial start"]) == 2
