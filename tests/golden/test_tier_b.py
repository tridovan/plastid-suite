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
                # MISA is not bundled (D-q: its zip states no licence), so the test says where
                # the copy it found is — the same `--set misa=` a user without it on PATH passes.
                "--set",
                f"misa={tools['misa.pl']}",
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
    assert record["version"] == "2.1"  # the pinned release, read out of the script's own header
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
    assert values["Total SSRs"] == "37"  # golden number
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
    """The golden Ts/Tv, measured by the toolkit's own path: **2.6917 = 323 Ts / 120 Tv** over the
    **79** shared CDS and **68,448** alignment columns, IR excluded.

    The scope is unchanged since 02/10 — CDS only · shared genes · MAFFT · IR excluded — but the
    *set* it ranges over was re-measured on 04/10: three loci (atpI, psaI, psbB) that
    NC_024936.1 spells differently are the same loci, so the shared set is 79, not 76, and the
    ratio is 2.6917. The old 2.686 = 317/118 over 66,066 columns was correct
    for a shared set counted by spelling, and stays in the log as that.

    This check is why the module's CDS scope exists at all — the golden number is *not* the
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
                # NC_024936.1's annotated IR is 99.80 %; the sequence's own repeat is the
                # explicit path for this comparison (D11-2).
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    summary = _table(tmp_path / "runs", "substitutions_summary")
    row = next(summary)
    assert row["Shared genes"] == "79"
    assert row["Alignment columns (bp)"] == "68448"
    assert (row["Transitions"], row["Transversions"]) == ("323", "120")
    assert row["Ts/Tv"] == "2.6917"

    # The twelve directed changes recombine into exactly those totals — this is the table the
    # manuscript's Figure 7 reads, and it is the golden number, measured on the golden pair.
    by_type = list(_table(tmp_path / "runs", "substitutions_by_type"))
    # twelve directed types for the one query genome — the reference is not a column
    assert len(by_type) == 12
    transitions = sum(int(r["Count"]) for r in by_type if r["Class"] == "Transition")
    transversions = sum(int(r["Count"]) for r in by_type if r["Class"] == "Transversion")
    assert (transitions, transversions) == (323, 120)


def test_the_naming_step_does_not_move_the_golden_ts_tv(
    tools: dict[str, str | None],
    blast_tools: dict[str, str | None],
    data_dir: Path,
    tmp_path: Path,
) -> None:
    """P3B-06: the golden number re-measured **under the new naming direction** — and unchanged.

    P3B-04 flipped the alias direction (five of nine loci now print the naming reference's
    spelling) and the three comparison modules read those canonical names, so 2.6917 was a
    *prediction* until it was re-run. This is the re-run, in the order the
    pipeline actually uses: `gene_normalize` writes the sample in vocabulary A, and the CDS scope
    compares what comes out — not the published spelling as it arrived.

    What makes the match a result rather than a coincidence is *why* it cannot move: the alias
    table maps every member of a locus to that locus whichever member is canonical, so the flip
    **relabels** the shared set instead of re-selecting it. It could only re-select if a record
    carried two *different* genes whose names fall in one locus — measured 08/10/2026 on both
    Iris records: 86 CDS → 79 canonical names each, and no locus reached by two spellings.
    The artifact is checked to be in vocabulary A below, so a run where the
    naming step quietly did nothing cannot pass this test.
    """
    require_tools(tools, "mafft")
    require_tools(blast_tools, "blastp", "makeblastdb")
    assert (
        main(
            [
                "run",
                "-m",
                "gene_normalize",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "normalized"),
                "--set",
                "mode=normalize",
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    artifacts = sorted((tmp_path / "normalized").rglob("NC_024936.1.normalized.gb"))
    assert len(artifacts) == 1, "the naming step wrote no artifact for the comparison record"

    # the input really changed: the three CDS the two records spell differently are in A now
    spellings = {gene.name for gene in read_genbank(artifacts[0]).genes if gene.kind == "CDS"}
    assert {"atpI", "psaI", "psbB"} <= spellings, "the artifact is not in vocabulary A"

    assert (
        main(
            [
                "run",
                "-m",
                "substitutions",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(artifacts[0]),
                "-o",
                str(tmp_path / "runs"),
                "--set",
                "scope=cds",
                # the normalized artifact keeps NC_024936.1's annotation, whose IR is 99.80 %
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    row = next(_table(tmp_path / "runs", "substitutions_summary"))
    assert row["Shared genes"] == "79"
    assert row["Alignment columns (bp)"] == "68448"
    assert (row["Transitions"], row["Transversions"]) == ("323", "120")
    assert row["Ts/Tv"] == "2.6917"


def test_the_whole_alignment_scope_is_a_different_number_and_says_so(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """The default scope is not wrong — it is *wider*, and its label has to keep it apart from
    the golden's. Measured 03/10/2026: 1.8416 over the whole aligned plastome."""
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
                "ir_source=sequence",
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

    There is no golden value for π, so this check holds the *relationships* — which is what
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
                "--set",
                "ir_source=sequence",
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


def test_repeats_reproduces_the_measured_golden_count(
    mummer_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """M2.12 on the golden — measured 10/10/2026 with `repeat-match -n 30` (MUMmer 4.0.1).

    **9** maximal exact matches at 30 bp or more: 4 dispersed, 5 palindromic, and the longest is
    the inverted-repeat pair itself — 26,214 bp between IRa and IRb, which is the match the board
    asks to be pinned rather than checked by eye. The phase's own trap was inheriting the 26/9
    survey's "9 matches ≥30 bp", measured on other data; this number is this record's, and
    the coordinates below are the proof it is not the survey's.
    """
    require_tools(mummer_tools, "repeat-match")
    assert (
        main(
            [
                "run",
                "-m",
                "repeats",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-o",
                str(tmp_path / "runs"),
            ]
        )
        == 0
    )
    summary = next(_table(tmp_path / "runs", "repeat_summary"))
    assert summary["Minimum match length (bp)"] == "30"
    assert summary["Matches"] == "9"
    assert (summary["Dispersed"], summary["Palindromic"]) == ("4", "5")
    assert summary["Longest match (bp)"] == "26214"
    assert summary["Longest match regions"] == "IRa–IRb"
    assert summary["Tool"].startswith("repeat-match (MUMmer 4.0.1)")

    repeats = {row["Repeat"]: row for row in _table(tmp_path / "runs", "repeats")}
    assert len(repeats) == 9
    ir_pair = repeats["R01"]
    assert ir_pair["Class"] == "palindromic"
    assert (ir_pair["Copy 1 start"], ir_pair["Copy 1 end"]) == ("83128", "109341")
    assert (ir_pair["Copy 1 region"], ir_pair["Copy 2 region"]) == ("IRa", "IRb")
    assert (ir_pair["Copy 2 start"], ir_pair["Copy 2 end"]) == ("127511", "153724")
    # every row is at or above the floor, and the two classes are the tool's own markings
    assert all(int(row["Length (bp)"]) >= 30 for row in repeats.values())
    assert {row["Class"] for row in repeats.values()} == {"dispersed", "palindromic"}
    assert sorted(int(row["Length (bp)"]) for row in repeats.values()) == [
        31,
        32,
        42,
        42,
        42,
        42,
        52,
        52,
        26_214,
    ]


def test_vmatch_agrees_with_repeat_match_on_the_golden(
    mummer_tools: dict[str, str | None],
    vmatch_tools: dict[str, str | None],
    data_dir: Path,
    tmp_path: Path,
) -> None:
    """D-f: the count is published because **two independent engines** give the same matches.

    Vmatch is academic-only, so the module never needs it — but a count that only one tool
    produces is a count nobody can check. Measured 10/10/2026: nine matches either way, the same
    coordinates, the same two classes (Vmatch prints 0-based starts and D/P; repeat-match prints
    1-based and an `r` — both land on the same 0-based half-open intervals), the IR pair among
    them. Converting the coordinates in one place is what makes the two comparable.
    """
    require_tools(mummer_tools, "repeat-match")
    require_tools(vmatch_tools, "mkvtree", "vmatch")
    from plastid.core.genbank import read_genbank
    from plastid.external.repeatmatch import run_repeat_match
    from plastid.external.vmatch import run_vmatch

    genome = read_genbank(data_dir / "NC_050833.1.gb")
    ctx = _ctx(tmp_path)
    mine = run_repeat_match(genome, ctx, min_length=30)
    theirs = run_vmatch(genome, ctx, min_length=30, directory=Path(vmatch_tools["vmatch"]).parent)

    assert len(mine.matches) == len(theirs.matches) == 9
    assert {(m.start1, m.end1, m.start2, m.end2, m.palindromic) for m in mine.matches} == {
        (m.start1, m.end1, m.start2, m.end2, m.palindromic) for m in theirs.matches
    }


def test_diversity_measures_pi_per_intron_and_per_spacer(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P2B-01 on the golden pair — measured 10/10/2026, never inherited.

    The numbers below are the golden's own: 17 introns (the `introns` module's set — the
    cross-check the board asks for) and 91 spacers (79 LSC · 12 SSC, 32,134 bp). The three IR
    intron pairs are the reason D-g excludes the IR from the spacer set: each pair measures the
    same locus twice, to the last decimal — rps12 0.000000 both copies, rpl2 0.001506 both,
    ndhB 0.000000 both, over identical usable-column counts.
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
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )

    summary = next(_table(tmp_path / "runs", "diversity_summary"))
    assert summary["Introns (reference)"] == "17"
    assert summary["Intergenic spacers (reference)"] == "91"

    introns = {
        row["Intron (gene.copy.n)"]: row for row in _table(tmp_path / "runs", "diversity_introns")
    }
    assert len(introns) == 17
    assert {row["Region"] for row in introns.values()} == {"LSC", "IRa", "SSC", "IRb"}
    for first, second in (
        ("rps12.1.2", "rps12.2.2"),
        ("rpl2.1.1", "rpl2.2.1"),
        ("ndhB.1.1", "ndhB.2.1"),
    ):
        assert introns[first]["π"] == introns[second]["π"], "an IR mirror pair is one locus"
        assert introns[first]["Usable columns"] == introns[second]["Usable columns"]
    assert introns["rpoC1.1.1"]["π"] == "0.004138"  # least diverse outside the IR
    assert introns["clpP.1.1"]["π"] == "0.012658"  # most diverse outside the IR

    spacers = list(_table(tmp_path / "runs", "diversity_igs"))  # a DictReader reads once
    assert len(spacers) == 91
    assert {row["Region"] for row in spacers} == {"LSC", "SSC"}
    assert sum(int(row["Length (bp)"]) for row in spacers) == 32_134
    widest = max(spacers, key=lambda row: float(row["π"]))
    assert (widest["Spacer (left–right)"], widest["π"]) == ("psbT-psbN", "0.082353")


def test_synteny_reproduces_the_measured_iris_pair_result(
    mummer_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P4-01 on real data: the golden pair through `nucmer --maxmatch → delta-filter -1 →
    show-coords`, measured 04/10/2026 — **7 blocks, %IDY 89.66–99.11 (mean 97.37)** over
    152,332 paired bases, coverage 99.09 %.

    The board carried "21 blocks, %IDY 83.7–98.2 (mean 91.5)" as this check. It could not be
    reproduced from the probe's own inputs under any nucmer setting (a sweep of -l/-c and every
    delta-filter mode is); the golden pair are congeners and align at
    ~98–99 %, while a cross-genus pair aligns at ~89 %. The measured numbers are asserted here
    and the finding is recorded — the check was changed by evidence, not the number tuned.
    """
    require_tools(mummer_tools, "nucmer", "delta-filter", "show-coords")
    directory = str(Path(mummer_tools["nucmer"]).parent)
    assert (
        main(
            [
                "run",
                "-m",
                "synteny",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
                "--dataset",
                "golden",
                "--set",
                f"mummer={directory}",
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "runs" / "golden").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)
    records = manifest["modules"][0]["tools"]
    assert [record["name"] for record in records] == ["nucmer", "delta-filter", "show-coords"]
    assert all(record["version"] == "4.0.1" for record in records)
    assert "--maxmatch" in records[0]["argv"]

    summary = next(_table(tmp_path / "runs", "synteny_genomes"))
    assert summary["Alignment blocks"] == "7"
    assert summary["Min %IDY (%)"] == "89.66"
    assert summary["Max %IDY (%)"] == "99.11"
    assert summary["Mean %IDY (%)"] == "97.37"
    assert summary["Aligned base pairs (bp)"] == "152332"
    assert summary["Coverage (%)"] == "99.09"

    blocks = list(_table(tmp_path / "runs", "synteny_blocks"))
    assert len(blocks) == 7
    assert all(row["Orientation"] == "forward" for row in blocks)

    # the series is per window with its size stated, and a hole is blank, never a zero
    series = list(_table(tmp_path / "runs", "synteny_identity"))
    assert series[0]["Window start (bp)"] == "1"
    assert series[0]["Window end (bp)"] == "600"
    assert (run_dir / "01_synteny" / "mummer_NC_024936.1.delta").is_file()
    assert (run_dir / "01_synteny" / "mummer_NC_024936.1.coords").is_file()


def test_kaks_reproduces_the_measured_iris_pair_result(
    tools: dict[str, str | None],
    paml_tools: dict[str, str | None],
    data_dir: Path,
    tmp_path: Path,
) -> None:
    """P4-02 on real data: every shared CDS aligned on its own with MAFFT and scored by PAML's
    yn00 (Yang & Nielsen 2000), measured 04/10/2026 — **79 shared CDS, 78 scored**, pooled
    Ka 0.003586 · Ks 0.016550 · Ka/Ks 0.2167 over the site-weighted genes.

    There is no published Ka/Ks for this pair (measured 03/10), so the run's own values are the
    record. The three checks that are not numbers:

    * *identical CDS* — atpI's two copies are byte-identical — give Ka = Ks = 0 and a **blank**
      Ka/Ks, never 0 and never 1;
    * *accD*'s MAFFT alignment breaks the reading frame (a gap after 1,445 bases), so it keeps a
      row with that reason instead of a number yn00 would have misread;
    * the three loci NC_024936.1 spells differently (atpL/psaL/bspB) are **in** this table under
      their canonical names — before the aliases were added they looked like genes the genome
      does not carry.
    """
    require_tools(tools, "mafft")
    require_tools(paml_tools, "yn00")
    directory = str(Path(paml_tools["yn00"]).parent)
    assert (
        main(
            [
                "run",
                "-m",
                "kaks",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
                "--dataset",
                "golden",
                "--set",
                f"paml={directory}",
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    run_dir = sorted((tmp_path / "runs" / "golden").iterdir())[0]
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    validate(manifest)
    names = [record["name"] for record in manifest["modules"][0]["tools"]]
    assert "mafft" in names and "yn00" in names
    yn00 = next(record for record in manifest["modules"][0]["tools"] if record["name"] == "yn00")
    assert yn00["version"] == "4.10.10"

    summary = next(_table(tmp_path / "runs", "kaks_summary"))
    assert summary["Genes considered"] == "79"
    assert summary["Genes shared with every annotated genome"] == "79"
    assert summary["Estimator"] == "Yang & Nielsen 2000 (yn00)"

    pooled = next(_table(tmp_path / "runs", "kaks_pooled"))
    assert pooled["Genes scored"] == "78"
    assert pooled["Genes without a value"] == "1"
    assert float(pooled["Pooled Ka"]) == pytest.approx(0.003586, rel=1e-3)
    assert float(pooled["Pooled Ks"]) == pytest.approx(0.016550, rel=1e-3)
    assert float(pooled["Pooled Ka/Ks"]) == pytest.approx(0.2167, rel=1e-3)

    genes = {row["Gene"]: row for row in _table(tmp_path / "runs", "kaks_genes")}
    assert len(genes) == 79  # every shared gene keeps a row, scored or not
    assert {"rbcL", "atpI", "psaI", "psbB", "accD"} <= set(genes)
    assert float(genes["rbcL"]["Ka (dN)"]) == pytest.approx(0.0019, abs=1e-4)
    assert float(genes["rbcL"]["Ks (dS)"]) == pytest.approx(0.0299, abs=1e-4)
    assert float(genes["rbcL"]["Ka/Ks"]) == pytest.approx(0.0635, abs=1e-3)
    assert genes["atpI"]["Status"] == "ks = 0"
    assert genes["atpI"]["Ka/Ks"] == ""  # undefined — never 0, never 1
    assert genes["accD"]["Status"] == "frameshift in the alignment"


def test_barcodes_reuses_diversitys_pi_window_for_window(
    tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P4-04's Done-when, as a cross-check between two modules rather than a re-implementation:
    the barcode module's window series and the diversity module's must be **the same numbers**
    for the same window — same alignment, same window/step, same core function.

    Measured 04/10/2026 on the golden pair: 773 windows, every π identical, and the candidate
    loci of ``data/barcode_loci.yaml`` checked against the record's own annotation — 8 of 10
    found. ``trnH-psbA`` is **not** one of them, and the row says exactly why: in both golden
    records rps19 sits between trnH-GUG and psbA, so the classic spacer is not a spacer there
    (the module's first version reported it as "not carried", which is a different claim).
    """
    require_tools(tools, "mafft")
    assert (
        main(
            [
                "run",
                "-m",
                "barcodes",
                "-m",
                "diversity",
                "-i",
                str(data_dir / "NC_050833.1.gb"),
                "-i",
                str(data_dir / "NC_024936.1.gb"),
                "-o",
                str(tmp_path / "runs"),
                "--dataset",
                "golden",
                "--set",
                "ir_source=sequence",
            ]
        )
        == 0
    )
    barcode_windows = list(_table(tmp_path / "runs", "barcodes_windows"))
    diversity_windows = list(_table(tmp_path / "runs", "diversity_windows"))
    assert len(barcode_windows) == len(diversity_windows) == 773
    assert [row["π"] for row in barcode_windows] == [row["π"] for row in diversity_windows]

    summary = next(_table(tmp_path / "runs", "barcodes_summary"))
    assert summary["Hypervariable regions"] == "26"
    assert summary["Candidate loci"] == "10"
    assert summary["Candidates found"] == "8"
    assert float(summary["Mean π"]) == pytest.approx(0.007332, rel=1e-3)
    assert float(summary["Region threshold"]) == pytest.approx(0.014157, rel=1e-3)

    regions = list(_table(tmp_path / "runs", "barcodes_regions"))
    assert regions[0]["Region"] == "R1"  # ranked by mean π
    assert regions[0]["Reference start (bp)"] != ""  # coordinates, not just alignment columns
    assert all(row["Span (bp)"] != "" for row in regions)

    candidates = {row["Locus"]: row for row in _table(tmp_path / "runs", "barcodes_candidates")}
    assert candidates["rbcL"]["Found"] == "yes"
    assert candidates["atpF-atpH"]["Rank by π"] == "1"  # the most variable candidate here
    assert candidates["trnH-psbA"]["Found"] == "no"
    assert "not adjacent" in candidates["trnH-psbA"]["Note"]
    assert candidates["ITS"]["Found"] == "no"
    assert "plastid" in candidates["ITS"]["Note"]
    # an intron-containing locus keeps both its parts: rps16's span is its two exons (252 bp)
    assert candidates["rps16"]["Span (bp)"] == "252"
    assert float(candidates["rps16"]["π"]) == 0.0  # measured: the exons are identical here


PANEL = Path(__file__).parents[2] / ".panel" / "oldenlandia-coleus"


def test_popgen_over_the_pinned_panel(tools: dict[str, str | None], tmp_path: Path) -> None:
    """P4-03's real run: the pinned stand-in panel, every statistic with its n.

    Measured 04/10/2026 — 11 records in two groups (Coleus 6 · Oldenlandia 5; three of them carry
    no CDS), 80 genes aligned across them. **n = 4 for most genes and the neutrality tests then
    exist; n = 3 and n = 1 rows keep a blank with the reason**, which is the point of the check:
    a statistic that needed more sequences must not appear as a number.

    The panel is gitignored (re-fetch with `plastid panel pin`) so this skips with a named
    reason when it is absent — and what it does read is verified against ``PINNED.json``'s
    sha256 first, so a drifted download cannot quietly change the numbers.
    """
    require_tools(tools, "mafft")
    if not (PANEL / "PINNED.json").is_file():
        pytest.skip(
            f"the pinned stand-in panel is not on this machine ({PANEL}) — re-fetch it with the "
            f"two 'plastid panel pin --groups' commands"
        )
    import hashlib
    import json as _json

    pinned = _json.loads((PANEL / "PINNED.json").read_text(encoding="utf-8"))
    paths = sorted(PANEL.glob("*.gb"))
    assert len(paths) == len(pinned["genomes"])
    for entry in pinned["genomes"]:
        digest = hashlib.sha256((PANEL / Path(entry["path"]).name).read_bytes()).hexdigest()
        assert digest == entry["sha256"], f"{entry['accession']} has drifted from its pin"

    assert (
        main(
            [
                "run",
                "--module",
                "popgen",
                *[part for path in paths for part in ("--input", str(path))],
                "--out",
                str(tmp_path / "runs"),
                "--dataset",
                "panel",
                "--groups",
                str(Path(__file__).parents[2] / "data" / "oldenlandia-coleus.yaml"),
            ]
        )
        == 0
    )
    summary = next(_table(tmp_path / "runs", "popgen_summary"))
    assert summary["Genomes"] == "11"
    assert summary["Groups"] == "2"
    assert summary["Group sizes"] == "Coleus: 6, Oldenlandia: 5"
    assert summary["Genes analysed"] == "80"

    groups = {row["Row"]: row for row in _table(tmp_path / "runs", "popgen_groups")}
    rbc = groups["rbcL#Oldenlandia"]
    assert rbc["Sequences"] == "4"
    assert rbc["Segregating sites"] == "41"
    assert float(rbc["π per site"]) == pytest.approx(0.01844071, rel=1e-4)
    assert float(rbc["Tajima's D"]) == pytest.approx(1.8434, abs=1e-3)

    thin = groups["ndhK#Oldenlandia"]
    assert thin["Sequences"] == "3"
    assert thin["Tajima's D"] == ""  # n = 3 cannot carry it
    assert "at least 4 sequences" in thin["Note"]

    alone = groups["ycf15#Oldenlandia"]
    assert alone["Sequences"] == "1"
    assert alone["Haplotype diversity"] == ""
    assert "nothing is estimable" in alone["Note"]

    fst = {row["Row"]: row for row in _table(tmp_path / "runs", "popgen_fst")}
    assert float(fst["rbcL#Coleus#Oldenlandia"]["Fst"]) == pytest.approx(0.1852, abs=1e-3)
    assert len(fst) == 80  # every analysed gene has its between-group row


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


def _ctx(tmp_path: Path):
    from plastid.run.context import RunContext
    from plastid.run.layout import RunLayout

    layout = RunLayout.create(tmp_path / "ctx", "tier-b", fingerprint="test")
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
    NCBI's form generates from Trang's details (``tests/data/PROVENANCE.md``). What is checked
    here is the pipeline: the tool exits 0, the ``.sqn`` exists, and every error NCBI reports
    belongs to the one class a duplicated inverted-repeat plastome produces (measured 2026-10-03
    against a published plastome submission) — a new class would fail this test,
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
