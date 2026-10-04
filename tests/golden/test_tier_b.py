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
    """PROJECT.md §6, measured by the toolkit's own path: **2.6917 = 323 Ts / 120 Tv** over the
    **79** shared CDS and **68,448** alignment columns, IR excluded.

    The scope is unchanged since 02/10 — CDS only · shared genes · MAFFT · IR excluded — but the
    *set* it ranges over was re-measured on 04/10: three loci (atpI, psaI, psbB) that
    NC_024936.1 spells differently are the same loci, so the shared set is 79, not 76, and the
    ratio is 2.6917 (PROJECT.md §8.27). The old 2.686 = 317/118 over 66,066 columns was correct
    for a shared set counted by spelling, and stays in the log as that.

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
    assert row["Shared genes"] == "79"
    assert row["Alignment columns (bp)"] == "68448"
    assert (row["Transitions"], row["Transversions"]) == ("323", "120")
    assert row["Ts/Tv"] == "2.6917"


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


def test_synteny_reproduces_the_measured_iris_pair_result(
    mummer_tools: dict[str, str | None], data_dir: Path, tmp_path: Path
) -> None:
    """P4-01 on real data: the golden pair through `nucmer --maxmatch → delta-filter -1 →
    show-coords`, measured 04/10/2026 — **7 blocks, %IDY 89.66–99.11 (mean 97.37)** over
    152,332 paired bases, coverage 99.09 %.

    The board carried "21 blocks, %IDY 83.7–98.2 (mean 91.5)" as this check. It could not be
    reproduced from the probe's own inputs under any nucmer setting (a sweep of -l/-c and every
    delta-filter mode is recorded in PROJECT.md §9); the golden pair are congeners and align at
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
      does not carry (PROJECT.md §8.27).
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
    assert summary["Genes shared with every genome"] == "79"
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

    The panel is gitignored (re-fetch with PROJECT.md §8.6) so this skips with a named reason when
    it is absent — and what it does read is verified against ``PINNED.json``'s sha256 first, so a
    drifted download cannot quietly change the numbers.
    """
    require_tools(tools, "mafft")
    if not (PANEL / "PINNED.json").is_file():
        pytest.skip(
            f"the pinned stand-in panel is not on this machine ({PANEL}) — re-fetch it with the "
            f"two 'plastid panel pin --groups' commands in PROJECT.md §8.6"
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
