"""The PREP-Cp-compatible predictor against its acceptance set — golden, tier B (MAFFT).

Three independent checks, kept apart (docs/tools/prep_cp.md, "Acceptance evidence"):

1. **The archived genuine PREP-Cp outputs.** The Internet Archive kept two results of the original
   server for the form's sample — the *Vitis vinifera* accD CDS of NC_007957.1 — at C = 0.2 (eight
   sites) and C = 0.8 (two). Positions, effects and scores must all come back; the alignment column
   depends on the reference and is not compared.
2. **Reference B′ is reproducible.** Rebuilding it from the pinned records, the sourced edit list
   and the pinned MAFFT gives the packaged file byte for byte, so its digest names one thing.
3. **Aegilops tauschii KJ614412.1** (Wang et al. 2017) — the paper's original PREP-Cp prediction
   (Table 2), its RNA-seq calls (Table 3) and its RT-PCR checks, compared site by site and never
   pooled. Every difference from Table 2 is listed below with the reference column that decides
   it, and the test re-reads that column from the run's own alignment, so an explanation cannot
   outlive the evidence it rests on. Nothing is tuned to reach 34.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

import pytest

from conftest import require_tools
from plastid.cli import main
from plastid.core.resources import data_file
from plastid.core.rna_editing import SUPPORTED_GENES, load_reference

pytestmark = pytest.mark.tier_b

ROOT = Path(__file__).resolve().parents[2]
RECORDS = ROOT / "tests" / "data" / "rna_editing"
# B′'s seven records are pinned in taxa.yaml and kept out of git (`plastid panel pin`).
B_PRIME_RECORDS = ROOT / ".panel" / "rna_editing"


def _run(tmp_path: Path, record: str, *sets: str) -> Path:
    out = tmp_path / "-".join((record, *sets))
    arguments = ["run", "-m", "rna_editing", "-i", str(RECORDS / f"{record}.gb"), "-o", str(out)]
    for item in sets:
        arguments += ["--set", item]
    assert main(arguments) == 0
    (run_dir,) = [path for path in out.glob("*/*") if (path / "manifest.json").is_file()]
    return run_dir


def _table(run_dir: Path, table: str) -> list[dict[str, str]]:
    lines = (run_dir / "01_rna_editing" / "tables" / f"{table}.csv").read_text().splitlines()
    return list(csv.DictReader(lines[2:]))


def _alignment(run_dir: Path, locus: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    name = ""
    path = run_dir / "01_rna_editing" / "alignments" / f"{locus.replace('#', '_')}.fasta"
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0]
        else:
            rows[name] = line
    return rows


def _column(run_dir: Path, locus: str, aa_position: int) -> dict[str, str]:
    """Taxon (genus) → reference residue in the column of the query's residue ``aa_position``."""
    rows = _alignment(run_dir, locus)
    query = rows.pop("query")
    seen = 0
    for column, residue in enumerate(query):
        if residue != "-":
            seen += 1
            if seen == aa_position:
                return {name.split("_")[0]: row[column] for name, row in rows.items()}
    raise AssertionError(f"{locus} has no residue {aa_position}")


@pytest.fixture(scope="module")
def vitis(tmp_path_factory: pytest.TempPathFactory, tools: Mapping[str, str | None]) -> dict:
    require_tools(dict(tools), "mafft")
    root = tmp_path_factory.mktemp("vitis")
    return {cutoff: _run(root, "NC_007957.1", f"cutoff={cutoff}") for cutoff in ("0.2", "0.8")}


def test_vitis_accd_reproduces_both_archived_prep_cp_outputs(vitis: dict) -> None:
    archived = list(
        csv.DictReader(
            (RECORDS / "prep_cp_archived_vitis_accD.tsv").read_text().splitlines(), delimiter="\t"
        )
    )
    for cutoff, run_dir in vitis.items():
        expected = [
            (int(row["nt_pos"]), int(row["aa_pos"]), row["effect"], row["score"])
            for row in archived
            if row["cutoff"] == cutoff
        ]
        found = [
            (
                int(row["CDS position (nt)"]),
                int(row["Amino-acid position"]),
                row["Effect"],
                row["Score"],
            )
            for row in _table(run_dir, "rna_editing_sites")
            if row["Locus"] == "accD"
        ]
        assert found == expected, cutoff
    query = {row["Locus"]: row for row in _table(vitis["0.2"], "rna_editing_query")}
    assert query["accD"]["CDS length (nt)"] == "1527"  # the form's sample, base for base


def test_the_c_0_8_result_is_the_c_0_2_sites_that_score_at_least_0_8(vitis: dict) -> None:
    low = {row["Site"]: row for row in _table(vitis["0.2"], "rna_editing_sites")}
    high = {row["Site"] for row in _table(vitis["0.8"], "rna_editing_sites")}
    assert high == {site for site, row in low.items() if float(row["Score"]) >= 0.8 - 1e-9}


def test_every_output_carries_the_label_and_the_reference_digest(vitis: dict) -> None:
    digest = load_reference().digest
    run_dir = vitis["0.8"]
    manifest = json.loads((run_dir / "manifest.json").read_text())
    (entry,) = [
        item for item in manifest["reference_sets"] if item["role"] == "rna_editing_reference"
    ]
    assert entry["digest"] == digest and len(entry["taxa"]) == 7
    assert manifest["params"]["rna_editing_reference_digest"] == digest
    for path in (run_dir / "01_rna_editing" / "tables").glob("*.csv"):
        scope = path.read_text().splitlines()[1]
        assert "predicted · PREP-Cp-compatible" in scope and digest in scope, path.name
    result = json.loads((run_dir / "01_rna_editing" / "result.json").read_text())
    assert all(digest in table["scope"]["label"] for table in result["tables"].values())
    assert not any(
        "original PREP-Cp output" in note and "never" not in note for note in result["notes"]
    )
    report = (run_dir / "01_rna_editing" / "report.md").read_text()
    assert "PREP-Cp-compatible" in report and digest in report
    sites = _table(run_dir, "rna_editing_sites")
    assert {row["Status"] for row in sites} == {"predicted · PREP-Cp-compatible"}
    assert {row["Reference digest"] for row in sites} == {digest}
    # The QC report reads the real manifest — tool records and reference sets as written (K-15).
    assert main(["qc", str(run_dir)]) == 0
    qc = (run_dir / "report" / "qc.md").read_text()
    assert "rna_editing_reference" in qc and digest in qc and "not the original PREP-Cp" in qc
    assert "| rna_editing | mafft |" in qc


@pytest.mark.panel
@pytest.mark.skipif(
    not B_PRIME_RECORDS.is_dir(),
    reason="the seven B′ records are not in this working tree (.panel/ is gitignored): "
    "plastid panel pin --accessions … --out .panel/rna_editing",
)
def test_reference_b_prime_rebuilds_byte_for_byte(
    tmp_path: Path, tools: Mapping[str, str | None]
) -> None:
    require_tools(dict(tools), "mafft")
    out = tmp_path / "reference_b_prime.json"
    code = main(
        [
            "editing-reference",
            "build",
            "--records",
            str(B_PRIME_RECORDS),
            "--out",
            str(out),
            "--runs",
            str(tmp_path / "runs"),
        ]
    )
    assert code == 0
    assert out.read_bytes() == data_file("rna_editing/reference_b_prime.json").read_bytes()


# --------------------------------------------------------------------------- Aegilops tauschii

# Table 2 sites this run does not reproduce. The column composition is the explanation: the residue
# each reference taxon holds where the query codon aligns, re-read from the run's alignment.
NOT_REPRODUCED = {
    # S 1/2: a grass-only insertion in B′'s alignment — five rows are gaps. The original's
    # alignment was adjusted by hand and never published, so this column cannot be compared.
    ("rpoC2", 2158, "rpoC2", 720): {
        "Arabidopsis": "-",
        "Nicotiana": "-",
        "Atropa": "-",
        "Zea": "S",
        "Oryza": "L",
        "Phalaenopsis": "-",
        "Pinus": "-",
    },
    # L 5/7. The S are rice (its list's source has no rpoC2 row) and Pinus (the 1996 survey did
    # not examine rpoC2); either as L would give 6/7 and the site.
    ("rpoC2", 4031, "rpoC2", 1344): {
        "Arabidopsis": "L",
        "Nicotiana": "L",
        "Atropa": "L",
        "Zea": "L",
        "Oryza": "S",
        "Phalaenopsis": "L",
        "Pinus": "S",
    },
}
# Sites this run predicts that Table 2 does not: both are single-residue columns, where one
# reference residue makes a score of 1/1.
NOT_IN_TABLE_2 = {
    # maize's petB CDS is annotated 19 codons longer at the N-terminus than the other taxa's
    ("petB", 20, "petB", 7): {"Zea": "V"},
    # inside the grass-only rpoC2 insertion
    ("rpoC2", 2123, "rpoC2", 708): {"Zea": "I"},
}


def _sets() -> dict[str, dict[tuple[str, int], dict[str, str]]]:
    sets: dict[str, dict[tuple[str, int], dict[str, str]]] = {}
    lines = (RECORDS / "wang2017_aegilops_sets.tsv").read_text().splitlines()
    for row in csv.DictReader(lines, delimiter="\t"):
        if row["cds_position"]:
            sets.setdefault(row["set"], {})[(row["gene"], int(row["cds_position"]))] = row
    return sets


def _assessable(row: Mapping[str, str]) -> bool:
    """A site the method can predict: one of the 35 genes, an amino-acid change, not a stop."""
    before, after = row["aa_change"].split("→")
    return row["gene"] in SUPPORTED_GENES and before != after and after.lower() != "stop"


@pytest.fixture(scope="module")
def aegilops(tmp_path_factory: pytest.TempPathFactory, tools: Mapping[str, str | None]) -> Path:
    require_tools(dict(tools), "mafft")
    return _run(tmp_path_factory.mktemp("aegilops"), "KJ614412.1")


def _predicted(run_dir: Path) -> set[tuple[str, int]]:
    """Predicted sites by gene and CDS position: the two IR copies of a gene are one site."""
    return {
        (row["Name (vocabulary A)"], int(row["CDS position (nt)"]))
        for row in _table(run_dir, "rna_editing_sites")
    }


def test_aegilops_against_the_original_prep_cp_prediction(aegilops: Path) -> None:
    """Table 2 — 34 sites by the original server at C = 0.8: 32 reproduced, each miss explained."""
    predicted = _predicted(aegilops)
    table2 = set(_sets()["table2_prepcp"])
    assert len(table2) == 34
    assert table2 - predicted == {key[:2] for key in NOT_REPRODUCED}
    assert predicted - table2 == {key[:2] for key in NOT_IN_TABLE_2}
    assert len(table2 & predicted) == 32
    for (gene, _, locus, aa_position), column in {**NOT_REPRODUCED, **NOT_IN_TABLE_2}.items():
        observed = _column(aegilops, locus, aa_position)
        assert {taxon: observed[taxon] for taxon in column} == column, (gene, aa_position)
        assert set(observed) - set(column) <= {
            taxon for taxon, residue in observed.items() if residue == "-"
        }, (gene, aa_position)


def test_aegilops_against_rna_seq(aegilops: Path) -> None:
    """Table 3 — 60 RNA-seq calls; only the ones the method could predict are compared."""
    predicted = _predicted(aegilops)
    rnaseq = _sets()["rnaseq"]
    assessable = {key for key, row in rnaseq.items() if _assessable(row)}
    assert len(assessable) == 26
    assert len(assessable & predicted) == 22
    # ndhG 347, petL 56 and rpoA 527 were not predicted by the original server either (Table 2).
    assert assessable - predicted == {("ndhG", 347), ("petL", 56), ("rpoA", 527), ("rpoC2", 4031)}


def test_aegilops_against_rt_pcr(aegilops: Path) -> None:
    """24 RT-PCR checks: 23 edited, 20 of them predictable; rpoB-617 was not edited."""
    predicted = _predicted(aegilops)
    sets = _sets()
    edited = {key for key, row in sets["rtpcr_edited"].items() if _assessable(row)}
    assert len(sets["rtpcr_edited"]) == 23 and len(edited) == 20
    # 17 of 20, as the original server's Table 2: the same three sites missed by both.
    assert edited - predicted == {("ndhG", 347), ("petL", 56), ("rpoA", 527)}
    # The original server predicted rpoB-617 too (Table 2): a conserved-residue false positive of
    # the method itself, not of this reference.
    assert set(sets["rtpcr_unedited"]) == {("rpoB", 617)}
    assert ("rpoB", 617) in predicted and ("rpoB", 617) in sets["table2_prepcp"]


def test_aegilops_coverage_states_every_gene(aegilops: Path) -> None:
    coverage = {row["Locus"]: row for row in _table(aegilops, "rna_editing_coverage")}
    assert coverage["accD"]["Status"] == "absent from the record"
    assert coverage["accD"]["Predicted sites"] == ""
    assert {coverage["ndhB#1"]["Status"], coverage["ndhB#2"]["Status"]} == {"evaluated"}
    not_supported = [row for row in coverage.values() if row["Status"] == "not supported"]
    assert not_supported and all(row["Predicted sites"] == "" for row in not_supported)


# ------------------------------------------------------------- Nicotiana sylvestris (unseen)

# Mower 2009 Table 3 tested PREP-Cp on N. sylvestris, a taxon outside its reference: TP 32, FP 2,
# FN 1 at C = 0.8 — the 33 sites of Sasaki et al. 2003's N. sylvestris column. The same 33,
# placed on the record (`sasaki2003_nsylvestris_sites.tsv`), against this run.
NSYLVESTRIS = B_PRIME_RECORDS / "NC_007500.1.gb"
NSYLVESTRIS_SHA256 = "786941e34e74a196fc4b5bb16885daa2cc2f621f9845d6c39ac80c315009c0e0"
# Known sites this run does not predict, with the column that decides each one.
NSYLVESTRIS_MISSED = {
    # I 4/7. The T are Atropa (rps2-1 was found in 2003, after the 2002 Atropa survey, so Atropa
    # was never assessed there) and the two grasses (their list's source has no rps2 row).
    ("rps2", 134, "rps2", 45): {
        "Arabidopsis": "I",
        "Nicotiana": "I",
        "Atropa": "T",
        "Zea": "T",
        "Oryza": "T",
        "Phalaenopsis": "I",
        "Pinus": "I",
    },
    # L 5/7 — the column of the Aegilops rpoC2 4031 site: S in rice and Pinus.
    ("rpoC2", 3731, "rpoC2", 1244): {
        "Arabidopsis": "L",
        "Nicotiana": "L",
        "Atropa": "L",
        "Zea": "L",
        "Oryza": "S",
        "Phalaenopsis": "L",
        "Pinus": "S",
    },
}


@pytest.mark.panel
@pytest.mark.skipif(
    not NSYLVESTRIS.is_file(),
    reason="NC_007500.1 is not in this working tree (.panel/ is gitignored): "
    "plastid panel pin --accessions NC_007500.1 --out .panel/rna_editing",
)
def test_nicotiana_sylvestris_unseen_as_in_mower_2009_table_3(
    tmp_path: Path, tools: Mapping[str, str | None]
) -> None:
    """31 of the 33 known sites predicted (Mower: 32), no prediction outside them (Mower: 2)."""
    require_tools(dict(tools), "mafft")
    assert hashlib.sha256(NSYLVESTRIS.read_bytes()).hexdigest() == NSYLVESTRIS_SHA256
    out = tmp_path / "nsylvestris"
    assert main(["run", "-m", "rna_editing", "-i", str(NSYLVESTRIS), "-o", str(out)]) == 0
    (run_dir,) = [path for path in out.glob("*/*") if (path / "manifest.json").is_file()]
    known = {
        (row["gene"], int(row["cds_position"]))
        for row in csv.DictReader(
            (RECORDS / "sasaki2003_nsylvestris_sites.tsv").read_text().splitlines(),
            delimiter="\t",
        )
    }
    predicted = _predicted(run_dir)
    assert len(known) == 33
    assert known - predicted == {key[:2] for key in NSYLVESTRIS_MISSED}
    assert predicted - known == set()
    for (_, _, locus, aa_position), column in NSYLVESTRIS_MISSED.items():
        assert _column(run_dir, locus, aa_position) == column, (locus, aa_position)
