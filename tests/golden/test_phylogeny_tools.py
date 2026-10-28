"""The phylogeny module against the real tools (M11.21) — golden, tier B.

These ran in the fast unit suite until 28/10/2026 and were ~5.5 minutes of every `make test`,
because six of them each run MAFFT and IQ-TREE. They are integration checks with the panel and
the pinned tools, so they belong here, where `make golden` runs them and a missing tool is named.
The pure parts of the module (declaration, refusals, partition text, MACSE wiring with stubs) stay
in `tests/unit/test_phylogeny_module.py`; this file is what actually calls the tools.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank

pytestmark = pytest.mark.tier_b

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation" / "source"


def _require_mafft() -> None:
    """The module aligns per gene with MAFFT before it ever reaches IQ-TREE."""
    from plastid.run.tools import which

    if which("mafft") is not None:
        return
    reason = "mafft is not on PATH — the phylogeny matrix is built per gene with MAFFT"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _iqtree_dir() -> Path:
    directory = os.environ.get("PTS_IQTREE_DIR")
    candidates = [Path(directory)] if directory else []
    candidates.append(Path("/root/toolbox/env/iqtree/bin"))
    for candidate in candidates:
        if (candidate / "iqtree2").is_file() or (candidate / "iqtree").is_file():
            return candidate
    reason = "iqtree2 is not installed; set PTS_IQTREE_DIR to the directory holding it"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def _panel_records() -> list[Path]:
    wanted = ("NC_008641.1", "NC_016708.2", "NC_028027.1", "NC_016736.1")
    paths = [PANEL / f"{accession}.gb" for accession in wanted]
    if not all(path.is_file() for path in paths):
        pytest.skip("the IR panel records are not in this tree (.panel/ is gitignored)")
    return paths


def _run_tree(tmp_path: Path, extras: list[str] | None = None) -> tuple[Path, dict]:
    from plastid.cli import main

    _require_mafft()
    records = _panel_records()
    argv = ["run", "-m", "phylogeny"]
    for record in records:
        argv += ["-i", str(record)]
    argv += ["-o", str(tmp_path), "--dataset", "tree", "--set", f"iqtree={_iqtree_dir()}"]
    for extra in extras or []:
        argv += ["--set", extra]
    code = main(argv)
    assert code == 0, "the phylogeny run failed — check the tool env"
    run_dir = next((tmp_path / "tree").iterdir())
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    return run_dir, manifest


def test_the_tree_names_every_taxon_once(tmp_path: Path) -> None:
    run_dir, _manifest = _run_tree(tmp_path)
    tree = next(run_dir.rglob("tree.treefile")).read_text(encoding="utf-8")
    for record in _panel_records():
        accession = read_genbank(record).accession
        assert tree.count(accession) == 1, f"{accession} appears {tree.count(accession)} times"
    assert tree.strip().endswith(";")


def test_the_tree_is_byte_identical_on_a_second_run(tmp_path: Path) -> None:
    """Two runs of one matrix give one tree, byte for byte.

    Measured 28/10/2026, and the first version of this test was **wrong**: it asserted the same
    thing while the module ran IQ-TREE with `-nt AUTO`, where two runs gave **different
    topologies and support values** (99.5 vs 99.7 on an internal node). IQ-TREE's parallel
    reductions are not associative. The fix was in the product, not the test: the module now
    passes `-nt 1` by default. `--set threads=N` still asks for speed at the cost of
    bit-reproducibility, and the module records that choice.
    """
    first, _ = _run_tree(tmp_path / "a")
    second, _ = _run_tree(tmp_path / "b")
    one = next(first.rglob("tree.treefile")).read_text(encoding="utf-8")
    two = next(second.rglob("tree.treefile")).read_text(encoding="utf-8")
    assert one == two
    assert "99.5/98" in one  # UFBoot/SH-aLRT present


def test_a_gene_present_in_one_genome_is_excluded_by_name(tmp_path: Path) -> None:
    run_dir, _ = _run_tree(tmp_path)
    table = next(run_dir.rglob("matrix_genes.csv")).read_text(encoding="utf-8")
    excluded = [line for line in table.splitlines() if line.startswith(("infA,", "ycf15,"))]
    assert excluded, "the genes carried by one genome should be listed"
    for line in excluded:
        assert ",no," in line
        assert "present in 1 genome" in line


def test_the_supermatrix_and_tree_are_artifacts_of_the_run(tmp_path: Path) -> None:
    """The alignment the tree came from is kept beside it, under the name it was written as."""
    run_dir, manifest = _run_tree(tmp_path)
    files = {entry["path"].split("/")[-1] for entry in manifest["files"]}
    assert {"supermatrix.fasta", "tree.treefile", "tree.iqtree"} <= files
    assert manifest["status"] == "ok"


def test_the_tool_and_its_argv_are_recorded(tmp_path: Path) -> None:
    _run_dir, manifest = _run_tree(tmp_path)
    tools = {t["name"]: t for module in manifest["modules"] for t in module.get("tools", [])}
    assert "iqtree2" in tools and tools["iqtree2"]["version"] == "2.3.6"
    argv = tools["iqtree2"]["argv"]
    assert argv[argv.index("-m") + 1] == "MFP"
    assert argv[argv.index("-bb") + 1] == "1000"
    assert argv[argv.index("-seed") + 1] == "1"
    assert "mafft" in tools


def test_a_gene_partition_file_reaches_iqtree_and_covers_every_included_gene(
    tmp_path: Path,
) -> None:
    """The NEXUS file is accepted by the real tool, and its charsets are the included genes.

    A partition that is written but never passed is a silent no-op; a file IQ-TREE rejects is a
    broken run. This is the only check that sees both, because NEXUS syntax is the tool's, not
    ours. The count is taken from ``matrix_genes`` — the table that says which genes went in —
    rather than hardcoded, so a change in the panel cannot make the two disagree silently.
    """
    run_dir, manifest = _run_tree(tmp_path, ["partition=gene"])

    partition = next(run_dir.rglob("supermatrix.partition.nex"))
    text = partition.read_text(encoding="utf-8")
    genes = next(run_dir.rglob("matrix_genes.csv")).read_text(encoding="utf-8")
    included = [line for line in genes.splitlines() if ",yes," in line]
    assert text.startswith("#nexus\nbegin sets;")
    assert text.count("  charset ") == len(included) > 1

    tools = {t["name"]: t for module in manifest["modules"] for t in module.get("tools", [])}
    argv = tools["iqtree2"]["argv"]
    assert argv[argv.index("-p") + 1] == "partition.nex"

    summary = next(run_dir.rglob("phylogeny_summary.csv")).read_text(encoding="utf-8")
    assert f"gene ({len(included)} partitions)" in summary
    assert next(run_dir.rglob("tree.treefile")).read_text(encoding="utf-8").strip().endswith(";")


def test_a_bootstrap_below_one_thousand_is_refused(tmp_path: Path) -> None:
    """A support value from fewer than 1000 replicates is not one this toolkit publishes."""
    from plastid.cli import main

    _require_mafft()
    records = _panel_records()
    argv = ["run", "-m", "phylogeny"]
    for record in records:
        argv += ["-i", str(record)]
    argv += [
        "-o",
        str(tmp_path),
        "--set",
        f"iqtree={_iqtree_dir()}",
        "--set",
        "ultrafast_bootstrap=100",
    ]
    assert main(argv) == 3


def test_an_intron_matrix_is_built_and_treed_by_the_real_tools(tmp_path: Path) -> None:
    """`matrix=introns` is CGAS module 14's separate-matrix route, on this toolkit's spans.

    The matrix must be intron sequences — its sites are the loci table's, not the CDS matrix's —
    and the IR gate runs for every genome before a span is read, because an intron boundary is a
    coordinate. A short panel keeps this cheap: introns are hundreds of bases, not a plastome.
    """
    run_dir, manifest = _run_tree(tmp_path, ["matrix=introns"])

    summary = next(run_dir.rglob("phylogeny_summary.csv")).read_text(encoding="utf-8")
    assert ",introns," in summary

    loci = next(run_dir.rglob("matrix_genes.csv")).read_text(encoding="utf-8")
    included = [line for line in loci.splitlines()[1:] if ",yes," in line]
    assert len(included) > 1
    sites = sum(int(line.split(",")[3]) for line in included)
    assert f",{sites}," in summary  # the supermatrix is the included loci, nothing else

    tree = next(run_dir.rglob("tree.treefile")).read_text(encoding="utf-8")
    assert tree.strip().endswith(";")
    assert manifest["status"] == "ok"
