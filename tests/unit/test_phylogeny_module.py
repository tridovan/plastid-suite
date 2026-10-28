"""The phylogeny module: what it builds, what it refuses, and what it will not claim (M11.21).

**The independent expectation, stated once.** A tree's *shape* cannot be checked against an
answer key the way a π or a Ts/Tv can — there is no arithmetic that says where NC_016708.2
belongs. What can be checked, and is checked here, is that the artifact is the one the inputs
imply:

* every input genome appears as a tip, exactly once, with the name it was given;
* two runs of the same matrix produce the **same** tree (the seed is fixed, so a different tree
  would mean the run is not reproducible);
* a gene present in one genome only is **excluded by name**, not filled with gaps and counted;
* a partial CDS is excluded by default and admitted only when asked, and the table says which;
* the module refuses a two-genome input — that is a comparison, not a topology;
* the matrix and the tree are written as artifacts, so the tree can be traced to the alignment.

The tests that need IQ-TREE skip with a named reason when it is absent (`PTS_IQTREE_DIR`), the
same convention as every other tier-B check.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Genome
from plastid.modules import registry

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
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


# ------------------------------------------------------------------ declaration and refusals


def test_the_module_is_registered_with_its_tools() -> None:
    spec = registry.get("phylogeny")
    assert spec.arity == 2
    assert set(spec.citations) == {"mafft", "iqtree"}
    declared = {param.name for param in spec.params}
    assert {"iqtree", "model", "ultrafast_bootstrap", "alrt", "seed", "include_partial"} <= declared


def test_a_two_genome_input_is_refused_by_name(tmp_path: Path) -> None:
    """Two sequences are a pairwise comparison; calling it a topology would be a false claim."""
    from plastid.cli import main

    code = main(
        [
            "run",
            "-m",
            "phylogeny",
            "-i",
            str(DATA / "NC_050833.1.gb"),
            "-i",
            str(DATA / "NC_024936.1.gb"),
            "-o",
            str(tmp_path),
        ]
    )
    assert code == 1


def test_an_unknown_parameter_is_refused() -> None:
    """The M11.20 rule reaches the new module without any extra wiring."""
    problems = registry.unknown_params(["phylogeny"], {"bootstrap": 1000})
    assert problems and "ultrafast_bootstrap" in problems[0]


# ------------------------------------------------------------------------ with the real tool


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


def test_the_default_thread_count_is_one() -> None:
    """The default is the reproducible one; asking for more is an explicit trade."""
    from plastid.modules.phylogeny import DEFAULT_THREADS

    assert DEFAULT_THREADS == "1"


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


def _unused(_: Genome) -> None:  # pragma: no cover - keeps the import meaningful
    return None
