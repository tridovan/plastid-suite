"""K-06 / M11.7 — the tree figure records the tree it actually drew.

Checking the caption alone would not see a root or a drop that never reached R. The script now
writes ``<tag>_used.treefile`` beside its images, so this check reads the transformed tree: the
tips that were dropped are gone, the outgroup root is where it was asked to be, and ``--root
none`` keeps the input's own root. ``argv`` is asserted too, because a parameter that is not
forwarded cannot reach either output.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from conftest import require_r
from plastid.external.rscript import missing_packages, run_figure
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

NEWICK = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.25,(E:0.1,F:0.12):0.3);"


def _context(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "tree", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="r_tree", module_order=1, ir=None)


def _require_r(ctx: RunContext) -> str:
    binary = shutil.which("Rscript")
    if binary is None:
        require_r("no Rscript on PATH (docs/tools/r.md)")
    missing = missing_packages(ctx, binary)
    if missing:
        require_r(f"R is missing {', '.join(missing)}")
    return binary


def _used_tree(out_dir: Path, tag: str):
    from Bio import Phylo

    path = out_dir / f"{tag}_used.treefile"
    assert path.is_file(), f"the tree as drawn was not written: {path}"
    return Phylo.read(path, "newick")


def test_the_dropped_tip_is_absent_from_the_tree_that_was_drawn(tmp_path: Path) -> None:
    ctx = _context(tmp_path)
    binary = _require_r(ctx)
    tree = tmp_path / "input.treefile"
    tree.write_text(NEWICK + "\n", encoding="utf-8")
    out = tmp_path / "keep"

    run_figure(
        ctx,
        "tree",
        input_file=tree,
        out_dir=out,
        extra=["--tag", "keep", "--root", "none", "--drop-tips", "C"],
        binary=binary,
        check_packages=False,
    )

    drawn = _used_tree(out, "keep")
    assert sorted(terminal.name for terminal in drawn.get_terminals()) == ["A", "B", "D", "E", "F"]
    argv = list(ctx.tools[-1].argv)
    assert argv[argv.index("--drop-tips") :][:2] == ["--drop-tips", "C"]
    assert argv[argv.index("--root") :][:2] == ["--root", "none"]


def test_the_outgroup_root_reaches_the_tree_and_the_caption_does_not_claim_otherwise(
    tmp_path: Path,
) -> None:
    ctx = _context(tmp_path)
    binary = _require_r(ctx)
    tree = tmp_path / "input.treefile"
    tree.write_text(NEWICK + "\n", encoding="utf-8")
    out = tmp_path / "og"

    run_figure(
        ctx,
        "tree",
        input_file=tree,
        out_dir=out,
        extra=["--tag", "og", "--root", "outgroup", "--outgroup", "E,F", "--outgroup-label", ""],
        binary=binary,
        check_packages=False,
    )

    drawn = _used_tree(out, "og")
    tips = sorted(terminal.name for terminal in drawn.get_terminals())
    assert tips == ["A", "B", "C", "D", "E", "F"]
    root_sides = [
        {terminal.name for terminal in clade.get_terminals()} for clade in drawn.root.clades
    ]
    assert {"E", "F"} in root_sides, f"the outgroup is not a root side: {root_sides}"


def test_no_parameter_means_no_drop_and_a_midpoint_root(tmp_path: Path) -> None:
    """The neutral default: every tip stays, and the tree is rooted without a chosen outgroup."""
    ctx = _context(tmp_path)
    binary = _require_r(ctx)
    tree = tmp_path / "input.treefile"
    tree.write_text(NEWICK + "\n", encoding="utf-8")
    out = tmp_path / "plain"

    run_figure(
        ctx,
        "tree",
        input_file=tree,
        out_dir=out,
        extra=["--tag", "plain"],
        binary=binary,
        check_packages=False,
    )

    drawn = _used_tree(out, "plain")
    tips = sorted(terminal.name for terminal in drawn.get_terminals())
    assert tips == ["A", "B", "C", "D", "E", "F"]
    assert len(drawn.root.clades) == 2, "a midpoint root has two sides"
    argv = ctx.tools[-1].argv
    assert "--drop-tips" not in argv and "--outgroup" not in argv


def _draw(tmp_path: Path, newick: str, tag: str, extra: list[str]):
    """The tree as drawn, or the ToolError R stopped with."""
    ctx = _context(tmp_path)
    binary = _require_r(ctx)
    tree = tmp_path / f"{tag}.input.treefile"
    tree.write_text(newick + "\n", encoding="utf-8")
    out = tmp_path / tag
    run_figure(
        ctx,
        "tree",
        input_file=tree,
        out_dir=out,
        extra=["--tag", tag, *extra],
        binary=binary,
        check_packages=False,
    )
    return _used_tree(out, tag)


def _root_sides(drawn) -> list[set[str]]:
    return [{terminal.name for terminal in clade.get_terminals()} for clade in drawn.root.clades]


def test_one_outgroup_tip_roots_the_tree(tmp_path: Path) -> None:
    """A single outgroup is the ordinary case in a small tree, and it is a valid root."""
    drawn = _draw(tmp_path, NEWICK, "one", ["--root", "outgroup", "--outgroup", "E"])
    assert {"E"} in _root_sides(drawn)


def test_naming_an_outgroup_is_enough_to_root_on_it(tmp_path: Path) -> None:
    """The caption and docs/tools/r.md say a named outgroup is the root.

    A,B and not E,F: the longest path here runs through the E,F stem, so a midpoint root would
    land on E,F and pass this test for the wrong reason.
    """
    drawn = _draw(tmp_path, NEWICK, "named", ["--outgroup", "A,B"])
    assert {"A", "B"} in _root_sides(drawn)


def test_an_outgroup_that_is_not_in_the_tree_stops_the_figure(tmp_path: Path) -> None:
    """A midpoint tree drawn under a name that matched nothing reads as the rooted one."""
    with pytest.raises(ToolError, match="not tips of this tree: Z"):
        _draw(tmp_path, NEWICK, "miss", ["--root", "outgroup", "--outgroup", "E,Z"])


def test_a_dropped_name_that_is_not_in_the_tree_stops_the_figure(tmp_path: Path) -> None:
    """Otherwise the figure says it left out a tip, and every tip is still there."""
    with pytest.raises(ToolError, match="not tips of this tree: Z"):
        _draw(tmp_path, NEWICK, "drop", ["--drop-tips", "Z"])


def test_rooting_on_an_outgroup_needs_one_named(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="--root outgroup needs --outgroup"):
        _draw(tmp_path, NEWICK, "unnamed", ["--root", "outgroup"])


def test_two_records_of_one_species_stay_two_tips(tmp_path: Path) -> None:
    """The script strips `_NC_...` from a tip, which here would leave one name for two tips."""
    newick = "((Iris_domestica_NC_1.1:0.1,Iris_domestica_NC_2.1:0.1):0.2,(C:0.1,D:0.1):0.25);"
    drawn = _draw(tmp_path, newick, "same", [])
    tips = sorted(terminal.name for terminal in drawn.get_terminals())
    assert tips[:2] == ["C", "D"]
    assert len(set(tips)) == 4, tips


def test_a_tip_is_named_as_the_tree_file_spells_it(tmp_path: Path) -> None:
    """The script draws `_` as a space, and a user copies a name from the tree file."""
    newick = "((NC_1.1:0.1,NC_2.1:0.1):0.2,(C:0.1,D:0.1):0.25,(E:0.1,F:0.12):0.3);"
    drawn = _draw(tmp_path, newick, "spelt", ["--outgroup", "NC_1.1,NC_2.1"])
    assert {"NC_1.1", "NC_2.1"} in _root_sides(drawn)
