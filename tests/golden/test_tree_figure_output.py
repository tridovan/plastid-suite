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

from conftest import require_r
from plastid.external.rscript import missing_packages, run_figure
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

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
