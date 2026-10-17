"""The R adapter: the loud paths, and the scripts arriving with the package.

Nothing here needs R — which is the point. Every way this adapter can fail has to be readable on a
machine that has no R at all, and the one test that draws a figure is ``tier_b``.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from plastid.external.rscript import (
    R_PACKAGES,
    R_SCRIPTS,
    figure_script,
    missing_packages,
    rscript_binary,
    run_figure,
)
from plastid.run.tools import ToolError


def test_an_explicit_path_that_is_not_a_file_names_the_doc(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="docs/tools/r.md"):
        rscript_binary(tmp_path / "absent" / "Rscript")


def test_an_unknown_figure_lists_the_ones_there_are() -> None:
    with pytest.raises(ToolError, match="no figure 'nope'"):
        figure_script("nope")


def test_every_figure_the_adapter_names_ships_with_the_package() -> None:
    """The scripts are package data: a missing one is a packaging bug, not a run-time surprise."""
    for name in R_SCRIPTS:
        path = figure_script(name)
        assert path.is_file(), f"{name}: {path} is not in the package"
        assert path.parent.name == "r_scripts"


def test_the_shared_helper_is_carried_and_makes_input_paths_absolute() -> None:
    """The one change to the bundle's scripts, and the reason it exists.

    ``khoi_tao()`` calls ``setwd(OUT_DIR)``, so a relative ``--inputdir`` would point into the
    output directory. Measured 17/10/2026: every script fails to find its own data that way.
    """
    helper = figure_script("rscu").parent / "00_ham_chung.R"
    text = helper.read_text(encoding="utf-8")
    assert "setwd(OUT_DIR)" in text
    assert "normalizePath(IN_DIR" in text
    assert text.index("normalizePath(IN_DIR") < text.index("setwd(OUT_DIR)")


@pytest.mark.tier_b
def test_the_nine_packages_are_all_present_here() -> None:
    r = os.environ.get("PTS_RSCRIPT")
    if r is None and not _on_path("Rscript"):
        pytest.skip(f"no Rscript: set PTS_RSCRIPT=/path/to/Rscript (needs {', '.join(R_PACKAGES)})")
    ctx = _StubContext()
    assert missing_packages(ctx, rscript_binary(r)) == ()


@pytest.mark.tier_b
def test_every_script_runs_and_writes_its_figure(tmp_path: Path) -> None:
    """The whole engine, over the bundle's own sample tables — all seven, each with a PNG.

    Measured 17/10/2026 on this machine: 22 files in about 37 s, and every script writes a PNG and
    a PDF. The bundle's author verified the same seven on R 4.1.2; this machine has 4.3.3, so the
    sweep is also the evidence that the figures do not move with the interpreter.
    """
    r = os.environ.get("PTS_RSCRIPT")
    if r is None and not _on_path("Rscript"):
        pytest.skip("no Rscript: set PTS_RSCRIPT=/path/to/Rscript")
    ctx = _StubContext()
    data = _sample_data()
    missing: list[str] = []
    for name in R_SCRIPTS:
        out = tmp_path / name
        try:
            figures = run_figure(ctx, name, input_dir=data, out_dir=out, binary=r)
        except ToolError as error:  # a machine with R but without the nine packages
            pytest.skip(str(error))
        if not any(path.suffix == ".png" for path in figures.iterdir()):
            missing.append(name)
    assert not missing, f"no PNG from: {', '.join(missing)}"


def _on_path(name: str) -> bool:
    from shutil import which

    return which(name) is not None


def _sample_data() -> Path:
    """The bundle's own sample tables, when they are on this machine.

    They are the team's material and gitignored, so this skips rather than fails: the point of the
    test is that the engine runs, not that a particular machine holds the paper's data.
    """
    root = (
        Path(__file__).parents[2]
        / "feedback"
        / "goi_tool_hinh_IR_2026-10-16"
        / "01_R_scripts"
        / "data_mau"
    )
    if not root.is_dir():
        pytest.skip(f"the bundle's sample tables are not here: {root}")
    return root


class _StubContext:
    """A run context's tool door, without a run — the adapter only needs ``ctx.tool``."""

    def __init__(self) -> None:
        self.tools: list[object] = []

    def tool(self, name: str, argv, **kwargs):  # type: ignore[no-untyped-def]
        import tempfile

        from plastid.run.tools import run_tool

        record = run_tool(name, argv, log_dir=tempfile.mkdtemp(prefix="rscript_test_"), **kwargs)
        self.tools.append(record)
        return record
