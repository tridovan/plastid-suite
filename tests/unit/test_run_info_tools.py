"""`00_run_info/run_info.md` must name **every** tool the run called, not the last module's.

**The defect this file records (K-08, found in M11.8, 28/10/2026).** A run is a list of modules;
each one gets its own `RunContext`, so each one's tool calls live on that context. `_finish`
rendered the Methods page from the *last* context's records only. A run of `genome_stats … ssr`
therefore printed "no external tools ran" while `misa 2.1` had run and been logged, and a run of
`ssr` alone printed the MISA row correctly — which is what made the mistake invisible in a
single-module test.

That page is the Materials & Methods record: it is the place a reader is told which programs
produced an annotation or a count, and a false "(no external tools ran)" is worse than a missing
row because it reads as a positive statement.

The per-module manifest entries were already right (`modules[*].tools`), so the check below reads
the page and the manifest together: the page must name every tool the manifest records.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from plastid.run.runinfo import render_run_info
from plastid.run.tools import ToolRecord

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "data" / "NC_050833.1.gb"


def _record(name: str, version: str) -> ToolRecord:
    return ToolRecord(
        name=name,
        kind="binary",
        version=version,
        version_source="pinned",
        path=f"/usr/bin/{name}",
        argv=(name,),
        cwd=".",
        rc=0,
        started_at="2026-10-28T10:00:00",
        duration_s=0.1,
        stdout_file="out.log",
        stderr_file="err.log",
        stdout_tail="",
    )


def _manifest(modules: list[dict]) -> dict:
    return {
        "run_id": "r1",
        "package_version": "1.0.0",
        "command": ["plastid", "run"],
        "started_at": "2026-10-28T10:00:00",
        "finite": True,
        "status": "ok",
        "modules": modules,
        "params": {},
        "inputs": [],
        "coordinate_system": {"name": "genbank", "rotation_offset": 0},
    }


def test_the_page_names_the_tools_of_every_module_not_only_the_last() -> None:
    """The page is built from the manifest's per-module records, so no context can hide one."""
    tools = [_record("misa", "2.1"), _record("repeat-match", "4.0.1")]
    manifest = _manifest(
        [
            {"id": "ssr", "version": "0.4.0", "status": "ok", "tools": [tools[0]]},
            {"id": "repeats", "version": "0.1.0", "status": "ok", "tools": [tools[1]]},
        ]
    )
    page = render_run_info(manifest)
    assert "| misa | 2.1 |" in page
    assert "| repeat-match | 4.0.1 |" in page
    assert "no external tools ran" not in page


def test_a_tool_called_by_two_modules_is_named_once() -> None:
    """MAFFT is called by five modules in one comparison run; five identical rows are noise."""
    tool = _record("mafft", "7.526")
    manifest = _manifest(
        [
            {"id": "substitutions", "version": "0.3.0", "status": "ok", "tools": [tool]},
            {"id": "diversity", "version": "0.3.0", "status": "ok", "tools": [tool]},
        ]
    )
    page = render_run_info(manifest)
    assert page.count("| mafft | 7.526 |") == 1


def test_a_run_with_no_tool_still_says_so() -> None:
    """The honest negative must survive: a pure-Python run is not a missing row."""
    manifest = _manifest([{"id": "genome_stats", "version": "0.1.0", "status": "ok", "tools": []}])
    assert "no external tools ran" in render_run_info(manifest)


def _require_misa() -> Path:
    """MISA, from ``PTS_MISA_DIR`` — skip with a named reason, or fail under PTS_REQUIRE_TOOLS."""
    directory = os.environ.get("PTS_MISA_DIR")
    candidate = Path(directory) / "misa.pl" if directory else None
    if candidate is not None and candidate.is_file():
        return candidate
    reason = "the end-to-end half needs misa.pl; set PTS_MISA_DIR to the directory holding it"
    if os.environ.get("PTS_REQUIRE_TOOLS") == "1":
        pytest.fail(reason)
    pytest.skip(reason)


def test_a_real_two_module_run_reports_its_tool_on_the_page(tmp_path: Path) -> None:
    """The reproduction, through the CLI, in the order that failed: a tool module, then a
    pure-Python one. The page must still name MISA after `genome_stats` has run last."""
    from plastid.cli import main

    script = _require_misa()
    code = main(
        [
            "run",
            "-m",
            "ssr",
            "-m",
            "genome_stats",
            "-i",
            str(GOLDEN),
            "-o",
            str(tmp_path / "runs"),
            "--dataset",
            "repro",
            "--set",
            f"misa={script}",
        ]
    )
    assert code == 0
    run_dir = next((tmp_path / "runs" / "repro").iterdir())
    page = (run_dir / "00_run_info" / "run_info.md").read_text(encoding="utf-8")
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    recorded = {tool["name"] for module in manifest["modules"] for tool in module.get("tools", [])}
    assert "misa" in recorded, (
        "the manifest did not record the call — check the fixture, not the page"
    )
    assert "| misa | 2.1 |" in page, (
        "the page says no tool ran while the manifest records misa — K-08"
    )
