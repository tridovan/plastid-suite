"""`plastid tools` is where a user learns what to do about a missing tool.

A bare "not found on PATH" leaves them to guess the package, the version and whether it matters.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.external.inventory import check_tools

TOOL_DIRECTORIES = (
    "PTS_ASSEMBLY_DIR",
    "PTS_BLAST_DIR",
    "PTS_IQTREE_DIR",
    "PTS_MISA_DIR",
    "PTS_MUMMER_DIR",
    "PTS_PAML_DIR",
)


@pytest.fixture
def bare_machine(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("PATH", str(tmp_path))
    for name in TOOL_DIRECTORIES:
        monkeypatch.delenv(name, raising=False)


@pytest.mark.usefixtures("bare_machine")
def test_a_missing_tool_is_named_with_the_way_to_get_it() -> None:
    statuses = check_tools()
    assert all(status.status == "missing" for status in statuses), statuses
    for status in statuses:
        assert "docs/" in status.message, status.message
        assert not status.message.endswith("not found on PATH"), status.message


@pytest.mark.usefixtures("bare_machine")
def test_macse_says_it_can_be_given_as_a_jar() -> None:
    macse = next(status for status in check_tools() if status.name == "macse")
    assert "--set macse=" in macse.message


@pytest.mark.usefixtures("bare_machine")
def test_a_tool_only_some_runs_ask_for_is_marked(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["tools", "--json"]) == 0
    tools = {tool["name"]: tool for tool in json.loads(capsys.readouterr().out)["tools"]}
    assert tools["macse"]["used_by"] == ["phylogeny (if asked)"]
    assert tools["table2asn"]["used_by"] == ["submission (if asked)"]
    assert "phylogeny" in tools["iqtree"]["used_by"]
