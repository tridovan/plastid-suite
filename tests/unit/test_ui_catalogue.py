"""Everything the command line can do has a place in the interface, read from the parser itself."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pytest

from plastid.cli import build_parser
from plastid.modules.registry import all_modules
from plastid.run.state import EXIT_STATES
from plastid.ui import catalogue, commands
from plastid.ui.catalogue import COMMANDS, FROZEN, PAGES, page_of_module


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    found = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    return dict(found[0].choices) if found else {}


def _leaves(parser: argparse.ArgumentParser) -> set[str]:
    leaves: set[str] = set()
    for name, child in _subparsers(parser).items():
        below = _subparsers(child)
        leaves |= {f"{name} {action}" for action in below} if below else {name}
    return leaves


def _uncovered(leaves: set[str], covered: dict[str, str]) -> set[str]:
    return leaves - set(covered) - set(FROZEN)


def test_every_command_has_a_page() -> None:
    assert _uncovered(_leaves(build_parser()), COMMANDS) == set()


def test_the_check_sees_a_command_with_no_page() -> None:
    assert _uncovered({"run", "a new command"}, {"run": "analyse"}) == {"a new command"}


def test_a_page_named_for_a_command_exists() -> None:
    assert set(COMMANDS.values()) <= set(PAGES)


def test_every_module_has_a_page_but_the_frozen_two() -> None:
    pages = {module: page_of_module(module) for module in all_modules()}
    assert {module for module, page in pages.items() if page is None} == {
        "mito_annotate",
        "mito_genes",
    }
    assert {page for page in pages.values() if page} <= set(PAGES)


def test_every_command_that_does_work_can_be_built() -> None:
    """`run` builds `batch` too; the lists are what the server reads, not what it starts."""
    not_started = {"ui", "modules", "tools", "runs", "batch"}
    assert set(commands.BUILDERS) == set(COMMANDS) - not_started


def test_the_page_names_every_way_a_run_can_end() -> None:
    """A status the page does not know would show as a raw word, or as nothing at all."""
    page = (Path(catalogue.__file__).parent / "assets" / "app.js").read_text(encoding="utf-8")
    known = re.search(r"const STATUS = \{(.*?)\};", page, re.DOTALL)
    assert known is not None
    for state in EXIT_STATES:
        assert f"{state.status}:" in known.group(1), state.status


@pytest.mark.parametrize("frozen", FROZEN)
def test_the_frozen_commands_stay_out_of_the_table(frozen: str) -> None:
    assert frozen not in COMMANDS
