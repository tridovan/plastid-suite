"""Smoke tests: the package imports and the CLI answers."""

from __future__ import annotations

import pytest

import plastid
from plastid.cli import main


def test_version_exposed() -> None:
    assert plastid.__version__
    assert plastid.__version__ != "0.0.0+unknown"  # the test env installs the package


def test_cli_version(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--version"]) == 0
    assert capsys.readouterr().out.startswith("plastid ")


def test_cli_no_args_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    assert main([]) == 0
    assert "usage: plastid" in capsys.readouterr().out


def test_cli_run_stub_not_implemented(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["run"]) == 1
    assert "not implemented" in capsys.readouterr().out
