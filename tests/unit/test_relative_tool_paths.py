"""A path typed at the prompt keeps working once a tool runs in a directory of its own.

Measured 01/11/2026: `plastid ir-junction --cpjsdraw relative/CPJSdraw.pl` failed with "Can't open
perl script". The file was there; the tool's working directory was not the user's.
"""

from __future__ import annotations

from pathlib import Path

from plastid.cli import _checked_sets, main
from plastid.external.cpjsdraw import cpjsdraw_script

CPJSDRAW = '#!/usr/bin/perl\nmy $version = "0.0.1";\n'


def test_a_relative_tool_file_is_absolute_before_a_module_reads_it(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.chdir(tmp_path)
    params = _checked_sets(["misa=tools/misa.pl"], ["ssr"])
    assert params == {"misa": str(tmp_path / "tools" / "misa.pl")}


def test_a_relative_tool_directory_is_absolute_before_a_module_reads_it(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.chdir(tmp_path)
    params = _checked_sets(["mummer=tools/bin"], ["repeats"])
    assert params == {"mummer": str(tmp_path / "tools" / "bin")}


def test_a_value_that_names_no_path_is_left_as_it_was(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert _checked_sets(["scope=cds"], ["substitutions"]) == {"scope": "cds"}


def test_a_relative_cpjsdraw_is_absolute_when_the_adapter_hands_it_on(
    tmp_path: Path, monkeypatch
) -> None:
    (tmp_path / "CPJSdraw.pl").write_text(CPJSDRAW, encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    assert cpjsdraw_script("CPJSdraw.pl") == tmp_path / "CPJSdraw.pl"


def test_figures_takes_a_relative_rscript(tmp_path: Path, data_dir: Path, monkeypatch) -> None:
    """The figure pass runs R inside a scratch directory, so the path has to survive the move."""
    out = tmp_path / "runs"
    assert (
        main(["run", "-m", "codon_usage", "-i", str(data_dir / "NC_050833.1.gb"), "-o", str(out)])
        == 0
    )
    run_dir = next(out.glob("*/*/manifest.json")).parent
    script = tmp_path / "r.sh"
    script.write_text(
        '#!/bin/sh\n[ "$1" = --version ] && echo "Rscript (R) version 4.3.3"\nexit 0\n',
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    assert main(["figures", str(run_dir), "--figure", "r_rscu", "--set", "rscript=r.sh"]) != 2
    logged = (run_dir / "00_run_info" / "logs").glob("*__rscript_rscu__stdout.log")
    assert list(logged), "a relative rscript did not survive the move into the scratch directory"
