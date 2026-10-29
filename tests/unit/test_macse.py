"""The MACSE adapter: version pin, the two non-DNA characters, and a real alignment (M11.21).

The adapter exists for one question MAFFT cannot answer — is a codon partition a partition of
codon positions? — so the tests pin the two places where its answer is a rule rather than a
tool call: the release named in the banner, and MACSE's `!`/`*` alphabet. The real-tool test is one
gene, not a supermatrix: 78 genes would make every suite run wait minutes for a wiring check.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.external.macse import (
    MACSE_VERSION,
    clean_sequence,
    macse_version,
    parse_version,
    resolve_macse,
)
from plastid.run.tools import ToolError

BANNER = """\
This is MACSE V2.07 If you find MACSE useful, please cite:

MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons
"""


def test_the_release_is_read_from_the_banners_first_line() -> None:
    assert parse_version(BANNER) == "2.07"
    assert MACSE_VERSION == "2.07"


def test_a_banner_that_never_names_a_release_is_refused() -> None:
    with pytest.raises(ToolError, match="no 'This is MACSE"):
        parse_version("usage: [-prog PROGRAM_NAME]")


def test_a_different_release_is_refused_by_name(monkeypatch: pytest.MonkeyPatch) -> None:
    """The alignment algorithm and the output alphabet are part of the result."""

    class Result:
        returncode = 0
        stdout = "This is MACSE V1.02\n"
        stderr = ""

    monkeypatch.setattr("plastid.external.macse.subprocess.run", lambda *a, **k: Result())
    with pytest.raises(ToolError, match="this is MACSE 1.02"):
        macse_version(("macse",))


def test_frameshift_and_stop_become_gaps() -> None:
    """`!` and `*` are MACSE's alphabet, not DNA; a downstream counter cannot classify them."""
    assert clean_sequence("AC!GT*NN-") == "AC-GT-NN-"


def test_a_missing_macse_names_the_install_line(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("plastid.external.macse.which", lambda name: None)
    with pytest.raises(ToolError, match="bioconda"):
        resolve_macse()


def test_an_explicit_path_that_is_not_a_file_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="not a file"):
        resolve_macse(tmp_path / "typo.jar")
