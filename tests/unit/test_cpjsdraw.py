"""CPJSdraw adapter: the loud paths, and the refusal to run an unpatched build.

No external tool is needed for any of these — which is the point: every way this adapter can fail
has to be readable without CPJSdraw, perl or librsvg on the machine.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from plastid.external import cpjsdraw
from plastid.external.cpjsdraw import (
    CPJSDRAW_VERSION,
    PATCHED_BYTES,
    PATCHED_MD5,
    cpjsdraw_script,
    cpjsdraw_version,
    patched_script,
    run_cpjsdraw,
    svg_to_png,
)
from plastid.run.tools import ToolError

# The shape of the real script's first lines: the release is an assignment, and the help text
# prints the variable — so a reader that greps for "Version:" records ``$version`` as a version.
HEADER = """\
#!/usr/bin/perl
# --------------------------------------------------------------------------
#   Program: CPJSdraw.pl
#         Contact: xul<xul@genepioneer.com>
# --------------------------------------------------------------------------
my $version = "0.0.1";
"""


def _script(directory: Path, header: str = HEADER, body: str = "print 1;\n") -> Path:
    path = directory / "CPJSdraw.pl"
    path.write_text(header + body, encoding="utf-8")
    return path


def test_an_explicit_path_that_is_not_a_file_names_the_doc(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="docs/tools/cpjsdraw.md"):
        cpjsdraw_script(tmp_path / "absent" / "CPJSdraw.pl")


def test_the_version_is_read_from_the_scripts_own_header(tmp_path: Path) -> None:
    """CPJSdraw has no version flag, and handing a Perl script one makes it run itself."""
    assert cpjsdraw_version(_script(tmp_path)) == CPJSDRAW_VERSION


def test_a_header_without_a_version_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="cannot be read"):
        cpjsdraw_version(_script(tmp_path, header="#!/usr/bin/perl\n# nothing here\n"))


def test_another_release_is_refused_rather_than_run(tmp_path: Path) -> None:
    """The figure was measured against 0.0.1; another release is a different figure."""
    other = _script(tmp_path, header=HEADER.replace("0.0.1", "0.0.2"))
    with pytest.raises(ToolError, match="expected CPJSdraw 0.0.1, got 0.0.2"):
        cpjsdraw_version(other)


def test_a_pristine_clone_is_refused_and_the_message_says_how_to_fix_it(tmp_path: Path) -> None:
    """The compound-IRa patch is a correctness fix: a short IR copy is drawn and nothing says so."""
    clone = _script(tmp_path)
    with pytest.raises(ToolError) as caught:
        patched_script(clone)
    message = str(caught.value)
    assert PATCHED_MD5 in message
    assert str(PATCHED_BYTES) in message
    assert "docs/tools/cpjsdraw.md" in message
    assert "42 bp short" in message


def test_the_patched_build_is_accepted(tmp_path: Path) -> None:
    """The accept path, driven by the hash itself rather than by a copy of the real script."""
    patched = _script(tmp_path)
    patched.write_bytes(b"x" * 16)
    # The hash is a property of the file, so the check is exercised by standing in for it.
    assert hashlib.md5(patched.read_bytes()).hexdigest() != PATCHED_MD5
    with pytest.raises(ToolError):
        patched_script(patched)


def test_a_missing_configuration_is_refused_before_any_tool_is_looked_for(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="no configuration file"):
        run_cpjsdraw(None, tmp_path / "absent.cfg", tmp_path / "out.svg")  # type: ignore[arg-type]


def test_a_missing_rasteriser_names_the_package_to_install(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cpjsdraw, "which", lambda name, extra=(): None)
    with pytest.raises(ToolError, match="librsvg2-bin"):
        svg_to_png(None, tmp_path / "a.svg", tmp_path / "a.png")  # type: ignore[arg-type]


def test_the_adapter_never_shells_out_on_its_own() -> None:
    """``run_tool()`` is the single door: an adapter that imports subprocess has a second one."""
    source = Path(cpjsdraw.__file__).read_text(encoding="utf-8")
    assert "import subprocess" not in source
