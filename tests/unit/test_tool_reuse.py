"""One MAFFT call per identical alignment in a run, and the reuse on the record (offline).

A stub stands in for MAFFT: it counts its calls and prints its input as the "alignment", so the
test checks the door, not the aligner.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from plastid.external import mafft
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROWS = [("ref", "acgtacgtac"), ("other", "acgaacgtac")]


@pytest.fixture
def stub_mafft(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A ``mafft`` on PATH that appends one line per call to the returned file."""
    calls = tmp_path / "calls.txt"
    binary = tmp_path / "bin" / "mafft"
    binary.parent.mkdir()
    binary.write_text(
        f'#!/bin/sh\necho call >> "{calls}"\nfor last; do :; done\ncat "$last"\n',
        encoding="utf-8",
    )
    binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{binary.parent}:{os.environ['PATH']}")
    return calls


def _context(layout: RunLayout, module_id: str, order: int, computed: dict) -> RunContext:
    return RunContext(
        layout=layout, params={}, module_id=module_id, module_order=order, computed=computed
    )


def test_a_second_step_reads_the_first_step_s_alignment_and_says_so(
    tmp_path: Path, stub_mafft: Path
) -> None:
    layout = RunLayout.create(tmp_path / "runs", "reuse", fingerprint="0123456789abcdef")
    layout.make_dirs()
    computed: dict = {}
    first = _context(layout, "substitutions", 1, computed)
    second = _context(layout, "diversity", 2, computed)

    aligned_first = mafft.align(ROWS, first)
    aligned_second = mafft.align(ROWS, second)

    assert stub_mafft.read_text().count("call") == 1
    assert aligned_second == aligned_first
    assert len(first.tools) == 1 and first.reused == []
    assert second.tools == []  # this step called no tool itself
    assert second.reused == [
        {
            "name": "mafft",
            "argv": list(first.tools[0].argv),
            "computed_by": "substitutions",
            "stdout_file": first.tools[0].stdout_file,
        }
    ]


def test_other_rows_or_another_run_call_the_aligner_again(tmp_path: Path, stub_mafft: Path) -> None:
    layout = RunLayout.create(tmp_path / "runs", "reuse", fingerprint="0123456789abcdef")
    layout.make_dirs()
    computed: dict = {}
    mafft.align(ROWS, _context(layout, "substitutions", 1, computed))
    mafft.align([ROWS[0], ("other", "acgtacgtaa")], _context(layout, "diversity", 2, computed))
    # A context of its own (another run) shares nothing.
    mafft.align(ROWS, RunContext(layout=layout, params={}, module_id="barcodes", module_order=3))

    assert stub_mafft.read_text().count("call") == 3
