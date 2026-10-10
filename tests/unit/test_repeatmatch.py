"""MUMmer's ``repeat-match`` — the two coordinate systems, and a parse that refuses to guess.

The output is three columns and an ``r`` suffix, and the suffix is the whole trap: it marks the
copy read on the reverse strand, and it sits on that copy's **start in the reverse-complement
frame** — not on the forward coordinate a table wants. The golden's own IR pair is the check:
``83128 153724r 26214`` must convert to ``[83127, 109341)`` and ``[127510, 153724)``, which are
exactly IRa and IRb (§8.58). Getting this wrong gives two intervals that still look plausible.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.core.model import Genome
from plastid.external.repeatmatch import (
    REPEAT_MATCH_VERSION,
    parse_repeat_match,
    repeat_match_version,
    run_repeat_match,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

# The golden's own output, measured 10/10/2026 with `repeat-match -n 30` (MUMmer 4.0.1) — the
# header is the tool's, verbatim.
GOLDEN_OUTPUT = """Long Exact Matches:
   Start1     Start2    Length
    88397     148434r       42
    88418     148455r       42
    88397      88418        42
   111914     111957        32
   111933     111982        31
    29417      29468r       52
   148393     148414        42
    83128     153724r    26214
     5918       5969r       52
"""


def _ctx(tmp_path: Path, module_id: str = "repeats") -> RunContext:
    layout = RunLayout.create(tmp_path / "ctx", "unit", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module_id, module_order=1)


def test_the_ir_pair_converts_to_the_two_ir_copies() -> None:
    """The one row that proves the reverse-coordinate rule: 153724r is the *start* in the RC
    frame, and its forward interval is IRb — not the 26 kb ending at 153,724 that reading it as a
    forward start would give."""
    matches = parse_repeat_match(GOLDEN_OUTPUT)
    ir_pair = next(match for match in matches if match.length == 26_214)
    assert ir_pair.to_1based() == (83_128, 109_341, 127_511, 153_724)
    assert (ir_pair.start1, ir_pair.end1) == (83_127, 109_341)  # IRa, 0-based half-open
    assert (ir_pair.start2, ir_pair.end2) == (127_510, 153_724)  # IRb
    assert ir_pair.palindromic


def test_every_golden_row_parses_and_the_copy_order_is_normalised() -> None:
    matches = parse_repeat_match(GOLDEN_OUTPUT)
    assert len(matches) == 9
    assert sum(match.palindromic for match in matches) == 5
    # a row printed with the reverse copy first comes back with copy 1 upstream, always
    assert all(match.start1 <= match.start2 for match in matches)
    shortest = min(matches, key=lambda match: match.length)
    assert (shortest.length, shortest.start1, shortest.end1) == (31, 111_932, 111_963)


def test_a_self_palindrome_has_both_copies_at_the_same_place() -> None:
    """Two of the golden's matches fold onto themselves — copy 1 *is* copy 2, and saying so is
    the honest reading of a hairpin rather than a table that looks broken."""
    matches = parse_repeat_match(GOLDEN_OUTPUT)
    folded = sorted(
        (match.start1, match.length) for match in matches if match.start1 == match.start2
    )
    assert folded == [(5_917, 52), (29_416, 52)]


def test_a_row_that_will_not_parse_stops_the_run_with_its_line_number() -> None:
    with pytest.raises(ValueError, match="line 3"):
        parse_repeat_match(
            "Long Exact Matches:\n   Start1     Start2    Length\n  88397   ????  42\n"
        )


def test_a_match_below_the_floor_is_refused_rather_than_published(tmp_path: Path) -> None:
    """The floor is what the count means, so a tool disagreeing with it is an error, not a row."""
    stub = _stub(tmp_path / "bin", "  1000  2000  20\n")
    with pytest.raises(ToolError, match="shorter than the 30 bp floor"):
        run_repeat_match(
            Genome(name="synthetic", sequence="ACGT" * 600),
            _ctx(tmp_path),
            min_length=30,
            directory=stub,
        )


def test_the_version_is_read_from_the_sibling_mummer(tmp_path: Path) -> None:
    """repeat-match has no version flag, so its release is its package's — and a mismatch stops."""
    good = _stub(tmp_path / "good", "  1000  2000  40\n", version=REPEAT_MATCH_VERSION)
    assert repeat_match_version(str(good / "repeat-match")) == REPEAT_MATCH_VERSION

    wrong = _stub(tmp_path / "wrong", "  1000  2000  40\n", version="3.23")
    with pytest.raises(ToolError, match="expected MUMmer 4.0.1"):
        repeat_match_version(str(wrong / "repeat-match"))


def test_the_adapter_runs_the_floor_it_was_given(tmp_path: Path) -> None:
    """The argv is the record: `-n 30` must reach the tool, or the count means something else."""
    stub = _stub(tmp_path / "bin", "  1000  2000  40\n")
    ctx = _ctx(tmp_path)
    result = run_repeat_match(
        Genome(name="synthetic", sequence="ACGT" * 600),
        ctx,
        min_length=30,
        directory=stub,
    )
    assert result.min_length_bp == 30
    assert result.version == REPEAT_MATCH_VERSION
    record = ctx.tools[0]  # the provenance the module's Materials & Methods reads
    assert record.name == "repeat-match" and record.rc == 0
    assert "-n" in record.argv and "30" in record.argv
    assert record.version == REPEAT_MATCH_VERSION and record.version_source == "pinned"


def _stub(
    directory: Path, body: str, *, version: str = REPEAT_MATCH_VERSION, echo_argv: bool = False
) -> Path:
    """A ``repeat-match`` that prints canned output, beside a ``mummer`` that answers --version."""
    directory.mkdir(parents=True, exist_ok=True)
    program = directory / "repeat-match"
    printed = 'echo "$@"\n' if echo_argv else ""
    program.write_text(
        "#!/bin/sh\nLong Exact Matches:\n   Start1     Start2    Length\n"
        + printed
        + f"printf '%s' '{body}'\n",
        encoding="utf-8",
    )
    program.chmod(program.stat().st_mode | stat.S_IXUSR)
    sibling = directory / "mummer"
    sibling.write_text(f'#!/bin/sh\necho "{version}"\n', encoding="utf-8")
    sibling.chmod(sibling.stat().st_mode | stat.S_IXUSR)
    return directory
