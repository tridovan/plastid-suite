"""Comments and docstrings state reasons, not history (CLAUDE.md, Code rules: Comments).

Two shapes are refused inside them. A date narrates when something happened rather than why a
rule exists. A planning code (a milestone such as M11.16, a decision such as D5-09, a finding such
as H-21, a task such as P3B-10) points at maintainer records this repository does not carry, so a
reader of the public source cannot follow it. Known-issue ids (K-07) stay: they name entries of
docs/known_issues.md, which ships with the code.
"""

from __future__ import annotations

import ast
import io
import re
import tokenize
from collections.abc import Iterator
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src" / "plastid"

PLANNING_CODE = re.compile(
    r"\b(?:M\d+[A-Z]?\.\d+[a-z]?|D\d+[A-Z]?-\d+[a-z]?|D-[a-z]\b|H-\d+|P\d+[A-Z]?-\d+)"
)
DATE = re.compile(r"\b(?:[0-3]\d/[01]\d/20\d\d|20\d\d-[01]\d-[0-3]\d)\b")


def _comments_and_docstrings(path: Path) -> Iterator[tuple[int, str]]:
    source = path.read_text(encoding="utf-8")
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            yield token.start[0], token.string
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                for offset, line in enumerate(first.value.value.splitlines()):
                    yield first.lineno + offset, line


def test_no_comment_or_docstring_carries_a_date_or_a_planning_code() -> None:
    found = [
        f"{path.relative_to(SRC.parent)}:{line}: {match.group(0)!r}"
        for path in sorted(SRC.rglob("*.py"))
        for line, text in _comments_and_docstrings(path)
        for match in (*PLANNING_CODE.finditer(text), *DATE.finditer(text))
    ]
    assert not found, "say the reason instead:\n" + "\n".join(found)


def test_the_scan_sees_both_shapes_and_spares_known_issue_ids(tmp_path: Path) -> None:
    """A negative control: the patterns bite on what they refuse and pass a K- id."""
    sample = tmp_path / "sample.py"
    sample.write_text(
        '"""Measured 09/10/2026 (M11.16)."""\n# D5-09 says so, and K-07 is a known issue\nx = 1\n',
        encoding="utf-8",
    )
    texts = [text for _, text in _comments_and_docstrings(sample)]
    codes = [match.group(0) for text in texts for match in PLANNING_CODE.finditer(text)]
    dates = [match.group(0) for text in texts for match in DATE.finditer(text)]
    assert sorted(codes) == ["D5-09", "M11.16"]
    assert dates == ["09/10/2026"]
