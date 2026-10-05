"""The evidence archive stays citable — every `§8.N` and `§9.X` the repo quotes must exist.

`PROJECT.md` §8 (the log of every real run) was split into `docs/log/` on 04/10/2026: 33 entries in
one 118 KB file, so reading §8.27 meant opening a 216 KB document. The split **kept every entry
number**, because **271 references** point at them by number — 27 in code docstrings, 36 in tests,
23 in `docs/`, 115 across the four session files, and 70 inside `PROJECT.md` itself. Renumbering was
never an option: a docstring cannot be told to re-derive "§8.27".

This test is what keeps that promise. It is the same shape as the rule it serves — a citation that
resolves to nothing is a claim nobody can check, and it must fail loudly rather than read as
authoritative prose.

Two numbering schemes, deliberately checked differently:

* **`§8.N`** is a section number. Each one must be a heading in exactly one file under `docs/log/`,
  and appear in `PROJECT.md` §8's index.
* **`§9.X`** — `X` a letter — is a subsection. It must be a `### X.` heading in `PROJECT.md`.
* **`§9.N`** — `N` a number — is an *item* inside §9, numbered continuously across its subsections,
  so it is matched against the numbered items themselves. That scheme is weaker by construction
  (inserting an item renumbers the ones after it) and this test says so rather than pretending
  otherwise.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]
PROJECT = ROOT / "PROJECT.md"
LOG_DIR = ROOT / "docs" / "log"

# Every file that may quote the archive: the session documents, docs/, and the source tree.
# This module is excluded because it *describes* the two numbering schemes — its docstring writes
# `§8.N` and `§9.X` as shapes, not as citations of an entry that has to exist.
SEARCHED = sorted(
    path
    for path in [
        *(ROOT.glob("*.md")),
        *(ROOT / "docs").rglob("*.md"),
        *(ROOT / "src").rglob("*.py"),
        *(ROOT / "tests").rglob("*.py"),
    ]
    if path != Path(__file__)
)

ENTRY = re.compile(r"^### (8\.\d+)\.", re.MULTILINE)
SUBSECTION = re.compile(r"^### ([A-Z])\.", re.MULTILINE)
ITEM = re.compile(r"^(\d+)\. \*\*", re.MULTILINE)

CITES_8 = re.compile(r"§\s*(8\.\d+)")  # the whole "8.27", so it compares with a heading key
CITES_9 = re.compile(r"§\s*9\.([A-Z0-9]+)")


def _log_entries() -> dict[str, str]:
    """Every `§8.N` heading and its file. One number in two files means the split went wrong."""
    found: dict[str, str] = {}
    for path in sorted(LOG_DIR.glob("*.md")):
        for number in ENTRY.findall(path.read_text(encoding="utf-8")):
            assert number not in found, (
                f"§{number} appears in both {found[number]} and {path.name} — one entry, one file"
            )
            found[number] = path.name
    return found


def _citations(pattern: re.Pattern[str], paths: list[Path] | None = None) -> dict[str, list[str]]:
    """Citation → the files that make it, so a failure can name where the broken reference is.

    ``paths`` is injectable so the check can be shown failing on a citation that resolves to
    nothing — a guard nobody has watched fail is not a guard.
    """
    seen: dict[str, list[str]] = {}
    for path in SEARCHED if paths is None else paths:
        text = path.read_text(encoding="utf-8")
        where = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        for hit in pattern.findall(text):
            seen.setdefault(hit, []).append(where)
    return seen


def test_the_log_entries_did_not_move_house_twice() -> None:
    """One entry per number, none lost — 33 from the 04/10 split plus phase 5's first entry."""
    entries = _log_entries()
    assert len(entries) == 34, f"expected 34 logged entries, found {len(entries)}"
    assert set(entries) >= {"8.1", "8.31", "8.34"}, "the first, the corpus and the newest entry"


def test_every_cited_log_entry_exists() -> None:
    """`§8.N` is a section number, so it must resolve to a heading — wherever that heading lives."""
    entries = _log_entries()
    broken = {
        number: files for number, files in _citations(CITES_8).items() if number not in entries
    }
    assert not broken, f"cited but not written: {broken}"


def test_project_md_indexes_every_entry() -> None:
    """The index in §8 is how a `§8.N` citation resolves for a reader; a gap makes it a dead end."""
    index = PROJECT.read_text(encoding="utf-8")
    entries = _log_entries()
    in_index = lambda number: bool(  # noqa: E731 — a one-line predicate reads better inline
        re.search(rf"^\| {re.escape(number)} \|", index, re.M)
    )
    missing = sorted(
        (number for number in entries if not in_index(number)),
        key=lambda n: int(n.split(".")[1]),
    )
    assert not missing, f"§8 entries missing from the index in PROJECT.md: {missing}"

    files = {LOG_DIR / name for name in entries.values()}
    for name in (path.name for path in files):
        assert f"docs/log/{name}" in index, f"{name} is not linked from the index"


def test_every_cited_error_subsection_exists() -> None:
    """`§9.A` … `§9.I` are subsections of `PROJECT.md` §9 and must still be there."""
    body = PROJECT.read_text(encoding="utf-8")
    start = body.index("## 9.")
    section = body[start : body.index("## 10.", start)]
    letters = set(SUBSECTION.findall(section))
    cited = {letter[:1] for letter, files in _citations(CITES_9).items() if letter[:1].isalpha()}
    assert cited <= letters, f"cited but not written: {sorted(cited - letters)} across §9"


def test_every_cited_error_item_exists() -> None:
    """`§9.14` is an *item* inside §9, numbered continuously — checked, and weak by construction."""
    body = PROJECT.read_text(encoding="utf-8")
    start = body.index("## 9.")
    items = set(ITEM.findall(body[start : body.index("## 10.", start)]))
    cited = {number for number in _citations(CITES_9) if number.isdigit()}
    missing = sorted(cited - items, key=int)
    assert not missing, (
        f"§9 items cited but not present: {missing} — §9's items are numbered by position, so an "
        "insertion in the middle renumbers every citation after it"
    )


@pytest.mark.parametrize(
    "name",
    ["phase-1.md", "phase-2-3.md", "phase-3a.md", "phase-4.md", "phase-4a.md", "phase-5.md"],
)
def test_each_log_file_says_which_numbers_it_holds(name: str) -> None:
    """A reader who lands in the file (not the index) must be able to tell what it is."""
    text = (LOG_DIR / name).read_text(encoding="utf-8")
    assert text.startswith("# Nhật ký đo lường —"), name
    assert "PROJECT.md` §8" in text, f"{name} does not say it is part of PROJECT.md §8"
    assert re.search(r"§8\.\d+–§?8\.\d+|§8\.\d+–8\.\d+", text) or "§8." in text


def test_a_dangling_citation_is_caught(tmp_path: Path) -> None:
    """The guard has to fail on a citation that resolves to nothing.

    The project's own rule for gates — "a gate that cannot fail loudly is not done" — applies to a
    documentation check too: one that only ever passes is indistinguishable from no check at all.
    """
    note = tmp_path / "note.md"
    note.write_text("see §8.99 and §9.44 for the measurement\n", encoding="utf-8")

    cited_8 = _citations(CITES_8, [note])
    assert cited_8 == {"8.99": [str(note)]}
    assert "8.99" not in _log_entries(), "a citation to an entry nobody wrote must not resolve"

    items = _citations(CITES_9, [note])
    assert "44" in items, "§9.44 is read as an item number, not a subsection"
