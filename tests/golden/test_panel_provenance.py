"""The private IR validation panel — provenance and integrity, checked before it is believed.

**Why this file exists.** `.panel/ir-validation/` is 170 records (~250 MB) that live outside git
on purpose, and every conclusion drawn from them depends on their being *the* records that were
pinned. Three things can go wrong silently and none of them is visible from a test that passes:

* a record's bytes change after the pin — a re-export, an editor, a partial copy;
* the aggregate `PINNED.json` stops describing a group directory (a group added by hand, a file
  dropped) so a run reads a record no pin ever covered;
* the panel is described as "golden". It is third-party data with its own source and date, and
  mixing it with the pinned fixture would make one input's numbers look like another's.

A missing panel is **not** a failure — it is gitignored by design and a fresh clone has none.
These checks skip with a named reason there, exactly as the tier-B tool checks do.

The measurements themselves (identity, gate verdict) are not re-derived here; they are in
`tools/validate_ir_panel.py` and in `measurements/`, and re-measuring them in a test would make a
slow, duplicated oracle. This file checks the *bytes and the bookkeeping* a number rests on.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "ir-validation"
SOURCE = PANEL / "source"
INDEX = PANEL / "PINNED.json"

pytestmark = pytest.mark.panel


def _require_panel() -> None:
    if not INDEX.is_file():
        reason = f"the IR validation panel is not in this tree ({PANEL} — .panel/ is gitignored)"
        if os.environ.get("PTS_REQUIRE_PANEL") == "1":
            pytest.fail(reason)
        pytest.skip(reason)


def _index() -> dict:
    return json.loads(INDEX.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def test_every_source_record_still_hashes_to_its_pin() -> None:
    """The pin's sha256 is the claim; a record that no longer hashes to it is a different record.

    Measured read-back, not inferred from the directory listing: a file present with the right
    name and the wrong bytes is exactly the failure a listing cannot see.
    """
    _require_panel()
    genomes = _index()["groups"]["source"]["genomes"]
    assert len(genomes) == 80, f"the source group is {len(genomes)} records, not the pinned 80"
    changed: list[str] = []
    missing: list[str] = []
    for record in genomes:
        path = ROOT / record["path"]
        if not path.is_file():
            missing.append(record["accession"])
            continue
        if _sha256(path) != record["sha256"]:
            changed.append(record["accession"])
    assert not missing, f"pinned source records missing from the tree: {missing}"
    assert not changed, (
        f"source records whose bytes no longer match the pin: {changed} — a result computed "
        f"from them is not a result from the pinned panel"
    )


def test_the_aggregate_index_describes_every_record_directory() -> None:
    """A record directory the aggregate index does not describe is data nobody pinned.

    The aggregate index has two shapes on purpose: `groups` holds the two directories
    `plastid panel pin --file` wrote (each with its own `PINNED.json`), and the flat `records`
    list covers all four record groups, including the generated `fasta` and `synthetic` sets
    that the pin command cannot read as records. A directory outside both is the failure —
    added by hand, or left behind after an index rebuild.
    """
    _require_panel()
    index = _index()
    described = {Path(record["path"]).parent.name for record in index["records"]}
    assert described == {"source", "annotation-variants", "fasta", "synthetic"}, (
        f"the aggregate records list describes {sorted(described)}"
    )
    assert set(index["groups"]) <= described, (
        f"a pinned group is absent from the records list: {set(index['groups']) - described}"
    )
    for directory in _record_directories():
        if directory.name in described:
            continue
        assert (directory / "PINNED.json").is_file(), (
            f"{directory}: not in the aggregate index and no PINNED.json of its own — a "
            f"record set with no provenance cannot be a validation group"
        )
    for name in described:
        directory = PANEL / name
        assert directory.is_dir(), f"the index describes {name}, which is not on disk"
        assert (directory / "PINNED.json").is_file() or (directory / "EXPECTED.tsv").is_file(), (
            f"{directory}: neither PINNED.json nor EXPECTED.tsv — a generated set needs its "
            f"expected answer beside it"
        )


def test_the_panel_is_not_called_golden_and_carries_its_source() -> None:
    """Provenance, as text a later reader can act on: the release, the date, the sha256s.

    The rule this holds is the project's own: a third-party panel is never merged with the golden
    pair in one table (M11.18), and its records carry source · date · digest.
    """
    _require_panel()
    readme = (PANEL / "README.md").read_text(encoding="utf-8")
    for expected in ("2026-09-03", "sha256", "ftp.ncbi.nlm.nih.gov"):
        assert expected in readme, f"the panel README no longer records {expected!r}"
    assert "golden" not in readme.lower().replace("the golden pair", ""), (
        "the panel README calls itself or its records golden — the panel is third-party "
        "validation data, and the golden pair is one pinned input"
    )
    manifest = (PANEL / "manifest.tsv").read_text(encoding="utf-8")
    header = manifest.splitlines()[0].split("\t")
    for column in ("accession", "source_url_or_command", "retrieved_date", "sha256"):
        assert column in header, f"manifest.tsv lost the {column!r} column"


def test_the_manifest_rows_match_the_records_on_disk() -> None:
    """170 rows, one per record — the number the roadmap quotes as the panel's size.

    Checked by name set rather than by counting files: the pinned 80 are `.gb`, the rotated set is
    `.fa`, and the synthetic set is `.gb` again, so a single-extension scan would silently drop a
    whole group while leaving a plausible total.
    """
    _require_panel()
    manifest = (PANEL / "manifest.tsv").read_text(encoding="utf-8").splitlines()
    header = manifest[0].split("\t")
    column = header.index("accession")
    rows = {line.split("\t")[column] for line in manifest[1:] if line.strip()}
    assert len(manifest) - 1 == 170, (
        f"manifest.tsv has {len(manifest) - 1} rows, not the pinned 170"
    )
    on_disk = {
        path.stem
        for directory in _record_directories()
        for path in directory.iterdir()
        if path.suffix in {".gb", ".fa", ".fasta"}
    }
    assert not rows - on_disk, f"rows whose record is not on disk: {sorted(rows - on_disk)}"
    # `pp078674/` is a record pinned on its own (its own PINNED.json, a candidate second copy
    # outside the 170-row design), so it is described by that pin rather than by manifest.tsv.
    separate = (PANEL / "pp078674").is_dir()
    unlisted = on_disk - rows - ({"PP078674.1"} if separate else set())
    assert not unlisted, f"records with no manifest row: {sorted(unlisted)}"


def _record_directories() -> list[Path]:
    """Every directory of the panel that holds records, whichever index describes it."""
    return sorted(
        path
        for path in PANEL.iterdir()
        if path.is_dir()
        and any(
            child.suffix in {".gb", ".fa", ".fasta"} for child in path.iterdir() if child.is_file()
        )
    )
