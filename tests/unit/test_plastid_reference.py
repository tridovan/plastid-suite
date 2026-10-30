"""The plastid reference set: the pin, its integrity, and the leave-one-out rule (P3B-01).

Two halves, tested separately on purpose. The *rule* is arithmetic over accessions and must hold
with nothing on disk; the *pin* is 28 MB of gitignored records, so the test that re-hashes it
skips with a named reason when the working tree does not carry it — and the committed manifest
(`data/plastid_reference.json`) is what keeps the design readable without the bytes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.core.panel import pin_files
from plastid.core.plastid_reference import (
    ReferenceSet,
    build_manifest,
    leave_one_out,
    read_manifest,
    verify_pin,
    write_manifest,
)

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "plastid-reference"
MANIFEST = ROOT / "src" / "plastid" / "data" / "plastid_reference.json"
CORPUS = ROOT / "tests" / "data" / "corpus"


def test_leave_one_out_excludes_every_target_from_its_own_reference_set() -> None:
    """The rule in one line: target t is annotated by the set minus t, and by nothing else."""
    sets = leave_one_out(["NC_3", "NC_1", "NC_2"])
    assert sets == {
        "NC_1": ("NC_2", "NC_3"),
        "NC_2": ("NC_1", "NC_3"),
        "NC_3": ("NC_1", "NC_2"),
    }
    for target, references in sets.items():
        assert target not in references
        assert len(references) == 2


def test_leave_one_out_refuses_a_repeated_accession() -> None:
    """Two files under one accession cannot be told apart by "minus itself" — so it stops."""
    with pytest.raises(ValueError, match="more than once"):
        leave_one_out(["NC_1", "NC_2", "NC_1"])


def test_a_target_outside_the_set_is_annotated_by_the_whole_set() -> None:
    reference_set = ReferenceSet(
        directory=".panel/x",
        family="Iridaceae",
        rule="same family first (D-b)",
        built_at="2026-10-08T00:00:00",
        records=("NC_1", "NC_2"),
        leave_one_out=leave_one_out(["NC_1", "NC_2"]),
    )
    assert reference_set.references_for("NC_2") == ("NC_1",)
    assert reference_set.references_for("PX_not_in_the_set") == ("NC_1", "NC_2")


def test_the_committed_manifest_is_a_leave_one_out_over_its_own_records() -> None:
    """The design, readable without the pin: 34 targets, each annotated by the other 33."""
    reference_set = read_manifest(MANIFEST)
    assert reference_set.family == "Iridaceae"
    assert len(reference_set.records) == 34
    assert len(set(reference_set.records)) == 34
    assert set(reference_set.leave_one_out) == set(reference_set.records)
    for target, references in reference_set.leave_one_out.items():
        assert target not in references, f"{target} would annotate itself"
        assert set(references) == set(reference_set.records) - {target}
    # and the design is not a hand-written table: it is what the rule produces from the records
    assert reference_set.leave_one_out == leave_one_out(reference_set.records)


def test_the_manifest_round_trips_through_its_file(tmp_path: Path) -> None:
    built = ReferenceSet(
        directory=".panel/x",
        family="Iridaceae",
        rule="same family first (D-b)",
        built_at="2026-10-08T00:00:00",
        records=("NC_1", "NC_2"),
        leave_one_out=leave_one_out(["NC_1", "NC_2"]),
    )
    written = write_manifest(built, tmp_path / "REFERENCE_SET.json")
    assert read_manifest(written) == built
    payload = json.loads(written.read_text(encoding="utf-8"))
    assert payload["leave_one_out"]["NC_1"] == ["NC_2"]  # a list in the file, tuples in memory


def _pinned_pair(tmp_path: Path) -> Path:
    """Two real GenBank records pinned into a throwaway set — the writer's own path, not a stub."""
    panel = tmp_path / "reference"
    panel.mkdir()
    sources = [
        CORPUS / "Iris_domestica_NC_050833.1.pga.gb",
        CORPUS / "Iris_hermona_NC_081931.1.published.gb",
    ]
    pin_files(sources, panel, note="test")
    return panel


def test_verify_pin_is_silent_on_a_pin_that_still_matches(tmp_path: Path) -> None:
    assert verify_pin(_pinned_pair(tmp_path)) == ()


def test_verify_pin_names_a_file_whose_bytes_changed(tmp_path: Path) -> None:
    """A changed byte is not a warning: the set is no longer the set that was pinned."""
    panel = _pinned_pair(tmp_path)
    victim = next(panel.glob("*.gb"))
    victim.write_bytes(victim.read_bytes() + b"\n")
    problems = verify_pin(panel)
    assert len(problems) == 1
    assert problems[0].startswith(f"{victim.name} — sha256 ")
    assert "the pin records" in problems[0]


def test_verify_pin_names_a_record_the_directory_lost(tmp_path: Path) -> None:
    panel = _pinned_pair(tmp_path)
    victim = next(panel.glob("*.gb"))
    victim.unlink()
    assert verify_pin(panel) == (f"{victim.name} — in PINNED.json but not in the directory",)


def test_verify_pin_works_from_another_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A pin made with a relative ``--out`` records a relative path; verify still finds it."""
    monkeypatch.chdir(tmp_path)
    pin_files(
        [CORPUS / "Iris_hermona_NC_081931.1.published.gb"],
        Path(".panel") / "reference",
        note="test",
    )
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    assert verify_pin(tmp_path / ".panel" / "reference") == ()


def test_verify_pin_names_a_record_the_pin_never_listed(tmp_path: Path) -> None:
    """An unlisted record is one nobody checked, and a run must not search it silently."""
    panel = _pinned_pair(tmp_path)
    (panel / "NC_999999.1.gb").write_text("LOCUS       NC_999999\n", encoding="utf-8")
    assert verify_pin(panel) == ("NC_999999.1.gb — in the directory but not in PINNED.json",)


def test_the_cli_reports_a_broken_pin_and_exits_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    panel = _pinned_pair(tmp_path)
    assert main(["reference", "verify", "--panel", str(panel)]) == 0
    again = main(["reference", "verify", "--panel", str(tmp_path / "nothing-here")])
    assert again == 1
    assert "is this a pinned reference set?" in capsys.readouterr().err


@pytest.mark.panel
@pytest.mark.skipif(
    not PANEL.is_dir(),
    reason="the pinned plastid reference set is not in this working tree (.panel/ is gitignored)",
)
def test_the_committed_manifest_is_what_the_pinned_set_builds(tmp_path: Path) -> None:
    """Every sha256 of the pin, and the committed manifest re-derived from it.

    This is P3B-01's Done-when: the set is what was pinned (34 records, no drift), and its
    leave-one-out manifest is not a transcription of the pin — it is what the rule produces.
    """
    assert verify_pin(PANEL) == ()
    built = build_manifest(PANEL, family="Iridaceae", rule=read_manifest(MANIFEST).rule)
    committed = read_manifest(MANIFEST)
    assert built.records == committed.records
    assert built.leave_one_out == committed.leave_one_out

    # the same conclusion through the CLI, which is the door a person actually uses
    assert main(["reference", "verify", "--panel", str(PANEL)]) == 0
    copied = tmp_path / "REFERENCE_SET.json"
    assert (
        main(
            [
                "reference",
                "loo",
                "--panel",
                str(PANEL),
                "--family",
                "Iridaceae",
                "--rule",
                committed.rule,
                "--out",
                str(copied),
            ]
        )
        == 0
    )
    assert read_manifest(copied).leave_one_out == committed.leave_one_out


@pytest.mark.panel
@pytest.mark.skipif(
    not PANEL.is_dir() or not (ROOT / ".panel" / "iridaceae-34").is_dir(),
    reason="the corpus panel is not in this working tree (.panel/ is gitignored)",
)
def test_the_reference_set_holds_the_corpus_records_byte_for_byte() -> None:
    """The reference bytes *are* the corpus bytes — re-pinning must not have changed them.

    Worth a test of its own: the answer keys (the published annotations) belong to the corpus
    files, so a reference set that quietly held a re-downloaded or edited copy would be a
    different set than the one the validation is measured against.
    """
    corpus = ROOT / ".panel" / "iridaceae-34"
    for record in json.loads((PANEL / "PINNED.json").read_text(encoding="utf-8"))["genomes"]:
        other = corpus / Path(record["path"]).name
        assert other.is_file(), f"{other} is missing — the corpus pin is how this set was checked"
        assert hashlib.sha256(other.read_bytes()).hexdigest() == record["sha256"], (
            f"{other.name} differs from the record pinned as the reference set"
        )
