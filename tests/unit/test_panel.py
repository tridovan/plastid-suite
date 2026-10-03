"""Reference Panel: tiering, proposal modes, and pinning (offline — the fetcher is injected)."""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plastid.core.panel import (
    IndexSnapshot,
    PanelEntry,
    entries_from_summaries,
    entry_from_summary,
    load_index,
    load_pinned,
    pin,
    propose,
    save_index,
    tier_of,
)

LINEAGE = "cellular organisms; Eukaryota; Lamiales; Lamiaceae; Coleus"

TARGET = PanelEntry(
    accession="NC_050833.1",
    organism="Iris domestica",
    species="Iris domestica",
    genus="Iris",
    family="Iridaceae",
    order="Asparagales",
    length=153_724,
    annotated=True,
)


def _entry(
    accession: str,
    species: str,
    genus: str,
    family: str,
    order: str,
    *,
    length: int = 153_000,
    annotated: bool = True,
    refseq: bool = True,
) -> PanelEntry:
    return PanelEntry(
        accession=accession,
        organism=species,
        species=species,
        genus=genus,
        family=family,
        order=order,
        length=length,
        annotated=annotated,
        refseq=refseq,
    )


def _index() -> IndexSnapshot:
    return IndexSnapshot(
        retrieved_at="2026-10-02",
        entries=(
            TARGET,
            _entry("NC_111111.1", "Iris domestica", "Iris", "Iridaceae", "Asparagales"),  # tier 0
            _entry("NC_222222.1", "Iris gatesii", "Iris", "Iridaceae", "Asparagales"),  # tier 1
            _entry("NC_333333.1", "Crocus sativus", "Crocus", "Iridaceae", "Asparagales"),  # tier 2
            _entry(
                "NC_444444.1", "Allium cepa", "Allium", "Amaryllidaceae", "Asparagales"
            ),  # tier 3
            _entry("NC_555555.1", "Oryza sativa", "Oryza", "Poaceae", "Poales"),  # unrelated
            _entry("NC_666666.1", "Iris dwarf", "Iris", "Iridaceae", "Asparagales", length=40_000),
            _entry(
                "NC_777777.1", "Iris plain", "Iris", "Iridaceae", "Asparagales", annotated=False
            ),
        ),
    )


def test_tier_of() -> None:
    entries = {entry.accession: entry for entry in _index().entries}
    assert tier_of(TARGET, entries["NC_111111.1"]) == 0
    assert tier_of(TARGET, entries["NC_222222.1"]) == 1
    assert tier_of(TARGET, entries["NC_333333.1"]) == 2
    assert tier_of(TARGET, entries["NC_444444.1"]) == 3
    assert tier_of(TARGET, entries["NC_555555.1"]) is None
    assert tier_of(TARGET, TARGET) is None  # never proposes the target itself


def test_quality_filter_drops_short_and_unannotated() -> None:
    proposal = propose(_index(), TARGET.accession, mode="A", level="genus")
    assert [candidate.entry.accession for candidate in proposal.candidates] == [
        "NC_111111.1",
        "NC_222222.1",
    ]
    assert any("dropped" in note for note in proposal.notes)


def test_mode_a_takes_one_level() -> None:
    proposal = propose(_index(), TARGET.accession, mode="A", level="family")
    assert proposal.mode == "A"
    assert [candidate.tier for candidate in proposal.candidates] == [0, 1, 2]
    assert any("mode A" in note for note in proposal.notes)


def test_mode_b_takes_representatives_per_tier() -> None:
    proposal = propose(_index(), TARGET.accession, mode="B", per_tier=1)
    assert [candidate.tier for candidate in proposal.candidates] == [0, 1, 2, 3]
    assert [candidate.entry.accession for candidate in proposal.candidates] == [
        "NC_111111.1",
        "NC_222222.1",
        "NC_333333.1",
        "NC_444444.1",
    ]


def test_unknown_mode_and_level_raise() -> None:
    with pytest.raises(ValueError, match="mode"):
        propose(_index(), TARGET.accession, mode="C")
    with pytest.raises(ValueError, match="level"):
        propose(_index(), TARGET.accession, mode="A", level="kingdom")
    with pytest.raises(KeyError, match="refresh"):
        propose(_index(), "NC_999999.9")


def test_index_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "index.json"
    save_index(path, _index())
    snapshot = load_index(path)
    assert snapshot.retrieved_at == "2026-10-02"
    assert len(snapshot.entries) == len(_index().entries)


def test_pin_writes_checksums_and_verifies_the_file(tmp_path: Path, data_dir: Path) -> None:
    calls: list[str] = []
    golden = (data_dir / "NC_050833.1.gb").read_text(encoding="utf-8")

    def fake_fetch(accession: str) -> str:
        calls.append(accession)
        return golden

    entries = tuple(entry for entry in _index().entries if entry.accession.startswith("NC_222"))
    pinned = pin(entries, fake_fetch, tmp_path, timestamp="2026-10-02T18:00:00")
    assert calls == ["NC_222222.1"]
    assert pinned[0].sha256 and pinned[0].md5 and pinned[0].bytes > 0
    assert (tmp_path / "NC_222222.1.gb").is_file()
    # the pin is verified by reading the file back with our own reader
    assert pinned[0].parse_error is None
    assert pinned[0].length == 153_724
    assert pinned[0].cds == 86
    assert (pinned[0].trna, pinned[0].rrna) == (38, 8)
    assert pinned[0].topology == "circular"
    # the IR outcome is measured at pin time — recorded as data, never enforced here
    assert (pinned[0].ir_length, pinned[0].ir_source) == (26_214, "annotation")
    assert pinned[0].ir_identity == pytest.approx(1.0)
    assert pinned[0].ir_error is None

    manifest = json.loads((tmp_path / "PINNED.json").read_text(encoding="utf-8"))
    assert manifest["pinned_at"] == "2026-10-02T18:00:00"
    assert manifest["genomes"][0]["accession"] == "NC_222222.1"

    # The offline path: later runs read PINNED.json and never call fetch again.
    again = load_pinned(tmp_path)
    assert again == pinned
    assert calls == ["NC_222222.1"]


def test_pin_records_a_parse_failure_instead_of_crashing(tmp_path: Path) -> None:
    entries = tuple(entry for entry in _index().entries if entry.accession.startswith("NC_222"))
    pinned = pin(entries, lambda accession: "this is not a GenBank file", tmp_path)
    assert pinned[0].parse_error is not None  # loud in the data, never silent
    assert pinned[0].length is None and pinned[0].cds is None
    assert pinned[0].ir_length is None and pinned[0].ir_error is None


def _genbank_text(record: SeqRecord) -> str:
    buffer = io.StringIO()
    SeqIO.write(record, buffer, "genbank")
    return buffer.getvalue()


def test_pin_records_a_failed_ir_gate_instead_of_raising(tmp_path: Path, data_dir: Path) -> None:
    """A diverged IR is data at pin time; only the run-time gate may stop an analysis."""
    record = next(
        record
        for record in SeqIO.parse(str(data_dir / "iridaceae_edge_cases.gb"), "genbank")
        if record.id == "NC_063920.1"
    )
    entries = (PanelEntry(accession="NC_063920.1", organism="Iris petrana"),)
    pinned = pin(entries, lambda accession: _genbank_text(record), tmp_path)

    assert pinned[0].parse_error is None
    assert pinned[0].ir_identity == pytest.approx(0.9263, abs=0.001)
    assert "IR gate failed" in (pinned[0].ir_error or "")


def test_pin_records_an_absent_ir(tmp_path: Path) -> None:
    """No IR pair at all — a mitogenome, a fragment — is pinned with the reason recorded."""
    record = SeqRecord(
        Seq("ATGC" * 1000),
        id="XX_TEST.1",
        name="XX_TEST.1",
        description="synthetic record with no inverted repeat",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    entries = (PanelEntry(accession="XX_TEST.1", organism="synthetic"),)
    pinned = pin(entries, lambda accession: _genbank_text(record), tmp_path)

    assert pinned[0].parse_error is None
    assert pinned[0].ir_length is None
    assert "IRDetectionError" in (pinned[0].ir_error or "")


def test_entries_from_summaries_records_problems_instead_of_dropping() -> None:
    records = [
        (
            {"accessionversion": "NC_050833.1", "organism": "Iris domestica", "slen": 153724},
            LINEAGE,
        ),
        ({"organism": "record without an accession"}, LINEAGE),
    ]
    entries, problems = entries_from_summaries(records)
    assert [entry.accession for entry in entries] == ["NC_050833.1"]
    assert len(problems) == 1 and "accession" in problems[0]


def test_entry_from_summary_never_assumes_genbank_annotation() -> None:
    """RefSeq is annotated by NCBI's build; a GenBank submission's annotation is unverified."""
    refseq = entry_from_summary({"accessionversion": "NC_050833.1", "organism": "Iris"}, LINEAGE)
    genbank = entry_from_summary({"accessionversion": "PX923121.1", "organism": "Coleus"}, LINEAGE)
    assert refseq.annotated is True
    assert genbank.annotated is False
