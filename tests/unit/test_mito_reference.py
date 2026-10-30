"""The mito reference set (P3A-02) — what it can find, and what it refuses to include.

The committed FASTA is an artifact, so the tests hold it against the vocabulary it must cover and
against the panel it was built from. The exclusions get synthetic records, because the reasons
they are excluded (an anonymous ORF, a product sentence, an rRNA of the wrong size) are exactly
the ones that would otherwise pass unnoticed.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from plastid.core.mito_reference import (
    RRNA_LENGTH_BANDS,
    build_proteins,
    build_rrna,
    read_reference,
    render_fasta,
)

DATA = Path(__file__).resolve().parents[2] / "src" / "plastid" / "data"
PANEL_MITO_REFERENCE = Path(__file__).resolve().parents[2] / ".panel" / "mito-reference"
NEEDS_PANEL = pytest.mark.skipif(
    not PANEL_MITO_REFERENCE.is_dir(),
    reason="the reference records are not on this machine — re-pin them (data/mito_reference.md)",
)


def _vocabulary() -> list[str]:
    payload = yaml.safe_load((DATA / "mito_gene_names.yaml").read_text(encoding="utf-8"))
    return [entry["canonical"] for entry in payload["loci"]]


def _record(path: Path, accession: str, features: list[SeqFeature]) -> Path:
    record = SeqRecord(Seq("ATG" * 300), id=accession, name=accession)
    record.annotations = {
        "molecule_type": "DNA",
        "topology": "circular",
        "organism": "Testus testus",
    }
    record.features = [
        SeqFeature(
            SimpleLocation(0, 8, strand=1),
            type="source",
            qualifiers={"organism": ["Testus testus"]},
        ),
        *features,
    ]
    SeqIO.write(record, str(path), "genbank")
    return path


def _cds(gene: str, translation: str | None, start: int = 0, end: int = 9) -> SeqFeature:
    qualifiers: dict[str, list[str]] = {"gene": [gene]}
    if translation is not None:
        qualifiers["translation"] = [translation]
    return SeqFeature(SimpleLocation(start, end, strand=1), type="CDS", qualifiers=qualifiers)


def test_the_committed_protein_set_covers_the_whole_vocabulary() -> None:
    """Done-when: every one of the 31 loci has at least one reference protein behind it."""
    entries = read_reference(DATA / "mito_reference_proteins.faa")
    covered = {entry.vocabulary for entry in entries if entry.vocabulary}
    assert covered == set(_vocabulary()), "a vocabulary locus no reference can find is unfindable"


def test_the_committed_rrna_set_holds_the_three_rrna_genes() -> None:
    entries = read_reference(DATA / "mito_reference_rrna.fa")
    assert {entry.vocabulary for entry in entries} == set(RRNA_LENGTH_BANDS)
    for entry in entries:
        low, high = RRNA_LENGTH_BANDS[entry.vocabulary or ""]
        assert low <= len(entry.sequence) <= high


def test_a_header_carries_the_record_a_hit_came_from() -> None:
    """ "Which record, which gene, can we name it, whose genome" — a hit has to be traceable."""
    entries = read_reference(DATA / "mito_reference_proteins.faa")
    entry = entries[0]
    assert entry.accession.startswith("NC_")
    assert entry.organism
    assert entry.sequence.isupper()
    # a gene the vocabulary cannot name is kept and marked, never dropped (measured: rpl2, rps12,
    # rps19 … are in the references and not in data/mito_gene_names.yaml)
    assert any(other.vocabulary is None for other in entries)


def test_an_anonymous_open_reading_frame_is_excluded_and_reported(tmp_path: Path) -> None:
    """Measured: Nicotiana's record annotates 121 of these as CDS. Keeping them would let an
    unnamed ORF outscore the real gene it sits next to."""
    path = _record(
        tmp_path / "ref.gb", "NC_999991.1", [_cds("orf101a", "MKK"), _cds("atp1", "MKKV")]
    )
    build = build_proteins([path], ["atp1"])
    assert [entry.source_gene for entry in build.entries] == ["atp1"]
    assert build.excluded == (("NC_999991.1", "orf101a", "anonymous open reading frame"),)


def test_a_feature_named_by_a_sentence_is_excluded(tmp_path: Path) -> None:
    path = _record(tmp_path / "ref.gb", "NC_999991.1", [_cds("hypothetical protein", "MKK")])
    build = build_proteins([path], ["atp1"])
    assert build.entries == ()
    assert build.excluded == (
        ("NC_999991.1", "hypothetical protein", "named by a product sentence, not a symbol"),
    )


def test_a_cds_without_a_published_translation_is_excluded(tmp_path: Path) -> None:
    """Self-translating is not an option: measured, plant mito CDS start with codons no table
    accepts as an initiator, so a record that publishes no protein is a record we cannot use."""
    path = _record(tmp_path / "ref.gb", "NC_999991.1", [_cds("atp1", None)])
    build = build_proteins([path], ["atp1"])
    assert build.entries == ()
    assert build.excluded[0][2] == "the record publishes no /translation"


def test_an_rrna_outside_its_length_band_is_excluded(tmp_path: Path) -> None:
    """Measured: Psychotria's record calls a 722 bp feature "5S" and an 1112 bp one "18S". A
    reference that believed those would search with the wrong sequence."""
    path = _record(
        tmp_path / "ref.gb",
        "NC_999991.1",
        [
            SeqFeature(
                SimpleLocation(0, 722, strand=1), type="rRNA", qualifiers={"gene": ["rrn5"]}
            ),
            SeqFeature(
                SimpleLocation(0, 119, strand=1), type="rRNA", qualifiers={"gene": ["rrn5"]}
            ),
        ],
    )
    build = build_rrna([path])
    assert [len(entry.sequence) for entry in build.entries] == [119]
    assert build.excluded == (("NC_999991.1", "rrn5(722 bp)", "length outside the rRNA band"),)


def test_a_reference_alias_maps_the_records_own_spelling(tmp_path: Path) -> None:
    """``ccmFC``/``ccmFN1`` (Arabidopsis) and ``mat-R`` (Nicotiana) are the same loci as ours."""
    path = _record(
        tmp_path / "ref.gb",
        "NC_999991.1",
        [_cds("ccmFC", "MKK"), _cds("mat-R", "MKK"), _cds("rps12", "MKK")],
    )
    build = build_proteins([path], ["ccmFc", "matR"])
    assert [entry.vocabulary for entry in build.entries] == ["ccmFc", "matR", None]
    assert [entry.source_gene for entry in build.entries] == ["ccmFC", "mat-R", "rps12"]


def test_render_then_read_gives_the_same_entries(tmp_path: Path) -> None:
    """The header is the only record of where a hit came from, so it has to survive the trip."""
    entries = read_reference(DATA / "mito_reference_proteins.faa")[:5]
    written = tmp_path / "reference.faa"
    written.write_text(render_fasta(entries), encoding="utf-8")
    assert read_reference(written) == entries


@pytest.mark.panel
@NEEDS_PANEL
def test_the_committed_set_is_what_the_pinned_panel_builds(tmp_path: Path) -> None:
    """The artifact is reproducible: rebuild it from the pinned records and it is the same file.

    This is the check that makes the FASTA evidence rather than a leftover — the panel's own
    sha256s are in ``PINNED.json``, so the whole chain is pinned.
    """
    paths = sorted(PANEL_MITO_REFERENCE.glob("*.gb"))
    assert len(paths) == 8
    build = build_proteins(paths, _vocabulary())
    assert render_fasta(build.entries) == (DATA / "mito_reference_proteins.faa").read_text(
        encoding="utf-8"
    )
    assert build.vocabulary_coverage == set(_vocabulary())
