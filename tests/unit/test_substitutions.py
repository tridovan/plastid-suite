"""The substitutions module's CDS scope with MAFFT stubbed — the panel's real shape.

Measured 04/10/2026: 3 of the 11 pinned stand-in plastomes carry **no annotated CDS**, and the
first version of this scope returned an empty table for the whole run. A record with no
annotation cannot share a CDS with anything, so it is excluded **and named** — the run continues
for the records that can be compared.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import plastid.modules.substitutions as substitutions
from plastid.core.align import Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

_GENE = "ATGGCTTCTAAACCCGGGTTTTAA"


@pytest.fixture
def stubbed(monkeypatch):
    def fake_align(rows, ctx):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    monkeypatch.setattr(substitutions, "align", fake_align)


def _genome(accession: str, genes: dict[str, str]) -> Genome:
    chunks, records, offset = [], [], 0
    for name, sequence in genes.items():
        chunks.append(sequence)
        records.append(
            Gene(name, "CDS", (Interval(offset, offset + len(sequence)),), transl_table=11)
        )
        offset += len(sequence)
    return Genome(
        name=f"synthetic {accession}",
        accession=accession,
        sequence="".join(chunks),
        genes=tuple(records),
        topology="linear",
    )


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "substitutions-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="substitutions", module_order=1)


def test_a_record_without_annotation_is_named_not_silently_dropped(tmp_path: Path, stubbed) -> None:
    reference = _genome("REF.1", {"rbcL": _GENE})
    query = _genome("QRY.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTA")})
    empty = Genome(name="unannotated", accession="EMPTY.1", sequence="ACGT" * 20)
    result = substitutions.run((reference, query, empty), {"scope": "cds"}, _ctx(tmp_path))
    summary = result.tables["substitutions_summary"].rows[0]
    assert summary["n_genomes"] == 2  # the empty record is not compared …
    assert any("EMPTY.1" in note for note in result.notes)  # … and it is named


def test_the_cds_scope_refuses_to_report_nothing_when_no_genome_is_annotated(
    tmp_path: Path, stubbed
) -> None:
    reference = _genome("REF.1", {"rbcL": _GENE})
    empty = Genome(name="unannotated", accession="EMPTY.1", sequence="ACGT" * 20)
    with pytest.raises(ValueError, match="no genome other than the reference"):
        substitutions.run((reference, empty), {"scope": "cds"}, _ctx(tmp_path))


def test_the_cds_scope_reports_each_query_genome_separately(tmp_path: Path, stubbed) -> None:
    """A manuscript table is per pair: Ts/Tv between the reference and *each* congener, not one
    number pooled over every query. Measured 04/10/2026 — the pooled number alone cannot fill
    that table, which is what P4-05 needs."""
    reference = _genome("REF.1", {"rbcL": _GENE})
    close = _genome("CLOSE.1", {"rbcL": _GENE.replace("GGGTTT", "GGGTTC")})  # one transition
    far = _genome("FAR.1", {"rbcL": _GENE.replace("GGGTTT", "GGGAAA")})  # three T→A
    result = substitutions.run((reference, close, far), {"scope": "cds"}, _ctx(tmp_path))
    by_genome = {row["genome"]: row for row in result.tables["substitutions_by_genome"].rows}
    assert set(by_genome) == {"CLOSE.1", "FAR.1"}
    assert (by_genome["CLOSE.1"]["transitions"], by_genome["CLOSE.1"]["transversions"]) == (1, 0)
    assert (by_genome["FAR.1"]["transitions"], by_genome["FAR.1"]["transversions"]) == (0, 3)
    assert by_genome["CLOSE.1"]["genes"] == 1
    # and the pooled summary is still the sum of the pairs
    summary = result.tables["substitutions_summary"].rows[0]
    assert summary["transitions"] == 1
    assert summary["transversions"] == 3
