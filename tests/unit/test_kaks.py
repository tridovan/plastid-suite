"""The Ka/Ks module with the two external boundaries stubbed (P4-02).

MAFFT and yn00 are replaced at the module's own imports, so what is under test is this module's
decisions: which genes get a row, what a row says when a gene cannot be scored, and how the
pooled row is weighted. The real numbers — a real alignment, a real yn00 run — are asserted in
``tests/golden/test_tier_b.py``.

The cases below are the ones the project's rules single out: identical CDS give Ka = Ks = 0 with
the ratio **undefined, never zero or one**; a gene with no value keeps a row and a named reason;
a gene one genome does not carry is listed rather than dropped.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import plastid.modules.kaks as kaks
from plastid.core.align import Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.external.paml import Yn00Result
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError


def _genome(
    accession: str, genes: dict[str, str], *, extra: dict[str, str] | None = None
) -> Genome:
    """A genome whose sequence is its genes laid end to end, with the coordinates to match."""
    parts = list(genes.items()) + list((extra or {}).items())
    chunks, records, offset = [], [], 0
    for name, sequence in parts:
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


_CDS = {
    "rbcL": "ATGGCTTCTAAACCCGGGTTTTAA",  # 24 bp
    "matK": "ATGAAACCCGGGTTTAAACCCAAATAG",
}


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "kaks-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="kaks", module_order=1)


@pytest.fixture
def stubbed(monkeypatch):
    """Replace MAFFT, the two PAML entry points, and the version probe."""
    calls: list[tuple[str, str]] = []

    def fake_align(rows, ctx):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    def fake_run_yn00(name_a, sequence_a, name_b, sequence_b, ctx, **kwargs):
        calls.append((name_a, sequence_b))
        differences = sum(1 for a, b in zip(sequence_a, sequence_b, strict=True) if a != b)
        if differences == 0:
            return Yn00Result(
                dN=0.0,
                dS=0.0,
                omega=None,
                kappa=None,
                t=0.0,
                sites_synonymous=100.0,
                sites_nonsynonymous=300.0,
                codons=len(sequence_a) // 3,
                raw="",
            )
        return Yn00Result(
            dN=0.01,
            dS=0.05,
            omega=0.2,
            kappa=2.0,
            t=0.06,
            sites_synonymous=100.0,
            sites_nonsynonymous=300.0,
            codons=len(sequence_a) // 3,
            raw="",
        )

    monkeypatch.setattr(kaks, "align", fake_align)
    monkeypatch.setattr(kaks, "run_yn00", fake_run_yn00)
    monkeypatch.setattr(kaks, "resolve_binary", lambda **kwargs: "/fake/yn00")
    monkeypatch.setattr(kaks, "yn00_version", lambda binary: "4.10.10")
    return calls


def test_identical_cds_report_zero_and_an_undefined_ratio(tmp_path: Path, stubbed) -> None:
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"]})
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    row = result.tables["kaks_genes"].rows[0]
    assert (row["ka"], row["ks"]) == (0.0, 0.0)
    assert row["ka_ks"] is None  # undefined, never 0 and never 1
    assert row["status"] == kaks.STATUS_KS_ZERO
    assert "undefined" in row["detail"]


def test_the_ratio_is_divided_only_when_ks_is_not_zero(tmp_path: Path, stubbed) -> None:
    # same 24 bases, one substitution, terminal stop kept so both sides trim to 21
    query_sequence = "ATGGCTTCTAAACCCGGGTATTAA"
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": query_sequence})
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    row = result.tables["kaks_genes"].rows[0]
    assert row["ka_ks"] == pytest.approx(0.2)
    assert row["status"] == kaks.STATUS_OK


def test_a_gene_one_genome_does_not_carry_keeps_a_row(tmp_path: Path, stubbed) -> None:
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"], "matK": _CDS["matK"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"], "psaA": _CDS["matK"]})
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    rows = {row["gene"]: row for row in result.tables["kaks_genes"].rows}
    assert set(rows) == {"rbcL", "matK", "psaA"}
    assert rows["matK"]["status"] == kaks.STATUS_ABSENT_GENOME
    assert "QRY.1" in rows["matK"]["detail"]
    assert rows["psaA"]["status"] == kaks.STATUS_ABSENT_REFERENCE
    assert any("psaA" in note for note in result.notes)


def test_a_cds_that_is_not_whole_codons_is_reported_not_scored(tmp_path: Path, stubbed) -> None:
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"][:-1]})  # 23 bases
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    row = result.tables["kaks_genes"].rows[0]
    assert row["status"] == kaks.STATUS_NOT_WHOLE_CODONS
    assert row["ka"] is None


def test_a_frameshifted_alignment_is_reported_not_scored(
    tmp_path: Path, stubbed, monkeypatch
) -> None:
    def shifting_align(rows, ctx):
        sequences = tuple(seq for _, seq in rows)
        return Alignment(
            names=tuple(name for name, _ in rows),
            sequences=(sequences[0], sequences[1][:4] + "-" + sequences[1][5:]),
        )

    monkeypatch.setattr(kaks, "align", shifting_align)
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"]})
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    row = result.tables["kaks_genes"].rows[0]
    assert row["status"] == kaks.STATUS_FRAMESHIFT
    assert "codon" in row["detail"]
    assert not stubbed  # yn00 was never reached — the frame check comes first


def test_a_pair_yn00_refuses_is_a_row_not_a_crash(tmp_path: Path, stubbed, monkeypatch) -> None:
    def refusing(*args, **kwargs):
        raise ToolError("stop codon TAA in seq # 2")

    monkeypatch.setattr(kaks, "run_yn00", refusing)
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"]})
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    row = result.tables["kaks_genes"].rows[0]
    assert row["status"] == kaks.STATUS_TOOL_FAILED
    assert "stop codon" in row["detail"]


def test_a_missing_estimator_stops_the_run(tmp_path: Path, monkeypatch) -> None:
    def missing(**kwargs):
        raise ToolError("yn00 not found — install the pinned PAML 4.10.10")

    monkeypatch.setattr(kaks, "resolve_binary", missing)
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"]})
    query = _genome("QRY.1", {"rbcL": _CDS["rbcL"]})
    with pytest.raises(ToolError, match="not found"):
        kaks.run((reference, query), {}, _ctx(tmp_path))


def test_the_pooled_row_weights_by_sites_not_by_genes(tmp_path: Path, stubbed) -> None:
    """Two scored genes: the stub gives every pair the same Ka/Ks, so the pooled ratio is that
    value — which is the point of weighting by the estimator's own site counts rather than
    averaging ratios (an average over ratios is a different quantity, and not one yn00 reports).
    """
    reference = _genome("REF.1", {"rbcL": _CDS["rbcL"], "matK": _CDS["matK"]})
    query = _genome(
        "QRY.1",
        {"rbcL": "ATGGCTTCTAAACCCGGGTATTAA", "matK": "ATGAAACCCGGGTTTAAACCCAAGTAG"},
    )
    result = kaks.run((reference, query), {}, _ctx(tmp_path))
    pooled = result.tables["kaks_pooled"].rows[0]
    assert pooled["n_computed"] == 2
    assert pooled["n_not_computed"] == 0
    assert pooled["pooled_ka"] == pytest.approx(0.01)
    assert pooled["pooled_ks"] == pytest.approx(0.05)
    assert pooled["pooled_ka_ks"] == pytest.approx(0.2)
