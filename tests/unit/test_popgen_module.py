"""The popgen module with MAFFT stubbed — groups, blanks, and the notes that explain them.

The estimators themselves are tested in ``test_popgen.py``; this file is about the wiring: which
rows exist, what a group with too few sequences looks like, and that a genome without a group is
refused rather than quietly dropped from every table.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import plastid.modules.popgen as popgen
from plastid.core.align import Alignment
from plastid.core.model import Gene, Genome, Interval
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

# One gene per genome; the two groups differ by one base in matK.
_RBCL = "ATGGCTTCTAAACCCGGGTTT"
_MATK_A = "ATGAAACCCGGGTTTAAACCCAAA"
_MATK_B = "ATGAAACCCGGGTTTAAACCCAAG"


@pytest.fixture
def stubbed(monkeypatch):
    def fake_align(rows, ctx):
        return Alignment(
            names=tuple(name for name, _ in rows), sequences=tuple(seq for _, seq in rows)
        )

    monkeypatch.setattr(popgen, "align", fake_align)


def _genome(accession: str, genes: dict[str, str]) -> Genome:
    chunks, records, offset = [], [], 0
    for name, sequence in genes.items():
        chunks.append(sequence)
        records.append(Gene(name, "CDS", (Interval(offset, offset + len(sequence)),)))
        offset += len(sequence)
    return Genome(
        name=f"synthetic {accession}",
        accession=accession,
        sequence="".join(chunks),
        genes=tuple(records),
        topology="linear",
    )


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "popgen-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="popgen", module_order=1)


def _sample() -> tuple[Genome, ...]:
    return (
        _genome("A1", {"rbcL": _RBCL, "matK": _MATK_A}),
        _genome("A2", {"rbcL": _RBCL, "matK": _MATK_A}),
        _genome("A3", {"rbcL": _RBCL, "matK": _MATK_A}),
        _genome("B1", {"rbcL": _RBCL, "matK": _MATK_B}),
        _genome("B2", {"rbcL": _RBCL, "matK": _MATK_B}),
    )


def test_each_gene_and_group_gets_a_row_with_its_n(tmp_path: Path, stubbed) -> None:
    groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}
    result = popgen.run(_sample(), {"groups": groups}, _ctx(tmp_path))
    rows = {row["row_id"]: row for row in result.tables["popgen_groups"].rows}
    assert set(rows) == {"rbcL#A", "rbcL#B", "matK#A", "matK#B"}
    assert rows["matK#A"]["n"] == 3
    assert rows["matK#B"]["n"] == 2
    # n = 3 and n = 2 cannot carry the neutrality tests: blank, with the reason in the note
    assert rows["matK#A"]["tajimas_d"] is None
    assert rows["matK#A"]["fu_li_d_star"] is None
    assert "4 sequences" in rows["matK#A"]["note"]
    # a constant gene inside a group is a measured zero, not a blank
    assert rows["rbcL#A"]["segregating"] == 0
    assert rows["rbcL#A"]["hd"] == 0
    assert rows["rbcL#A"]["pi"] == 0


def test_fst_is_computed_between_groups_per_gene(tmp_path: Path, stubbed) -> None:
    groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B"}
    result = popgen.run(_sample(), {"groups": groups}, _ctx(tmp_path))
    fst = {row["row_id"]: row for row in result.tables["popgen_fst"].rows}
    assert set(fst) == {"rbcL#A#B", "matK#A#B"}
    # rbcL is the same sequence everywhere: H_T = 0 → nothing to partition → blank, said why
    assert fst["rbcL#A#B"]["fst"] is None
    assert "not computable" in fst["rbcL#A#B"]["note"]
    # matK differs between the groups with no within-group variation → Fst = 1
    assert fst["matK#A#B"]["fst"] == pytest.approx(1.0)
    assert (fst["matK#A#B"]["n_a"], fst["matK#A#B"]["n_b"]) == (3, 2)


def test_a_genome_without_a_group_is_refused(tmp_path: Path, stubbed) -> None:
    groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B"}  # B2 missing
    with pytest.raises(ValueError, match="B2"):
        popgen.run(_sample(), {"groups": groups}, _ctx(tmp_path))


def test_no_groups_at_all_is_refused_with_the_fix_named(tmp_path: Path, stubbed) -> None:
    with pytest.raises(ValueError, match="--groups"):
        popgen.run(_sample(), {}, _ctx(tmp_path))


def test_a_record_with_no_cds_is_named_in_the_notes(tmp_path: Path, stubbed) -> None:
    genomes = (*_sample(), Genome(name="empty", accession="EMPTY.1", sequence="ACGT" * 10))
    groups = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "EMPTY.1": "B"}
    result = popgen.run(genomes, {"groups": groups}, _ctx(tmp_path))
    assert any("EMPTY.1" in note and "no CDS" in note for note in result.notes)
    summary = result.tables["popgen_summary"].rows[0]
    assert summary["n_genomes"] == 6
    assert summary["group_sizes"] == "A: 3, B: 3"
