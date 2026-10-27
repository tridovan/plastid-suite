"""The ``introns`` module's counting rules (M11.15) — CDS and tRNA are different facts.

The golden's numbers are pinned in ``tests/golden/test_tier_a.py``. This file pins the rules a
small record can show without any tool: the kind is published per row, a missing feature class is
"not measured" rather than zero, and tRNA gaps the derivation cannot classify are named in the
notes.
"""

from __future__ import annotations

from pathlib import Path

from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, Interval
from plastid.modules import introns
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

LENGTH = 200
IR = IRResult(
    ira=Interval(50, 100),
    irb=Interval(150, 200),
    identity=1.0,
    mismatches=0,
    compared_length=50,
    exact_core=50,
    mirror_constant=249,
    source="manual",
)


def _gene(name: str, kind: str, parts: tuple[tuple[int, int], ...]) -> Gene:
    return Gene(name=name, kind=kind, parts=tuple(Interval(*part) for part in parts))


def _context(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "introns", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="introns", module_order=1, ir=IR)


def _genome(*genes: Gene) -> Genome:
    return Genome(name="synthetic", sequence="ACGT" * (LENGTH // 4), genes=genes)


def test_cds_and_trna_introns_are_counted_and_published_with_their_kind(tmp_path: Path) -> None:
    genome = _genome(
        _gene("psaA", "CDS", ((1, 5), (8, 12))),  # CDS intron [5,8)
        _gene("trnK", "tRNA", ((14, 18), (22, 26))),  # tRNA intron [18,22)
    )

    result = introns.run(genome, {}, _context(tmp_path))

    summary = result.tables["intron_summary"].rows[0]
    assert (summary["n_cds_introns"], summary["n_trna_introns"], summary["n_introns"]) == (1, 1, 2)
    rows = result.tables["introns"].rows
    assert {(row["kind"], row["gene"]) for row in rows} == {("CDS", "psaA"), ("tRNA", "trnK")}
    assert summary["assessment"] == "CDS and tRNA introns measured"


def test_a_record_without_trna_reports_not_measured_not_zero(tmp_path: Path) -> None:
    genome = _genome(_gene("psaA", "CDS", ((1, 5), (8, 12))))

    result = introns.run(genome, {}, _context(tmp_path))

    summary = result.tables["intron_summary"].rows[0]
    assert summary["n_cds_introns"] == 1
    assert summary["n_trna_introns"] == "not measured"
    assert "tRNA introns not measured" in summary["assessment"]
    assert any("not measured — the record annotates no tRNA" in note for note in result.notes)


def test_a_record_without_cds_reports_not_measured_not_zero(tmp_path: Path) -> None:
    genome = _genome(_gene("trnK", "tRNA", ((14, 18), (22, 26))))

    result = introns.run(genome, {}, _context(tmp_path))

    summary = result.tables["intron_summary"].rows[0]
    assert summary["n_cds_introns"] == "not measured"
    assert summary["n_trna_introns"] == 1


def test_an_annotation_with_no_genes_measures_nothing_and_says_so(tmp_path: Path) -> None:
    """Do not infer "no introns" from an annotation that carries no CDS or tRNA at all."""
    genome = _genome()

    result = introns.run(genome, {}, _context(tmp_path))

    summary = result.tables["intron_summary"].rows[0]
    assert summary["n_cds_introns"] == "not measured"
    assert summary["n_trna_introns"] == "not measured"
    assert summary["n_introns"] == 0
    assert summary["assessment"].startswith("nothing measured")


def test_an_unassessed_trna_gap_is_named_in_the_notes(tmp_path: Path) -> None:
    """A gap across a region boundary is not an intron; it is reported, never dropped."""
    genome = _genome(_gene("trnZ", "tRNA", ((45, 49), (105, 109))))  # LSC → SSC

    result = introns.run(genome, {}, _context(tmp_path))

    assert result.tables["intron_summary"].rows[0]["n_trna_introns"] == 0
    assert any("cross-region gap" in note for note in result.notes)
    assert any("not assessed" in note for note in result.notes)
