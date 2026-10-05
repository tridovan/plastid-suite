"""The `stats` module: a real two-group comparison, from genomes to an FDR table.

The genomes here are built by hand — four short CDS over a repeating sequence — so the numbers
the tests assert are arithmetic, not a fixture nobody can re-derive.
"""

from __future__ import annotations

import importlib.util
import re

import pytest

from plastid.core.model import Gene, Genome, Interval
from plastid.modules import stats as stats_module
from plastid.result.result import Result

requires_stats = pytest.mark.skipif(
    importlib.util.find_spec("scipy") is None,
    reason="scipy is an optional dependency (pip install 'plastid-suite[stats]')",
)

# Four CDS over a fixed repeat: ATG + a codon run + TAA. The codons differ between the two
# groups on purpose, so the comparison has something to find.
_CDS = {
    "rbcL": ("ATGGCTGCTGCTTAA", 15),
    "psbA": ("ATGGCTGCAGCTTAA", 15),
    "matK": ("ATGTTTAAACCTTAA", 15),
}


def _genome(accession: str) -> Genome:
    """One short genome: the three CDS laid out head to tail after a filler sequence."""
    filler = "ACGT" * 3
    sequence = filler
    genes: list[Gene] = []
    cursor = len(filler)
    for name, (coding, _length) in _CDS.items():
        genes.append(Gene(name=name, kind="CDS", parts=(Interval(cursor, cursor + len(coding)),)))
        sequence += coding
        cursor += len(coding)
    return Genome(
        name=f"Testus {accession}",
        sequence=sequence,
        accession=accession,
        genes=tuple(genes),
    )


def _genomes() -> list[Genome]:
    return [_genome(name) for name in ("A1", "A2", "B1", "B2")]


def _groups() -> dict[str, str]:
    return {"A1": "Alpha", "A2": "Alpha", "B1": "Beta", "B2": "Beta"}


@requires_stats
def test_the_summary_counts_what_the_comparison_actually_used() -> None:
    result = stats_module.run(_genomes(), {"groups": _groups()}, None)  # type: ignore[arg-type]
    summary = result.tables["stats_summary"]
    assert summary.value("n_genomes") == 4
    assert summary.value("n_groups") == 2
    assert summary.value("groups") == "Alpha, Beta"
    assert summary.value("n_loci") == 3
    assert summary.value("n_genomes_without_cds") == 0
    assert summary.value("n_observations") == 12  # 3 loci × 4 genomes


@requires_stats
def test_a_locus_is_one_observation_per_genome_even_with_two_ir_copies() -> None:
    genome = _genomes()[0]
    doubled = Genome(
        name=genome.name,
        sequence=genome.sequence,
        accession=genome.accession,
        genes=genome.genes + (Gene(name="rbcL", kind="CDS", parts=(Interval(12, 27),)),),
    )
    result = stats_module.run([doubled, doubled], {"groups": {"A1": "Alpha"}}, None)  # type: ignore[arg-type]
    summary = result.tables["stats_summary"]
    assert summary.value("n_collapsed_cds") == 2  # one per genome, not silently dropped
    assert summary.value("n_observations") == 6


@requires_stats
def test_a_genome_with_no_cds_is_named_rather_than_counted_as_a_contribution() -> None:
    """Measured on the stand-in panel: 3 of its 11 records carry no annotation at all."""
    bare = Genome(name="Bare record", sequence="ACGT" * 10, accession="BARE1")
    result = stats_module.run([*_genomes(), bare], {"groups": {**_groups(), "BARE1": "Beta"}}, None)  # type: ignore[arg-type]
    summary = result.tables["stats_summary"]
    assert summary.value("n_genomes") == 5
    assert summary.value("n_genomes_without_cds") == 1
    assert summary.value("genomes_without_cds") == "BARE1"
    assert any("contribute no CDS" in note for note in result.notes)


@requires_stats
def test_a_group_of_one_is_not_computable_in_the_table_not_missing_from_it() -> None:
    result = stats_module.run(
        # Beta has one genome, Alpha two: Beta cannot be ranked against anything, and that must
        # be a row saying so rather than a p-value nobody could reproduce.
        _genomes()[:3],
        {"groups": {"A1": "Alpha", "A2": "Alpha", "B1": "Beta"}},
        None,
    )  # type: ignore[arg-type]
    tests = result.tables["stats_tests"]
    not_computable = [row for row in tests.rows if row["status"] != "ok"]
    assert not_computable, "a group with one genome must not quietly produce a p-value"
    row = not_computable[0]
    assert row["p"] is None and row["p_bh"] is None
    assert "at least 2 observations per group" in row["note"]
    assert tests.value(f"{row['test_id']}.statistic") is None


@requires_stats
def test_the_fdr_table_carries_every_test_with_both_corrections() -> None:
    result = stats_module.run(_genomes(), {"groups": _groups()}, None)  # type: ignore[arg-type]
    tests = result.tables["stats_tests"]
    assert len(tests.rows) == 3 * len(stats_module.MEASURES)
    for row in tests.rows:
        assert row["p_bonferroni"] is None or row["p_bonferroni"] >= row["p"] - 1e-12
        assert row["p_bh"] is None or row["p_bh"] >= row["p"] - 1e-12
    assert result.tables["stats_correlations"].rows


@requires_stats
def test_every_number_in_a_generated_sentence_is_a_cell_of_the_summary() -> None:
    """D5-05's mechanism: the interpretation may only say what the table already carries."""
    result = stats_module.run(_genomes(), {"groups": _groups()}, None)  # type: ignore[arg-type]
    summary = result.tables["stats_summary"]
    rendered = {summary.format(f"{column.name}") for column in summary.columns}
    assert result.notes, "no sentence was generated"
    for sentence in result.notes:
        for number in re.findall(r"\d+(?:\.\d+)?", sentence):
            assert any(number in value for value in rendered), (
                f"{number!r} in {sentence!r} is not carried by any summary cell"
            )


def test_a_run_without_groups_says_how_to_give_them() -> None:
    with pytest.raises(ValueError, match="--groups"):
        stats_module.run(_genomes(), {}, None)  # type: ignore[arg-type]


def test_an_unknown_measure_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="unknown measure"):
        stats_module.run(_genomes(), {"groups": _groups(), "measures": ["nope"]}, None)  # type: ignore[arg-type]


@requires_stats
def test_the_module_renders_a_report_from_its_own_template() -> None:
    from plastid.report.citations import load_citations
    from plastid.report.engine import find_bare_digits, find_refs, load_template, render

    result: Result = stats_module.run(_genomes(), {"groups": _groups()}, None)  # type: ignore[arg-type]
    template = load_template("stats.en.md")
    assert find_bare_digits(template) == []
    for ref in find_refs(template):
        result.value(ref)
    text = render(template, result, citations=load_citations())
    assert "Alpha, Beta" in text
    assert "Draft — needs human review." in text
