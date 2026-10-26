"""Where PTS already draws the figure, the R version must show the same numbers.

A disagreement is a defect in one of them, and the test says which table each read.
Both figures read the *same* table — that is the whole point of the R pass being a second source —
so any difference is one of the two paths reading a cell wrong, and the assertion names the table
under test rather than leaving a picture to be eyeballed.

**π is the exception, and it is stated rather than faked.** PTS's ``pi_window`` draws π along the
alignment in sliding windows (``diversity.diversity_windows``); the manuscript's Figure 6 draws π
per **feature** (``diversity_genes``/``_introns``/``_igs``). They are two different pictures of two
different tables, so there is no number they must agree on — and the last test here pins that
disjointness so the reason cannot quietly change.
"""

from __future__ import annotations

import csv
import io

import pytest
from tests.unit.test_figures_r_tables import _codon_usage_result, _ssr_result

from plastid.result.figures import r_rscu, r_ssr, rscu_heatmap, ssr_distribution
from plastid.result.figures.r_common import AMINO_ACID_NAMES


def _row(text: str) -> dict[str, str]:
    header, row = list(csv.reader(io.StringIO(text)))[:2]
    return dict(zip(header, row, strict=True))


def test_the_rscu_figure_and_the_matplotlib_one_show_the_same_codons_and_values() -> None:
    """Both read ``codon_usage.rscu``: same codons, same order, same RSCU — to the cell's own
    precision.
    """
    results = {"codon_usage": _codon_usage_result()}
    drawn = rscu_heatmap.series(results)
    written = r_rscu.frame(results)
    assert [entry[1] for entry in drawn.codons] == [entry[1] for entry in written.codons], (
        "the two figures list different codons — one of them read codon_usage.rscu wrongly"
    )
    for (amino_acid, codon, _count, rscu), (name, written_codon, text) in zip(
        drawn.codons, written.codons, strict=True
    ):
        assert written_codon == codon
        # the drawn figure keeps the raw float; the written file the column's 4-decimal text
        assert (text == "") is (rscu is None), f"{codon}: blank in one figure, a value in the other"
        if text:
            assert float(text) == pytest.approx(rscu, abs=5e-5), f"{codon}: RSCU differs"
        assert AMINO_ACID_NAMES[amino_acid] == name, f"{codon}: family label differs"


def test_the_ssr_figure_and_the_matplotlib_one_show_the_same_counts() -> None:
    """Both read the ``ssr`` tables — regions, units and the total must match cell for cell."""
    results = {"ssr": _ssr_result()}
    drawn = ssr_distribution.series(results)
    written = _row(r_ssr.frame(results).inputs()[0][1])
    by_region = {region.region: region for region in drawn.regions}

    assert written["Total_SSRs"] == drawn.total
    assert written["Species"] == drawn.accession
    for region in ("LSC", "SSC"):
        assert written[region] == str(by_region[region].total), f"{region} total differs"
    # the paper draws one IR bar; PTS reports the two copies apart — the sum is the agreement
    assert written["IR"] == str(by_region["IRa"].total + by_region["IRb"].total)
    for label, unit in (("Mono", "mono"), ("Di", "di"), ("Tri", "tri")):
        drawn_total = sum(
            count
            for region in drawn.regions
            for name, count, _text in region.segments
            if name == unit and count is not None
        )
        assert written[label] == str(drawn_total), f"{label} ({unit}) differs"


def test_the_two_diversity_figures_read_disjoint_tables_and_that_is_why_they_cannot_agree() -> None:
    """PTS draws π per window; the paper draws π per feature. Stated, not papered over."""
    from plastid.result.figures import pi_window, r_diversity

    window_tables = {table for _module, table in pi_window.SPEC.requires}
    feature_tables = {table for _module, table in r_diversity.SPEC.requires}
    assert window_tables == {"diversity_windows", "diversity_summary"}
    # What one draws and the other does not — this is the whole reason they cannot be compared.
    assert window_tables - feature_tables == {"diversity_windows"}
    assert feature_tables - window_tables == {
        "diversity_genes",
        "diversity_introns",
        "diversity_igs",
    }
