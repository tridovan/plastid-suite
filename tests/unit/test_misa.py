"""MISA adapter: parsing and the loud missing-tool path (no external tools needed)."""

from __future__ import annotations

import pytest

from plastid.core.model import Genome
from plastid.external.misa import MisaResult, SsrFeature, _parse_gff, _parse_statistics, run_misa
from plastid.run.tools import ToolError

STATISTICS = """\
RESULTS OF MICROSATELLITE SEARCH
================================

Total number of sequences examined:              1
Total number of identified SSRs:                 37
Number of SSRs present in compound formation:    3

Distribution to different repeat type classes
---------------------------------------------

Unit size\tNumber of SSRs
1\t31
2\t4
3\t2
"""

GFF = """\
##gff-version 3
##sequence-region NC_050833.1 1 153724
NC_050833.1\tMISA\tregion\t1\t153724\t.\t.\t.\tID=NC_050833.1.1
NC_050833.1\tMISA\tmonomeric_repeat\t4562\t4572\t.\t.\t.\tNote=monomeric_repeat,(A)11;ID=NC_050833.1.2
NC_050833.1\tMISA\tdinucleotide_repeat\t48659\t48679\t.\t.\t.\tNote=dinucleotide_repeat,(TA)6;ID=NC_050833.1.3
"""


def test_statistics_parsing() -> None:
    total, by_unit, compound = _parse_statistics(STATISTICS)
    assert total == 37
    assert by_unit == {1: 31, 2: 4, 3: 2}
    assert compound == 3


def test_gff_parsing() -> None:
    features = _parse_gff(GFF)
    assert len(features) == 2  # the region line is not an SSR
    assert features[0] == SsrFeature(unit_size=1, motif="(A)11", start=4562, end=4572)
    assert features[1].unit_label == "di"
    assert features[1].motif == "(TA)6"


def test_missing_perl_raises_a_named_tool_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("plastid.external.misa.which", lambda name: None)
    genome = Genome(name="synthetic", sequence="ACGT" * 25)
    with pytest.raises(ToolError, match="perl"):
        run_misa(genome, ctx=None)  # type: ignore[arg-type] — perl is checked before ctx is used


def test_misa_result_shape() -> None:
    result = MisaResult(total=37, by_unit={1: 31}, compound=3, features=())
    assert result.total == 37
