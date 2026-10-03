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

# Verbatim line shapes from a real MISA 2.1 run (golden genome, 2026-10-03). The earlier fixture
# invented a `dinucleotide_repeat` type MISA never writes — which is exactly how the parser's
# silent loss of every non-mono SSR went unnoticed until P2-06 reconciled the counts.
GFF = """\
##gff-version 3
##sequence-region NC_050833.1 1 153724
NC_050833.1\tMISA\tregion\t1\t153724\t.\t.\t.\tID=NC_050833.1.1
NC_050833.1\tMISA\tmonomeric_repeat\t4562\t4572\t.\t.\t.\tNote=monomeric_repeat,(A)11;ID=NC_050833.1.2
NC_050833.1\tMISA\tmonomeric_repeat\t17767\t17776\t.\t.\t.\tNote=compound_repeat,monomeric_repeat,(T)10;ID=NC_050833.1.9
NC_050833.1\tMISA\tmicrosatellite\t58875\t58886\t.\t.\t.\tNote=dinucleotide_repeat_microsatellite_feature,(AT)6;ID=NC_050833.1.23
NC_050833.1\tMISA\tmicrosatellite\t126493\t126513\t.\t.\t.\tNote=trinucleotide_repeat_microsatellite_feature,(CTT)7;ID=NC_050833.1.36
NC_050833.1\tMISA\trepeat_region\t48659\t48679\t.\t.\t.\tNote=repeat_region,dinucleotide_repeat_microsatellite_feature,(TA)6,monomeric_repeat,(A)10;ID=NC_050833.1.20
"""


def test_statistics_parsing() -> None:
    total, by_unit, compound = _parse_statistics(STATISTICS)
    assert total == 37
    assert by_unit == {1: 31, 2: 4, 3: 2}
    assert compound == 3


def test_gff_parsing_reads_every_vocabulary_misa_writes() -> None:
    """``microsatellite`` carries the class in its Note; a ``region`` line is not an SSR."""
    features = _parse_gff(GFF)
    labels = [feature.unit_label for feature in features]
    assert labels == ["mono", "mono", "di", "tri", "di", "mono"]
    assert features[0] == SsrFeature(unit_size=1, motif="(A)11", start=4562, end=4572)
    assert features[1].motif == "(T)10"  # a compound_repeat-prefixed monomeric line
    assert features[2].motif == "(AT)6"
    assert features[3].unit_label == "tri"


def test_a_compound_cluster_yields_one_feature_per_component() -> None:
    """MISA writes a cluster as one line and its statistics count each component separately."""
    cluster = [feature for feature in _parse_gff(GFF) if feature.compound]
    assert [(feature.unit_label, feature.motif) for feature in cluster] == [
        ("di", "(TA)6"),
        ("mono", "(A)10"),
    ]
    # the GFF does not separate them, so both components carry the cluster's own span
    assert {(feature.start, feature.end) for feature in cluster} == {(48659, 48679)}


def test_missing_perl_raises_a_named_tool_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("plastid.external.misa.which", lambda name: None)
    genome = Genome(name="synthetic", sequence="ACGT" * 25)
    with pytest.raises(ToolError, match="perl"):
        run_misa(genome, ctx=None)  # type: ignore[arg-type] — perl is checked before ctx is used


def test_misa_result_shape() -> None:
    result = MisaResult(total=37, by_unit={1: 31}, compound=3, features=())
    assert result.total == 37
