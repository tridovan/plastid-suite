"""MISA adapter: parsing and the loud missing-tool path (no external tools needed)."""

from __future__ import annotations

from pathlib import Path

import pytest

from plastid.core.model import Genome
from plastid.external.misa import (
    MisaResult,
    SsrFeature,
    _parse_gff,
    _parse_statistics,
    misa_ini_text,
    misa_script,
    misa_version,
    parse_unit_repeats,
    run_misa,
)
from plastid.run.params import ParamValueError
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
    """A ``repeat_region`` line is one cluster, and each of its components becomes a feature.

    The *other* encoding MISA uses — a ``compound_repeat`` note on the component's own line — is
    the subject of ``test_a_compound_component_is_flagged_on_its_own_line_too``; this test is
    about the cluster line, so it selects by position.
    """
    cluster = [
        feature for feature in _parse_gff(GFF) if feature.compound and feature.start == 48659
    ]
    assert [(feature.unit_label, feature.motif) for feature in cluster] == [
        ("di", "(TA)6"),
        ("mono", "(A)10"),
    ]
    # the GFF does not separate them, so both components carry the cluster's own span
    assert {(feature.start, feature.end) for feature in cluster} == {(48659, 48679)}


def test_a_compound_component_is_flagged_on_its_own_line_too() -> None:
    """MISA marks a cluster's components two ways, and one of them is the component's own line.

    Measured 10/10/2026 on the Crataegus answer-key corpus (P8B-01): MISA 2.1 writes *some*
    clusters as a single ``repeat_region`` line listing their components, and others as one line
    per component carrying a ``compound_repeat`` prefix — and a file can carry both. The parser
    understood only the first, so every row of the second kind read "no": on PX413282 the
    per-SSR table flagged **0** of its **16** compound components while MISA's own statistics
    said 8, and the reconciliation did not catch it because it compares counts, not flags. The
    golden genome hid the same fault — its table flagged 2 of its 6 components.
    """
    labelled = {feature.start: feature for feature in _parse_gff(GFF)}
    assert labelled[17767].compound is True  # compound_repeat,monomeric_repeat,(T)10
    assert labelled[4562].compound is False  # a plain SSR, no cluster
    assert labelled[48659].compound is True  # the repeat_region line's components — unchanged


def test_missing_perl_raises_a_named_tool_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("plastid.external.misa.which", lambda name: None)
    genome = Genome(name="synthetic", sequence="ACGT" * 25)
    with pytest.raises(ToolError, match="perl"):
        run_misa(genome, ctx=None)  # type: ignore[arg-type] — perl is checked before ctx is used


def test_misa_result_shape() -> None:
    result = MisaResult(total=37, by_unit={1: 31}, compound=3, features=(), unit_repeats={1: 10})
    assert result.total == 37


# --- MISA is not bundled (D-q, 11/10/2026) ------------------------------------------------
#
# The script left the package because its source zip states no licence and redistribution is the
# one act a missing licence forbids (docs/licensing.md). What replaces "a file is missing" is a
# message naming the download line, and a release check — MISA has no version flag, so the only
# place it says which release it is, is its own header.

HEADER = """\
#!/usr/bin/perl -w

###_______________________________________________________________________________
###
###Program name: misa.pl
###Release date: 25/08/20 (version 2.1)
###
"""


def test_misa_script_names_the_download_line_when_it_is_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The negative case that matters: no `misa.pl` anywhere says how to get one."""
    monkeypatch.setenv("PATH", str(tmp_path))  # a PATH with nothing on it
    with pytest.raises(ToolError, match="docs/tools/misa.md"):
        misa_script()


def test_misa_script_takes_an_explicit_path(tmp_path: Path) -> None:
    """`--set misa=…` is the route for a machine that keeps its tools off PATH."""
    script = tmp_path / "misa.pl"
    script.write_text(HEADER, encoding="utf-8")
    assert misa_script(script) == script


def test_misa_script_rejects_an_explicit_path_that_is_not_a_file(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="not a file"):
        misa_script(tmp_path / "typo.pl")


def test_misa_version_is_read_from_the_scripts_own_header(tmp_path: Path) -> None:
    script = tmp_path / "misa.pl"
    script.write_text(HEADER, encoding="utf-8")
    assert misa_version(script) == "2.1"


def test_misa_version_refuses_another_release(tmp_path: Path) -> None:
    """A different release is a different caller — the golden numbers were measured with 2.1."""
    script = tmp_path / "misa.pl"
    script.write_text(HEADER.replace("version 2.1", "version 2.2"), encoding="utf-8")
    with pytest.raises(ToolError, match="expected MISA 2.1"):
        misa_version(script)


def test_misa_version_refuses_a_script_that_never_says(tmp_path: Path) -> None:
    script = tmp_path / "misa.pl"
    script.write_text("#!/usr/bin/perl -w\nprint 'hello';\n", encoding="utf-8")
    with pytest.raises(ToolError, match="cannot be read"):
        misa_version(script)


def test_the_reported_thresholds_are_the_ones_in_the_shipped_ini() -> None:
    """The summary's threshold string and the file MISA actually reads must not drift (M11.17)."""
    from plastid.core.resources import data_file

    line = next(
        line
        for line in data_file("misa.ini").read_text(encoding="utf-8").splitlines()
        if line.startswith("definition(unit_size,min_repeats)")
    )
    expected = " ".join(f"({value.replace('-', '/')})" for value in line.split(":", 1)[1].split())
    from plastid.external.misa import threshold_text

    assert threshold_text() == expected
    assert expected == "(1/10) (2/6) (3/5) (4/5) (5/5) (6/5)"


# --- the user's thresholds reach the ini (M11.20) -------------------------------------------
#
# Before this, `--set unit_repeats=…` was recorded in the manifest and MISA still ran the packaged
# ini: the summary read the shipped numbers, so the manifest stated something false. The fix is
# not a looser parser — it is that the ini handed to MISA is built from the value, and the result
# reports the thresholds read back out of that ini.


def test_unit_repeats_are_parsed_in_unit_order() -> None:
    assert parse_unit_repeats(" 2-8  1-12 3-7 ") == {1: 12, 2: 8, 3: 7}


@pytest.mark.parametrize(
    ("text", "reason"),
    [
        ("", "at least one"),
        ("1-12 nope", "not '<unit>-<min_repeats>'"),
        ("7-5", "outside 1–6"),
        ("0-5", "outside 1–6"),
        ("1-0", "at least 1"),
        ("1-12 1-8", "named twice"),
    ],
)
def test_a_malformed_unit_repeats_value_is_refused_by_name(text: str, reason: str) -> None:
    """A bad value is a usage error with the accepted shape, never a silently dropped token."""
    with pytest.raises(ParamValueError, match=reason):
        parse_unit_repeats(text)


def test_the_generated_ini_changes_only_the_definition_line() -> None:
    """The threshold choice must not move MISA's other settings."""
    from plastid.core.resources import data_file

    shipped = data_file("misa.ini").read_text(encoding="utf-8")
    generated = misa_ini_text({1: 12, 2: 8, 3: 7, 4: 6, 5: 5, 6: 5})

    shipped_lines = shipped.splitlines()
    generated_lines = generated.splitlines()
    assert generated_lines[0].startswith("definition(unit_size,min_repeats):")
    assert generated_lines[0].split(":", 1)[1].split() == [
        "1-12",
        "2-8",
        "3-7",
        "4-6",
        "5-5",
        "6-5",
    ]
    assert generated_lines[1:] == shipped_lines[1:]
    assert generated.endswith("\n") is shipped.endswith("\n")


def test_no_choice_returns_the_shipped_ini_unchanged() -> None:
    from plastid.core.resources import data_file

    assert misa_ini_text() == data_file("misa.ini").read_text(encoding="utf-8")
