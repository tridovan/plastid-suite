"""The PAML/yn00 adapter — the ctl it writes, the parser, and the refusals (P4-02).

The output fixture is verbatim from the measured run of 2026-10-04 (the golden pair's *rbcL*,
MAFFT-aligned, terminal stop trimmed). The tool is otherwise stubbed the way the BLAST adapter's
tests stub theirs: what is under test is this module's contract with yn00, not PAML's maths — the
real numbers are asserted in ``tests/golden/test_tier_b.py``.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from plastid.external.paml import (
    PAML_VERSION,
    parse_yn00_output,
    read_matrix,
    run_yn00,
    yn00_version,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolError

# The Yang & Nielsen section of the measured rbcL run — the row is what the parser reads.
YN_OUTPUT = """YN00          input.dat

ns =   2\tls = 486

Codon usage for each species
----------------------------------------------------------------------

(A) Nei-Gojobori (1986) method

Nei & Gojobori 1986. dN/dS (dN, dS)

dom
gat                  0.0558 (0.0018 0.0323)

(B) Yang & Nielsen (2000) method

Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous \
substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43

(equal weighting of pathways)

seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE

   2    1   378.2  1079.8   0.0274 16.9092  0.0620 0.0019 +- 0.0013  0.0299 +- 0.0091


(C) LWL85, LPB93 & LWLm methods

2 (gat) vs. 1 (dom)
"""

DN_MATRIX = """     2
dom
gat                    0.0019
"""

DS_MATRIX = """     2
dom
gat                    0.0299
"""

IDENTICAL_OUTPUT = YN_OUTPUT.replace(
    "   2    1   378.2  1079.8   0.0274 16.9092  0.0620 0.0019 +- 0.0013  0.0299 +- 0.0091",
    "   2    1   378.2  1079.8   0.0000  0.0000  -nan  -0.0000 +- 0.0000  0.0000 +- 0.0000",
)


def _ctx(tmp_path: Path) -> RunContext:
    layout = RunLayout.create(tmp_path / "runs", "paml-test", fingerprint="test")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id="kaks", module_order=1)


def _stub(directory: Path, *, output: str = YN_OUTPUT, dn: str = DN_MATRIX, ds: str = DS_MATRIX):
    """A stand-in yn00: it answers the banner probe and writes the canned result files."""
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / "yn00"
    script.write_text(
        "#!/bin/sh\n"
        'if [ ! -f yn00.ctl ]; then echo "YN00 in paml version ' + PAML_VERSION + ', 29 Jan 2026";'
        ' echo "error: ctl file open error." >&2; exit 1; fi\n'
        f"cat > yn00.out <<'OUT'\n{output}OUT\n"
        f"cat > 2YN.dN <<'DN'\n{dn}DN\n"
        f"cat > 2YN.dS <<'DS'\n{ds}DS\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return directory


def test_the_parser_reads_the_yang_and_nielsen_row() -> None:
    result = parse_yn00_output(YN_OUTPUT)
    assert result.dN == 0.0019
    assert result.dS == 0.0299
    assert result.omega == 0.0620
    assert result.kappa == 16.9092
    assert result.t == 0.0274
    assert result.sites_synonymous == 378.2
    assert result.sites_nonsynonymous == 1079.8
    assert result.codons == 486


def test_nan_and_negative_zero_are_read_as_what_they_mean() -> None:
    """yn00 prints ``-nan`` for dN/dS when dS is zero, and ``-0.0000`` for an identical pair.

    A nan that survived into a table would be a number nobody can read; a minus sign on a zero is
    a "less than zero" in a figure.
    """
    result = parse_yn00_output(IDENTICAL_OUTPUT)
    assert result.omega is None
    assert result.dN == 0.0 and result.dS == 0.0
    assert str(result.dN) == "0.0"  # not -0.0


def test_a_missing_yang_and_nielsen_section_is_refused() -> None:
    with pytest.raises(ToolError, match="Yang & Nielsen"):
        parse_yn00_output("YN00  input.dat\nns = 2 ls = 30\n(no sections here)\n")


def test_read_matrix_checks_the_names_it_was_given(tmp_path: Path) -> None:
    path = tmp_path / "2YN.dN"
    path.write_text(DN_MATRIX, encoding="utf-8")
    assert read_matrix(path, "dom", "gat") == 0.0019
    with pytest.raises(ToolError):
        read_matrix(path, "other", "gat")
    with pytest.raises(ToolError):
        read_matrix(path, "dom", "other")
    with pytest.raises(ToolError):
        read_matrix(tmp_path / "absent.dN", "dom", "gat")


def test_the_version_probe_reads_the_no_ctl_banner(tmp_path: Path) -> None:
    stub = _stub(tmp_path / "tools")
    assert yn00_version(str(stub / "yn00")) == PAML_VERSION


def test_a_different_paml_is_refused_by_name(tmp_path: Path) -> None:
    directory = tmp_path / "tools"
    directory.mkdir()
    script = directory / "yn00"
    script.write_text(
        "#!/bin/sh\necho 'YN00 in paml version 4.9.6, Jan 2020'\nexit 1\n", encoding="utf-8"
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    with pytest.raises(ToolError, match="4.9.6"):
        yn00_version(str(script))


def test_run_yn00_writes_the_ctl_and_cross_checks_the_matrices(tmp_path: Path) -> None:
    stub = _stub(tmp_path / "tools")
    ctx = _ctx(tmp_path)
    result = run_yn00(
        "dom",
        "ATGGCTTCTAAACCCGGGTTT",
        "gat",
        "ATGGCTTCTAAACCCGGGTTT",
        ctx,
        directory=stub,
    )
    assert (result.dN, result.dS) == (0.0019, 0.0299)
    assert [record.name for record in ctx.tools] == ["yn00"]
    assert ctx.tools[0].version == PAML_VERSION


def test_run_yn00_refuses_when_the_tool_disagrees_with_itself(tmp_path: Path) -> None:
    """The matrices and the printed table are two code paths inside PAML; this is the check."""
    stub = _stub(tmp_path / "tools", dn=DN_MATRIX.replace("0.0019", "0.0910"))
    with pytest.raises(ToolError, match="disagree"):
        run_yn00(
            "dom",
            "ATGGCTTCTAAACCCGGGTTT",
            "gat",
            "ATGGCTTCTAAACCCGGGTTT",
            _ctx(tmp_path),
            directory=stub,
        )


def test_run_yn00_refuses_inputs_it_cannot_score(tmp_path: Path) -> None:
    stub = _stub(tmp_path / "tools")
    ctx = _ctx(tmp_path)
    with pytest.raises(ValueError, match="whole number of codons"):
        run_yn00("dom", "ATGGCTTC", "gat", "ATGGCTTC", ctx, directory=stub)
    with pytest.raises(ValueError, match="alignment"):
        run_yn00("dom", "ATGGCT", "gat", "ATGGCTT", ctx, directory=stub)
    with pytest.raises(ValueError, match="whitespace"):
        run_yn00("dom estica", "ATGGCT", "gat", "ATGGCT", ctx, directory=stub)
    with pytest.raises(ValueError, match="itself"):
        run_yn00("dom", "ATGGCT", "dom", "ATGGCT", ctx, directory=stub)


def test_an_empty_result_file_is_refused_rather_than_read(tmp_path: Path) -> None:
    stub = _stub(tmp_path / "tools", output="")
    with pytest.raises(ToolError, match="empty"):
        run_yn00(
            "dom",
            "ATGGCTTCTAAACCCGGGTTT",
            "gat",
            "ATGGCTTCTAAACCCGGGTTT",
            _ctx(tmp_path),
            directory=stub,
        )
