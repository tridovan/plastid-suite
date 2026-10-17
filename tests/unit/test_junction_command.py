"""`plastid ir-junction`: the table, the gate's verdict in it, and the loud paths.

The figure itself needs CPJSdraw, which this toolkit does not ship (no licence —
`docs/tools/cpjsdraw.md`), so everything that can be checked without it is checked here, and the
one test that draws a real figure is ``tier_b`` and skips with a named reason.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from plastid.cli import EXIT_OK, EXIT_TOOL, main
from plastid.core.genbank import read_genbank
from plastid.core.model import Genome, Interval
from plastid.core.sequence import reverse_complement
from plastid.run.junction import junction_table

# A record whose two IR copies are exact reverse complements: IRb at 5,001-10,000 and IRa at
# 12,001-17,000 in a 20,000 bp circle. The filler is not a repeat, so nothing else can be mistaken
# for a copy.
_COPY = ("ACGTTGCA" * 625)[:5_000]
_IR_SEQUENCE = "A" * 5_000 + _COPY + "A" * 2_000 + reverse_complement(_COPY) + "A" * 3_000

ALIGNED = Genome(
    name="Aligned example",
    accession="XX_100001.1",
    sequence=_IR_SEQUENCE,
    repeat_regions=(Interval(5_000, 10_000), Interval(12_000, 17_000)),
)

# The same sequence with the second copy **shifted two bases**, its length unchanged. The mirror
# constant `measure_ir` derives from the two boundaries is then two bases away from the one that
# pairs the copies, so every base is compared against the wrong partner and the identity collapses
# — the shape measured on three of the bundle's eight records (H-25). Note what does *not* cause
# it: lengthening a copy at its start leaves the constant right and the gate green.
MISALIGNED = Genome(
    name="Misaligned example",
    accession="XX_100002.1",
    sequence=_IR_SEQUENCE,
    repeat_regions=(Interval(5_000, 10_000), Interval(11_998, 16_998)),
)


def _write(tmp_path: Path, genome: Genome, filename: str) -> Path:
    """A record on disk — the command reads files, so the test writes one."""
    from plastid.core.genbank import to_genbank

    path = tmp_path / filename
    path.write_text(to_genbank(genome), encoding="utf-8")
    return path


def test_the_table_carries_a_row_per_record_and_english_headers() -> None:
    from plastid.run.junction import junction_rows

    text = junction_table(junction_rows([(Path("a.gb"), ALIGNED), (Path("b.gb"), MISALIGNED)]))
    header, *rows = text.strip().splitlines()
    assert header.split(",")[:4] == ["species", "accession", "genome_bp", "LSC"]
    assert "ir_gate" in header and "JLA" in header
    assert len(rows) == 2
    assert "Aligned example" in rows[0] and "Misaligned example" in rows[1]


def test_the_gate_verdict_is_in_the_table_not_left_to_be_rediscovered() -> None:
    """A record the gate refuses is still drawn — and the table says so, with the number."""
    from plastid.run.junction import junction_rows

    aligned, misaligned = junction_rows([(Path("a.gb"), ALIGNED), (Path("b.gb"), MISALIGNED)])
    assert aligned.gate == "passed"
    assert aligned.identity == pytest.approx(1.0)
    assert misaligned.gate == "failed"
    assert misaligned.identity is not None and misaligned.identity < 0.5


def test_a_junction_names_the_genes_around_it() -> None:
    """The column is what lets a reader check the drawing: nearest first, with the distance."""
    from plastid.run.junction import _flanking

    genome = Genome(
        name="genes",
        accession="XX_100003.1",
        sequence="A" * 20_000,
        genes=(
            _gene("rps19", 4_800, 5_400),  # straddles the junction at 5,001
            _gene("rpl2", 5_600, 6_500),  # 599 bp after it
            _gene("trnH", 3_000, 3_100),  # 1,900 bp before it
            _gene("far", 1_000, 1_100),  # out of range
        ),
    )
    text = _flanking(genome, 5_001)
    assert "rps19 (contains 0 bp)" in text
    assert "rpl2 (after 599 bp)" in text
    assert "trnH (before 1901 bp)" in text
    assert "far" not in text


def _gene(name: str, start_1b: int, end_1b: int):
    from plastid.core.model import Gene

    return Gene(name=name, kind="CDS", parts=(Interval(start_1b - 1, end_1b),))


def test_the_command_stops_naming_a_missing_cpjsdraw(tmp_path: Path) -> None:
    """Exit 3 and the tool's name — never a traceback, never a figure that is quietly absent."""
    record = _write(tmp_path, ALIGNED, "aligned.gb")
    code = main(
        [
            "ir-junction",
            "-i",
            str(record),
            "-o",
            str(tmp_path / "runs"),
            "--cpjsdraw",
            str(tmp_path / "not-a-tool.pl"),
        ]
    )
    assert code == EXIT_TOOL


def test_the_command_refuses_a_record_it_cannot_cut_into_four_regions(tmp_path: Path) -> None:
    """A record with no IR annotation is refused by name, not drawn empty."""
    bare = Genome(name="no IR", accession="XX_100004.1", sequence="ACGT" * 5_000)
    record = _write(tmp_path, bare, "bare.gb")
    assert main(["ir-junction", "-i", str(record), "-o", str(tmp_path / "runs")]) == 1


def test_the_command_refuses_a_missing_input(tmp_path: Path) -> None:
    assert main(["ir-junction", "-i", str(tmp_path / "absent.gb"), "-o", str(tmp_path / "r")]) == 1


@pytest.mark.tier_b
def test_a_real_run_draws_the_figure(tmp_path: Path, data_dir: Path) -> None:
    """The whole path, when a patched CPJSdraw is on this machine.

    Pointed at with ``PTS_CPJSDRAW`` — the script is never bundled, so a machine that has it keeps
    it somewhere of its own choosing, exactly like ``PTS_MISA_DIR``.

    The record is the vendored golden, not a synthetic one: CPJSdraw refuses a file with no
    ``gene`` feature ("has no gene info(no primary gene tag)"), and ``to_genbank`` writes CDS, tRNA
    and rRNA features only — measured 17/10/2026, written down in ``docs/tools/cpjsdraw.md``.
    """
    script = os.environ.get("PTS_CPJSDRAW")
    if not script or not Path(script).is_file():
        pytest.skip("no patched CPJSdraw: set PTS_CPJSDRAW=/path/to/CPJSdraw.pl")
    record = data_dir / "NC_050833.1.gb"
    out = tmp_path / "runs"
    code = main(
        ["ir-junction", "-i", str(record), "-o", str(out), "--cpjsdraw", script, "--width", "800"]
    )
    assert code == EXIT_OK
    figures = next(out.glob("*/*/figures"))
    assert (figures / "ir_junction.svg").is_file()
    assert (figures / "ir_junction.png").is_file()
    table = (figures / "ir_junction.csv").read_text(encoding="utf-8")
    assert "NC_050833" in table and "passed" in table
    assert read_genbank(record).accession == "NC_050833.1"
