"""M11.17 panel cases: the same sequence under two annotations.

The 30 `annotation-variants/` records are PGA's annotation of the **same sequence** as
`source/` (the panel README states it; the first test checks it). That gives an independent
expectation for every module that reads only the sequence: its result must not move when the
annotation does. A module whose output *does* move is annotation-dependent by definition, and the
case says so rather than calling the difference an error.

Three modules are covered here — `genome_stats` (length/GC), `ssr` (MISA) and `repeats`
(repeat-match). `ssr`/`repeats` need `misa.pl`/`repeat-match`; the tests skip with the tool's name
when absent, as every tier-B check does.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import require_tools
from plastid.core.genbank import read_genbank
from plastid.core.ir import detect_ir, enforce_ir_gate
from plastid.modules import genome_stats, repeats, ssr
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / ".panel" / "ir-validation" / "source"
VARIANTS = ROOT / ".panel" / "ir-validation" / "annotation-variants"
ACCESSIONS = ("NC_008641.1", "NC_016708.2", "NC_028027.1")

pytestmark = [
    pytest.mark.tier_b,
    pytest.mark.skipif(
        not SOURCE.is_dir() or not VARIANTS.is_dir(),
        reason="the IR validation panels are not in this tree (.panel/ is gitignored)",
    ),
]


def _ctx(tmp_path: Path, module: str) -> RunContext:
    layout = RunLayout.create(tmp_path, module, fingerprint="panel")
    layout.make_dirs()
    return RunContext(layout=layout, params={}, module_id=module, module_order=1)


def _gate(genome) -> None:
    return enforce_ir_gate(detect_ir(genome, source="sequence"))


def test_the_pair_is_one_sequence_with_two_annotations() -> None:
    """The premise of every case below: only the annotation differs."""
    for accession in ACCESSIONS:
        source = read_genbank(SOURCE / f"{accession}.gb")
        variant = read_genbank(VARIANTS / f"{accession}.gb")
        assert source.sequence == variant.sequence, f"{accession}: the sequence differs"
        assert len(source) == len(variant)


def test_genome_stats_length_and_gc_do_not_depend_on_the_annotation(tmp_path: Path) -> None:
    """Length and GC are properties of the sequence; an annotation change cannot move them."""
    for accession in ACCESSIONS:
        values = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "genome_stats")
            ctx.ir = _gate(genome)
            row = genome_stats.run(genome, {}, ctx).tables["genome_stats"].rows[0]
            values.append((row["length_bp"], row["gc_pct"]))
        assert values[0] == values[1], f"{accession}: {values}"


def test_ssr_totals_do_not_depend_on_the_annotation(
    tools: dict[str, str | None], tmp_path: Path
) -> None:
    """MISA reads the sequence; the IR is the sequence's own, so every SSR row must be identical."""
    require_tools(tools, "perl", "misa.pl")
    for accession in ACCESSIONS:
        tables = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "ssr")
            ctx.ir = _gate(genome)
            tables.append(ssr.run(genome, {"misa": tools["misa.pl"]}, ctx).tables)
        source_summary = tables[0]["ssr_summary"].rows[0]
        variant_summary = tables[1]["ssr_summary"].rows[0]
        for key in ("total", "mono", "di", "tri", "compound", "thresholds"):
            assert source_summary[key] == variant_summary[key], f"{accession}.{key}"
        rows = [
            [
                (row["motif"], row["unit"], row["start_1b"], row["end_1b"], row["compound"])
                for row in table["ssrs"].rows
            ]
            for table in tables
        ]
        assert rows[0] == rows[1], f"{accession}: SSR positions moved with the annotation"


def test_exact_repeats_do_not_depend_on_the_annotation(
    mummer_tools: dict[str, str | None], tmp_path: Path
) -> None:
    """repeat-match reads the sequence; the same sequence gives the same pairs."""
    require_tools(mummer_tools, "repeat-match")
    directory = str(Path(mummer_tools["repeat-match"]).parent)
    for accession in ACCESSIONS:
        tables = []
        for folder in (SOURCE, VARIANTS):
            genome = read_genbank(folder / f"{accession}.gb")
            ctx = _ctx(tmp_path / folder.name / accession, "repeats")
            ctx.ir = _gate(genome)
            tables.append(repeats.run(genome, {"repeats": directory}, ctx).tables)
        source_rows = [
            (
                row["class"],
                row["length_bp"],
                row["start1_1b"],
                row["end1_1b"],
                row["start2_1b"],
                row["end2_1b"],
            )
            for row in tables[0]["repeats"].rows
        ]
        variant_rows = [
            (
                row["class"],
                row["length_bp"],
                row["start1_1b"],
                row["end1_1b"],
                row["start2_1b"],
                row["end2_1b"],
            )
            for row in tables[1]["repeats"].rows
        ]
        assert source_rows == variant_rows, f"{accession}: exact repeats moved with the annotation"
