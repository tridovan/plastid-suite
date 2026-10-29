"""Reproductions of the review's defects, with independently derived expected answers (M11.2).

Each test names the issue in [docs/known_issues.md](../../docs/known_issues.md) and the milestone
that owned it. K-01, K-03 and K-04 were fixed by M11.3/M11.4 and their markers were removed in the
same change; **K-06 was fixed by the M11.7 forwarding work** and its marker was removed with the
fix. Their expected values come from what the record was constructed to be — the golden's
annotated 26,214 bp copies, the synthetic fixture's LSC 1,000 / IR 6,000 / SSC 1,000 — never from
the reader or detector under test.

K-02 and K-05 are exercised where their fixes live: ``test_run.py`` (cache) and
``test_normalize.py`` (rotation). They are not duplicated here.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.ir import (
    IRDetectionError,
    IRGateError,
    detect_ir,
    enforce_ir_gate,
    ir_from_annotation,
)
from plastid.core.junctions import junction_regions
from plastid.core.model import Genome, Interval
from plastid.core.normalize import rotate_genome
from plastid.result.export import export_json
from plastid.result.result import Column, Result, Scope, Table

DATA = Path(__file__).parents[1] / "data"
IR_LENGTH = 26_214
GOLDEN_LENGTH = 153_724


@pytest.fixture(scope="module")
def golden() -> Genome:
    return read_genbank(DATA / "NC_050833.1.gb")


@pytest.fixture(scope="module")
def shapes() -> Genome:
    """The vendored synthetic record: LSC 1,000 · IR 6,000 · SSC 1,000 · IR 6,000."""
    return read_genbank(DATA / "synthetic_ir_shapes.gb")


# --------------------------------------------------------------------------------------
# K-01 — IR extent and identity depend on record shape (M11.3/M11.4)
# --------------------------------------------------------------------------------------


def _without_annotations(genome: Genome) -> Genome:
    return replace(genome, repeat_regions=())


def test_k01_a_single_snp_does_not_cut_the_de_novo_ir(golden: Genome) -> None:
    """One base differs in the second copy; the repeat is still the whole 26,214 bp pair."""
    position = 140_617
    replacement = "A" if golden.sequence[position] != "A" else "C"
    sequence = golden.sequence[:position] + replacement + golden.sequence[position + 1 :]
    mutated = _without_annotations(replace(golden, sequence=sequence))

    ir = detect_ir(mutated)

    assert ir.length == IR_LENGTH
    assert ir.mismatches == 1
    assert ir.identity == pytest.approx((IR_LENGTH - 1) / IR_LENGTH)


def test_k01_rotation_does_not_change_the_de_novo_extent(golden: Genome) -> None:
    """Rotating right by 5,000 bp is the same circle, so the IR is still 26,214 bp."""
    rotated = _without_annotations(rotate_genome(golden, len(golden) - 5_000))

    ir = detect_ir(rotated)

    assert ir.length == IR_LENGTH
    assert ir.identity == pytest.approx(1.0)
    assert ir.gate_passed


def test_k01_a_one_base_deletion_is_aligned_as_a_gap(golden: Genome) -> None:
    """The deletion shortens the second copy by one base; the other 26,213 bases still pair."""
    position = 140_617
    sequence = golden.sequence[:position] + golden.sequence[position + 1 :]
    deleted = replace(
        golden,
        sequence=sequence,
        # the annotation is adjusted to the shorter second copy, as a real record would be
        repeat_regions=(Interval(83_127, 109_341), Interval(127_510, 153_723)),
    )

    ir = ir_from_annotation(deleted)

    assert ir is not None
    assert ir.compared_length == IR_LENGTH
    assert ir.matched == IR_LENGTH - 1
    assert ir.mismatches == 0
    assert ir.internal_gaps == 1
    assert ir.identity == pytest.approx((IR_LENGTH - 1) / IR_LENGTH)
    assert ir.gate_passed


def test_k01_a_joined_copy_keeps_its_whole_length(golden: Genome) -> None:
    """Rotating one base splits a copy into a 26,213 bp run and a 1 bp piece; the copy is both."""
    rotated = rotate_genome(golden, len(golden) - 1)
    assert sum(region.length for region in rotated.repeat_regions) == 2 * IR_LENGTH

    ir = ir_from_annotation(rotated)

    assert ir is not None
    assert ir.ira.length == IR_LENGTH
    assert ir.irb.length == IR_LENGTH
    assert ir.identity == pytest.approx(1.0)


# --------------------------------------------------------------------------------------
# K-01 — shape coverage on the small synthetic record (M11.2)
# --------------------------------------------------------------------------------------


def test_the_synthetic_fixture_is_the_pair_it_was_built_as(shapes: Genome) -> None:
    """Pass: the construction itself is the oracle — LSC 1,000 · IR 6,000 · SSC 1,000."""
    ir = ir_from_annotation(shapes)
    assert ir is not None
    assert (ir.ira.start, ir.ira.end, ir.ira.length) == (1_000, 7_000, 6_000)
    assert (ir.irb.start, ir.irb.end, ir.irb.length) == (8_000, 14_000, 6_000)
    assert ir.identity == pytest.approx(1.0)
    assert ir.compared_length == 6_000


@pytest.mark.parametrize("shift", [500, 1_000])
def test_an_origin_in_a_single_copy_or_at_an_ir_boundary_keeps_the_pair(
    shapes: Genome, shift: int
) -> None:
    """Pass: these rotations leave both annotated copies whole."""
    ir = ir_from_annotation(rotate_genome(shapes, shift))
    assert ir is not None
    assert ir.length == 6_000
    assert ir.identity == pytest.approx(1.0)


def test_an_origin_inside_an_ir_copy_keeps_the_pair(shapes: Genome) -> None:
    """The molecule is unchanged; both copies are still 6,000 bp from an origin in their middle."""
    rotated = rotate_genome(shapes, 4_000)

    ir = ir_from_annotation(rotated)

    assert ir is not None
    assert ir.length == 6_000
    assert ir.identity == pytest.approx(1.0)


def test_a_joined_copy_with_a_small_and_a_large_part_keeps_its_length(shapes: Genome) -> None:
    """The join is one 6,000 bp copy: a 5,999 bp run plus the 1 bp piece at the record's end."""
    rotated = rotate_genome(shapes, 8_001)

    ir = ir_from_annotation(rotated)

    assert ir is not None
    assert ir.length == 6_000
    assert ir.identity == pytest.approx(1.0)


def test_the_annotation_strand_does_not_change_the_measurement(shapes: Genome) -> None:
    """Pass: a repeat_region on strand -1 describes the same bases as the same span on +1."""
    flipped = replace(
        shapes,
        repeat_regions=(
            Interval(1_000, 7_000, strand=-1),
            Interval(8_000, 14_000, strand=-1),
        ),
    )
    ir = ir_from_annotation(flipped)
    assert ir is not None
    assert ir.length == 6_000
    assert ir.identity == pytest.approx(1.0)


@pytest.mark.parametrize("offset", [0, 3_000, 5_999])
def test_a_snp_is_counted_at_the_start_middle_or_end(shapes: Genome, offset: int) -> None:
    """Pass: one known mismatch in the second copy, wherever it sits, is one mismatch."""
    position = 8_000 + offset
    replacement = "A" if shapes.sequence[position] != "A" else "C"
    sequence = shapes.sequence[:position] + replacement + shapes.sequence[position + 1 :]

    ir = ir_from_annotation(replace(shapes, sequence=sequence))

    assert ir is not None
    assert ir.mismatches == 1
    assert ir.compared_length == 6_000
    assert ir.identity == pytest.approx(5_999 / 6_000)


@pytest.mark.parametrize("offset", [0, 5_999])
def test_an_indel_at_an_end_is_aligned_as_a_gap(shapes: Genome, offset: int) -> None:
    """Pass (the trimmed copy aligns): the deletion is a gap, not 25 % identity."""
    position = 8_000 + offset
    sequence = shapes.sequence[:position] + shapes.sequence[position + 1 :]
    deleted = replace(
        shapes,
        sequence=sequence,
        repeat_regions=(Interval(1_000, 7_000), Interval(8_000, 13_999)),
    )

    ir = ir_from_annotation(deleted)

    assert ir is not None
    assert ir.compared_length == 6_000
    assert ir.matched == 5_999
    assert ir.mismatches == 0
    assert ir.terminal_gaps == 1
    assert ir.identity == pytest.approx(5_999 / 6_000)
    assert ir.gate_passed


def test_an_indel_in_the_middle_is_aligned_as_a_gap(shapes: Genome) -> None:
    """The other 5,999 bases pair; a middle deletion is one internal gap, not 2,251 mismatches."""
    position = 8_000 + 3_000
    sequence = shapes.sequence[:position] + shapes.sequence[position + 1 :]
    deleted = replace(
        shapes,
        sequence=sequence,
        repeat_regions=(Interval(1_000, 7_000), Interval(8_000, 13_999)),
    )

    ir = ir_from_annotation(deleted)

    assert ir is not None
    assert ir.compared_length == 6_000
    assert ir.matched == 5_999
    assert ir.mismatches == 0
    assert ir.internal_gaps == 1
    assert ir.terminal_gaps == 0
    assert ir.identity == pytest.approx(5_999 / 6_000)
    assert ir.gate_passed


def test_an_annotation_that_disagrees_with_the_sequence_fails_the_gate(shapes: Genome) -> None:
    """Pass: measured identity follows the sequence, not the annotation's claim."""
    sequence = list(shapes.sequence)
    for position in range(8_000, 14_000, 60):  # 100 substitutions
        original = sequence[position]
        sequence[position] = "A" if original != "A" else "C"
    disagreeing = replace(shapes, sequence="".join(sequence))

    ir = ir_from_annotation(disagreeing)

    assert ir is not None
    assert ir.mismatches == 100
    assert ir.identity == pytest.approx(5_900 / 6_000)
    assert not ir.gate_passed
    with pytest.raises(IRGateError):
        enforce_ir_gate(ir)


def test_a_genome_with_no_valid_pair_is_refused() -> None:
    """Pass: no annotation and no de-novo repeat is an error, never a silent empty result."""
    genome = Genome("no-ir", "ACGT" * 3_500)
    with pytest.raises(IRDetectionError, match="no valid IR pair"):
        detect_ir(genome)


def test_the_sequence_source_is_explicit_and_records_the_annotation(
    data_dir: Path,
) -> None:
    """D11-2: annotation-first refuses; the sequence path is reached only by asking for it.

    NC_024936.1's annotated copies are 26,220/26,186 bp and align at 99.80 %, so the default path
    refuses the run. ``ir_source=sequence`` uses the sequence's own repeat and records the
    annotation it did not use, so the disagreement is visible in the manifest.
    """
    genome = read_genbank(data_dir / "NC_024936.1.gb")

    with pytest.raises(IRGateError, match="ir_source=sequence"):
        enforce_ir_gate(detect_ir(genome))

    chosen = detect_ir(genome, source="sequence")
    assert chosen.source == "denovo_diag"
    assert chosen.annotated_lengths == (26_220, 26_186)
    assert chosen.annotation_note is None
    assert chosen.identity >= 0.999
    enforce_ir_gate(chosen)


def test_a_malformed_annotation_is_named_and_only_sequence_source_bypasses_it(
    golden: Genome,
) -> None:
    """One/three/overlapping copies are refused by name, not silently re-read (D11-2)."""
    three = replace(
        golden,
        repeat_regions=golden.repeat_regions + (Interval(5_000, 30_000),),
    )
    with pytest.raises(IRDetectionError, match="expected exactly two"):
        detect_ir(three)

    chosen = detect_ir(three, source="sequence")
    assert chosen.source == "denovo_diag"
    assert chosen.length == IR_LENGTH
    assert chosen.annotated_lengths == ()  # the malformed annotation supplied no usable lengths
    assert chosen.annotation_note is not None
    assert "expected exactly two" in chosen.annotation_note


# --------------------------------------------------------------------------------------
# K-03 — the junction reader cannot serve as an independent oracle (M11.3)
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("shift", [83_127, 127_510])
def test_k03_a_copy_starting_at_position_one_is_one_copy(golden: Genome, shift: int) -> None:
    """The same two copies and gaps, with one copy beginning at GenBank position 1."""
    regions = junction_regions(rotate_genome(golden, shift))

    assert [(region.name, region.length) for region in regions.regions] == [
        ("LSC", 83_127),
        ("IRb", 26_214),
        ("SSC", 18_169),
        ("IRa", 26_214),
    ]
    assert regions.tiles_genome


@pytest.mark.parametrize("shift", [88_127, 96_234])
def test_k03_a_joined_copy_with_two_large_parts_is_one_copy(golden: Genome, shift: int) -> None:
    """Both pieces are at least 5 kb; together they are still one 26,214 bp copy."""
    rotated = rotate_genome(golden, shift)
    assert len(rotated.repeat_regions) == 2  # one joined copy, not three parts

    regions = junction_regions(rotated)

    assert [(region.name, region.length) for region in regions.regions] == [
        ("LSC", 83_127),
        ("IRb", 26_214),
        ("SSC", 18_169),
        ("IRa", 26_214),
    ]


# --------------------------------------------------------------------------------------
# K-04 — unknown bases are counted as matching evidence (M11.4)
# --------------------------------------------------------------------------------------


def test_k04_unknown_bases_are_not_matching_evidence(golden: Genome) -> None:
    """Every annotated IR base is unknown: ambiguous columns, no observed matches, no PASS.

    The denominator stays the full alignment (26,214 ambiguous columns); unknowns are not dropped
    from it, so a tiny known fragment cannot look fully supported.
    """
    sequence = list(golden.sequence)
    for start, end in ((83_127, 109_341), (127_510, 153_724)):
        sequence[start:end] = ["N"] * (end - start)
    unknown = replace(golden, sequence="".join(sequence))

    ir = ir_from_annotation(unknown)

    assert ir is not None
    assert ir.compared_length == 26_214
    assert ir.matched == 0
    assert ir.known_columns == 0
    assert ir.ambiguous_columns == 26_214
    assert ir.identity == 0.0
    assert not ir.gate_passed


# --------------------------------------------------------------------------------------
# K-06 — the CLI did not forward tree options to R (fixed by M11.7's forwarding work)
# --------------------------------------------------------------------------------------
#
# The marker was removed with the fix: `_render_r_figures` now passes the tree script's
# `--set root/outgroup/drop-tips` (and its drawing options) through to R, and the script no longer
# defaults to the manuscript's tips or outgroup. The check below holds the argv, which is what the
# defect was about; the output tree is the R script's own concern and is exercised by the R figure
# pass on a machine with R.


def test_k06_the_tree_figure_forwards_its_options(tmp_path: Path) -> None:
    """``--set root=none --set outgroup=... --set drop-tips=...`` must reach the script argv."""
    from plastid.cli import _render_r_figures

    run_dir = tmp_path / "dataset" / "run1"
    (run_dir / "stub").mkdir(parents=True)
    result = Result(
        module="stub",
        module_version="0.1.0",
        tables={
            "stub": Table(
                id="stub",
                title="stub",
                columns=(Column("value", "Value"),),
                rows=({"value": 1},),
                scope=Scope(label="synthetic"),
            )
        },
    )
    (run_dir / "stub" / "result.json").write_text(export_json(result), encoding="utf-8")

    tree = tmp_path / "tree.treefile"
    tree.write_text("((A:1,B:1):1,C:1);\n", encoding="utf-8")
    rscript = tmp_path / "Rscript"
    rscript.write_text(
        "#!/bin/sh\n"
        'case "$1" in\n'
        "  --version|-version|-v|-V) echo 'R 4.1.2'; exit 0;;\n"
        "  -e) exit 0;;\n"
        "esac\n"
        'out=""\n'
        'prev=""\n'
        'for a in "$@"; do\n'
        '  if [ "$prev" = "--outdir" ]; then out="$a"; fi\n'
        '  prev="$a"\n'
        "done\n"
        'mkdir -p "$out/Figures"\n'
        'echo png > "$out/Figures/tree.png"\n',
        encoding="utf-8",
    )
    rscript.chmod(0o755)

    args = SimpleNamespace(
        figure=["r_tree"],
        sets=[
            f"tree={tree}",
            f"rscript={rscript}",
            "root=none",
            "outgroup=Outgroup one,Outgroup two",
            "drop-tips=Tip one",
        ],
    )
    outcome, ctx = _render_r_figures(run_dir, sets=args.sets, figures=args.figure)

    assert outcome.files, "the fake Rscript should have produced a figure"
    argv = next(record.argv for record in ctx.tools if record.name == "rscript:tree")
    assert "--root" in argv and argv[argv.index("--root") + 1] == "none"
    assert "--outgroup" in argv and argv[argv.index("--outgroup") + 1] == (
        "Outgroup one,Outgroup two"
    )
    assert "--drop-tips" in argv and argv[argv.index("--drop-tips") + 1] == "Tip one"


# The K-02/K-05 paths live in their own files; this assertion keeps the pointer honest.
def test_the_fixed_regressions_are_owned_by_their_module_tests() -> None:
    root = Path(__file__).parent
    run_tests = (root / "test_run.py").read_text(encoding="utf-8")
    normalize_tests = (root / "test_normalize.py").read_text(encoding="utf-8")
    assert "test_a_changed_translation_table_is_not_served_the_cached_orf_result" in run_tests
    assert "test_a_changed_tool_version_invalidates_a_cached_payload" in run_tests
    assert "test_rotation_keeps_partial_flags_on_the_endpoints_they_describe" in normalize_tests
    assert "test_read_rotate_write_read_keeps_the_annotation_model" in normalize_tests
