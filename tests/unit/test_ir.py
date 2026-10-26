"""IR detection and the mandatory gate — the highest-risk component.

The measured failure modes this file guards (each one stays forever):
  * naive detection was silently wrong on 2/9 published plastomes;
  * identity measured on the exact-match core hides divergence (22,998 bp core at 100 %
    while the true 26,162 bp copies are only 92.62 % identical);
  * an off-by-one pairing constant collapses identity to ~29 % — it must fail, loudly.
"""

from __future__ import annotations

import random

import pytest

from plastid.core.ir import (
    IRGateError,
    detect_ir,
    enforce_ir_gate,
    ir_de_novo,
    ir_from_annotation,
    measure_ir,
    mirror_constant_for,
    pair_identity,
)
from plastid.core.model import Genome, Interval
from plastid.core.sequence import reverse_complement


def test_annotation_ir_on_the_golden_genome(golden_genome: Genome) -> None:
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    assert ir.source == "annotation"
    assert ir.length == 26_214
    assert ir.identity == 1.0
    assert ir.mismatches == 0
    assert ir.mirror_constant == 236_850  # = a1 + b2 - 1, the measured pairing constant
    enforce_ir_gate(ir)  # must not raise


def test_off_by_one_pairing_constant_collapses_identity(golden_genome: Genome) -> None:
    ir = ir_from_annotation(golden_genome)
    assert ir is not None
    shifted_identity, _, _, _ = pair_identity(
        golden_genome.sequence, ir.ira, ir.mirror_constant + 1
    )
    assert shifted_identity < 0.5  # measured ~29 %; never silently accepted


def _unequal_pair(
    d: int, *, exact: bool = True, length: int = 600, extra_at: str = "outer"
) -> tuple[str, Interval, Interval]:
    """A plastome with an inverted repeat whose **first copy is ``d`` bases longer**.

    ``d`` may be negative — the second copy longer, which the annotation path also tolerates.
    ``extra_at="outer"`` puts the extra bases at the pair's outer ends (the first copy's start,
    or the second copy's end); ``"inner"`` puts them at the inner ends. ``exact=False`` scrambles
    the second copy so the pair is not a repeat at all — the negative case the gate must refuse.
    """
    rng = random.Random(11)
    core = "".join(rng.choices("ACGT", k=length))
    first_copy = core
    second_copy = reverse_complement(core)
    if d > 0:
        extra = "".join(rng.choices("ACGT", k=d))
        first_copy = extra + core if extra_at == "outer" else core + extra
    elif d < 0:
        extra = "".join(rng.choices("ACGT", k=-d))
        rev = reverse_complement(core)
        second_copy = rev + extra if extra_at == "outer" else extra + rev
    if not exact:
        second_copy = "".join(rng.choices("ACGT", k=len(second_copy)))
    # LSC · first copy · SSC · second copy · a tail, so neither copy touches an end
    sequence = "ACGT" * 40 + first_copy + "TTAA" * 40 + second_copy + "GGCC" * 40
    first_start = 160
    second_start = first_start + len(first_copy) + 160
    return (
        sequence,
        Interval(first_start, first_start + len(first_copy), 1),
        Interval(second_start, second_start + len(second_copy), -1),
    )


@pytest.mark.parametrize("d", [1, 2, 4])
def test_a_length_difference_is_measured_at_the_end_it_sits_at(d: int) -> None:
    """H-25: the mirror constant assumed equal copies, so *d* bases of excess made an exact
    repeat measure as unrelated (~25 %). Measured 17/10/2026 — Chasmanthe NC_072562 (*d* = 1),
    Iris hermona NC_081931 (*d* = 2) and Trimezia NC_072662 (*d* = 4) were all refused."""
    sequence, ira, irb = _unequal_pair(d)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == 1.0, f"d={d}: an exact repeat measured {result.identity:.4f}"
    assert result.mismatches == 0
    assert result.compared_length == 600  # the shorter copy, unchanged
    assert result.excess_at == "outer"


@pytest.mark.parametrize("d", [-1, -3])
def test_the_second_copy_may_be_the_longer_one(d: int) -> None:
    """The annotation path tolerates a difference either way, so the rule must too.

    Trimming the *longer* copy is what makes this symmetric; shifting only the first copy's span
    would read bases outside the copy and report a number that means nothing.
    """
    sequence, ira, irb = _unequal_pair(d)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == 1.0, f"d={d}: an exact repeat measured {result.identity:.4f}"
    assert result.excess_at == "outer"
    assert result.compared_length == 600


def test_a_pair_that_is_not_a_repeat_still_measures_as_unrelated() -> None:
    """The two alignments a length difference implies must not become a search for a match.

    A pair whose copies are simply different sequences scores ~25 % under *both* alignments, so
    neither can carry it over the gate.
    """
    sequence, ira, irb = _unequal_pair(3, exact=False)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity < 0.5
    with pytest.raises(IRGateError):
        enforce_ir_gate(result)


def test_the_excess_is_taken_at_the_inner_end_when_that_is_where_it_sits() -> None:
    sequence, ira, irb = _unequal_pair(4, extra_at="inner")
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == 1.0
    assert result.excess_at == "inner"


def test_a_difference_with_no_end_named_is_refused_rather_than_defaulted() -> None:
    """An `IRResult` says "none" only when the copies are equal. Anything else must name an end —
    a silent default would make the mirror constant depend on a choice nobody made."""
    _, ira, irb = _unequal_pair(3)
    with pytest.raises(ValueError, match="must say where the difference sits"):
        mirror_constant_for(ira, irb, "none")


def test_equal_copies_are_measured_exactly_as_before() -> None:
    """*d* = 0 is the golden genome's case: a correct fix is a no-op there."""
    sequence, ira, irb = _unequal_pair(0)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == 1.0
    assert result.excess_at == "none"  # no excess to place: the copies are the same length
    assert result.mirror_constant == ira.start + irb.end - 1  # the measured pairing rule


def test_identity_is_measured_on_full_copies_not_the_exact_core(
    edge_genomes: dict[str, Genome],
) -> None:
    # NC_063920: a perfect 22,998 bp exact core, but the full 26,162 bp copies are 92.62 %
    # identical. The gate must use the full copies.
    ir = ir_from_annotation(edge_genomes["NC_063920.1"])
    assert ir is not None
    assert ir.exact_core == 22_998
    # 1,929 mismatches over the 26,161 compared bases (the shorter copy; the two differ by
    # 1 bp). The original measurement rounds this to 92.63 %.
    assert ir.compared_length == 26_161
    assert ir.mismatches == 1_929
    assert ir.identity == pytest.approx(1 - 1_929 / 26_161)
    with pytest.raises(IRGateError):
        enforce_ir_gate(ir)


def test_gate_flags_the_second_diverged_genome(edge_genomes: dict[str, Genome]) -> None:
    # Re-measured 2026-10-02 with the structural pairing (a1 + b2 - 1): 98.87 %, not the
    # 53.65 %. Verdict unchanged: below 99.9 % → gate fails.
    ir = ir_from_annotation(edge_genomes["NC_063921.1"])
    assert ir is not None
    assert ir.identity == pytest.approx(0.988695, abs=1e-5)
    with pytest.raises(IRGateError):
        enforce_ir_gate(ir)


def test_short_repeat_regions_are_ignored(edge_genomes: dict[str, Genome]) -> None:
    # NC_063921 annotates five repeat_regions, three of them spurious (1 bp, 393 bp, 2 bp).
    # Structural filtering (length >= 5 kb, exactly two candidates) must pick the real pair.
    genome = edge_genomes["NC_063921.1"]
    assert len(genome.repeat_regions) == 5
    ir = ir_from_annotation(genome)
    assert ir is not None
    assert ir.length == 26_183


def test_copies_of_unequal_length_are_measured_over_the_shorter(
    edge_genomes: dict[str, Genome],
) -> None:
    # NC_063920 annotates 26,162 bp and 26,161 bp copies — a 1 bp difference.
    genome = edge_genomes["NC_063920.1"]
    ir = ir_from_annotation(genome)
    assert ir is not None
    assert ir.compared_length == 26_161
    assert ir.length == 26_161
    # …and its excess sits at the **inner** end: 92.63 % here against 36.46 % at the outer end.
    # This is the measurement that says H-25's fix must weigh both ends rather than assume one —
    # "take it at the outer end" would have turned the 92.63 % finding into 36.46 % and taken
    # the corpus's known-divergent case with it.
    assert ir.excess_at == "inner"
    assert ir.identity == pytest.approx(1 - 1_929 / 26_161)


def test_de_novo_reproduces_the_two_known_failures(edge_genomes: dict[str, Genome]) -> None:
    # Measured: the validated detector matches GenBank on 7/9 annotated genomes and
    # fails on exactly these two — one loudly, one silently. Both must stay visible.
    loud = ir_de_novo(edge_genomes["NC_063921.1"])
    assert loud is None  # no candidate passes the validity window — loud failure

    silent = ir_de_novo(edge_genomes["NC_063920.1"])
    assert silent is not None
    assert silent.length == 22_998  # the exact-match core, shorter than the true 26,161


def test_detect_ir_records_the_de_novo_disagreement(edge_genomes: dict[str, Genome]) -> None:
    ir = detect_ir(edge_genomes["NC_063920.1"])
    assert ir.source == "annotation"  # annotation wins ...
    assert ir.length == 26_161
    assert ir.denovo_length == 22_998
    assert ir.denovo_delta == 3_163  # ... but the disagreement is stored, never dropped


def test_de_novo_agrees_with_annotation_on_the_golden(golden_genome: Genome) -> None:
    ir = detect_ir(golden_genome)
    assert ir.denovo_length == 26_214
    assert ir.denovo_delta == 0


def test_de_novo_on_a_genome_without_any_annotation(edge_genomes: dict[str, Genome]) -> None:
    genome = edge_genomes["NC_085783.1"]
    assert genome.repeat_regions == ()
    ir = detect_ir(genome)
    assert ir.source == "denovo_diag"
    assert ir.length == 21_148
    assert ir.identity == 1.0
    enforce_ir_gate(ir)


def test_measure_ir_rejects_intervals_that_are_not_a_pair() -> None:
    from plastid.core.model import Interval

    with pytest.raises(ValueError, match="disjoint"):
        measure_ir("ACGT" * 10, Interval(0, 10), Interval(5, 15), source="manual")
