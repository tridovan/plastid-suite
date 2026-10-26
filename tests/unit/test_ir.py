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
    align_copies,
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
    Iris hermona NC_081931 (*d* = 2) and Trimezia NC_072662 (*d* = 4) were all refused.

    The end-to-end alignment (M11.4) counts the *d* overhang bases as terminal gap columns: every
    known base matches, but the identity denominator includes the overhang. On a 600 bp
    constructed pair that is 600/(600+d); at plastome scale the same overhang is far below 0.1 %.
    """
    sequence, ira, irb = _unequal_pair(d)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == pytest.approx(600 / (600 + d)), f"d={d}"
    assert result.mismatches == 0
    assert result.terminal_gaps == d
    assert result.internal_gaps == 0
    assert result.compared_length == 600 + d
    assert result.excess_at == "outer"


@pytest.mark.parametrize("d", [-1, -3])
def test_the_second_copy_may_be_the_longer_one(d: int) -> None:
    """The annotation path tolerates a difference either way, so the rule must too.

    Counting from either copy is symmetric: the overhang is *d* terminal gap columns at the same
    outer end, whichever copy carries the extra bases.
    """
    sequence, ira, irb = _unequal_pair(d)
    result = measure_ir(sequence, ira, irb, source="annotation")
    assert result.identity == pytest.approx(600 / (600 + abs(d)))
    assert result.excess_at == "outer"
    assert result.terminal_gaps == abs(d)
    assert result.compared_length == 600 + abs(d)


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
    assert result.identity == pytest.approx(600 / 604)
    assert result.excess_at == "inner"
    assert result.terminal_gaps == 4


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
    """NC_063920: the exact core is 22,998 bp; the annotated copies are 26,162/26,161 bp.

    Re-measured 26/10/2026 with the end-to-end gap-aware alignment (M11.4): the old ungapped
    mirror pairing reported 92.63 % because a 1 bp indel shifted every base after it. The full
    copies differ by 10 mismatches and 5 internal gaps over 26,164 columns — 99.94 %, so the gate
    passes. The core is still shorter than the copies, which is the point here: the measurement
    is made on the copies, not on the exact core.
    """
    ir = ir_from_annotation(edge_genomes["NC_063920.1"])
    assert ir is not None
    assert ir.exact_core == 22_998
    assert ir.length == 26_161
    assert ir.exact_core < ir.length
    assert ir.matched == 26_149
    assert ir.mismatches == 10
    assert ir.internal_gaps == 5
    assert ir.compared_length == 26_164
    assert ir.identity == pytest.approx(26_149 / 26_164)
    enforce_ir_gate(ir)  # the old 92.63 % was an ungapped-alignment artefact, not divergence


def test_gate_flags_the_second_diverged_genome(edge_genomes: dict[str, Genome]) -> None:
    """NC_063921's copies are 99.16 % identical: below 99.9 %, so the gate still refuses."""
    ir = ir_from_annotation(edge_genomes["NC_063921.1"])
    assert ir is not None
    assert ir.matched == 25_999
    assert ir.mismatches == 147
    assert ir.compared_length == 26_219
    assert ir.identity == pytest.approx(25_999 / 26_219)
    assert ir.identity < 0.999
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


def test_copies_of_unequal_length_are_measured_end_to_end(
    edge_genomes: dict[str, Genome],
) -> None:
    """NC_063920 annotates 26,162 bp and 26,161 bp copies — a 1 bp difference.

    The gap-aware alignment finds no terminal overhang: the length difference is absorbed by five
    internal gaps, so ``excess_at`` is ``"none"`` and the denominator counts all 26,164 columns.
    """
    genome = edge_genomes["NC_063920.1"]
    ir = ir_from_annotation(genome)
    assert ir is not None
    assert ir.annotated_lengths == (26_162, 26_161)
    assert ir.compared_length == 26_164
    assert ir.length == 26_161
    assert ir.excess_at == "none"
    assert ir.terminal_gaps == 0
    assert ir.internal_gaps == 5
    assert ir.identity == pytest.approx(26_149 / 26_164)


def test_de_novo_finds_shorter_repeats_inside_the_two_diverged_records(
    edge_genomes: dict[str, Genome],
) -> None:
    """The sequence path must not return the full divergent pair, and must not hide it either.

    Re-measured 26/10/2026 (M11.4): NC_063921 contains a perfect 16,197 bp reverse-complement
    repeat inside the 26,183 bp annotated pair; NC_063920 contains a 23,446 bp repeat
    (99.945 %) whose exact core is 22,998 bp. Both are shorter than the annotated copies, so the
    annotation path is what measures the full pair; the de-novo lengths are cross-checks.
    """
    short = ir_de_novo(edge_genomes["NC_063921.1"])
    assert short is not None
    assert short.length == 16_197
    assert (short.matched, short.mismatches, short.ambiguous_columns) == (16_196, 0, 1)
    assert short.identity == pytest.approx(16_196 / 16_197)
    assert short.length < 26_183

    shorter = ir_de_novo(edge_genomes["NC_063920.1"])
    assert shorter is not None
    assert shorter.length == 23_446
    assert shorter.exact_core == 22_998
    assert (shorter.matched, shorter.mismatches, shorter.internal_gaps) == (23_435, 9, 4)
    assert shorter.length < 26_161


def test_detect_ir_records_the_de_novo_disagreement(edge_genomes: dict[str, Genome]) -> None:
    ir = detect_ir(edge_genomes["NC_063920.1"])
    assert ir.source == "annotation"  # annotation wins ...
    assert ir.length == 26_161
    assert ir.denovo_length == 23_446
    assert ir.denovo_delta == 2_715  # ... but the disagreement is stored, never dropped


def test_de_novo_agrees_with_annotation_on_the_golden(golden_genome: Genome) -> None:
    ir = detect_ir(golden_genome)
    assert ir.denovo_length == 26_214
    assert ir.denovo_delta == 0


def test_de_novo_on_a_genome_without_any_annotation(edge_genomes: dict[str, Genome]) -> None:
    """NC_085783.1 has no repeat_region; the origin sits inside one copy.

    The old detector returned a 21,148 bp truncated run. The chained detector finds 26,184 bp;
    on the rotated sequence the two copies are exact reverse complements over that whole extent,
    and extending by one base breaks the equality — independent evidence that this is the full
    copy, not an over-extension.
    """
    genome = edge_genomes["NC_085783.1"]
    assert genome.repeat_regions == ()
    ir = detect_ir(genome)
    assert ir.source == "denovo_diag"
    assert ir.length == 26_184
    assert ir.identity == pytest.approx(1.0)
    enforce_ir_gate(ir)


def test_overlapping_exact_runs_at_an_indel_are_trimmed() -> None:
    """The panel's unequal-copy shape, constructed so the two exact runs overlap by the motif.

    ``b`` carries one extra copy of ``motif``; the exact run after the insertion starts inside the
    run before it. Requiring non-overlapping anchors dropped one run and mis-measured the rest
    (86.7 % on NC_065154.1). The construction gives the answer directly: every base of ``a``
    matches, and the ten inserted bases are gap columns — 2,000/2,010.
    """
    rng = random.Random(20261026)
    motif = "ACGTACGTAC"
    x = "".join(rng.choice("ACGT") for _ in range(990)) + motif
    y = "".join(rng.choice("ACGT") for _ in range(1_000))
    a = x + y
    b = x + motif + y

    stats = align_copies(a, b)

    assert stats.matched == len(a)
    assert stats.mismatches == 0
    assert stats.internal_gaps + stats.terminal_gaps == len(motif)
    assert stats.columns == len(a) + len(motif)
    assert stats.identity == pytest.approx(len(a) / (len(a) + len(motif)))


def test_a_repetitive_pair_does_not_explode_the_anchor_index() -> None:
    """A homopolymer has no informative 20-mer; the measurement must fall back, not run for hours.

    The 32-occurrence cap removes low-complexity words from the index. The fallback still counts
    observed A/C/G/T columns, so an identical homopolymer pair is 100 % and an all-``N`` pair has
    no known evidence at all.
    """
    repeated = align_copies("A" * 5_000, "A" * 5_000)
    assert repeated.method == "ungapped_fallback"
    assert (repeated.identity, repeated.matched, repeated.known) == (1.0, 5_000, 5_000)

    unknown = align_copies("N" * 5_000, "N" * 5_000)
    assert unknown.method == "ungapped_fallback"
    assert (unknown.matched, unknown.known, unknown.ambiguous) == (0, 0, 5_000)


def test_measure_ir_rejects_intervals_that_are_not_a_pair() -> None:
    from plastid.core.model import Interval

    with pytest.raises(ValueError, match="disjoint"):
        measure_ir("ACGT" * 10, Interval(0, 10), Interval(5, 15), source="manual")
