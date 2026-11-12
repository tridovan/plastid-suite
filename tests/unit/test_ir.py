"""IR detection and the mandatory gate — the highest-risk component.

The measured failure modes this file guards (each one stays forever):
  * naive detection was silently wrong on 2/9 published plastomes;
  * identity measured on the exact-match core hides divergence (22,998 bp core at 100 %
    while the true 26,162 bp copies are only 92.62 % identical);
  * an off-by-one pairing constant collapses identity to ~29 % — it must fail, loudly.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import replace

import pytest

from plastid.core import ir
from plastid.core.ir import (
    IRGateError,
    align_copies,
    detect_ir,
    enforce_ir_gate,
    identity_percent,
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
    """The sequence path must not manufacture a passing IR out of a divergent pair.

    Re-measured 28/10/2026 with the overlap-aware chain (K-07): the two records now behave
    differently, and the difference is the record, not the detector.

    * ``NC_063921`` diverges too far inside the pair for the chain to cross it, so the longest pair
      it can build is the perfect 16,197 bp reverse-complement repeat inside the annotated copies.
    * ``NC_063920`` carries five internal gaps and ten mismatches across the whole annotated pair,
      so the chain now reaches the pair's own boundaries — 26,161 bp at 99.9427 % over 26,164
      columns, the same interval the record annotates, which clears the 99.9 % gate. The old
      23,446 bp value was the chain stopping at the first gap; the de-novo path is a cross-check,
      and here it corroborates the annotation instead of truncating it.

    Independent evidence for the 26,161 bp reading: re-measuring the unedited pair over the
    annotation's own interval gives the same 10 mismatches and 5 gaps (26,149 / 26,164 columns,
    0.999427).
    """
    short = ir_de_novo(edge_genomes["NC_063921.1"])
    assert short is not None
    assert short.length == 16_197
    assert (short.matched, short.mismatches, short.ambiguous_columns) == (16_196, 0, 1)
    assert short.identity == pytest.approx(16_196 / 16_197)
    assert short.length < 26_183

    full = ir_de_novo(edge_genomes["NC_063920.1"])
    assert full is not None
    assert full.length == 26_161
    assert full.exact_core == 22_998  # the exact core is still not what the identity measures
    assert (full.matched, full.mismatches, full.internal_gaps) == (26_149, 10, 5)
    # columns = matched + mismatched + gaps = 26,164, the identity denominator
    assert full.identity == pytest.approx(26_149 / 26_164)
    assert full.gate_passed  # five internal gaps and ten mismatches still clear the 99.9 % gate


def test_detect_ir_records_the_de_novo_disagreement(
    edge_genomes: dict[str, Genome],
) -> None:
    """The annotation wins, and a de-novo result that disagrees with it stays visible.

    ``NC_063921`` is the record where the two paths really disagree (16,197 bp inside a 26,183 bp
    annotated pair), so it is the one that can hold this contract. ``NC_063920`` is not: with the
    overlap-aware chain the sequence path reaches the same interval the annotation names, and its
    delta is 0 — dropping that record from the test is the honest consequence, not a weakened
    check.
    """
    ir = detect_ir(edge_genomes["NC_063921.1"])
    assert ir.source == "annotation"  # annotation wins ...
    assert ir.length == 26_183
    assert ir.denovo_length == 16_197
    assert ir.denovo_delta == 9_986  # ... but the disagreement is stored, never dropped

    agrees = detect_ir(edge_genomes["NC_063920.1"])
    assert agrees.source == "annotation" and agrees.length == 26_161
    assert agrees.denovo_length == 26_161 and agrees.denovo_delta == 0


def test_skipping_the_cross_check_changes_nothing_but_the_recorded_de_novo_length(
    edge_genomes: dict[str, Genome], golden_genome: Genome
) -> None:
    """``cross_check=False`` spares the scan beside an annotation that was read, and only that.

    The copies, their measurement and the gate's verdict are the annotation's either way, on a
    record where the two paths agree and on the one where they disagree by 9,986 bp; only the
    recorded ``denovo_length`` is left unmeasured.
    """
    for genome in (golden_genome, edge_genomes["NC_063920.1"], edge_genomes["NC_063921.1"]):
        checked = detect_ir(genome)
        spared = detect_ir(genome, cross_check=False)
        assert spared.denovo_length is None and spared.denovo_delta is None
        assert spared == replace(checked, denovo_length=None)
        assert spared.gate_passed == checked.gate_passed


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


def test_a_refused_identity_never_reads_as_the_threshold(golden_genome: Genome) -> None:
    """KX548222.1 measures 99.8964 %. To two places that printed 99.90 % beside its refusal."""
    assert identity_percent(0.9989644) == 99.89
    assert identity_percent(0.99801746) == 99.80
    assert identity_percent(0.999) == 99.90
    assert identity_percent(1.0) == 100.0
    below = replace(detect_ir(golden_genome), identity=0.9989644)
    with pytest.raises(IRGateError, match=r"99\.89 % identity"):
        enforce_ir_gate(below)


# ------------------------------------------------ the vectorized scan equals the word-by-word walk
#
# The two functions below are the scan as it was written word by word, kept here as the oracle:
# the vectorized scan in core/ir.py is a speed-up of exactly this code and must return the same
# chains, run for run, including the order of diagonals with equal hit counts.


def _reference_index(sequence: str, kmer: int) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for position in range(len(sequence) - kmer + 1):
        word = sequence[position : position + kmer]
        if ir._KNOWN.issuperset(word):
            index.setdefault(word, []).append(position)
    return {
        word: positions
        for word, positions in index.items()
        if len(positions) <= ir._MAX_KMER_OCCURRENCES
    }


def _reference_anchor_chain(a: str, b: str) -> list[tuple[int, int, int]]:
    kmer = ir.IR_ALIGN_KMER
    index = _reference_index(b, kmer)

    hits_by_diagonal: dict[int, list[tuple[int, int]]] = {}
    for i in range(len(a) - kmer + 1):
        word = a[i : i + kmer]
        if not ir._KNOWN.issuperset(word):
            continue
        for j in index.get(word, ()):
            diagonal = j - i
            if abs(diagonal) <= ir.IR_ALIGN_MAX_DIAGONAL:
                hits_by_diagonal.setdefault(diagonal, []).append((i, j))

    anchors: list[tuple[int, int, int]] = []
    for hits in hits_by_diagonal.values():
        hits.sort()
        start_i = start_j = previous_i = previous_j = None
        for i, j in hits:
            if start_i is None:
                start_i, start_j = i, j
            elif (i, j) != (previous_i + 1, previous_j + 1):
                length = previous_i - start_i + kmer
                if length >= kmer:
                    anchors.append((start_i, start_j, length))
                start_i, start_j = i, j
            previous_i, previous_j = i, j
        if start_i is not None:
            length = previous_i - start_i + kmer
            if length >= kmer:
                anchors.append((start_i, start_j, length))

    anchors.sort()
    if not anchors:
        return []
    n = len(anchors)
    best = [anchor[2] for anchor in anchors]
    previous: list[int] = [-1] * n
    for i, (ai, bi, li) in enumerate(anchors):
        for j in range(i):
            aj, bj, _ = anchors[j]
            if aj < ai and bj < bi and best[j] + li > best[i]:
                best[i] = best[j] + li
                previous[i] = j
    end = max(range(n), key=lambda idx: best[idx])
    chain: list[tuple[int, int, int]] = []
    while end != -1:
        chain.append(anchors[end])
        end = previous[end]
    chain.reverse()
    kept: list[tuple[int, int, int]] = []
    for anchor in chain:
        if kept:
            last = kept[-1]
            if (
                anchor[0] + anchor[2] <= last[0] + last[2]
                and anchor[1] + anchor[2] <= last[1] + last[2]
            ):
                continue
        kept.append(anchor)
    return kept


def _reference_mirror_chains(sequence: str, mirror: str) -> list[tuple[int, int, int]]:
    length = len(sequence)
    index = _reference_index(sequence, ir._KMER)

    diagonals: Counter[int] = Counter()
    hits_by_diagonal: dict[int, list[tuple[int, int]]] = {}
    for position in range(length - ir._KMER + 1):
        word = mirror[position : position + ir._KMER]
        if not ir._KNOWN.issuperset(word):
            continue
        for mate in index.get(word, ()):
            diagonal = mate - position
            hits_by_diagonal.setdefault(diagonal, []).append((mate, position))
            diagonals[diagonal] += 1
    if not diagonals:
        return []

    runs: list[tuple[int, int, int]] = []
    for diagonal, _ in diagonals.most_common(ir._DIAGONALS_CONSIDERED):
        hits = sorted(hits_by_diagonal[diagonal])
        start = previous = None
        for position, mate in hits:
            if start is None:
                start = previous = (position, mate)
            elif (position, mate) == (previous[0] + 1, previous[1] + 1):
                previous = (position, mate)
            else:
                run_length = previous[0] - start[0] + ir._KMER
                if run_length >= ir._KMER:
                    runs.append((start[0], start[1], run_length))
                start = previous = (position, mate)
        if start is not None:
            run_length = previous[0] - start[0] + ir._KMER
            if run_length >= ir._KMER:
                runs.append((start[0], start[1], run_length))

    runs.sort()
    if not runs:
        return []
    count = len(runs)
    score = [run[2] for run in runs]
    covered_a = [run[2] for run in runs]
    covered_b = [run[2] for run in runs]
    previous_run: list[int] = [-1] * count
    for i, (ai, bi, li) in enumerate(runs):
        for j in range(i):
            aj, bj, lj = runs[j]
            a_gap = ai - (aj + lj)
            b_gap = bi - (bj + lj)
            if a_gap <= -li or a_gap <= -lj or b_gap <= -li or b_gap <= -lj:
                continue
            if a_gap > ir.IR_ALIGN_GAP_BP or b_gap > ir.IR_ALIGN_GAP_BP:
                continue
            new_a = covered_a[j] + li - max(0, -a_gap)
            new_b = covered_b[j] + li - max(0, -b_gap)
            new_score = min(new_a, new_b)
            if new_score > score[i]:
                score[i] = new_score
                covered_a[i], covered_b[i] = new_a, new_b
                previous_run[i] = j
    end = max(range(count), key=lambda index: score[index])
    chain: list[tuple[int, int, int]] = []
    while end != -1:
        chain.append(runs[end])
        end = previous_run[end]
    chain.reverse()
    return chain


def _scan_shapes() -> list[tuple[str, str]]:
    """Sequences that exercise every branch of the scan: SNPs, indels, N, a word seen > 32 times."""
    generator = random.Random(20261112)

    def random_bases(length: int) -> str:
        return "".join(generator.choice("ACGT") for _ in range(length))

    def mutate(sequence: str, snps: int, indels: int) -> str:
        bases = list(sequence)
        for _ in range(snps):
            index = generator.randrange(len(bases))
            bases[index] = generator.choice("ACGT".replace(bases[index], ""))
        for _ in range(indels):
            index = generator.randrange(1, len(bases) - 1)
            if generator.random() < 0.5:
                del bases[index]
            else:
                bases.insert(index, generator.choice("ACGT"))
        return "".join(bases)

    repeat = random_bases(6_000)
    plastome = (
        random_bases(9_000)
        + repeat
        + random_bases(2_500)
        + reverse_complement(mutate(repeat, 6, 3))
    )
    with_unknowns = (
        plastome[:3_000] + "N" * 40 + plastome[3_040:12_000] + "nnRY" + plastome[12_004:]
    )
    homopolymer = random_bases(4_000) + "A" * 90 + repeat[:3_000] + "T" * 90 + random_bases(1_500)
    homopolymer += reverse_complement(repeat[:3_000])
    lowercase = plastome[:8_000] + plastome[8_000:9_500].lower() + plastome[9_500:]
    return [
        ("plastome-like", plastome),
        ("unknown bases", with_unknowns),
        ("homopolymers", homopolymer),
        ("lower-case run", lowercase),
        ("shorter than a word", "ACGTACGTAC"),
        ("exactly one word", random_bases(ir.IR_ALIGN_KMER)),
        ("empty", ""),
    ]


@pytest.mark.parametrize("shape", _scan_shapes(), ids=lambda shape: shape[0])
def test_the_vectorized_scan_returns_the_walk_s_chains(shape: tuple[str, str]) -> None:
    _name, sequence = shape
    mirror = reverse_complement(sequence)
    assert ir._mirror_chains(sequence, mirror) == _reference_mirror_chains(sequence, mirror)
    half = len(sequence) // 2
    first, second = sequence[:half], reverse_complement(sequence[half:])
    assert ir._anchor_chain(first, second) == _reference_anchor_chain(first, second)


def test_the_vectorized_scan_returns_the_walk_s_chains_on_the_golden(golden_genome: Genome) -> None:
    sequence = golden_genome.sequence
    mirror = reverse_complement(sequence)
    assert ir._mirror_chains(sequence, mirror) == _reference_mirror_chains(sequence, mirror)
    annotated = ir_from_annotation(golden_genome)
    assert annotated is not None
    copy_a = sequence[annotated.ira.start : annotated.ira.end]
    copy_b = reverse_complement(sequence[annotated.irb.start : annotated.irb.end])
    assert ir._anchor_chain(copy_a, copy_b) == _reference_anchor_chain(copy_a, copy_b)
