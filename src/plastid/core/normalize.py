"""Canonical origin and the four-region split.

Decision 2026-10-02: per-genome tables keep **GenBank
coordinates**; cross-genome comparisons use the **normalized frame**; every table states
which frame it uses, and the rotation offset is recorded in the run info.

Normalized frame: starts at the first base of the LSC; region order is
LSC → IRa → SSC → IRb. ``rotation_offset`` is the GenBank-frame position that becomes
the new position 0.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

from plastid.core.ir import IRResult, mirror_constant_for
from plastid.core.model import Gene, Genome, Interval, RepeatRegion

__all__ = [
    "NormalizedFrame",
    "Region",
    "lsc_origin",
    "normalize",
    "regions_from_ir",
    "rotate_gene",
    "rotate_genes",
    "rotate_genome",
]


@dataclass(frozen=True, slots=True)
class Region:
    """A named region in the frame the caller states; never wraps after normalization."""

    name: str  # "LSC" | "IRa" | "SSC" | "IRb"
    start: int  # 0-based half-open
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class NormalizedFrame:
    """The rotated genome plus everything needed to map back to GenBank coordinates."""

    genome: Genome
    ir: IRResult
    regions: tuple[Region, ...]
    rotation_offset: int  # 0-based GenBank position that is now position 0

    def region_at(self, position: int) -> Region | None:
        """The region holding a position given in the frame this frame was **built from**.

        Region names live in the normalized frame; a record's features are in the frame it was
        read in. Comparing the two directly is right only while nothing moved — the golden's
        offset is 0, which is how such a comparison survives unnoticed until a record that starts
        anywhere else (tests/unit/test_spans.py rotates the golden and checks the labels hold).
        """
        offset = (position - self.rotation_offset) % len(self.genome)
        for region in self.regions:
            if region.start <= offset < region.end:
                return region
        return None


def lsc_origin(genome_length: int, ira: Interval, irb: Interval) -> int:
    """0-based position of the LSC's first base: the start of the longer single-copy gap."""
    first, second = sorted((ira, irb), key=lambda interval: interval.start)
    gap_between = second.start - first.end
    gap_wrapping = (genome_length - second.end) + first.start
    return second.end if gap_wrapping > gap_between else first.end


def regions_from_ir(genome: Genome, ir: IRResult) -> tuple[Region, ...]:
    """Split a **normalized** genome into LSC → IRa → SSC → IRb."""
    first, second = sorted((ir.ira, ir.irb), key=lambda interval: interval.start)
    if first.start <= 0 or second.end != len(genome):
        raise ValueError("genome is not in the normalized frame — call normalize() first")
    return (
        Region("LSC", 0, first.start),
        Region("IRa", first.start, first.end),
        Region("SSC", first.end, second.start),
        Region("IRb", second.start, second.end),
    )


def normalize(genome: Genome, ir: IRResult) -> NormalizedFrame:
    """Rotate the genome so it starts at the LSC's first base; genes and IR move with it."""
    length = len(genome)
    # The IR may be expressed in a rotated frame (detector retry); the origin it gives is
    # in *that* frame, while genome features are in the genome frame. Keep both.
    origin_in_ir_frame = lsc_origin(length, ir.ira, ir.irb)
    origin = (origin_in_ir_frame + ir.frame_offset) % length
    rotated = rotate_genome(genome, origin)
    shifted_ira = _shift_interval(ir.ira, -origin_in_ir_frame, length)
    shifted_irb = _shift_interval(ir.irb, -origin_in_ir_frame, length)
    if len(shifted_ira) != 1 or len(shifted_irb) != 1:
        raise ValueError("rotation would split an IR copy — refusing rather than corrupting it")
    # After normalization: IRa = the copy adjacent to the LSC's 3' end, IRb = the copy
    # ending at the origin — reorder so the names match the region labels.
    first, second = sorted((shifted_ira[0], shifted_irb[0]), key=lambda iv: iv.start)
    # `mirror_constant` is provenance for mapping a position between copies. That mapping is
    # defined only for equal copies; unequal copies need the alignment-based mapping of M11.11.
    # Normalization itself (origin, regions, sequence) does not need the constant, so it records 0
    # rather than refusing a record the audit must still be able to describe.
    constant = (
        mirror_constant_for(first, second, "none")
        if first.length == second.length and ir.excess_at == "none"
        else 0
    )
    normalized_ir = replace(
        ir,
        ira=first,
        irb=second,
        mirror_constant=constant,
        frame_offset=0,
        # The alignment's cursor tables describe the frame it was measured in; this frame moved
        # both copies. A consumer that needs the mapping must rebuild it from the frame it uses,
        # not read a stale one (M11.11).
        mapping=None,
    )
    shifted_repeat_regions = tuple(
        _shift_region(region, -origin, length) for region in genome.repeat_regions
    )
    rotated = replace(rotated, repeat_regions=shifted_repeat_regions)
    return NormalizedFrame(
        genome=rotated,
        ir=normalized_ir,
        regions=regions_from_ir(rotated, normalized_ir),
        rotation_offset=origin,
    )


def rotate_genome(genome: Genome, shift: int) -> Genome:
    """Return a copy of ``genome`` whose sequence and features start ``shift`` bases later."""
    length = len(genome)
    shift %= length
    if shift == 0:
        return genome
    sequence = genome.sequence[shift:] + genome.sequence[:shift]
    genes = rotate_genes(genome.genes, shift, length)
    repeat_regions = tuple(
        _shift_region(region, -shift, length) for region in genome.repeat_regions
    )
    return replace(genome, sequence=sequence, genes=genes, repeat_regions=repeat_regions)


def rotate_genes(genes: Sequence[Gene], shift: int, length: int) -> tuple[Gene, ...]:
    """Features expressed in the frame ``rotate_genome`` produces — a part may split at the origin.

    Public because a module can hold genes without holding the genome that carries them: the
    annotator receives its calls as :class:`Gene` objects and must move them into the IR's frame
    before mapping a copy (M11.11), then move the result back. One implementation, so the two
    directions cannot disagree about transcript order or partial flags.
    """
    shift %= length
    if shift == 0:
        return tuple(genes)
    return tuple(rotate_gene(gene, shift, length) for gene in genes)


def rotate_gene(gene: Gene, shift: int, length: int) -> Gene:
    """One feature shifted into the frame ``rotate_genome`` produces."""
    shift %= length
    if shift == 0:
        return gene
    parts = tuple(
        shifted for part in gene.parts for shifted in _shift_interval(part, -shift, length)
    )
    return replace(gene, parts=parts)


def _shift_region(region: RepeatRegion, delta: int, length: int) -> RepeatRegion:
    """Shift one repeat feature, keeping its parts together as one copy.

    The parts are shifted independently (a part may itself cross the new origin and split), but
    they are collected back into one :class:`RepeatRegion`; flattening them would turn one copy
    into several candidates again (M11.3). Order is not semantic for a repeat_region — an
    intergenic region has no transcript — so the parts are left in the order they came out.
    """
    return RepeatRegion(
        tuple(shifted for part in region.parts for shifted in _shift_interval(part, delta, length))
    )


def _shift_interval(interval: Interval, delta: int, length: int) -> tuple[Interval, ...]:
    """Shift an interval by ``delta`` (mod length); a result crossing the origin splits in two.

    Three pieces of meaning travel with the interval, and each has bitten:

    * **Partial flags belong to the coordinates they describe**, not to the low/high end of the
      new linear record. They therefore stay on the original endpoints. The cut the rotation
      itself introduces is artificial and must not be marked partial — that would turn "the
      record's coordinate is incomplete" into a claim this toolkit made up (K-05).
    * **The fragments are returned in transcript order.** On the plus strand that is high then
      low (the transcript runs through the genomic end before continuing at position 1); on the
      minus strand it is low then high, because the transcript reads the plus strand backwards.
      Returning a fixed order would silently change a minus-strand CDS's translation across a
      rotation, which is the kind of quiet error the IR gate exists to prevent.
    * ``strand`` is kept on both fragments; a split does not make a feature mixed-strand.
    """
    start = (interval.start + delta) % length
    end = start + interval.length
    if end <= length:
        return (
            Interval(
                start,
                end,
                interval.strand,
                partial_start=interval.partial_start,
                partial_end=interval.partial_end,
            ),
        )
    high = Interval(
        start,
        length,
        interval.strand,
        partial_start=interval.partial_start,
        # the split at `length` is where the new linear record ends, not an incomplete coordinate
        partial_end=False,
    )
    low = Interval(
        0,
        end - length,
        interval.strand,
        # the join at position 1 is artificial for the same reason
        partial_start=False,
        partial_end=interval.partial_end,
    )
    return (high, low) if interval.strand > 0 else (low, high)
