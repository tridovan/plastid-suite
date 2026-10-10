"""Canonical origin and the four-region split.

Decision 2026-10-02 (PROJECT.md §3, closing §10.6): per-genome tables keep **GenBank
coordinates**; cross-genome comparisons use the **normalized frame**; every table states
which frame it uses, and the rotation offset is recorded in the run info.

Normalized frame: starts at the first base of the LSC; region order is
LSC → IRa → SSC → IRb. ``rotation_offset`` is the GenBank-frame position that becomes
the new position 0.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from plastid.core.ir import IRResult
from plastid.core.model import Gene, Genome, Interval

__all__ = [
    "NormalizedFrame",
    "Region",
    "lsc_origin",
    "normalize",
    "regions_from_ir",
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
    normalized_ir = replace(
        ir,
        ira=first,
        irb=second,
        mirror_constant=first.start + second.end - 1,
        frame_offset=0,
    )
    shifted_repeat_regions: list[Interval] = []
    for region in genome.repeat_regions:
        shifted_repeat_regions.extend(_shift_interval(region, -origin, length))
    rotated = replace(rotated, repeat_regions=tuple(shifted_repeat_regions))
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
    genes = tuple(_shift_gene(gene, -shift, length) for gene in genome.genes)
    repeat_regions = tuple(
        shifted
        for region in genome.repeat_regions
        for shifted in _shift_interval(region, -shift, length)
    )
    return replace(genome, sequence=sequence, genes=genes, repeat_regions=repeat_regions)


def _shift_gene(gene: Gene, delta: int, length: int) -> Gene:
    parts = tuple(
        shifted for part in gene.parts for shifted in _shift_interval(part, delta, length)
    )
    return replace(gene, parts=parts)


def _shift_interval(interval: Interval, delta: int, length: int) -> tuple[Interval, ...]:
    """Shift an interval by ``delta`` (mod length); a result crossing the origin splits in two."""
    start = (interval.start + delta) % length
    end = start + interval.length
    if end <= length:
        return (Interval(start, end, interval.strand),)
    return (
        Interval(start, length, interval.strand),
        Interval(0, end - length, interval.strand),
    )
