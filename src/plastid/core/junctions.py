"""The four regions a plastome is cut into, and the four junctions between them.

**Frame.** Everything here is the record's own GenBank frame, 1-based inclusive — the frame a
published junction table prints and the frame CPJSdraw draws in. That is deliberately *not* the
normalized frame ``structure`` exports (LSC → IRa → SSC → IRb, origin moved to the LSC/IRa
boundary): a junction is a statement about the deposited coordinates, and rotating the record would
move every one of them.

**Why this is not the detector.** :mod:`plastid.core.ir` answers "are these two copies the same
inverted repeat, well enough to analyse positions?" — and the IR gate is the guard that acts on the
answer. This module answers a different question: "where does the record *say* the four regions
are?" It reads the record's own ``repeat_region`` annotation as written. The two answers can
disagree, and when they do the disagreement is reported rather than resolved here (H-25).

**A copy that wraps the origin.** A ``repeat_region`` spanning position 1 is written in GenBank as a
``join()``, and the reader hands its pieces over one by one — so one copy can arrive as two
intervals: the main run, plus a piece at the start of the record. Both belong to that copy. Dropping
the small piece shortens the copy, and a copy one base short is a copy whose bases no longer pair
with the other's: measured 17/10/2026 on Chasmanthe NC_072562, dropping a 1 bp piece turned an
identity of 1.0000 into 0.2924. So the piece is kept, and the region is reported as wrapping.
"""

from __future__ import annotations

from dataclasses import dataclass

from plastid.core.model import Genome

__all__ = [
    "JunctionError",
    "JunctionRegions",
    "RegionSpan",
    "junction_regions",
]

# The two copies must be at least this long to be the inverted repeat rather than a stray
# ``repeat_region`` — the same floor the detector uses, so the two readers agree on what a copy is.
_MIN_COPY_BP = 5_000
_ORDER = ("LSC", "IRb", "SSC", "IRa")


class JunctionError(Exception):
    """The record's own annotation does not cut into four regions — say which way it fails."""


@dataclass(frozen=True, slots=True)
class RegionSpan:
    """One region of the record, in the record's own frame (1-based inclusive).

    ``wrap_start_1b``/``wrap_end_1b`` are set when the region runs across the origin: the run
    starts after ``start_1b``-``end_1b`` and continues at position 1. ``length`` counts both parts.
    """

    name: str
    start_1b: int
    end_1b: int
    length: int
    wrap_end_1b: int | None = None

    @property
    def wraps(self) -> bool:
        return self.wrap_end_1b is not None

    @property
    def text(self) -> str:
        """CPJSdraw's spelling: ``start-end``, plus ``;1-x`` for the part across the origin."""
        if self.wrap_end_1b is None:
            return f"{self.start_1b}-{self.end_1b}"
        return f"{self.start_1b}-{self.end_1b};1-{self.wrap_end_1b}"


@dataclass(frozen=True, slots=True)
class JunctionRegions:
    """The four regions and the four junction positions, in the record's own frame.

    ``jlb``/``jsb``/``jsa``/``jla`` are the base the junction *follows*, so the first base of the
    region the name leads into: JLB (LSC→IRb) is IRb's first base, JSB (IRb→SSC) the base after
    IRb's last, JSA (SSC→IRa) the base before IRa's first, JLA (IRa→LSC) IRa's first.
    """

    lsc: RegionSpan
    irb: RegionSpan
    ssc: RegionSpan
    ira: RegionSpan
    jlb: int
    jsb: int
    jsa: int
    jla: int
    genome_bp: int

    def region(self, name: str) -> RegionSpan:
        if name not in _ORDER:
            raise KeyError(f"no region {name!r} (have {_ORDER})")
        return {"LSC": self.lsc, "IRb": self.irb, "SSC": self.ssc, "IRa": self.ira}[name]

    @property
    def regions(self) -> tuple[RegionSpan, ...]:
        return (self.lsc, self.irb, self.ssc, self.ira)

    @property
    def total_bp(self) -> int:
        return sum(region.length for region in self.regions)

    @property
    def tiles_genome(self) -> bool:
        """Whether the four regions account for every base.

        Two disjoint copies always tile the circle, so this is a property a caller *checks and
        reports*, not a branch inside :func:`junction_regions` — the derivation cannot leave a gap,
        and a guard that can never fire is a guard nobody has tested.
        """
        return self.total_bp == self.genome_bp

    @property
    def config_line(self) -> str:
        """The one line CPJSdraw reads: ``LSC:…,IRb:…,SSC:…,IRa:…``, in GenBank coordinates."""
        return ",".join(f"{region.name}:{region.text}" for region in self.regions)


def junction_regions(genome: Genome) -> JunctionRegions:
    """Cut one record into LSC/IRb/SSC/IRa from its own ``repeat_region`` annotation.

    Each ``repeat_region`` feature is one copy, read whole: a copy the record writes as
    ``join(m..L,1..k)`` is one circular arc, not two candidates (M11.3, K-03). The two copies are
    told apart by geometry, not by qualifier text: the larger of the two gaps between them is the
    LSC, and the copy that follows the LSC is IRb.
    """
    length = len(genome)
    if length == 0:
        raise JunctionError(f"{genome.accession or genome.name}: empty sequence")

    copies = [region for region in genome.repeat_regions if region.length >= _MIN_COPY_BP]
    if len(copies) != 2:
        raise JunctionError(
            f"{genome.accession or genome.name}: expected two repeat_region copies of at least "
            f"{_MIN_COPY_BP} bp, found {len(copies)} "
            f"({[(r.start + 1, r.end) for r in genome.repeat_regions]}) — the junction figure is "
            f"drawn from the record's own annotation and cannot be drawn without it"
        )
    arcs: list[tuple[int, int]] = []
    for region in copies:
        try:
            arcs.append(region.arc(length))
        except ValueError as error:
            raise JunctionError(
                f"{genome.accession or genome.name}: repeat_region "
                f"{region.start + 1}..{region.end} is not one circular copy — {error}"
            ) from None
    first, second = sorted(arcs, key=lambda arc: arc[0])
    if _arcs_overlap(first, second, length):
        raise JunctionError(
            f"{genome.accession or genome.name}: the two IR copies overlap "
            f"({first[0] + 1}..{first[0] + first[1]} and {second[0] + 1}..{second[0] + second[1]})"
        )

    # The larger gap is the LSC; the copy after it is IRb.
    between = (second[0] - (first[0] + first[1])) % length
    across = length - first[1] - second[1] - between
    if between >= across:
        irb_arc, ira_arc = second, first
    else:
        irb_arc, ira_arc = first, second

    ira_span = _arc_span("IRa", ira_arc, length)
    irb_span = _arc_span("IRb", irb_arc, length)

    # LSC runs from the end of IRa to the start of IRb; SSC from the end of IRb to the start of
    # IRa. Either copy may wrap the origin — which of them does depends on where the record's
    # origin happens to fall — so "the base after this region" is asked of the region, not of its
    # main run, and the length is computed circularly.
    lsc_span = _span_between("LSC", _after(ira_span, length), irb_span.start_1b - 1, length)
    ssc_span = _span_between("SSC", _after(irb_span, length), ira_span.start_1b - 1, length)

    return JunctionRegions(
        lsc=lsc_span,
        irb=irb_span,
        ssc=ssc_span,
        ira=ira_span,
        jlb=irb_span.start_1b,
        jsb=_after(irb_span, length),
        jsa=ira_span.start_1b - 1,
        jla=ira_span.start_1b,
        genome_bp=length,
    )


def _arcs_overlap(first: tuple[int, int], second: tuple[int, int], length: int) -> bool:
    (first_start, first_length), (second_start, second_length) = first, second
    return (second_start - first_start) % length < first_length or (
        first_start - second_start
    ) % length < second_length


def _after(region: RegionSpan, length: int) -> int:
    """The base that follows a region — position 1 when the region ends at the record's end."""
    last = region.wrap_end_1b if region.wrap_end_1b is not None else region.end_1b
    return 1 if last >= length else last + 1


def _arc_span(name: str, arc: tuple[int, int], length: int) -> RegionSpan:
    """One whole copy as a region span; a copy across the origin keeps both pieces."""
    start, region_length = arc
    if start + region_length <= length:
        return RegionSpan(
            name=name, start_1b=start + 1, end_1b=start + region_length, length=region_length
        )
    return RegionSpan(
        name=name,
        start_1b=start + 1,
        end_1b=length,
        length=region_length,
        wrap_end_1b=start + region_length - length,
    )


def _span_between(name: str, start_1b: int, end_1b: int, length: int) -> RegionSpan:
    """A region given by its first and last base; ``end < start`` runs across the origin.

    ``end_1b == 0`` means the region runs to the record's end and the next base is 1, so there is
    no wrapped tail to print: a copy at position 1 is not a piece wrapping onto itself.
    """
    if end_1b == 0:
        return RegionSpan(name=name, start_1b=start_1b, end_1b=length, length=length - start_1b + 1)
    if end_1b < start_1b:
        return RegionSpan(
            name=name,
            start_1b=start_1b,
            end_1b=length,
            length=length - start_1b + 1 + end_1b,
            wrap_end_1b=end_1b,
        )
    return RegionSpan(name=name, start_1b=start_1b, end_1b=end_1b, length=end_1b - start_1b + 1)
