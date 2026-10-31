"""Reading an assembly graph without trusting the assembler that wrote it.

An assembler's "complete" makes two claims: that a path closes, and that it is the plastome. Reads
mapped back test the first (``core/assembly.py``). They cannot test the second. Two arrangements of
the same pieces use the same links and so hold the same k-mers, and a read shorter than a repeat
fits both. Only the graph says how many arrangements there are.

Three readings, from the graph and the written paths alone:

* **Is each path a walk of its graph?** Every step must be a link, and the walk must spell the path
  as a circle. A directory where it does not is not one assembly, and is refused.
* **Is the graph closed?** One component, no free end, every segment and every link used.
* **How many structures close through it?** Every circular walk that uses each segment as often as
  the path does is enumerated here, whatever the assembler wrote. One structure, or two a single
  inversion apart, which is what recombination between inverted copies gives, is resolved. More is
  a choice no short read can make, and the assembler's first path is not an answer to it: on
  *Picea abies* the published order is the third of the 15 written.

Depth does not decide anything here. It is reported per segment, because on 3 of 47 closed published
runs a segment's depth rounds to another copy number than the path gives it.

A step is ``(segment, strand)``, strand ``+`` or ``-``. A link joins two steps and overlaps them by
a number of bases, and the same link read from the other strand joins the two steps flipped and in
reverse order. A walk is circular: its last step links to its first.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from plastid.core.assembly import (
    STATUS_COMPLETE,
    AssembledPath,
    AssemblyError,
    same_circular,
)
from plastid.core.sequence import reverse_complement

__all__ = [
    "GRAPH_NOT_MEASURED",
    "GRAPH_RESOLVED",
    "GRAPH_UNRESOLVED",
    "AssemblyGraph",
    "GraphReading",
    "Segment",
    "SegmentReading",
    "are_isomers",
    "closed_walks",
    "graph_verdict",
    "longest_repeat",
    "parse_gfa",
    "parse_walk",
    "read_selected_graph",
    "spell",
]

GRAPH_RESOLVED = "resolved"
GRAPH_UNRESOLVED = "unresolved"
GRAPH_NOT_MEASURED = "not measured"

Step = tuple[str, str]
Walk = tuple[Step, ...]
Link = tuple[Step, Step]

_FLIP = {"+": "-", "-": "+"}
_CIRCULAR = "(circular)"
_OVERLAP = re.compile(r"(\d+)M")
# Search steps before the count of structures is given up. A graph with this many structures is not
# resolved whatever the exact number; the largest of the 47 published closed graphs takes far fewer.
_WALK_BUDGET = 2_000_000


@dataclass(frozen=True, slots=True)
class Segment:
    """One node of the graph. ``depth`` is the assembler's own, in k-mers per base, when stated."""

    name: str
    sequence: str
    depth: float | None


@dataclass(frozen=True, slots=True)
class AssemblyGraph:
    """Segments, and for each step the steps that follow it with their overlap in bp.

    A link is held from both ends, so walking never has to ask which strand it was written on.
    """

    segments: Mapping[str, Segment]
    following: Mapping[Step, Mapping[Step, int]]

    def bases(self, step: Step) -> str:
        """The step's sequence on its strand."""
        sequence = self.segments[step[0]].sequence
        return sequence if step[1] == "+" else reverse_complement(sequence)

    @property
    def links(self) -> frozenset[Link]:
        """Every link once, whichever strand it was written on."""
        return frozenset(
            _link(step, other) for step, targets in self.following.items() for other in targets
        )

    @property
    def free_ends(self) -> int:
        """Segment ends with no link: a closed molecule has none, in bases, of either end."""
        return sum(
            1 for name in self.segments for strand in "+-" if not self.following.get((name, strand))
        )

    @property
    def components(self) -> int:
        """Pieces the graph falls into when strands are ignored."""
        parent = {name: name for name in self.segments}

        def root(name: str) -> str:
            while parent[name] != name:
                parent[name] = parent[parent[name]]
                name = parent[name]
            return name

        for (name, _), targets in self.following.items():
            for other, _ in targets:
                parent[root(name)] = root(other)
        return len({root(name) for name in self.segments})


@dataclass(frozen=True, slots=True)
class SegmentReading:
    """One segment against the first path: how often the path uses it, and its depth against that.

    ``depth_ratio`` is the segment's depth over one copy's depth: the path's own mean depth per
    copy on a closed graph, the length-weighted median on an open one, which has no path to count.
    """

    name: str
    length: int
    depth: float | None
    copies: int | None
    depth_ratio: float | None


@dataclass(frozen=True, slots=True)
class GraphReading:
    """What one selected graph says about the paths written from it."""

    name: str
    call: str
    reasons: tuple[str, ...]
    segments: int
    links: int
    components: int
    free_ends: int
    paths_written: int
    structures: int | None  # circular walks that close; None when open or past the budget
    isomers: bool | None  # two structures one inversion apart; None when not two
    repeat: tuple[str, int] | None  # the first path's longest repeat: its kind, and its bp
    rows: tuple[SegmentReading, ...]


def _flip(step: Step) -> Step:
    return step[0], _FLIP[step[1]]


def _link(first: Step, second: Step) -> Link:
    return min((first, second), (_flip(second), _flip(first)))


def parse_gfa(lines: Iterable[str]) -> AssemblyGraph:
    """The segments and links of a GFA 1 text. Depth is ``RC`` over the segment's length."""
    segments: dict[str, Segment] = {}
    links: list[tuple[Step, Step, int]] = []
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "S":
            if len(fields) < 3 or fields[2] == "*":
                raise AssemblyError(f"a GFA segment with no sequence: {line[:80]!r}")
            tags = {tag[:2]: tag[5:] for tag in fields[3:]}
            count = tags.get("RC")
            sequence = fields[2].upper()
            depth = int(count) / len(sequence) if count is not None else None
            segments[fields[1]] = Segment(fields[1], sequence, depth)
        elif fields[0] == "L":
            if len(fields) < 6 or fields[2] not in _FLIP or fields[4] not in _FLIP:
                raise AssemblyError(f"not a GFA link: {line[:80]!r}")
            match = _OVERLAP.fullmatch(fields[5])
            if match is None:
                raise AssemblyError(f"a GFA link whose overlap is not plain: {line[:80]!r}")
            links.append(((fields[1], fields[2]), (fields[3], fields[4]), int(match.group(1))))
    following: dict[Step, dict[Step, int]] = {}
    for first, second, overlap in links:
        for name in (first[0], second[0]):
            if name not in segments:
                raise AssemblyError(f"a GFA link names a segment it does not hold: {name!r}")
        following.setdefault(first, {})[second] = overlap
        following.setdefault(_flip(second), {})[_flip(first)] = overlap
    return AssemblyGraph(segments=segments, following=following)


def parse_walk(title: str) -> Walk | None:
    """The circular walk a closed path's title states, or None when the title states none.

    GetOrganelle names a closed path ``a+,b-,c+(circular)``.
    """
    if not title.endswith(_CIRCULAR):
        return None
    steps: list[Step] = []
    for item in title.removesuffix(_CIRCULAR).split(","):
        if len(item) < 2 or item[-1] not in _FLIP:
            return None
        steps.append((item[:-1], item[-1]))
    return tuple(steps)


def spell(graph: AssemblyGraph, walk: Walk) -> str:
    """The sequence a circular walk spells; a step that is not a link is refused."""
    pieces: list[str] = []
    for step, after in zip(walk, (*walk[1:], walk[0]), strict=True):
        if step[0] not in graph.segments:
            raise AssemblyError(f"the walk names segment {step[0]!r}, which the graph lacks")
        overlap = graph.following.get(step, {}).get(after)
        if overlap is None:
            raise AssemblyError(f"the walk steps {step} to {after}, which is not a link")
        bases = graph.bases(step)
        pieces.append(bases[: len(bases) - overlap])
    return "".join(pieces)


def _canonical(walk: Walk) -> Walk:
    """One name for a circle: the smallest of its rotations, on either strand."""
    backward = tuple(_flip(step) for step in reversed(walk))
    return min(
        sequence[start:] + sequence[:start]
        for sequence in (walk, backward)
        for start in range(len(walk))
    )


def closed_walks(
    graph: AssemblyGraph,
    copies: Mapping[str, int],
    *,
    budget: int = _WALK_BUDGET,
) -> frozenset[Walk] | None:
    """Every circular walk using each segment exactly ``copies`` times; None past ``budget`` steps.

    Any such circle can be rotated to start at one chosen segment and read on the strand that has
    it ``+``, so one start covers them all. The fewest copies give the narrowest search.
    """
    start = min(copies, key=lambda name: (copies[name], name))
    size = sum(copies.values())
    left = dict(copies)
    left[start] -= 1
    found: set[Walk] = set()
    walk: list[Step] = [(start, "+")]
    searched = 0

    def descend() -> bool:
        nonlocal searched
        searched += 1
        if searched > budget:
            return False
        last = walk[-1]
        if len(walk) == size:
            if walk[0] in graph.following.get(last, {}):
                found.add(_canonical(tuple(walk)))
            return True
        for after in graph.following.get(last, {}):
            if left.get(after[0], 0) > 0:
                left[after[0]] -= 1
                walk.append(after)
                finished = descend()
                walk.pop()
                left[after[0]] += 1
                if not finished:
                    return False
        return True

    return frozenset(found) if descend() else None


def _used_links(walk: Walk) -> frozenset[Link]:
    return frozenset(
        _link(step, after) for step, after in zip(walk, (*walk[1:], walk[0]), strict=True)
    )


def are_isomers(first: Walk, second: Walk) -> bool:
    """Whether two circles are the flip-flop pair: one inversion apart, over the same links.

    Recombination between two inverted copies of a repeat turns over what lies between them, so
    the molecule exists in both forms and neither is an error. The two then cross the same joins,
    which is what lets the reads of one be read for both.
    """
    if _used_links(first) != _used_links(second):
        return False
    target = _canonical(second)
    size = len(first)
    for start in range(size):
        turned = first[start:] + first[:start]
        for cut in range(1, size):
            reversed_part = tuple(_flip(step) for step in reversed(turned[:cut]))
            if _canonical(reversed_part + turned[cut:]) == target:
                return True
    return False


def _spelled_length(graph: AssemblyGraph, walk: Sequence[Step]) -> int:
    """Bases a linear stretch of a walk spells: its segments, less what they share."""
    total = sum(len(graph.segments[name].sequence) for name, _ in walk)
    return total - sum(
        graph.following[step][after] for step, after in zip(walk, walk[1:], strict=False)
    )


def longest_repeat(graph: AssemblyGraph, walk: Walk) -> tuple[str, int] | None:
    """The longest stretch a circular walk holds twice: ``("inverted" | "direct", bp)``.

    This is a reading of the repeat independent of the sequence: where the graph puts its ends is
    where the two copies stop being one segment.
    """
    size = len(walk)
    best: tuple[str, int] | None = None
    for first in range(size):
        for second in range(size):
            if first == second:
                continue
            gap = (second - first) % size
            # Inverted: read forward from one copy and backward from the other, until they meet.
            stretch: list[Step] = []
            while 2 * len(stretch) <= gap - 1:
                step = walk[(first + len(stretch)) % size]
                if step != _flip(walk[(second - len(stretch)) % size]):
                    break
                stretch.append(step)
            if stretch:
                length = _spelled_length(graph, stretch)
                if best is None or length > best[1]:
                    best = ("inverted", length)
            # Direct: both copies forward, up to where the first would run into the second.
            stretch = []
            while len(stretch) < gap:
                step = walk[(first + len(stretch)) % size]
                if step != walk[(second + len(stretch)) % size]:
                    break
                stretch.append(step)
            if stretch:
                length = _spelled_length(graph, stretch)
                if best is None or length > best[1]:
                    best = ("direct", length)
    return best


def _median_depth(graph: AssemblyGraph) -> float | None:
    """The length-weighted median depth: a whole segment's vote is its length."""
    known = sorted(
        (segment.depth, len(segment.sequence))
        for segment in graph.segments.values()
        if segment.depth is not None
    )
    half = sum(width for _, width in known) / 2
    running = 0
    for depth, width in known:
        running += width
        if running >= half:
            return depth
    return None


def _rows(graph: AssemblyGraph, copies: Mapping[str, int]) -> tuple[SegmentReading, ...]:
    if copies:
        mass = sum(
            (graph.segments[name].depth or 0.0) * len(graph.segments[name].sequence)
            for name in copies
        )
        width = sum(count * len(graph.segments[name].sequence) for name, count in copies.items())
        one_copy = mass / width if mass else None
    else:
        one_copy = _median_depth(graph)
    return tuple(
        SegmentReading(
            name=segment.name,
            length=len(segment.sequence),
            depth=segment.depth,
            copies=copies.get(segment.name) if copies else None,
            depth_ratio=(
                segment.depth / one_copy if segment.depth is not None and one_copy else None
            ),
        )
        for segment in sorted(graph.segments.values(), key=lambda item: -len(item.sequence))
    )


def _walks(graph: AssemblyGraph, paths: Sequence[AssembledPath]) -> dict[str, Walk] | None:
    """Each path's walk, checked against the graph; None when a path does not state one."""
    walks: dict[str, Walk] = {}
    for path in paths:
        walk = parse_walk(path.title)
        if walk is None:
            return None
        spelled = spell(graph, walk)
        if len(path.contigs) != 1 or not same_circular(spelled, path.contigs[0]):
            raise AssemblyError(
                f"path {path.name} is not the sequence its walk spells through this graph: "
                "the path file and the graph are not from one assembly"
            )
        walks[path.name] = walk
    return walks


def read_selected_graph(
    name: str,
    graph: AssemblyGraph,
    paths: Sequence[AssembledPath],
    *,
    status: str,
) -> GraphReading:
    """Read one selected graph and the paths written from it.

    A path that its graph contradicts is refused, because that directory is not one assembly. Any
    other way the graph fails to settle the molecule is a named reason, and the call is unresolved.
    """
    if status == STATUS_COMPLETE and not paths:
        raise AssemblyError("no path was written from this graph, and the result is closed")
    reasons: list[str] = []
    components, free_ends = graph.components, graph.free_ends
    if components != 1 or free_ends:
        reasons.append(f"{components} component(s) and {free_ends} free end(s): it does not close")
    walks = _walks(graph, paths) if status == STATUS_COMPLETE else {}
    structures: int | None = None
    isomers: bool | None = None
    repeat: tuple[str, int] | None = None
    copies: Counter[str] = Counter()
    if status == STATUS_COMPLETE and walks is None:
        return GraphReading(
            name=name,
            call=GRAPH_NOT_MEASURED,
            reasons=("a path does not state its walk through the graph",),
            segments=len(graph.segments),
            links=len(graph.links),
            components=components,
            free_ends=free_ends,
            paths_written=len(paths),
            structures=None,
            isomers=None,
            repeat=None,
            rows=_rows(graph, {}),
        )
    if walks:
        ordered = list(walks.values())
        copies = Counter(step[0] for step in ordered[0])
        repeat = longest_repeat(graph, ordered[0])
        unused = set(graph.segments) - {step[0] for walk in ordered for step in walk}
        if unused:
            reasons.append(f"{len(unused)} segment(s) no path uses")
        used = frozenset().union(*(_used_links(walk) for walk in ordered))
        if used != graph.links:
            reasons.append(f"{len(graph.links - used)} link(s) no path uses")
        written = {_canonical(walk) for walk in ordered}
        allowed: set[Walk] = set()
        counted = True
        for vector in {
            tuple(sorted(Counter(step[0] for step in walk).items())) for walk in ordered
        }:
            found = closed_walks(graph, dict(vector))
            if found is None:
                counted = False
                break
            allowed |= found
        if not counted:
            reasons.append("more structures close through it than can be counted here")
        else:
            structures = len(allowed)
            if written != allowed:
                reasons.append(
                    f"{structures} structure(s) close through it and "
                    f"{len(written & allowed)} of them were written"
                )
            if len(allowed) == 2:
                first, second = sorted(allowed)
                isomers = are_isomers(first, second)
                if not isomers:
                    reasons.append("two structures close through it, and they are not isomers")
            elif len(allowed) > 2:
                reasons.append(
                    f"{structures} structures close through it, and reads shorter than its "
                    "repeats cannot choose among them"
                )
    elif status != STATUS_COMPLETE and not reasons:
        reasons.append(f"the assembler's result is {status}")
    return GraphReading(
        name=name,
        call=GRAPH_UNRESOLVED if reasons else GRAPH_RESOLVED,
        reasons=tuple(reasons),
        segments=len(graph.segments),
        links=len(graph.links),
        components=components,
        free_ends=free_ends,
        paths_written=len(paths),
        structures=structures,
        isomers=isomers,
        repeat=repeat,
        rows=_rows(graph, copies),
    )


def graph_verdict(readings: Sequence[GraphReading]) -> tuple[str, tuple[str, ...]]:
    """The call for a whole assembly: one graph, resolved, or the reasons it is not."""
    if not readings:
        return GRAPH_NOT_MEASURED, ("no selected graph to read",)
    if any(reading.call == GRAPH_NOT_MEASURED for reading in readings):
        reasons = tuple(reason for reading in readings for reason in reading.reasons)
        return GRAPH_NOT_MEASURED, reasons
    reasons = tuple(reason for reading in readings for reason in reading.reasons)
    if len(readings) > 1:
        reasons = (f"{len(readings)} graphs were selected for one sample", *reasons)
    return (GRAPH_UNRESOLVED if reasons else GRAPH_RESOLVED), reasons
