"""Reading an assembly graph: each shape a real selected graph takes, built by hand.

The expected results here come from the construction of each graph, not from an assembler: a
quadripartite circle has two closed walks by arithmetic, and four pieces between direct repeats have
3! circular orders. The real graphs are in `tests/golden/test_assembly_graph_panel.py`.
"""

from __future__ import annotations

import random

import pytest

from plastid.core.assembly import (
    STATUS_COMPLETE,
    STATUS_FRAGMENTED,
    STATUS_GAPPED,
    AssembledPath,
    AssemblyError,
)
from plastid.core.assembly_graph import (
    GRAPH_NOT_MEASURED,
    GRAPH_RESOLVED,
    GRAPH_UNRESOLVED,
    AssemblyGraph,
    are_isomers,
    closed_walks,
    graph_verdict,
    longest_repeat,
    parse_gfa,
    parse_walk,
    read_selected_graph,
    spell,
)
from plastid.core.sequence import reverse_complement

_K = 5
_SOURCE = random.Random(20261031)


def _bases(length: int) -> str:
    return "".join(_SOURCE.choice("ACGT") for _ in range(length))


def _gfa(
    segments: dict[str, str],
    links: list[tuple[str, str, str, str]],
    *,
    depths: dict[str, float] | None = None,
) -> AssemblyGraph:
    lines = []
    for name, bases in segments.items():
        depth = (depths or {}).get(name, 100.0)
        lines.append(f"S\t{name}\t{bases}\tLN:i:{len(bases)}\tRC:i:{int(depth * len(bases))}\n")
    lines.extend(f"L\t{a}\t{sa}\t{b}\t{sb}\t{_K}M\n" for a, sa, b, sb in links)
    return parse_gfa(lines)


class Quadripartite:
    """LSC, IR, SSC, reversed IR: a segment holds the k bases it shares with the next."""

    def __init__(self) -> None:
        self.lsc, self.ir, self.ssc = _bases(300), _bases(120), _bases(60)
        mirror = reverse_complement(self.ir)
        self.sequence = self.lsc + self.ir + self.ssc + mirror
        self.isomer = self.lsc + self.ir + reverse_complement(self.ssc) + mirror
        self.graph = _gfa(
            {
                "L": mirror[-_K:] + self.lsc + self.ir[:_K],
                "I": self.ir,
                "S": self.ir[-_K:] + self.ssc + mirror[:_K],
            },
            [
                ("L", "+", "I", "+"),
                ("I", "+", "S", "+"),
                ("S", "+", "I", "-"),
                ("I", "-", "L", "+"),
            ],
            depths={"L": 80.0, "I": 160.0, "S": 80.0},
        )

    def paths(self) -> list[AssembledPath]:
        return [
            AssembledPath("p.1", (self.sequence,), True, "L+,I+,S+,I-(circular)"),
            AssembledPath("p.2", (self.isomer,), True, "L+,I+,S-,I-(circular)"),
        ]


class DirectRepeats:
    """``count`` pieces between direct copies of one repeat, as on *Selaginella* or a conifer."""

    def __init__(self, count: int) -> None:
        self.repeat = _bases(100)
        self.pieces = {f"p{index}": _bases(200) for index in range(count)}
        segments = {
            name: self.repeat[-_K:] + bases + self.repeat[:_K]
            for name, bases in self.pieces.items()
        }
        segments["r"] = self.repeat
        links = [(name, "+", "r", "+") for name in self.pieces]
        links += [("r", "+", name, "+") for name in self.pieces]
        self.graph = _gfa(segments, links)

    def path(self, order: tuple[str, ...], number: int = 1) -> AssembledPath:
        sequence = "".join(self.pieces[name] + self.repeat for name in order)
        title = ",".join(f"{name}+,r+" for name in order)
        return AssembledPath(f"p.{number}", (sequence,), True, f"{title}(circular)")


@pytest.fixture
def plastome() -> Quadripartite:
    return Quadripartite()


def test_a_gfa_holds_each_link_from_both_ends(plastome: Quadripartite) -> None:
    graph = plastome.graph
    assert set(graph.segments) == {"L", "I", "S"}
    assert graph.following[("L", "+")] == {("I", "+"): _K}
    # The same link read from the other strand, which no line of the file states.
    assert graph.following[("I", "-")][("L", "-")] == _K
    assert len(graph.links) == 4


def test_depth_is_the_read_count_tag_over_the_length(plastome: Quadripartite) -> None:
    assert plastome.graph.segments["I"].depth == pytest.approx(160.0)


def test_a_segment_without_the_tag_has_no_depth_rather_than_none_at_all() -> None:
    graph = parse_gfa(["S\tx\tACGT\n"])
    assert graph.segments["x"].depth is None


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("S\tx\t*\n", "no sequence"),
        ("L\tx\t+\ty\t-\t5M\n", "does not hold"),
        ("L\tx\t?\ty\t-\t5M\n", "not a GFA link"),
        ("L\tx\t+\tx\t+\t3M2I\n", "not plain"),
    ],
)
def test_a_line_that_cannot_be_read_is_refused(line: str, message: str) -> None:
    with pytest.raises(AssemblyError, match=message):
        parse_gfa(["S\tx\tACGT\n", line])


def test_a_walk_is_read_from_a_closed_paths_title() -> None:
    assert parse_walk("a+,b_c-,d+(circular)") == (("a", "+"), ("b_c", "-"), ("d", "+"))


@pytest.mark.parametrize("title", ["a+,b-", "scaffold_1--17+", "a,b(circular)", "(circular)"])
def test_a_title_that_states_no_closed_walk_is_none(title: str) -> None:
    assert parse_walk(title) is None


def test_a_walk_spells_the_path_as_a_circle(plastome: Quadripartite) -> None:
    from plastid.core.assembly import same_circular

    for path, sequence in zip(plastome.paths(), (plastome.sequence, plastome.isomer), strict=True):
        walk = parse_walk(path.title)
        assert walk is not None
        assert same_circular(spell(plastome.graph, walk), sequence)


def test_spelling_refuses_a_step_that_is_not_a_link(plastome: Quadripartite) -> None:
    with pytest.raises(AssemblyError, match="not a link"):
        spell(plastome.graph, (("L", "+"), ("S", "+")))


def test_spelling_refuses_a_segment_the_graph_lacks(plastome: Quadripartite) -> None:
    with pytest.raises(AssemblyError, match="lacks"):
        spell(plastome.graph, (("Z", "+"), ("L", "+")))


def test_a_closed_graph_has_one_component_and_no_free_end(plastome: Quadripartite) -> None:
    assert plastome.graph.components == 1
    assert plastome.graph.free_ends == 0


def test_a_broken_graph_counts_its_pieces_and_their_ends() -> None:
    graph = _gfa({"a": _bases(50), "b": _bases(50), "c": _bases(50)}, [("a", "+", "b", "+")])
    assert graph.components == 2
    assert graph.free_ends == 4


def test_one_segment_closing_on_itself_has_no_free_end() -> None:
    # Cytinus: a 19 kb plastome with no repeat is one segment linked to itself.
    bases = _bases(200)
    graph = _gfa({"x": bases + bases[:_K]}, [("x", "+", "x", "+")])
    assert graph.components == 1
    assert graph.free_ends == 0
    assert spell(graph, (("x", "+"),)) == bases
    assert closed_walks(graph, {"x": 1}) == {(("x", "+"),)}


def test_a_quadripartite_circle_has_two_structures(plastome: Quadripartite) -> None:
    assert len(closed_walks(plastome.graph, {"L": 1, "I": 2, "S": 1}) or ()) == 2


def test_two_direct_repeats_give_one_structure() -> None:
    shape = DirectRepeats(2)
    assert len(closed_walks(shape.graph, {"p0": 1, "p1": 1, "r": 2}) or ()) == 1


def test_four_pieces_between_direct_repeats_give_six_structures() -> None:
    shape = DirectRepeats(4)
    copies = {"p0": 1, "p1": 1, "p2": 1, "p3": 1, "r": 4}
    assert len(closed_walks(shape.graph, copies) or ()) == 6


def test_counting_structures_stops_at_its_budget_rather_than_guessing() -> None:
    shape = DirectRepeats(4)
    copies = {"p0": 1, "p1": 1, "p2": 1, "p3": 1, "r": 4}
    assert closed_walks(shape.graph, copies, budget=3) is None


def test_the_flip_flop_pair_are_isomers(plastome: Quadripartite) -> None:
    first, second = plastome.paths()
    assert are_isomers(parse_walk(first.title), parse_walk(second.title))


def test_two_orders_of_direct_repeats_are_not_isomers() -> None:
    shape = DirectRepeats(3)
    first = parse_walk(shape.path(("p0", "p1", "p2")).title)
    second = parse_walk(shape.path(("p0", "p2", "p1")).title)
    assert not are_isomers(first, second)


def test_one_inversion_over_other_links_is_not_an_isomer() -> None:
    # The second circle turns ``b`` over between flanks that are not copies of one repeat, so it
    # crosses joins the first does not: reads of one say nothing about the other.
    first = (("a", "+"), ("b", "+"), ("c", "+"))
    second = (("a", "+"), ("b", "-"), ("c", "+"))
    assert not are_isomers(first, second)


def test_the_longest_repeat_of_a_plastome_is_its_inverted_repeat(plastome: Quadripartite) -> None:
    walk = parse_walk(plastome.paths()[0].title)
    assert longest_repeat(plastome.graph, walk) == ("inverted", len(plastome.ir))


def test_the_longest_repeat_of_direct_copies_is_direct() -> None:
    shape = DirectRepeats(2)
    walk = parse_walk(shape.path(("p0", "p1")).title)
    assert longest_repeat(shape.graph, walk) == ("direct", len(shape.repeat))


def test_a_repeat_of_several_segments_is_read_as_one_stretch() -> None:
    # The inverted repeat in two segments, which a variant between the copies' ends would give.
    one, two, lsc, ssc = _bases(60), _bases(70), _bases(200), _bases(90)
    first = one + two[:_K]
    second = two
    mirror = reverse_complement(one + two)
    graph = _gfa(
        {
            "L": mirror[-_K:] + lsc + one[:_K],
            "A": first,
            "B": second,
            "S": two[-_K:] + ssc + mirror[:_K],
        },
        [
            ("L", "+", "A", "+"),
            ("A", "+", "B", "+"),
            ("B", "+", "S", "+"),
            ("S", "+", "B", "-"),
            ("B", "-", "A", "-"),
            ("A", "-", "L", "+"),
        ],
    )
    walk = (("L", "+"), ("A", "+"), ("B", "+"), ("S", "+"), ("B", "-"), ("A", "-"))
    assert longest_repeat(graph, walk) == ("inverted", len(one) + len(two))


def test_a_walk_with_no_repeated_step_has_no_repeat() -> None:
    bases = _bases(200)
    graph = _gfa({"x": bases + bases[:_K]}, [("x", "+", "x", "+")])
    assert longest_repeat(graph, (("x", "+"),)) is None


def test_a_plastome_is_resolved(plastome: Quadripartite) -> None:
    reading = read_selected_graph(
        "graph1", plastome.graph, plastome.paths(), status=STATUS_COMPLETE
    )
    assert reading.call == GRAPH_RESOLVED
    assert reading.reasons == ()
    assert reading.structures == 2
    assert reading.isomers is True
    assert reading.repeat == ("inverted", len(plastome.ir))
    assert (reading.segments, reading.links, reading.components, reading.free_ends) == (3, 4, 1, 0)


def test_depth_is_read_against_the_paths_own_copies(plastome: Quadripartite) -> None:
    reading = read_selected_graph(
        "graph1", plastome.graph, plastome.paths(), status=STATUS_COMPLETE
    )
    rows = {row.name: row for row in reading.rows}
    assert rows["I"].copies == 2
    assert rows["L"].copies == 1
    assert rows["I"].depth_ratio == pytest.approx(2.0, abs=0.01)
    assert rows["L"].depth_ratio == pytest.approx(1.0, abs=0.01)


def test_depth_that_disagrees_with_the_copies_decides_nothing(plastome: Quadripartite) -> None:
    # Measured: on 3 of 47 published closed runs a segment's depth rounds to another copy number.
    uneven = Quadripartite()
    graph = _gfa(
        {name: segment.sequence for name, segment in uneven.graph.segments.items()},
        [("L", "+", "I", "+"), ("I", "+", "S", "+"), ("S", "+", "I", "-"), ("I", "-", "L", "+")],
        depths={"L": 80.0, "I": 110.0, "S": 80.0},
    )
    reading = read_selected_graph("graph1", graph, uneven.paths(), status=STATUS_COMPLETE)
    assert reading.call == GRAPH_RESOLVED


def test_one_structure_is_resolved_and_has_no_isomers() -> None:
    shape = DirectRepeats(2)
    reading = read_selected_graph(
        "graph1", shape.graph, [shape.path(("p0", "p1"))], status=STATUS_COMPLETE
    )
    assert reading.call == GRAPH_RESOLVED
    assert reading.structures == 1
    assert reading.isomers is None
    assert reading.repeat == ("direct", len(shape.repeat))


def test_many_structures_are_unresolved() -> None:
    shape = DirectRepeats(4)
    orders = [("p0", "p1", "p2", "p3"), ("p0", "p1", "p3", "p2"), ("p0", "p2", "p1", "p3")]
    orders += [("p0", "p2", "p3", "p1"), ("p0", "p3", "p1", "p2"), ("p0", "p3", "p2", "p1")]
    paths = [shape.path(order, number) for number, order in enumerate(orders, start=1)]
    reading = read_selected_graph("graph1", shape.graph, paths, status=STATUS_COMPLETE)
    assert reading.call == GRAPH_UNRESOLVED
    assert reading.structures == 6
    assert any("6 structures close through it" in reason for reason in reading.reasons)


def test_two_structures_that_are_not_isomers_are_unresolved() -> None:
    shape = DirectRepeats(3)
    paths = [shape.path(("p0", "p1", "p2"), 1), shape.path(("p0", "p2", "p1"), 2)]
    reading = read_selected_graph("graph1", shape.graph, paths, status=STATUS_COMPLETE)
    assert reading.call == GRAPH_UNRESOLVED
    assert reading.isomers is False


def test_a_structure_the_assembler_did_not_write_is_unresolved(plastome: Quadripartite) -> None:
    # An isomer pair with one written: the same state as a run cut off by --max-paths-num.
    reading = read_selected_graph(
        "graph1", plastome.graph, plastome.paths()[:1], status=STATUS_COMPLETE
    )
    assert reading.call == GRAPH_UNRESOLVED
    assert any("2 structure(s) close through it and 1" in reason for reason in reading.reasons)


def test_a_segment_no_path_uses_is_unresolved(plastome: Quadripartite) -> None:
    graph = _gfa(
        {
            **{name: segment.sequence for name, segment in plastome.graph.segments.items()},
            "extra": _bases(40),
        },
        [("L", "+", "I", "+"), ("I", "+", "S", "+"), ("S", "+", "I", "-"), ("I", "-", "L", "+")],
    )
    reading = read_selected_graph("graph1", graph, plastome.paths(), status=STATUS_COMPLETE)
    assert reading.call == GRAPH_UNRESOLVED
    assert any("1 segment(s) no path uses" in reason for reason in reading.reasons)


def test_a_link_no_path_uses_is_unresolved(plastome: Quadripartite) -> None:
    graph = _gfa(
        {name: segment.sequence for name, segment in plastome.graph.segments.items()},
        [
            ("L", "+", "I", "+"),
            ("I", "+", "S", "+"),
            ("S", "+", "I", "-"),
            ("I", "-", "L", "+"),
            ("L", "+", "S", "+"),
        ],
    )
    reading = read_selected_graph("graph1", graph, plastome.paths(), status=STATUS_COMPLETE)
    assert reading.call == GRAPH_UNRESOLVED
    assert any("1 link(s) no path uses" in reason for reason in reading.reasons)


def test_a_path_its_graph_cannot_spell_is_refused(plastome: Quadripartite) -> None:
    other = _bases(len(plastome.sequence))
    paths = [AssembledPath("p.1", (other,), True, "L+,I+,S+,I-(circular)")]
    with pytest.raises(AssemblyError, match="not from one assembly"):
        read_selected_graph("graph1", plastome.graph, paths, status=STATUS_COMPLETE)


def test_a_closed_result_with_no_path_from_its_graph_is_refused(plastome: Quadripartite) -> None:
    with pytest.raises(AssemblyError, match="no path was written"):
        read_selected_graph("graph1", plastome.graph, [], status=STATUS_COMPLETE)


def test_a_path_that_states_no_walk_is_not_measured(plastome: Quadripartite) -> None:
    paths = [AssembledPath("p.1", (plastome.sequence,), True, "edited by hand")]
    reading = read_selected_graph("graph1", plastome.graph, paths, status=STATUS_COMPLETE)
    assert reading.call == GRAPH_NOT_MEASURED


def test_an_open_graph_is_unresolved_and_counts_no_structure() -> None:
    graph = _gfa({"a": _bases(50), "b": _bases(70)}, [])
    paths = [AssembledPath("p.1", ("A" * 50, "C" * 70), False, "scaffold_1--a+")]
    reading = read_selected_graph("graph1", graph, paths, status=STATUS_FRAGMENTED)
    assert reading.call == GRAPH_UNRESOLVED
    assert reading.structures is None
    assert "2 component(s) and 4 free end(s)" in reading.reasons[0]


def test_an_open_graphs_depth_is_read_against_its_longest_pieces() -> None:
    # A repeat piece is short and deep: the length-weighted median takes the depth of one copy.
    graph = _gfa(
        {"long": _bases(500), "repeat": _bases(100)},
        [],
        depths={"long": 40.0, "repeat": 80.0},
    )
    reading = read_selected_graph("graph1", graph, [], status=STATUS_FRAGMENTED)
    rows = {row.name: row for row in reading.rows}
    assert rows["repeat"].depth_ratio == pytest.approx(2.0, abs=0.01)
    assert rows["long"].copies is None


def test_a_closed_graph_behind_a_path_with_gaps_is_still_unresolved(
    plastome: Quadripartite,
) -> None:
    reading = read_selected_graph("graph1", plastome.graph, [], status=STATUS_GAPPED)
    assert reading.call == GRAPH_UNRESOLVED
    assert "gapped" in reading.reasons[0]


def test_no_graph_at_all_is_not_measured() -> None:
    assert graph_verdict(()) == (GRAPH_NOT_MEASURED, ("no selected graph to read",))


def test_one_resolved_graph_is_the_verdict(plastome: Quadripartite) -> None:
    reading = read_selected_graph(
        "graph1", plastome.graph, plastome.paths(), status=STATUS_COMPLETE
    )
    assert graph_verdict((reading,)) == (GRAPH_RESOLVED, ())


def test_two_resolved_graphs_are_still_not_one_molecule(plastome: Quadripartite) -> None:
    reading = read_selected_graph(
        "graph1", plastome.graph, plastome.paths(), status=STATUS_COMPLETE
    )
    call, reasons = graph_verdict((reading, reading))
    assert call == GRAPH_UNRESOLVED
    assert "2 graphs were selected" in reasons[0]


def test_one_graph_not_measured_is_the_verdict(plastome: Quadripartite) -> None:
    paths = [AssembledPath("p.1", (plastome.sequence,), True, "edited")]
    reading = read_selected_graph("graph1", plastome.graph, paths, status=STATUS_COMPLETE)
    assert graph_verdict((reading,))[0] == GRAPH_NOT_MEASURED
