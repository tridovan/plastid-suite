"""The graph reading on real graphs: every run deposited for Jin et al. 2020 — panel tier.

`.panel/assembly/jin2020/` holds the selected graph, the paths and the log of the 50 runs with the
published settings (`auto--`) and of the 8 re-run with settings chosen per sample (`customized--`),
as the authors deposited them (GetOrganelle 1.6.2; `tools/fetch_assembly_panel.py`).

No tool runs here. The answers are not the reader's own: which runs closed comes from the authors'
logs, and the order of the pieces along a molecule from the published records, read with exact
anchors. Reads are not involved at all, which is the point: every structure of one graph crosses
the same joins, so what separates them cannot be something a short read sees.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
from Bio import SeqIO

from plastid.core.assembly import STATUS_COMPLETE, AssembledPath, AssemblyError
from plastid.core.assembly_graph import (
    GRAPH_RESOLVED,
    GRAPH_UNRESOLVED,
    graph_verdict,
    parse_gfa,
    parse_walk,
    read_selected_graph,
)
from plastid.core.ir import ir_de_novo
from plastid.core.model import Genome
from plastid.core.sequence import reverse_complement
from plastid.external.getorganelle import GetOrganelleOutput, read_output

pytestmark = pytest.mark.panel

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "assembly" / "jin2020"
RECORDS = ROOT / ".panel" / "assembly" / "records"

# The log of this directory is from the run before the authors repeated its last step by hand, so
# it says 9 contigs beside a closed path. The directory contradicts itself and is refused for that;
# its graph and path are read directly below.
REDONE = "customized--ERR964904-plastome"

SEVERAL_STRUCTURES = {
    "auto--ERR268390-plastome": 15,  # Picea abies
    "auto--SRR5028199-plastome": 15,  # Picea sitchensis
    "auto--SRR5602601-plastome": 12,  # Jasminum tortuosum
}
THE_THREE_LEFT_OPEN = ("ERR2206741", "SRR2057084", "SRR6478596")  # Ginkgo, Musa, Salvinia


@pytest.fixture(scope="module")
def outputs() -> dict[str, GetOrganelleOutput]:
    if not PANEL.is_dir():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {PANEL}")
    read = {}
    for root in sorted(PANEL.iterdir()):
        if root.name != REDONE:
            read[root.name] = read_output(root)
    return read


def _closed(outputs: dict[str, GetOrganelleOutput]) -> dict[str, GetOrganelleOutput]:
    return {name: out for name, out in outputs.items() if out.status == STATUS_COMPLETE}


def test_the_published_settings_closed_39_of_the_50(outputs) -> None:
    auto = {name: out for name, out in outputs.items() if name.startswith("auto--")}
    assert len(auto) == 50
    assert sum(out.status == STATUS_COMPLETE for out in auto.values()) == 39


def test_settings_chosen_per_sample_closed_eight_more_and_three_stayed_open(outputs) -> None:
    customized = {name for name in outputs if name.startswith("customized--")} | {REDONE}
    assert len(customized) == 8
    assert all(outputs[name].status == STATUS_COMPLETE for name in customized - {REDONE})
    rescued = {name.split("--")[1] for name in customized}
    open_runs = {
        name.split("--")[1]
        for name, out in outputs.items()
        if name.startswith("auto--") and out.status != STATUS_COMPLETE
    }
    assert sorted(item.removesuffix("-plastome") for item in open_runs - rescued) == sorted(
        THE_THREE_LEFT_OPEN
    )


def test_a_directory_whose_log_predates_its_paths_is_refused() -> None:
    if not (PANEL / REDONE).is_dir():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {PANEL / REDONE}")
    with pytest.raises(AssemblyError, match="the log says"):
        read_output(PANEL / REDONE)


def test_one_segment_closing_on_itself_is_a_resolved_plastome() -> None:
    # Cytinus hypocistis, 19,400 bp: no repeat at all, and the smallest graph a plastome can have.
    root = PANEL / REDONE
    if not root.is_dir():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {root}")
    [gfa] = root.glob("*.selected_graph.gfa")
    [file] = root.glob("*.path_sequence.fasta")
    record = next(SeqIO.parse(file, "fasta"))
    path = AssembledPath("path", (str(record.seq).upper(),), True, record.description)
    graph = parse_gfa(gfa.read_text(encoding="utf-8").splitlines())
    reading = read_selected_graph("graph1", graph, [path], status=STATUS_COMPLETE)
    assert reading.call == GRAPH_RESOLVED
    assert reading.structures == 1
    assert path.length == 19_400


def test_every_closed_run_states_its_walks_and_its_graph_spells_them(outputs) -> None:
    # read_output refuses a path its graph cannot spell, so reaching here is most of the check.
    closed = _closed(outputs)
    assert len(closed) == 46
    for name, out in closed.items():
        [reading] = out.graph_readings
        assert reading.call in (GRAPH_RESOLVED, GRAPH_UNRESOLVED), name
        assert reading.paths_written == len(out.paths), name


def test_every_closed_graph_is_closed_and_wholly_used(outputs) -> None:
    for name, out in _closed(outputs).items():
        [reading] = out.graph_readings
        assert (reading.components, reading.free_ends) == (1, 0), name
        assert not any("no path uses" in reason for reason in reading.reasons), name


def test_the_assembler_wrote_exactly_the_structures_the_graph_allows(outputs) -> None:
    for name, out in _closed(outputs).items():
        [reading] = out.graph_readings
        assert reading.structures == reading.paths_written, name


def test_the_graph_leaves_three_closed_runs_unresolved(outputs) -> None:
    calls = {name: graph_verdict(out.graph_readings)[0] for name, out in _closed(outputs).items()}
    assert Counter(calls.values()) == {GRAPH_RESOLVED: 43, GRAPH_UNRESOLVED: 3}
    for name, structures in SEVERAL_STRUCTURES.items():
        assert calls[name] == GRAPH_UNRESOLVED
        [reading] = outputs[name].graph_readings
        assert reading.structures == structures


def test_every_pair_is_a_flip_flop_pair(outputs) -> None:
    pairs = [
        out.graph_readings[0]
        for out in _closed(outputs).values()
        if out.graph_readings[0].structures == 2
    ]
    assert len(pairs) == 42
    assert all(reading.isomers for reading in pairs)


def test_every_open_run_is_unresolved_and_says_why(outputs) -> None:
    open_runs = {
        name: out
        for name, out in outputs.items()
        if out.status != STATUS_COMPLETE and name != REDONE
    }
    assert len(open_runs) == 11
    for name, out in open_runs.items():
        call, reasons = graph_verdict(out.graph_readings)
        assert call == GRAPH_UNRESOLVED, name
        assert "free end(s)" in reasons[0], name


def _order_along(record: str, segments: dict[str, str], names: list[str]) -> list[tuple[str, str]]:
    """Single-copy segments in the order the record has them, by an exact anchor of each.

    An anchor is 31 bp from a segment's middle that occurs once on the record, on one strand. The
    result is rotated and turned to start at the first name on ``+``, so two orders compare.
    """
    word_length = 31
    backward = reverse_complement(record)
    placed: list[tuple[int, str, str]] = []
    for name in names:
        sequence = segments[name]
        middle = len(sequence) // 2
        for offset in range(0, len(sequence) // 2 - word_length, 37):
            word = sequence[middle + offset : middle + offset + word_length]
            forward, reverse = record.count(word), backward.count(word)
            if forward == 1 and reverse == 0:
                placed.append((record.index(word), name, "+"))
                break
            if reverse == 1 and forward == 0:
                placed.append((len(record) - backward.index(word), name, "-"))
                break
    placed.sort()
    return [(name, strand) for _, name, strand in placed]


def _same_circle(first: list[tuple[str, str]], second: list[tuple[str, str]]) -> bool:
    flip = {"+": "-", "-": "+"}
    turned = [(name, flip[strand]) for name, strand in reversed(second)]
    doubled, doubled_turned = second + second, turned + turned
    size = len(first)
    return any(
        first == source[start : start + size]
        for source in (doubled, doubled_turned)
        for start in range(size)
    )


@pytest.mark.parametrize(
    ("run", "accession", "path_number"),
    [("ERR268390", "HF937082.1", 3), ("SRR5028199", "KU215903.2", 1)],
)
def test_the_published_order_is_one_of_the_structures_and_not_always_the_first(
    outputs, run: str, accession: str, path_number: int
) -> None:
    file = RECORDS / f"{accession}.gb"
    if not file.exists():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {file}")
    out = outputs[f"auto--{run}-plastome"]
    record = str(SeqIO.read(file, "genbank").seq).upper()
    [gfa] = PANEL.joinpath(f"auto--{run}-plastome").glob("*.selected_graph.gfa")
    graph = parse_gfa(gfa.read_text(encoding="utf-8").splitlines())
    segments = {name: segment.sequence for name, segment in graph.segments.items()}
    walks = {path.name.split(".")[-1]: parse_walk(path.title) for path in out.paths}
    copies = Counter(step[0] for step in walks["1"])
    single = [name for name, count in copies.items() if count == 1]
    along = _order_along(record, segments, single)
    assert len(along) == len(single)
    matching = [
        number
        for number, walk in walks.items()
        if _same_circle([step for step in walk if step[0] in single], along)
    ]
    assert matching == [str(path_number)]


def test_the_repeat_in_the_graph_is_the_repeat_in_the_sequence(outputs) -> None:
    agreed = 0
    for name, out in _closed(outputs).items():
        [reading] = out.graph_readings
        sequence = out.paths[0].contigs[0]
        detected = ir_de_novo(Genome(name=name, sequence=sequence, topology="circular"))
        if detected is None:
            continue
        assert reading.repeat == ("inverted", detected.length), name
        agreed += 1
    assert agreed == 41


def test_the_graph_sees_an_inverted_repeat_shorter_than_the_sequence_reader_looks_for(
    outputs,
) -> None:
    # Mankyua chejuensis: the sequence reader's validated window starts at 10 kb.
    out = outputs["auto--SRR7630500-plastome"]
    [reading] = out.graph_readings
    assert reading.repeat == ("inverted", 9_756)
    sequence = out.paths[0].contigs[0]
    assert ir_de_novo(Genome(name="Mankyua", sequence=sequence, topology="circular")) is None


def test_depth_is_reported_and_is_not_a_gate(outputs) -> None:
    # A rule that failed these three runs would refuse correct, resolved plastomes.
    departures = {
        name.split("--")[1].removesuffix("-plastome")
        for name, out in _closed(outputs).items()
        for reading in out.graph_readings
        for row in reading.rows
        if row.copies is not None and round(row.depth_ratio) != row.copies
    }
    assert departures == {"SRR1145775", "SRR5028199", "SRR7630500"}
    assert graph_verdict(outputs["auto--SRR7630500-plastome"].graph_readings)[0] == GRAPH_RESOLVED
