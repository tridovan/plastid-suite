"""The assembly step on long reads, against an answer from other data — panel tier, real tools.

*Arabidopsis thaliana* Col-0, one study (`PRJEB77203`): HiFi reads `ERR13987671` go through Oatk,
and the first 1.5 GB of the nanopore run `ERR13987449` ask what HiFi reads are too short to answer.
The answer is `AP000423.1`, assembled from BACs and PCR with Sanger reads in 1999, so it owes
nothing to these reads. `DRR173054`, stored for `test_assembly_panel.py`, is the accession in short
reads.

`.panel/assembly/results/ERR13987671/` is Oatk's own output of this step's default command.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from collections import Counter
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from conftest import require_tools
from plastid.cli import main
from plastid.core.assembly import compare_to_reference, same_circular
from plastid.core.assembly_graph import closed_walks, parse_gfa, parse_walk, spell
from plastid.core.sequence import reverse_complement
from plastid.external.minimap2 import map_long, resolve_minimap2
from plastid.external.oatk import read_output

pytestmark = [pytest.mark.tier_b, pytest.mark.panel]

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "assembly"
HIFI = PANEL / "reads" / "ERR13987671.fastq.gz"
NANOPORE = PANEL / "reads" / "ERR13987449.fastq.gz"
SHORT = tuple(PANEL / "reads" / f"DRR173054_{mate}.fastq.gz" for mate in (1, 2))
RESULT = PANEL / "results" / "ERR13987671"
RECORD = PANEL / "records" / "AP000423.1.gb"
SAMPLE = "Arabidopsis thaliana Col-0"


def _needed() -> list[str]:
    missing = [str(path) for path in (HIFI, RESULT, RECORD) if not path.exists()]
    if missing:
        pytest.skip(f"not in this tree (.panel/ is gitignored): {', '.join(missing)}")
    directory = os.environ.get("PTS_ASSEMBLY_DIR")
    tools = {"minimap2": str(Path(directory) / "minimap2") if directory else None}
    require_tools(tools, "minimap2")
    return ["--set", f"tools={directory}"]


def _table(run: Path, name: str) -> list[dict]:
    payload = json.loads((run / "01_assembly" / "result.json").read_text(encoding="utf-8"))
    return payload["tables"][name]["rows"] if name in payload["tables"] else []


def _run_dir(out: Path) -> Path:
    [manifest] = list(out.glob("*/*/manifest.json"))
    return manifest.parent


@pytest.fixture(scope="module")
def imported(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, int]:
    """Stored output through the import backend: the evaluation of a local run, in seconds."""
    sets = _needed()
    root = tmp_path_factory.mktemp("hifi")
    reads = ["--hifi", str(HIFI), "--sample", SAMPLE]
    assert main(["assemble", "prepare", *reads, "-o", str(root / "p"), *sets]) == 0
    [prepared] = list((root / "p").glob("*/*/01_assembly/prepared.json"))
    code = main(
        [
            "assemble",
            "run",
            "--backend",
            "import",
            "--prepared",
            str(prepared),
            "--result",
            str(RESULT),
            "--source-name",
            "panel: Oatk 1.0 on the development machine",
            "--retrieved-at",
            "2026-10-31",
            "--reference",
            str(RECORD),
            "-o",
            str(root / "runs"),
            *sets,
        ]
    )
    return _run_dir(root / "runs"), code


def test_the_hifi_assembly_is_two_bases_from_a_sanger_record(imported) -> None:
    run, code = imported
    assert code == 0
    [summary] = _table(run, "assembly_summary")
    assert summary["assembler"] == "Oatk 1.0"
    assert summary["status"] == "complete"
    assert summary["length_bp"] == 154_480
    [compared, _] = _table(run, "reference_comparison")
    assert (compared["mismatches"], compared["gap_columns"]) == (0, 2)
    places = [row["after_reference_base"] for row in _table(run, "reference_differences")]
    assert places == [28_685, 67_176]


def test_every_join_is_crossed_by_many_reads(imported) -> None:
    run, _ = imported
    [summary] = _table(run, "assembly_summary")
    assert summary["read_support"] == "supported"
    assert summary["fewest_spanning"] > 1_000
    assert summary["mean_depth"] == pytest.approx(1453, abs=15)
    assert _table(run, "weak_regions") == []


def test_the_graph_is_resolved_with_the_isomer_oatk_left_out(imported) -> None:
    run, _ = imported
    [summary] = _table(run, "assembly_summary")
    assert summary["graph_call"] == "resolved"
    assert summary["n_paths"] == 2
    [graph] = _table(run, "assembly_graph")
    assert (graph["structures"], graph["isomers"]) == (2, "yes")
    # The segment is 25,987 bp. The repeat is carried out to where the copies part.
    assert graph["repeat"] == f"{summary['ir_bp']} bp, inverted"
    assert summary["ir_bp"] == 26_264


def test_the_repeat_and_the_regions_are_those_of_the_short_read_assembly(imported) -> None:
    run, _ = imported
    [summary] = _table(run, "assembly_summary")
    short = sorted((PANEL / "results" / "DRR173054").glob("*.path_sequence.fasta"))
    sequences = [str(next(SeqIO.parse(path, "fasta")).seq).upper() for path in short]
    hifi = read_output(RESULT).paths
    # One of the two technologies' isomers against the other's: one base apart.
    compared = [compare_to_reference(hifi[0].contigs[0], sequence) for sequence in sequences]
    near = next(item for item in compared if item.identity is not None)
    assert (near.mismatches, near.gap_columns) == (0, 1)
    assert summary["length_bp"] == len(sequences[0]) + 1


def _run_lengths(files: tuple[Path, ...], left: str, right: str) -> Counter:
    """How long the run between two exact flanks is, read by read, on either strand."""
    forward = re.compile(left + "(.{0,40}?)" + right)
    backward = re.compile(reverse_complement(right) + "(.{0,40}?)" + reverse_complement(left))
    lengths: Counter = Counter()
    for file in files:
        with gzip.open(file, "rt") as handle:
            for _, sequence, _ in FastqGeneralIterator(handle):
                for hit in (forward.search(sequence), backward.search(sequence)):
                    if hit:
                        lengths[len(hit.group(1))] += 1
    return lengths


def test_two_technologies_read_14_t_where_the_record_has_13() -> None:
    """The base short reads left unexplained in M7.1: HiFi reads, with no Illumina chemistry, agree.

    Whether 1999's 13 was an error or a different stock of Col-0 these reads cannot say.
    """
    _needed()
    left, right = "TTAAGACGAGAAATTGAAATC", "GATAGTAAAATTTTATTCAAA"
    record = str(SeqIO.read(RECORD, "genbank").seq).upper()
    assert f"{left}{'T' * 13}{right}" in record
    for files in ((HIFI,), SHORT):
        lengths = _run_lengths(files, left, right)
        assert lengths.most_common(1)[0][0] == 14
        assert lengths[14] > 5 * lengths[13]


def test_the_other_difference_is_a_run_hifi_reads_do_not_settle() -> None:
    """The assembly's 14 G against the record's 13: an assembler with no polish keeps a majority.

    HiFi reads split, and nothing else here reads that base better. It is reported, not explained.
    """
    _needed()
    left, right = "GTAGAATAAATTAGAAAAGGT", "AAAGGTCATAAAAACGAGGGT"
    lengths = _run_lengths((HIFI,), left, right)
    [(first, most), (second, next_most)] = lengths.most_common(2)
    assert {first, second} == {13, 14}
    assert next_most > 0.7 * most


def _local(tmp_path: Path, *extra: str) -> str:
    sets = _needed()
    database = os.environ.get("PTS_OATK_DB")
    if not database:
        pytest.skip("set PTS_OATK_DB to run the assembler")
    command = ["assemble", "run", "--hifi", str(HIFI), "--sample", SAMPLE]
    code = main([*command, "-o", str(tmp_path), *sets, "--set", f"oatk_db={database}", *extra])
    assert code == 0
    fasta = _run_dir(tmp_path) / "01_assembly" / "assembly.fasta"
    return str(next(SeqIO.parse(fasta, "fasta")).seq).upper()


def _stored() -> str:
    return read_output(RESULT).paths[0].contigs[0]


def test_a_local_run_gives_the_stored_assembly(tmp_path: Path) -> None:
    assert same_circular(_local(tmp_path), _stored())


def test_a_tenth_of_the_reads_give_the_same_plastome(tmp_path: Path) -> None:
    assert same_circular(_local(tmp_path, "--set", "max_bases=150000000"), _stored())


def test_a_much_higher_coverage_threshold_gives_the_same_plastome(tmp_path: Path) -> None:
    assert same_circular(_local(tmp_path, "--set", "min_coverage=400"), _stored())


def test_single_nanopore_reads_cross_the_repeat_on_both_isomers_alike(tmp_path: Path) -> None:
    """What was only ever inferred from a graph, seen in single molecules.

    Both isomers are spelled from the graph and mapped together. A read counts when one alignment
    covers a whole repeat copy and 500 bp beyond each end, so it states which flanks that copy
    joins: it fits one isomer end to end, and the other only up to the repeat.
    """
    _needed()
    if not NANOPORE.exists():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {NANOPORE}")
    [gfa] = RESULT.glob("*.pltd.gfa")
    graph = parse_gfa(gfa.read_text(encoding="utf-8").splitlines())
    walk = parse_walk(read_output(RESULT).paths[0].title)
    copies = Counter(step[0] for step in walk)
    repeat = next(name for name, count in copies.items() if count == 2)
    margin, reach = 500, 60_000
    targets = tmp_path / "isomers.fasta"
    places: dict[str, tuple[int, list[tuple[int, int]]]] = {}
    with targets.open("w", encoding="ascii") as handle:
        for index, candidate in enumerate(sorted(closed_walks(graph, copies)), start=1):
            sequence = spell(graph, candidate)
            handle.write(f">isomer{index}\n{sequence + sequence[:reach]}\n")
            spans, position = [], 0
            for step, after in zip(candidate, (*candidate[1:], candidate[0]), strict=True):
                length = len(graph.segments[step[0]].sequence)
                if step[0] == repeat:
                    spans.append((position, position + length))
                position += length - graph.following[step][after]
            places[f"isomer{index}"] = (len(sequence), spans)

    class _Context:
        """The one thing ``map_long`` asks of a run, for a measurement that is not one."""

        def tool(self, name: str, argv: list[str], **kwargs: object) -> None:
            from plastid.run.tools import run_tool

            run_tool(name, argv, log_dir=tmp_path / "logs", **kwargs)

    binary = resolve_minimap2(os.environ["PTS_ASSEMBLY_DIR"])
    sam = map_long(
        targets,
        NANOPORE,
        preset="map-ont",
        out=tmp_path / "nanopore.sam",
        threads=4,
        ctx=_Context(),
        binary=binary,
    )
    across: Counter = Counter()
    with sam.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 6)
            if int(fields[1]) & 0x904:
                continue
            start = int(fields[3]) - 1
            span = sum(int(count) for count, kind in re.findall(r"(\d+)([=XDN])", fields[5]))
            length, spans = places[fields[2]]
            for low, high in spans:
                for shift in (0, length):
                    if start <= low + shift - margin and start + span >= high + shift + margin:
                        across[fields[2]] += 1
    assert across["isomer1"] == pytest.approx(305, abs=25)
    assert across["isomer2"] == pytest.approx(307, abs=25)
