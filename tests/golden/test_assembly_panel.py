"""The assembly step against published read sets — panel tier, real tools.

The answers are not this toolkit's. Each read set is one of the 50 of Jin et al. 2020 (Genome
Biology 21:241), whose Additional file 2 gives the assembled size, the published plastome, the
places they differ, and the mapped depth of each under Bowtie2. A second engine with the same
definition is what the depth comparison is: minimap2 here, Bowtie2 there, on the same sequence.

`.panel/assembly/` holds the reads, the published records and GetOrganelle's output for each set.
Evaluation runs on those outputs through the import backend, in minutes; two sets are assembled from
reads here, because that is the path a user takes. `results_galaxy/` holds one real run made on
usegalaxy.eu, which is the import route as a user meets it.
"""

from __future__ import annotations

import gzip
import json
import os
import statistics
from dataclasses import dataclass
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from conftest import require_tools
from plastid.cli import main
from plastid.core.assembly import same_circular
from plastid.core.sequence import reverse_complement

pytestmark = [pytest.mark.tier_b, pytest.mark.panel]

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / ".panel" / "assembly"


@dataclass(frozen=True)
class ReadSet:
    run: str
    sample: str
    record: str


MAGNOLIA = ReadSet("SRR5602589", "Magnolia officinalis", "KY085916.1")
LAURUS = ReadSet("SRR5602602", "Laurus nobilis", "KY085912.1")
SELAGINELLA = ReadSet("SRR2037123", "Selaginella kraussiana", "MH549643.1")
MUSA = ReadSet("SRR2057084", "Musa balbisiana", "KT595228.1")
ARABIDOPSIS = ReadSet("DRR173054", "Arabidopsis thaliana Col-0", "AP000423.1")


def _tool_directory() -> str | None:
    return os.environ.get("PTS_ASSEMBLY_DIR")


def _require(read_set: ReadSet) -> dict[str, Path]:
    files = {
        "reads_1": PANEL / "reads" / f"{read_set.run}_1.fastq.gz",
        "reads_2": PANEL / "reads" / f"{read_set.run}_2.fastq.gz",
        "record": PANEL / "records" / f"{read_set.record}.gb",
        "result": PANEL / "results" / read_set.run,
    }
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        pytest.skip(f"not in this tree (.panel/ is gitignored): {', '.join(missing)}")
    directory = _tool_directory()
    tools = {
        "minimap2": str(Path(directory) / "minimap2") if directory else None,
    }
    if directory is None:
        from plastid.run.tools import which

        tools["minimap2"] = which("minimap2")
    require_tools(tools, "minimap2")
    return files


def _sets() -> list[str]:
    directory = _tool_directory()
    return ["--set", f"tools={directory}"] if directory else []


def _table(run: Path, name: str) -> list[dict]:
    payload = json.loads((run / "01_assembly" / "result.json").read_text(encoding="utf-8"))
    return payload["tables"][name]["rows"] if name in payload["tables"] else []


def _imported(
    read_set: ReadSet,
    out: Path,
    *,
    result: Path | None = None,
    source: str = "panel: GetOrganelle 1.7.7.1 on the development machine",
) -> tuple[Path, int]:
    """Stored assembler output through the import backend, untrimmed as Jin ran it."""
    files = _require(read_set)
    prepared_root = out / "prepared"
    code = main(
        [
            "assemble",
            "prepare",
            "-1",
            str(files["reads_1"]),
            "-2",
            str(files["reads_2"]),
            "--sample",
            read_set.sample,
            "-o",
            str(prepared_root),
            "--set",
            "trim=no",
            *_sets(),
        ]
    )
    assert code == 0
    [prepared] = list(prepared_root.glob("*/*/01_assembly/prepared.json"))
    imported = out / "imported"
    code = main(
        [
            "assemble",
            "run",
            "--backend",
            "import",
            "--prepared",
            str(prepared),
            "--result",
            str(result if result is not None else files["result"]),
            "--source-name",
            source,
            "--retrieved-at",
            "2026-10-31",
            "--reference",
            str(files["record"]),
            "-o",
            str(imported),
            "--set",
            "trim=no",
            *_sets(),
        ]
    )
    [manifest] = list(imported.glob("*/*/manifest.json"))
    return manifest.parent, code


@pytest.fixture(scope="module")
def magnolia(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, int]:
    return _imported(MAGNOLIA, tmp_path_factory.mktemp("magnolia"))


@pytest.fixture(scope="module")
def laurus(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, int]:
    return _imported(LAURUS, tmp_path_factory.mktemp("laurus"))


@pytest.fixture(scope="module")
def musa(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, int]:
    return _imported(MUSA, tmp_path_factory.mktemp("musa"))


def test_magnolia_is_the_published_molecule_and_every_join_is_supported(magnolia) -> None:
    run, code = magnolia
    assert code == 0
    [summary] = _table(run, "assembly_summary")
    # Jin Table S1: GetOrganelle, NOVOPlasty and the record all give 160,136 bp, identical.
    assert summary["length_bp"] == 160_136
    assert summary["status"] == "complete"
    assert summary["read_support"] == "supported"
    assert any(row["identical"] == "yes" for row in _table(run, "reference_comparison"))
    assert {row["call"] for row in _table(run, "junction_support")} == {"supported"}
    assert _table(run, "weak_regions") == []


def test_the_magnolia_graph_is_resolved_and_reads_the_repeat_the_sequence_does(magnolia) -> None:
    run, _ = magnolia
    [summary] = _table(run, "assembly_summary")
    assert summary["graph_call"] == "resolved"
    [graph] = _table(run, "assembly_graph")
    assert (graph["segments"], graph["links"]) == (3, 4)
    assert (graph["components"], graph["free_ends"]) == (1, 0)
    assert (graph["structures"], graph["isomers"]) == (2, "yes")
    assert graph["repeat"] == f"{summary['ir_bp']} bp, inverted"
    rows = {row["copies"]: row for row in _table(run, "graph_segments")}
    assert rows[2]["depth_ratio"] == pytest.approx(2.0, abs=0.1)
    assert rows[1]["depth_ratio"] == pytest.approx(1.0, abs=0.1)


def test_the_magnolia_depth_agrees_with_the_published_bowtie2_depth(magnolia) -> None:
    run, _ = magnolia
    [summary] = _table(run, "assembly_summary")
    # Jin Table S2 on this sequence: 151.30 ± 29.22. A second aligner, the same definition.
    assert summary["mean_depth"] == pytest.approx(151.30, rel=0.01)
    assert summary["sd_depth"] == pytest.approx(29.22, rel=0.02)


def _circles(directory: Path) -> dict[str, str]:
    return {
        path.name: str(next(SeqIO.parse(path, "fasta")).seq).upper()
        for path in sorted(directory.glob("*.path_sequence.fasta"))
    }


def test_a_real_run_on_galaxy_comes_back_as_the_local_assembly(tmp_path: Path) -> None:
    """The import route as a user meets it: a server's output, on other versions of its tools."""
    galaxy = PANEL / "results_galaxy" / MAGNOLIA.run
    if not galaxy.is_dir():
        pytest.skip(f"not in this tree (.panel/ is gitignored): {galaxy}")
    run, code = _imported(MAGNOLIA, tmp_path, result=galaxy, source="usegalaxy.eu")
    assert code == 0
    # Both paths, base for base, although the server's assembler ran SPAdes 3.15.5 and not 4.3.0.
    assert _circles(galaxy) == _circles(PANEL / "results" / MAGNOLIA.run)
    [summary] = _table(run, "assembly_summary")
    assert summary["read_support"] == "supported"
    assert summary["where"] == "imported: usegalaxy.eu, retrieved 2026-10-31"
    assert "SPAdes 3.15.5" in summary["assembler_programs"]
    assert any(row["identical"] == "yes" for row in _table(run, "reference_comparison"))
    described = json.loads((run / "01_assembly" / "external_result.json").read_text("utf-8"))
    assert described["source"]["name"] == "usegalaxy.eu"
    assert described["tool"]["version"] == "1.7.7.1"
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    assert [item["format"] for item in manifest["artifacts"]] == ["fasta"]


def test_laurus_is_the_published_reassembly(laurus) -> None:
    run, code = laurus
    assert code == 0
    [summary] = _table(run, "assembly_summary")
    # Jin Table S1: reassembled 153,040 bp against 152,750 bp published.
    assert summary["length_bp"] == 153_040
    assert summary["read_support"] == "supported"


def test_laurus_differs_from_its_record_where_the_paper_says(laurus) -> None:
    run, _ = laurus
    differences = _table(run, "reference_differences")
    assert len(differences) == 2
    by_size = {row["gap_columns"]: row for row in differences}
    # Jin Table S1: "1 sites, 429 bp". The 139 bp is the record's own tail, at its origin.
    assert set(by_size) == {429, 139}
    assert by_size[429]["after_reference_base"] == 113_581
    assert by_size[429]["path_bp"] == 429
    assert by_size[139]["reference_bp"] == 139
    assert sum(row["mismatches"] for row in differences) == 0


def test_the_reads_side_with_the_assembly_at_both_differences(laurus) -> None:
    run, _ = laurus
    rows = {row["sequence"]: row for row in _table(run, "reference_read_support")}
    assert rows["assembly"]["read_support"] == "supported"
    assert rows["KY085912.1"]["read_support"] == "undetermined"
    assert rows["assembly"]["fewest_spanning"] > 50
    assert rows["KY085912.1"]["fewest_spanning"] == 0
    # Jin Table S2: 221.60 ± 31.43 reassembled, 221.65 ± 37.88 published. The mean cannot tell the
    # two apart and the spread can, which is the paper's own reading.
    assert rows["assembly"]["mean_depth"] == pytest.approx(221.60, rel=0.01)
    assert rows["KY085912.1"]["mean_depth"] == pytest.approx(221.65, rel=0.01)
    assert rows["assembly"]["sd_depth"] == pytest.approx(31.43, rel=0.03)
    assert rows["KY085912.1"]["sd_depth"] == pytest.approx(37.88, rel=0.03)
    assert rows["assembly"]["sd_depth"] < rows["KY085912.1"]["sd_depth"]

    record = {row["junction"]: row for row in _table(run, "reference_junctions")}
    assert record["JSB"]["call"] == "undetermined"
    assert record["JSB"]["after_base"] == 113_581
    assert record["origin"]["call"] == "undetermined"
    assert record["origin"]["spanning"] == 0
    # The reads are silent at exactly the two places the sequences differ, and nowhere else.
    [gap, origin] = _table(run, "reference_weak_regions")
    assert gap["first_after_base"] <= 113_581 <= gap["last_after_base"]
    assert origin["first_after_base"] > origin["last_after_base"]  # it wraps the origin


def test_musa_stays_open_and_passes_nothing_downstream(musa) -> None:
    run, code = musa
    assert code == 1
    [summary] = _table(run, "assembly_summary")
    # Jin Table S1: no closed path, 14 fragments with 1.6.2; 15 scaffolds with 1.7.7.1.
    assert summary["status"] == "fragmented"
    assert summary["read_support"] == "not measured"
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["artifacts"] == []
    assert manifest["status"] == "incomplete"


def test_the_musa_graph_is_unresolved_and_kept_for_a_person_to_read(musa) -> None:
    run, _ = musa
    [summary] = _table(run, "assembly_summary")
    assert summary["graph_call"] == "unresolved"
    assert "15 component(s) and 30 free end(s)" in summary["graph_reasons"]
    [graph] = _table(run, "assembly_graph")
    assert (graph["segments"], graph["links"]) == (15, 0)
    assert graph["structures"] == "not measured"
    kept = {path.name for path in (run / "01_assembly").iterdir()}
    assert "assembler.log" in kept
    assert "scaffolds.fasta" in kept
    assert "embplant_pt.K91.contigs.graph1.selected_graph.gfa" in kept
    assert "extended_K91.assembly_graph.fastg.extend-embplant_pt-embplant_mt.fastg" in kept


# The repeat's second copy starts here on KT595228.1, and only one copy may be indexed for its words
# to be unique. The library's insert is the log's 508.8 ± 13.7 bp.
_MUSA_ONE_REPEAT_COPY = 134_400
_MUSA_PAIRS = 2_000_000
_WORD = 31


def _musa_edges(names: tuple[str, ...]) -> set[str]:
    """The 31-mers of the named edges of the wider graph, on both strands."""
    root = PANEL / "results" / MUSA.run
    [fastg] = root.glob("extended_*.extend-*.fastg")
    words: set[str] = set()
    for record in SeqIO.parse(fastg, "fasta"):
        identifier = record.id.split("_")[1]
        if identifier in names and "'" not in record.id.split(":")[0]:
            sequence = str(record.seq).upper()
            for text in (sequence, reverse_complement(sequence)):
                words.update(text[i : i + _WORD] for i in range(len(text) - _WORD + 1))
    return words


def test_the_musa_arrays_are_longer_than_its_reads_and_than_its_record() -> None:
    """Why it stays open, measured from the reads with no aligner and no assembler.

    Three tandem arrays are loops of the graph. A read made wholly of an array's words can only
    exist if the array is at least a read long, so no read crosses it. Pairs, each mate placed by
    its first 31 bases, say how much longer: where the sample holds D bases more than the record,
    a pair spans D fewer on the record, and this library's insert is tight.
    """
    files = _require(MUSA)
    record = str(SeqIO.read(files["record"], "genbank").seq).upper()
    # The first array as the record has it: 84 bp, shorter than the 100 bp reads.
    array = record[61_250:61_334]
    assert all(array[i] == array[i + 24] for i in range(len(array) - 24))

    arrays = {
        "24 bp unit": _musa_edges(("377", "6827")),
        "48 bp unit": _musa_edges(("4261", "6631")),
        "self-loop": _musa_edges(("5",)),
    }
    indexed = record[:_MUSA_ONE_REPEAT_COPY]
    seen: set[str] = set()
    forward: dict[str, int] = {}
    for position in range(len(indexed) - _WORD + 1):
        word = indexed[position : position + _WORD]
        if word in seen:
            forward.pop(word, None)
        else:
            seen.add(word)
            forward[word] = position
    backward = {reverse_complement(word): position for word, position in forward.items()}

    wholly_inside = dict.fromkeys(arrays, 0)
    placed: list[tuple[int, int]] = []
    with gzip.open(files["reads_1"], "rt") as one, gzip.open(files["reads_2"], "rt") as two:
        for index, ((_, first, _), (_, second, _)) in enumerate(
            zip(FastqGeneralIterator(one), FastqGeneralIterator(two), strict=True)
        ):
            if index == _MUSA_PAIRS:
                break
            for name, words in arrays.items():
                for read in (first, second):
                    if read[:_WORD] in words and all(
                        read[i : i + _WORD] in words for i in range(0, len(read) - _WORD + 1, 7)
                    ):
                        wholly_inside[name] += 1
            head_1, head_2 = first[:_WORD], second[:_WORD]
            start, end = forward.get(head_1), backward.get(head_2)
            if start is None or end is None:
                start, end = forward.get(head_2), backward.get(head_1)
            if start is None or end is None:
                continue
            size = end + _WORD - start
            if 150 <= size <= 1_500:
                placed.append(((start + end + _WORD) // 2, size))

    assert all(count >= 100 for count in wholly_inside.values()), wholly_inside
    library = statistics.median(size for _, size in placed)
    assert library == pytest.approx(509, abs=2)
    for (low, high), longer in {
        (61_100, 61_500): 75,  # the 24 bp array
        (90_650, 90_800): 167,  # the self-loop
        (118_800, 119_050): 268,  # the 48 bp array
    }.items():
        spans = [size for middle, size in placed if low <= middle < high]
        assert library - statistics.median(spans) == pytest.approx(longer, abs=8)


def test_the_musa_record_has_joins_the_reads_do_not_cross(musa) -> None:
    run, _ = musa
    rows = {row["sequence"]: row for row in _table(run, "reference_read_support")}
    assert rows["KT595228.1"]["read_support"] == "undetermined"
    assert rows["KT595228.1"]["fewest_spanning"] == 0
    assert rows["KT595228.1"]["weak_regions"] == 8
    # Jin Table S2 on this sequence: 2205.78 ± 470.46.
    assert rows["KT595228.1"]["mean_depth"] == pytest.approx(2205.78, rel=0.02)


def test_selaginella_two_direct_repeats_and_no_inverted_repeat(tmp_path: Path) -> None:
    run, code = _imported(SELAGINELLA, tmp_path)
    assert code == 0
    [summary] = _table(run, "assembly_summary")
    # Jin Table S1 and Mower et al. 2019: 129,971 bp, identical to MH549643, repeats in direct
    # orientation. NOVOPlasty lost 60,647 bp between them; a closed path is not evidence by itself.
    assert summary["length_bp"] == 129_971
    assert summary["ir_bp"] == "no inverted repeat found"
    assert summary["assembler_repeat"] == "14597 bp, DRs detected"
    assert summary["read_support"] == "supported"
    # One structure: the graph reads the two copies as direct too, independently of the log.
    assert summary["graph_call"] == "resolved"
    [graph] = _table(run, "assembly_graph")
    assert (graph["structures"], graph["isomers"]) == (1, "not applicable")
    assert graph["repeat"] == "14597 bp, direct"
    assert any(row["identical"] == "yes" for row in _table(run, "reference_comparison"))
    assert [row["junction"] for row in _table(run, "junction_support")] == ["origin"]


def test_the_arabidopsis_assembly_is_one_base_from_a_sanger_record(tmp_path: Path) -> None:
    """An answer from other data and another technology: Sato et al. 1999, BAC and PCR, Sanger."""
    run, code = _imported(ARABIDOPSIS, tmp_path)
    assert code == 0
    [summary] = _table(run, "assembly_summary")
    assert summary["read_support"] == "supported"
    assert summary["ir_bp"] == 26_264
    [difference] = _table(run, "reference_differences")
    # One T in the run at 28,673-28,685: 13 in the record, 14 assembled, and 308 reads against 50
    # for 14. Whether that is an old error or a difference between stocks is not decided here.
    assert difference["after_reference_base"] == 28_685
    assert (difference["reference_bp"], difference["path_bp"]) == (0, 1)
    rows = {row["sequence"]: row for row in _table(run, "reference_read_support")}
    assert {row["read_support"] for row in rows.values()} == {"supported"}


def test_the_default_route_through_quality_control_keeps_the_molecule(tmp_path: Path) -> None:
    """What a user runs, fastp included, on reads unlike the MiSeq set it was first tried on."""
    files = _require(SELAGINELLA)
    directory, database = _tool_directory(), os.environ.get("PTS_GETORGANELLE_DB")
    if directory is None or database is None:
        pytest.skip("set PTS_ASSEMBLY_DIR and PTS_GETORGANELLE_DB to run the assembler here")
    code = main(
        [
            "assemble",
            "run",
            "-1",
            str(files["reads_1"]),
            "-2",
            str(files["reads_2"]),
            "--sample",
            SELAGINELLA.sample,
            "-o",
            str(tmp_path),
            "--reference",
            str(files["record"]),
            "--set",
            f"database={database}",
            *_sets(),
        ]
    )
    assert code == 0
    [manifest] = list(tmp_path.glob("*/*/manifest.json"))
    run = manifest.parent
    [qc] = _table(run, "read_qc")
    assert 0 < qc["reads_after"] < qc["reads_before"]
    assert any(row["identical"] == "yes" for row in _table(run, "reference_comparison"))
    [summary] = _table(run, "assembly_summary")
    assert summary["read_support"] == "supported"
    described = json.loads(manifest.read_text(encoding="utf-8"))
    assert [item["format"] for item in described["artifacts"]] == ["fasta"]
    assert [tool["name"] for tool in described["modules"][0]["tools"]][:2] == [
        "fastp",
        "getorganelle",
    ]


def _isomers(run: Path) -> list[str]:
    sequences: list[str] = []
    text = (run / "01_assembly" / "assembly.paths.fasta").read_text(encoding="utf-8")
    for record in text.split(">")[1:]:
        sequences.append("".join(record.splitlines()[1:]))
    return sequences


def test_a_local_run_from_reads_gives_the_imported_assembly(tmp_path: Path, magnolia) -> None:
    """The path a user takes: the two backends must give one assembly, or name the difference.

    The two isomers are compared as circles, because their order and origin are the assembler's and
    are not part of the molecule.
    """
    files = _require(MAGNOLIA)
    directory, database = _tool_directory(), os.environ.get("PTS_GETORGANELLE_DB")
    if directory is None or database is None:
        pytest.skip("set PTS_ASSEMBLY_DIR and PTS_GETORGANELLE_DB to run the assembler here")
    code = main(
        [
            "assemble",
            "run",
            "--backend",
            "local",
            "-1",
            str(files["reads_1"]),
            "-2",
            str(files["reads_2"]),
            "--sample",
            MAGNOLIA.sample,
            "-o",
            str(tmp_path),
            "--reference",
            str(files["record"]),
            "--set",
            "trim=no",
            "--set",
            f"database={database}",
            *_sets(),
        ]
    )
    assert code == 0
    [manifest] = list(tmp_path.glob("*/*/manifest.json"))
    local = manifest.parent
    imported, _ = magnolia
    first, second = _isomers(local), _isomers(imported)
    assert len(first) == len(second) == 2
    assert all(any(same_circular(one, other) for other in second) for one in first)
    keep = ("length_bp", "status", "read_support", "ir_bp", "lsc_bp", "ssc_bp")
    [local_row], [imported_row] = (
        _table(local, "assembly_summary"),
        _table(imported, "assembly_summary"),
    )
    assert {key: local_row[key] for key in keep} == {key: imported_row[key] for key in keep}
    assert any(row["identical"] == "yes" for row in _table(local, "reference_comparison"))
