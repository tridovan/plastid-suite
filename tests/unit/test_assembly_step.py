"""The assembly step through the command line, with its tools stubbed at their boundary.

The stubs write what the real tools write: a GetOrganelle output directory, a fastp report, a SAM.
What is under test is what the step decides from them, which a real assembly cannot exercise on
demand: when a sequence becomes an artifact for the next step, when the run stops and why, and
that the two backends give one result. The real tools run in `tests/golden/test_assembly_panel.py`.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import random
import stat
from collections.abc import Sequence
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.core.assembly import AssemblyError
from plastid.core.sequence import reverse_complement
from plastid.external import fastp, getorganelle, minimap2
from plastid.modules import assembly
from plastid.run.tools import ToolError

_LSC, _IR, _SSC = 20_000, 11_000, 5_000
_EXTENSION = 1_000
_READ = 150
_K = 127


def _sequence() -> str:
    """A circle laid out LSC, IRb, SSC, IRa, with a real inverted repeat."""
    source = random.Random(20261031)

    def random_bases(length: int) -> str:
        return "".join(source.choice("ACGT") for _ in range(length))

    lsc, ir, ssc = random_bases(_LSC), random_bases(_IR), random_bases(_SSC)
    return lsc + ir + ssc + reverse_complement(ir)


SEQUENCE = _sequence()
LENGTH = len(SEQUENCE)
ISOMER = (
    SEQUENCE[:_LSC]
    + SEQUENCE[_LSC : _LSC + _IR]
    + reverse_complement(SEQUENCE[_LSC + _IR : _LSC + _IR + _SSC])
    + SEQUENCE[_LSC + _IR + _SSC :]
)

_LOG_HEAD = """GetOrganelle v1.7.7.1

DEPENDENCIES: Bowtie2 2.5.5; SPAdes 4.3.0; Blast 2.5.0; Bandage 0.9.0
SEED  DB: embplant_pt 0.0.1; embplant_mt 0.0.1
LABEL DB: embplant_pt 0.0.1; embplant_mt 0.0.1
2026-10-31 09:35:18,285 - INFO: Reads used = 295+295
"""
_LOG_COMPLETE = (
    _LOG_HEAD
    + "2026-10-31 09:54:12,311 - INFO: Average embplant_pt base-coverage = 141.562\n"
    + "2026-10-31 09:54:12,500 - INFO: Detecting large repeats (>1000 bp) in PATH1 with "
    + f"IRs detected, Total:LSC:SSC:Repeat(bp) = {LENGTH}:{_LSC}:{_SSC}:{_IR}\n"
    + "2026-10-31 09:54:12,777 - INFO: Result status of embplant_pt: circular genome\n"
)
_LOG_SCAFFOLDS = _LOG_HEAD + "2026-10-31 INFO: Result status of embplant_pt: 3 scaffold(s)\n"
_LOG_NOTHING = _LOG_HEAD + "2026-10-31 INFO: Extracting embplant_pt from the assemblies finished.\n"

_FASTP_REPORT = {
    "summary": {
        "before_filtering": {
            "total_reads": 600,
            "total_bases": 90_000,
            "q30_rate": 0.9,
            "read1_mean_length": 150,
        },
        "after_filtering": {
            "total_reads": 590,
            "total_bases": 86_000,
            "q30_rate": 0.93,
            "read1_mean_length": 146,
        },
    },
    "adapter_cutting": {"adapter_trimmed_reads": 12},
}


def _sam(*, skip_between: tuple[int, int] | None = None, depth: int = 12) -> str:
    """Reads of 150 bp every 25 bp on the extended sequence, each placed ``depth`` times.

    ``skip_between`` leaves out every read overlapping a span, which is what a join the reads do
    not cross looks like.
    """
    lines = [f"@SQ\tSN:target\tLN:{LENGTH + _EXTENSION}\n"]
    for start in range(0, LENGTH + _EXTENSION - _READ, 25):
        if skip_between and start < skip_between[1] and start + _READ > skip_between[0]:
            continue
        for _ in range(depth):
            lines.append(f"r\t0\ttarget\t{start + 1}\t60\t{_READ}=\t*\t0\t0\t*\t*\n")
    return "".join(lines)


def _script(path: Path, body: str) -> None:
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _argument_reader(names: Sequence[str]) -> str:
    """Shell that sets ``$name`` for each ``--name value`` it finds, dashes as underscores."""
    cases = "".join(f'    {flag}) {variable}="$2";;\n' for flag, variable in names)
    return f'while [ $# -gt 0 ]; do\n  case "$1" in\n{cases}  esac\n  shift\ndone\n'


def _tools(directory: Path, *, sam: str) -> Path:
    """Stand-ins for the five programs, and the canned output of the assembler."""
    directory.mkdir(parents=True, exist_ok=True)
    canned = directory / "canned"
    canned.mkdir(exist_ok=True)
    _script(
        directory / "get_organelle_from_reads.py",
        "case \"$1\" in --version) echo 'GetOrganelle v1.7.7.1'; exit 0;; esac\n"
        + _argument_reader([("-o", "out")])
        + f'mkdir -p "$out"\ncp -r "{canned}/." "$out"\n',
    )
    _script(
        directory / "spades.py",
        "case \"$1\" in --version) echo 'SPAdes genome assembler v4.3.0'; exit 0;; esac\n",
    )
    _script(directory / "bowtie2", "echo 'bowtie2 version 2.5.5'\n")
    _script(directory / "blastn", "echo 'blastn: 2.5.0+'\n")
    report = directory / "fastp.report.json"
    report.write_text(json.dumps(_FASTP_REPORT), encoding="utf-8")
    _script(
        directory / "fastp",
        "case \"$1\" in --version) echo 'fastp 1.4.0'; exit 0;; esac\n"
        + _argument_reader(
            [
                ("--in1", "in1"),
                ("--in2", "in2"),
                ("--out1", "out1"),
                ("--out2", "out2"),
                ("--json", "report"),
                ("--html", "html"),
            ]
        )
        + f'cp "$in1" "$out1"\ncp "$in2" "$out2"\ncp "{report}" "$report"\necho x > "$html"\n',
    )
    canned_sam = directory / "canned.sam"
    canned_sam.write_text(sam, encoding="utf-8")
    _script(
        directory / "minimap2",
        "case \"$1\" in --version) echo '2.31-r1302'; exit 0;; esac\n"
        + _argument_reader([("-o", "out")])
        + f'cp "{canned_sam}" "$out"\n',
    )
    return directory


def _gfa(segments: dict[str, str], links: Sequence[tuple[str, str, str, str]]) -> str:
    """A selected graph as GetOrganelle writes it, each link overlapping by k."""
    lines = [
        f"S\t{name}\t{bases}\tLN:i:{len(bases)}\tRC:i:{100 * len(bases)}"
        for name, bases in segments.items()
    ]
    lines.extend(
        f"L\t{a}\t{a_strand}\t{b}\t{b_strand}\t{_K}M" for a, a_strand, b, b_strand in links
    )
    return "\n".join(lines) + "\n"


def _selected_graph() -> str:
    """The graph ``SEQUENCE`` was read from: a segment holds the k bases it shares with the next."""
    lsc = SEQUENCE[:_LSC]
    ir = SEQUENCE[_LSC : _LSC + _IR]
    ssc = SEQUENCE[_LSC + _IR : _LSC + _IR + _SSC]
    mirror = SEQUENCE[_LSC + _IR + _SSC :]
    return _gfa(
        {"lsc": mirror[-_K:] + lsc + ir[:_K], "ir": ir, "ssc": ir[-_K:] + ssc + mirror[:_K]},
        [
            ("lsc", "+", "ir", "+"),
            ("ir", "+", "ssc", "+"),
            ("ssc", "+", "ir", "-"),
            ("ir", "-", "lsc", "+"),
        ],
    )


def _canned_complete(tools: Path) -> None:
    canned = tools / "canned"
    (canned / "get_org.log.txt").write_text(_LOG_COMPLETE, encoding="utf-8")
    prefix = "embplant_pt.K127.complete.graph1"
    (canned / f"{prefix}.1.path_sequence.fasta").write_text(
        f">lsc+,ir+,ssc+,ir-(circular)\n{SEQUENCE}\n", encoding="utf-8"
    )
    (canned / f"{prefix}.2.path_sequence.fasta").write_text(
        f">lsc+,ir+,ssc-,ir-(circular)\n{ISOMER}\n", encoding="utf-8"
    )
    (canned / f"{prefix}.selected_graph.gfa").write_text(_selected_graph(), encoding="utf-8")


def _canned_many_structures(tools: Path) -> None:
    """Four pieces between four direct copies of one repeat: six circles over one set of joins.

    The six have the length of ``SEQUENCE``, so the reads fit, and the reads cannot choose: every
    join of every circle is crossed. This is the shape of *Picea abies* in the published comparison.
    """
    source = random.Random(20261101)

    def random_bases(length: int) -> str:
        return "".join(source.choice("ACGT") for _ in range(length))

    repeat = random_bases(600)
    size = (LENGTH - 4 * len(repeat)) // 4
    pieces = {name: random_bases(size) for name in "abc"}
    pieces["d"] = random_bases(LENGTH - 4 * len(repeat) - 3 * size)
    segments = {name: repeat[-_K:] + bases + repeat[:_K] for name, bases in pieces.items()}
    segments["r"] = repeat
    links = [(name, "+", "r", "+") for name in pieces] + [("r", "+", name, "+") for name in pieces]
    canned = tools / "canned"
    (canned / "get_org.log.txt").write_text(_LOG_COMPLETE, encoding="utf-8")
    prefix = "embplant_pt.K127.complete.graph1"
    for number, order in enumerate(itertools.permutations("bcd"), start=1):
        walk = ",".join(f"{name}+,r+" for name in ("a", *order))
        sequence = "".join(pieces[name] + repeat for name in ("a", *order))
        (canned / f"{prefix}.{number}.path_sequence.fasta").write_text(
            f">{walk}(circular)\n{sequence}\n", encoding="utf-8"
        )
    (canned / f"{prefix}.selected_graph.gfa").write_text(_gfa(segments, links), encoding="utf-8")


def _canned_scaffolds(tools: Path) -> None:
    canned = tools / "canned"
    (canned / "get_org.log.txt").write_text(_LOG_SCAFFOLDS, encoding="utf-8")
    (canned / "embplant_pt.K91.scaffolds.graph1.1.path_sequence.fasta").write_text(
        f">scaffold_1\n{SEQUENCE[:9000]}\n>scaffold_2\n{SEQUENCE[9000:15000]}\n"
        f">scaffold_3\n{SEQUENCE[15000:20000]}\n",
        encoding="utf-8",
    )
    (canned / "embplant_pt.K91.contigs.graph1.selected_graph.gfa").write_text(
        "S\t1\tACGT\n", encoding="utf-8"
    )


def _database(root: Path) -> Path:
    for kind in ("SeedDatabase", "LabelDatabase"):
        folder = root / kind
        folder.mkdir(parents=True)
        lines = []
        for name in ("embplant_pt", "embplant_mt"):
            (folder / f"{name}.fasta").write_text(">s\nACGT\n", encoding="utf-8")
            digest = hashlib.sha256(b">s\nACGT\n").hexdigest()
            lines.append(f"{name}\t0.0.1\t{digest}\n")
        (folder / "VERSION").write_text("".join(lines), encoding="utf-8")
    return root


class Bench:
    """Reads, a database and stub tools in one directory, and the command line over them."""

    def __init__(self, root: Path, *, sam: str | None = None) -> None:
        self.root = root
        self.tools = _tools(root / "tools", sam=sam if sam is not None else _sam())
        self.database = _database(root / "database")
        self.reads_1 = root / "r_1.fastq.gz"
        self.reads_2 = root / "r_2.fastq.gz"
        self.reads_1.write_bytes(b"reads one")
        self.reads_2.write_bytes(b"reads two")
        self.out = root / "runs"

    @property
    def sets(self) -> list[str]:
        return ["--set", f"tools={self.tools}", "--set", f"database={self.database}"]

    @property
    def reads(self) -> list[str]:
        return ["-1", str(self.reads_1), "-2", str(self.reads_2), "--sample", "Sample one"]

    def run(self, *extra: str) -> int:
        return main(["assemble", "run", *self.reads, "-o", str(self.out), *self.sets, *extra])

    def run_dir(self) -> Path:
        [run] = list(self.out.glob("*/*/manifest.json"))
        return run.parent

    def table(self, run_dir: Path, name: str) -> list[dict]:
        payload = json.loads((run_dir / "01_assembly" / "result.json").read_text("utf-8"))
        return payload["tables"][name]["rows"]

    def manifest(self, run_dir: Path) -> dict:
        return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


@pytest.fixture
def bench(tmp_path: Path) -> Bench:
    return Bench(tmp_path)


def test_a_closed_supported_assembly_is_the_artifact_the_next_step_reads(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    run = bench.run_dir()
    manifest = bench.manifest(run)
    assert manifest["status"] == "ok"
    [artifact] = manifest["artifacts"]
    assert artifact["format"] == "fasta"
    assert artifact["topology"] == "circular"
    assert artifact["sample_id"] == "Sample one"
    text = (run / artifact["path"]).read_text(encoding="utf-8")
    assert text.splitlines()[0] == ">Sample_one"
    assert "".join(text.splitlines()[1:]) == SEQUENCE

    [summary] = bench.table(run, "assembly_summary")
    assert summary["status"] == "complete"
    assert summary["read_support"] == "supported"
    assert summary["graph_call"] == "resolved"
    assert summary["graph_reasons"] == "nothing"
    assert summary["backend"] == "local"
    assert summary["assembler"] == "GetOrganelle 1.7.7.1"
    assert int(summary["length_bp"]) == LENGTH
    assert int(summary["ir_bp"]) == _IR
    assert int(summary["lsc_bp"]) == _LSC
    assert int(summary["ssc_bp"]) == _SSC
    assert int(summary["n_paths"]) == 2


def test_both_isomers_are_kept_and_the_junctions_are_read_on_the_first(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    run = bench.run_dir()
    assert (run / "01_assembly" / "assembly.paths.fasta").is_file()
    rows = {row["junction"]: row for row in bench.table(run, "junction_support")}
    assert set(rows) == {"origin", "JLB", "JSB", "JSA", "JLA"}
    assert {row["call"] for row in rows.values()} == {"supported"}
    # JLB follows the LSC, whose last base is 20000 in GenBank coordinates.
    assert int(rows["JLB"]["after_base"]) == _LSC
    assert int(rows["JSB"]["after_base"]) == _LSC + _IR
    assert int(rows["JSA"]["after_base"]) == _LSC + _IR + _SSC
    # The boundary at the origin follows the last base of the circle.
    assert int(rows["origin"]["after_base"]) == LENGTH


def test_two_readings_of_the_repeat_that_disagree_are_named(bench: Bench) -> None:
    _canned_complete(bench.tools)
    log = bench.tools / "canned" / "get_org.log.txt"
    log.write_text(log.read_text().replace(f":{_SSC}:{_IR}\n", f":{_SSC}:{_IR // 2}\n"), "utf-8")
    assert bench.run() == 0
    payload = json.loads((bench.run_dir() / "01_assembly" / "result.json").read_text("utf-8"))
    assert any(f"{_IR} bp here and {_IR // 2} bp" in note for note in payload["notes"])


def test_the_tool_record_and_the_database_are_in_the_manifest(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    manifest = bench.manifest(bench.run_dir())
    module = manifest["modules"][0]
    assert [tool["name"] for tool in module["tools"]] == ["fastp", "getorganelle", "minimap2"]
    [seed] = [item for item in manifest["reference_sets"] if item["role"] == "assembly_seed"]
    assert {entry["database"] for entry in seed["databases"]} == {
        "SeedDatabase/embplant_pt",
        "SeedDatabase/embplant_mt",
        "LabelDatabase/embplant_pt",
        "LabelDatabase/embplant_mt",
    }
    argv = module["tools"][1]["argv"]
    assert argv[argv.index("-k") + 1] == "21,55,91,127"
    assert argv[argv.index("-R") + 1] == "10"
    assert argv[argv.index("--config-dir") + 1] == str(bench.database)
    assert any(part.startswith("-m ") for part in argv)
    assert len(manifest["inputs"]) == 2


def test_the_map_back_reads_the_raw_reads_and_the_assembler_the_trimmed_ones(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    tools = {tool["name"]: tool for tool in bench.manifest(bench.run_dir())["modules"][0]["tools"]}
    assembler = tools["getorganelle"]["argv"]
    assert "reads_1.fastq.gz" in Path(assembler[assembler.index("-1") + 1]).name
    assert "work" in Path(assembler[assembler.index("-1") + 1]).parts
    assert tools["minimap2"]["argv"][-2:] == [str(bench.reads_1), str(bench.reads_2)]


def test_the_working_files_do_not_stay_unless_asked(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    assert not (bench.run_dir() / "01_assembly" / "work").exists()


def test_comparison_names_the_isomer_that_is_the_reference(bench: Bench) -> None:
    _canned_complete(bench.tools)
    reference = bench.root / "reference.fasta"
    # The reference is the second isomer, on another origin and the other strand.
    rotated = reverse_complement(ISOMER[777:] + ISOMER[:777])
    reference.write_text(f">ref\n{rotated}\n", encoding="utf-8")
    assert bench.run("--reference", str(reference)) == 0
    rows = {row["path"]: row for row in bench.table(bench.run_dir(), "reference_comparison")}
    second = "embplant_pt.K127.complete.graph1.2"
    assert rows[second]["identical"] == "yes"
    assert rows["embplant_pt.K127.complete.graph1.1"]["identical"] == "no"


def test_a_disagreement_with_the_reference_measures_both_sides(bench: Bench) -> None:
    _canned_complete(bench.tools)
    reference = bench.root / "reference.fasta"
    # One base away from both paths, so no path equals it. The stub SAM reads the same for both
    # sequences, so the two measurements must come out equal; the point is that both exist.
    changed = SEQUENCE[:5_000] + ("A" if SEQUENCE[5_000] != "A" else "C") + SEQUENCE[5_001:]
    reference.write_text(f">ref\n{changed}\n", encoding="utf-8")
    assert bench.run("--reference", str(reference)) == 0
    run = bench.run_dir()
    rows = bench.table(run, "reference_read_support")
    assert [row["sequence"] for row in rows] == ["assembly", "ref"]
    assert {row["read_support"] for row in rows} == {"supported"}
    assert {row["junction"] for row in bench.table(run, "reference_junctions")} >= {"origin", "JLB"}
    minimap2_calls = [
        tool for tool in bench.manifest(run)["modules"][0]["tools"] if tool["name"] == "minimap2"
    ]
    assert len(minimap2_calls) == 2


def test_an_assembly_identical_to_the_reference_is_not_measured_twice(bench: Bench) -> None:
    _canned_complete(bench.tools)
    reference = bench.root / "reference.fasta"
    reference.write_text(f">ref\n{SEQUENCE}\n", encoding="utf-8")
    assert bench.run("--reference", str(reference)) == 0
    payload = json.loads((bench.run_dir() / "01_assembly" / "result.json").read_text("utf-8"))
    assert "reference_read_support" not in payload["tables"]


def test_a_join_the_reads_do_not_cross_is_undetermined_and_no_artifact(tmp_path: Path) -> None:
    bench = Bench(tmp_path, sam=_sam(skip_between=(9_000, 9_100)))
    _canned_complete(bench.tools)
    assert bench.run() == 1
    run = bench.run_dir()
    manifest = bench.manifest(run)
    assert manifest["status"] == "incomplete"
    assert manifest["artifacts"] == []
    [summary] = bench.table(run, "assembly_summary")
    assert summary["status"] == "complete"
    assert summary["read_support"] == "undetermined"
    assert int(summary["fewest_spanning"]) == 0
    [region] = bench.table(run, "weak_regions")
    assert int(region["first_after_base"]) <= 9_000 <= int(region["last_after_base"])
    # The sequence is kept for a person to read; it is only not declared for the next step.
    assert (run / "01_assembly" / "assembly.fasta").is_file()


def test_a_junction_with_too_few_reads_is_undetermined(tmp_path: Path) -> None:
    bench = Bench(tmp_path, sam=_sam(depth=1))
    _canned_complete(bench.tools)
    assert bench.run("--set", "min_spanning=100000") == 1
    [summary] = bench.table(bench.run_dir(), "assembly_summary")
    assert summary["read_support"] == "undetermined"


def test_a_fragmented_assembly_is_named_and_its_support_is_not_measured(bench: Bench) -> None:
    _canned_scaffolds(bench.tools)
    assert bench.run() == 1
    run = bench.run_dir()
    manifest = bench.manifest(run)
    assert manifest["status"] == "incomplete"
    assert manifest["artifacts"] == []
    [summary] = bench.table(run, "assembly_summary")
    assert summary["status"] == "fragmented"
    assert summary["assembler_says"] == "3 scaffold(s)"
    assert summary["read_support"] == "not measured"
    assert summary["mean_depth"] == "not measured"
    assert not (run / "01_assembly" / "assembly.fasta").exists()
    scaffolds = (run / "01_assembly" / "scaffolds.fasta").read_text(encoding="utf-8")
    assert scaffolds.count(">") == 3
    assert ">Sample_one_2\n" in scaffolds


def test_an_open_assembly_still_measures_the_reference_it_was_given(bench: Bench) -> None:
    _canned_scaffolds(bench.tools)
    reference = bench.root / "reference.fasta"
    reference.write_text(f">ref\n{SEQUENCE}\n", encoding="utf-8")
    assert bench.run("--reference", str(reference)) == 1
    rows = {row["sequence"]: row for row in bench.table(bench.run_dir(), "reference_read_support")}
    assert rows["assembly"]["read_support"] == "not measured"
    assert rows["assembly"]["mean_depth"] == "not measured"
    assert rows["ref"]["read_support"] == "supported"


def test_no_result_is_stated_and_is_not_a_measured_zero(bench: Bench) -> None:
    (bench.tools / "canned" / "get_org.log.txt").write_text(_LOG_NOTHING, encoding="utf-8")
    assert bench.run() == 1
    run = bench.run_dir()
    [summary] = bench.table(run, "assembly_summary")
    assert summary["status"] == "no_result"
    assert summary["length_bp"] == "not measured"
    assert int(summary["n_paths"]) == 0


def test_a_missing_tool_is_named_with_the_route_and_starts_no_run(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    (bench.tools / "spades.py").unlink()
    assert bench.run() == 3
    error = capsys.readouterr().err
    assert "cannot run here" in error
    assert "spades.py not found" in error
    assert "plastid assemble prepare" in error
    assert "--backend import" in error
    assert "nothing was sent" in error
    assert not bench.out.exists()


def test_a_failure_the_assembler_exits_zero_on_is_named_with_the_route(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    # GetOrganelle catches every failure and exits 0, so only its log says SPAdes died.
    (bench.tools / "canned" / "get_org.log.txt").write_text(
        _LOG_HEAD + "2026-10-31 09:48:11,361 - ERROR: Assembling failed.\n", encoding="utf-8"
    )
    assert bench.run() == 3
    error = capsys.readouterr().err
    assert "Assembling failed" in error
    assert "plastid assemble prepare" in error
    manifest = bench.manifest(bench.run_dir())
    assert manifest["status"] == "failed"
    assert manifest["artifacts"] == []


def test_a_machine_below_the_floor_stops_auto_but_not_local(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    assert bench.run("--set", "min_ram_gb=100000") == 3
    assert "the machine has" in capsys.readouterr().err
    assert bench.run("--set", "min_ram_gb=100000", "--backend", "local") == 0


def test_a_disk_too_small_for_the_input_stops_auto(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    assert bench.run("--set", "disk_factor=100000000000000000000") == 3
    assert "of disk" in capsys.readouterr().err


def test_a_misspelled_parameter_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    assert bench.run("--set", "rounds_max=3") == 2
    assert "assemble does not take it" in capsys.readouterr().err
    assert not bench.out.exists()


def test_a_label_reaches_every_table_scope(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run("--label", "third-party reads") == 0
    path = bench.run_dir() / "01_assembly" / "tables" / "assembly_summary.csv"
    head = path.read_text(encoding="utf-8")
    assert "third-party reads" in head.split("\n# table", 1)[0].split("\n", 2)[1]


# --- the import backend --------------------------------------------------------------------------


def _prepare(bench: Bench) -> Path:
    code = main(["assemble", "prepare", *bench.reads, "-o", str(bench.root / "p"), *bench.sets])
    assert code == 0
    [prepared] = list((bench.root / "p").glob("*/*/01_assembly/prepared.json"))
    return prepared


def _import(bench: Bench, prepared: Path, result: Path, *extra: str) -> int:
    return main(
        [
            "assemble",
            "run",
            "--backend",
            "import",
            "--prepared",
            str(prepared),
            "--result",
            str(result),
            "--source-name",
            "usegalaxy.eu",
            "--retrieved-at",
            "2026-10-31",
            "-o",
            str(bench.out),
            *bench.sets,
            *extra,
        ]
    )


def test_prepare_hashes_what_it_hands_over_and_names_no_upload(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    prepared = _prepare(bench)
    payload = json.loads(prepared.read_text(encoding="utf-8"))
    for entry in payload["input_reads"]:
        file = prepared.parent / entry["path"]
        assert hashlib.sha256(file.read_bytes()).hexdigest() == entry["sha256"]
    assert payload["tool"]["version"] == getorganelle.GETORGANELLE_VERSION
    assert payload["tool"]["command"][0] == "get_organelle_from_reads.py"
    assert "--backend import" in capsys.readouterr().out
    assert not (prepared.parent / "work").exists()


def test_both_backends_give_one_result_and_differ_only_in_provenance(bench: Bench) -> None:
    _canned_complete(bench.tools)
    assert bench.run() == 0
    local = bench.run_dir()
    prepared = _prepare(bench)
    # A second output root, so the imported run is the only one found there.
    bench.out = bench.root / "imported"
    assert _import(bench, prepared, bench.tools / "canned") == 0
    imported = bench.run_dir()

    for name in ("junction_support", "assembly_paths", "read_qc"):
        assert bench.table(local, name) == bench.table(imported, name)
    first, second = (
        bench.table(local, "assembly_summary")[0],
        bench.table(imported, "assembly_summary")[0],
    )
    differing = {key for key in first if first[key] != second[key]}
    assert differing == {"backend", "where"}
    assert second["backend"] == "import"
    assert second["where"] == "imported: usegalaxy.eu, retrieved 2026-10-31"
    for file in ("assembly.fasta", "assembly.paths.fasta"):
        assert (local / "01_assembly" / file).read_bytes() == (
            imported / "01_assembly" / file
        ).read_bytes()


def test_an_import_records_the_external_result_contract(bench: Bench) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 0
    run = bench.run_dir()
    described = json.loads((run / "01_assembly" / "external_result.json").read_text("utf-8"))
    assert described["tool"]["name"] == "getorganelle"
    assert described["tool"]["version"] == "1.7.7.1"
    assert described["source"] == {"name": "usegalaxy.eu", "retrieved_at": "2026-10-31"}
    assert described["status"] == "ok"
    assert described["count"] == 2
    assert described["input"]["path"] == str(prepared)


def test_an_import_says_which_programs_it_ran_when_they_differ_from_the_pins(bench: Bench) -> None:
    # usegalaxy.eu ran SPAdes 3.15.5 under the same GetOrganelle: a real import, measured.
    _canned_complete(bench.tools)
    log = bench.tools / "canned" / "get_org.log.txt"
    log.write_text(log.read_text().replace("SPAdes 4.3.0", "SPAdes 3.15.5"), encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 0
    run = bench.run_dir()
    [summary] = bench.table(run, "assembly_summary")
    assert "SPAdes 3.15.5" in summary["assembler_programs"]
    assert summary["seed_database"] == "embplant_pt 0.0.1; embplant_mt 0.0.1"
    payload = json.loads((run / "01_assembly" / "result.json").read_text("utf-8"))
    assert any("SPAdes 3.15.5; the local route pins 4.3.0" in note for note in payload["notes"])
    assert "SPAdes 3.15.5" in (run / "01_assembly" / "report.md").read_text("utf-8")


def test_an_import_on_the_pinned_programs_carries_no_such_note(bench: Bench) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 0
    payload = json.loads((bench.run_dir() / "01_assembly" / "result.json").read_text("utf-8"))
    assert not any("the local route pins" in note for note in payload["notes"])


def test_a_path_its_graph_cannot_spell_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    # A person assembling a result directory by hand can pair a path with another run's graph.
    _canned_complete(bench.tools)
    path = bench.tools / "canned" / "embplant_pt.K127.complete.graph1.1.path_sequence.fasta"
    changed = "A" if SEQUENCE[5_000] != "A" else "C"
    edited = SEQUENCE[:5_000] + changed + SEQUENCE[5_001:]
    path.write_text(f">lsc+,ir+,ssc+,ir-(circular)\n{edited}\n", encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "not the sequence its walk spells" in capsys.readouterr().err


def test_a_walk_over_a_link_the_graph_lacks_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    path = bench.tools / "canned" / "embplant_pt.K127.complete.graph1.1.path_sequence.fasta"
    path.write_text(f">lsc+,ssc+,ir+,ir-(circular)\n{SEQUENCE}\n", encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "not a link" in capsys.readouterr().err


def test_a_graph_of_many_structures_is_not_passed_on_the_reads(bench: Bench) -> None:
    _canned_many_structures(bench.tools)
    assert bench.run() == 1
    run = bench.run_dir()
    [summary] = bench.table(run, "assembly_summary")
    assert summary["status"] == "complete"
    assert summary["read_support"] == "supported"  # the reads have nothing to say here
    assert summary["graph_call"] == "unresolved"
    assert "6 structures close through it" in summary["graph_reasons"]
    manifest = bench.manifest(run)
    assert manifest["status"] == "incomplete"
    assert manifest["artifacts"] == []
    assert manifest["modules"][0]["graph"] == "unresolved"
    [graph] = bench.table(run, "assembly_graph")
    assert int(graph["structures"]) == 6
    assert int(graph["paths_written"]) == 6
    # Every structure stays for a person to read; none is offered as the answer.
    text = (run / "01_assembly" / "assembly.paths.fasta").read_text(encoding="utf-8")
    assert text.count(">") == 6


def test_an_import_without_its_selected_graph_gives_no_artifact(bench: Bench) -> None:
    _canned_complete(bench.tools)
    (bench.tools / "canned" / "embplant_pt.K127.complete.graph1.selected_graph.gfa").unlink()
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 1
    run = bench.run_dir()
    [summary] = bench.table(run, "assembly_summary")
    assert summary["read_support"] == "supported"
    assert summary["graph_call"] == "not measured"
    assert summary["graph_reasons"] == "no selected graph to read"
    assert bench.manifest(run)["artifacts"] == []


def test_a_path_that_does_not_state_its_walk_is_not_measured(bench: Bench) -> None:
    _canned_complete(bench.tools)
    canned = bench.tools / "canned"
    path = canned / "embplant_pt.K127.complete.graph1.1.path_sequence.fasta"
    path.write_text(f">renamed by hand\n{SEQUENCE}\n", encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, canned) == 1
    run = bench.run_dir()
    [summary] = bench.table(run, "assembly_summary")
    assert summary["graph_call"] == "not measured"
    assert "does not state its walk" in summary["graph_reasons"]
    assert bench.manifest(run)["artifacts"] == []


def test_paths_of_a_graph_that_is_not_there_are_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    canned = bench.tools / "canned"
    extra = canned / "embplant_pt.K127.complete.graph2.1.path_sequence.fasta"
    extra.write_text(f">lsc+,ir+,ssc+,ir-(circular)\n{SEQUENCE}\n", encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, canned) == 1
    assert "no selected graph for its paths: graph" in capsys.readouterr().err


def test_a_graph_no_path_came_from_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    canned = bench.tools / "canned"
    source = canned / "embplant_pt.K127.complete.graph1.selected_graph.gfa"
    source.rename(canned / "embplant_pt.K127.complete.graph2.selected_graph.gfa")
    prepared = _prepare(bench)
    assert _import(bench, prepared, canned) == 1
    assert "no selected graph for its paths: graph" in capsys.readouterr().err


def test_two_selected_graphs_are_not_one_molecule(bench: Bench) -> None:
    _canned_complete(bench.tools)
    canned = bench.tools / "canned"
    for name in (".selected_graph.gfa", ".1.path_sequence.fasta", ".2.path_sequence.fasta"):
        source = canned / f"embplant_pt.K127.complete.graph1{name}"
        (canned / f"embplant_pt.K127.complete.graph2{name}").write_bytes(source.read_bytes())
    assert bench.run() == 1
    [summary] = bench.table(bench.run_dir(), "assembly_summary")
    assert summary["graph_call"] == "unresolved"
    assert "2 graphs were selected" in summary["graph_reasons"]


def test_an_import_maps_the_raw_reads_too(bench: Bench) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 0
    [tool] = bench.manifest(bench.run_dir())["modules"][0]["tools"]
    assert tool["name"] == "minimap2"
    assert tool["argv"][-2:] == [str(bench.reads_1.resolve()), str(bench.reads_2.resolve())]


def test_an_import_of_changed_raw_reads_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    bench.reads_1.write_bytes(b"other reads")
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "not the input the imported result was asked for" in capsys.readouterr().err


def test_an_import_whose_raw_reads_are_gone_says_so(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    bench.reads_2.unlink()
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "are gone" in capsys.readouterr().err


def test_an_import_needs_minimap2_and_says_so(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    (bench.tools / "minimap2").unlink()
    assert _import(bench, prepared, bench.tools / "canned") == 3
    assert "minimap2 not found" in capsys.readouterr().err


def test_prepared_reads_outside_the_run_are_named_absolutely(
    bench: Bench, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(bench.root)
    relative = ["-1", "r_1.fastq.gz", "-2", "r_2.fastq.gz", "--sample", "s", "-o", "p"]
    assert main(["assemble", "prepare", *relative, "--set", "trim=no", *bench.sets]) == 0
    [prepared] = list((bench.root / "p").glob("*/*/01_assembly/prepared.json"))
    monkeypatch.chdir(bench.root.parent)
    payload = assembly.load_prepared(prepared)
    assert Path(payload["input_reads"][0]["path"]).is_absolute()


def test_a_failing_quality_control_stops_with_the_route(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _script(
        bench.tools / "fastp",
        "case \"$1\" in --version) echo 'fastp 1.4.0'; exit 0;; esac\nexit 9\n",
    )
    assert bench.run() == 3
    error = capsys.readouterr().err
    assert "fastp exited 9" in error
    assert "plastid assemble prepare" in error


def test_an_import_of_changed_reads_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    payload = json.loads(prepared.read_text(encoding="utf-8"))
    (prepared.parent / payload["input_reads"][0]["path"]).write_bytes(b"other reads")
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "not the input the imported result was asked for" in capsys.readouterr().err


def test_an_import_with_more_reads_than_prepared_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    log = bench.tools / "canned" / "get_org.log.txt"
    log.write_text(log.read_text().replace("295+295", "3000+3000"), encoding="utf-8")
    prepared = _prepare(bench)
    assert _import(bench, prepared, bench.tools / "canned") == 1
    assert "it did not run on the prepared input" in capsys.readouterr().err


@pytest.mark.parametrize("missing", ["--prepared", "--result", "--source-name", "--retrieved-at"])
def test_an_import_without_its_provenance_is_refused(
    bench: Bench, missing: str, capsys: pytest.CaptureFixture[str]
) -> None:
    values = {
        "--prepared": "p.json",
        "--result": "dir",
        "--source-name": "x",
        "--retrieved-at": "2026-10-31",
    }
    del values[missing]
    argv = ["assemble", "run", "--backend", "import", "-o", str(bench.out)]
    for flag, value in values.items():
        argv += [flag, value]
    assert main(argv) == 2
    assert missing in capsys.readouterr().err


def test_an_import_with_a_bad_date_is_refused_by_the_contract(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _canned_complete(bench.tools)
    prepared = _prepare(bench)
    code = main(
        [
            "assemble",
            "run",
            "--backend",
            "import",
            "--prepared",
            str(prepared),
            "--result",
            str(bench.tools / "canned"),
            "--source-name",
            "usegalaxy.eu",
            "--retrieved-at",
            "last week",
            "-o",
            str(bench.out),
        ]
    )
    assert code == 1
    assert "not an ISO date" in capsys.readouterr().err


# --- the adapters --------------------------------------------------------------------------------


def test_the_log_states_the_result_the_versions_and_the_database() -> None:
    facts = getorganelle.parse_log(_LOG_COMPLETE)
    assert facts.status_text == "circular genome"
    assert facts.version == "1.7.7.1"
    assert facts.dependencies["SPAdes"] == "4.3.0"
    assert facts.seed_db["embplant_pt"] == "0.0.1"
    assert facts.base_coverage == 141.562
    assert facts.reads_used == 590
    assert facts.architecture == ("IRs detected", LENGTH, _LSC, _SSC, _IR)


def test_a_log_and_path_files_that_disagree_are_refused(bench: Bench, tmp_path: Path) -> None:
    _canned_complete(bench.tools)
    log = bench.tools / "canned" / "get_org.log.txt"
    log.write_text(_LOG_SCAFFOLDS, encoding="utf-8")
    with pytest.raises(AssemblyError, match="cannot be read from this directory"):
        getorganelle.read_output(bench.tools / "canned")


def test_a_directory_without_a_log_is_not_an_assembly(tmp_path: Path) -> None:
    with pytest.raises(AssemblyError, match="not a GetOrganelle output directory"):
        getorganelle.read_output(tmp_path)


def test_a_status_this_step_does_not_know_is_refused(tmp_path: Path) -> None:
    (tmp_path / "get_org.log.txt").write_text(
        "INFO: Result status of embplant_pt: something new\n", encoding="utf-8"
    )
    with pytest.raises(AssemblyError, match="does not know"):
        getorganelle.read_output(tmp_path)


def test_a_wrong_assembler_version_is_refused(bench: Bench) -> None:
    path = bench.tools / "get_organelle_from_reads.py"
    path.write_text(path.read_text().replace("1.7.7.1", "1.6.2"), encoding="utf-8")
    with pytest.raises(ToolError, match="expected GetOrganelle"):
        getorganelle.resolve_tools(bench.tools)


def test_a_wrong_spades_version_is_refused(bench: Bench) -> None:
    path = bench.tools / "spades.py"
    path.write_text(path.read_text().replace("4.3.0", "3.15.5"), encoding="utf-8")
    with pytest.raises(ToolError, match="expected SPAdes"):
        getorganelle.resolve_tools(bench.tools)


def test_a_database_that_is_not_its_recorded_bytes_is_refused(bench: Bench) -> None:
    (bench.database / "SeedDatabase" / "embplant_pt.fasta").write_text(">x\nA\n", encoding="utf-8")
    with pytest.raises(ToolError, match="is not the"):
        getorganelle.resolve_database(bench.database)


def test_a_missing_database_says_how_to_fetch_it(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="get_organelle_config.py -a embplant_pt,embplant_mt"):
        getorganelle.resolve_database(tmp_path)


def test_a_database_without_the_mitochondrial_labels_is_refused(bench: Bench) -> None:
    version = bench.database / "LabelDatabase" / "VERSION"
    version.write_text(version.read_text().splitlines()[0] + "\n", encoding="utf-8")
    with pytest.raises(ToolError, match="no embplant_mt entry"):
        getorganelle.resolve_database(bench.database)


def test_a_wrong_fastp_and_a_wrong_minimap2_are_refused(bench: Bench) -> None:
    for name, old, new in (("fastp", "1.4.0", "0.23.4"), ("minimap2", "2.31-r1302", "2.28-r1209")):
        path = bench.tools / name
        path.write_text(path.read_text().replace(old, new), encoding="utf-8")
    with pytest.raises(ToolError, match="expected fastp"):
        fastp.fastp_version(fastp.resolve_fastp(bench.tools))
    with pytest.raises(ToolError, match="expected minimap2"):
        minimap2.minimap2_version(minimap2.resolve_minimap2(bench.tools))


def test_a_fastp_report_without_its_counts_is_refused() -> None:
    with pytest.raises(ToolError, match="not a fastp JSON report"):
        fastp.parse_report('{"summary": {}}', source="report.json")


def test_kmers_may_be_written_as_a_list_on_the_command_line() -> None:
    assert assembly.resolve_params({"kmers": "21,45,65"})["kmers"] == (21, 45, 65)
    assert assembly.resolve_params({})["kmers"] == (21, 55, 91, 127)
