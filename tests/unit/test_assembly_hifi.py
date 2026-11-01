"""The assembly step on HiFi reads, through the command line, with Oatk stubbed at its boundary.

The stubs write what the real tools write: Oatk's contig and graph with its uneven overlaps, and the
`[E::...]` line it ends on when it stops. What is under test is what the step decides from them. The
real tools run in `tests/golden/test_assembly_long_panel.py`.
"""

from __future__ import annotations

import gzip
import json
import random
import stat
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.core.sequence import reverse_complement

_LSC, _IR, _SSC = 20_000, 11_000, 5_000
_READ = 12_000
_EXTENSION = 100_000


def _sequence() -> tuple[str, str, str, str]:
    source = random.Random(20261031)

    def random_bases(length: int) -> str:
        return "".join(source.choice("ACGT") for _ in range(length))

    # Ends known to differ, so the repeat is exactly the length it was built with.
    lsc = "A" + random_bases(_LSC - 2) + "A"
    ir = random_bases(_IR)
    ssc = "A" + random_bases(_SSC - 2) + "A"
    return lsc + ir + ssc + reverse_complement(ir), lsc, ir, ssc


SEQUENCE, LSC, IR, SSC = _sequence()
LENGTH = len(SEQUENCE)
ISOMER = LSC + IR + reverse_complement(SSC) + reverse_complement(IR)


def _gfa() -> str:
    """Oatk's graph of ``SEQUENCE``: overlaps that differ link by link, and coverage as ``SC``."""
    mirror = reverse_complement(IR)
    u1 = mirror[-900:] + LSC + IR[:1100]
    u2 = IR
    u3 = IR[-700:] + SSC + mirror[:1300]
    return (
        "H\tVN:Z:1.0\n"
        f"S\tu1\t{u1}\tLN:i:{len(u1)}\tKC:i:{100 * len(u1)}\tSC:f:100.000\n"
        f"S\tu2\t{u2}\tLN:i:{len(u2)}\tKC:i:{200 * len(u2)}\tSC:f:200.000\n"
        f"S\tu3\t{u3}\tLN:i:{len(u3)}\tKC:i:{100 * len(u3)}\tSC:f:100.000\n"
        "L\tu1\t+\tu2\t+\t1100M\tEC:i:90\n"
        "L\tu2\t+\tu3\t+\t700M\tEC:i:90\n"
        "L\tu3\t+\tu2\t-\t1300M\tEC:i:90\n"
        "L\tu2\t-\tu1\t+\t900M\tEC:i:90\n"
    )


def _contig(path: str, sequence: str, *, circular: bool = True, name: str = "ctg000001c") -> str:
    flag = "true" if circular else "false"
    fields = f"length={len(sequence)} wlength=1.0 nv=4 circular={flag} path={path}"
    return f">{name}\t{fields}\n{sequence}\n"


def _sam() -> str:
    """Reads of 12 kb every 200 bp on the extended sequence, each placed a few times."""
    lines = [f"@SQ\tSN:target\tLN:{LENGTH + LENGTH}\n"]
    for start in range(0, LENGTH + LENGTH - _READ, 200):
        for _ in range(4):
            lines.append(f"r\t0\ttarget\t{start + 1}\t60\t{_READ}=\t*\t0\t0\t*\t*\n")
    return "".join(lines)


def _script(path: Path, body: str) -> None:
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


class Bench:
    """Reads, a profile database and stub tools in one directory, and the command line over them."""

    def __init__(self, root: Path, *, oatk_body: str | None = None) -> None:
        self.root = root
        self.tools = root / "tools"
        self.tools.mkdir()
        self.canned = self.tools / "canned"
        self.canned.mkdir()
        self.canned_files(_contig("u1+,u2+,u3+,u2-", SEQUENCE), _gfa())
        success = (
            f'cp "{self.canned}"/oatk.* "$(dirname "$out")"/ ; '
            'echo "[M::syncasm] collected syncmers from 590 target sequence(s)" >&2\n'
            'echo "[M::main] Version: 1.0" >&2\n'
        )
        _script(
            self.tools / "oatk",
            'case "$1" in --version) echo 1.0; exit 0;; esac\n'
            # As the real tool does, from the directory it runs in.
            "for last; do :; done\n"
            '[ -f "$last" ] || { echo "[E::kseq] fail to open" >&2; exit 1; }\n'
            'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out="$2"; shift; done\n'
            + (oatk_body if oatk_body is not None else success),
        )
        _script(self.tools / "nhmmscan", 'echo "# HMMER 3.4 (Aug 2023); http://hmmer.org/"\n')
        sam = self.tools / "canned.sam"
        sam.write_text(_sam(), encoding="utf-8")
        _script(
            self.tools / "minimap2",
            "case \"$1\" in --version) echo '2.31-r1302'; exit 0;; esac\n"
            'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out="$2"; shift; done\n'
            f'cp "{sam}" "$out"\n',
        )
        self.database = root / "database"
        self.database.mkdir()
        for name in ("embryophyta_pltd.fam", "embryophyta_mito.fam"):
            for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
                (self.database / f"{name}{suffix}").write_bytes(b"profile")
        self.reads = root / "hifi.fastq.gz"
        self.reads.write_bytes(gzip.compress(f"@r\n{'A' * 50}\n+\n{'I' * 50}\n".encode()))
        self.out = root / "runs"

    def canned_files(self, contigs: str, gfa: str) -> None:
        (self.canned / "oatk.pltd.ctg.fasta").write_text(contigs, encoding="utf-8")
        (self.canned / "oatk.pltd.gfa").write_text(gfa, encoding="utf-8")
        (self.canned / "oatk.utg.final.gfa").write_text(gfa, encoding="utf-8")

    @property
    def sets(self) -> list[str]:
        return ["--set", f"tools={self.tools}", "--set", f"oatk_db={self.database}"]

    @property
    def given(self) -> list[str]:
        return ["--hifi", str(self.reads), "--sample", "Sample one"]

    def run(self, *extra: str) -> int:
        return main(["assemble", "run", *self.given, "-o", str(self.out), *self.sets, *extra])

    def run_dir(self) -> Path:
        [run] = list(self.out.glob("*/*/manifest.json"))
        return run.parent

    def table(self, name: str) -> list[dict]:
        payload = json.loads((self.run_dir() / "01_assembly" / "result.json").read_text("utf-8"))
        return payload["tables"][name]["rows"]

    def manifest(self) -> dict:
        return json.loads((self.run_dir() / "manifest.json").read_text(encoding="utf-8"))


@pytest.fixture
def bench(tmp_path: Path) -> Bench:
    return Bench(tmp_path)


def test_hifi_reads_make_the_artifact_the_next_step_reads(bench: Bench) -> None:
    assert bench.run() == 0
    [summary] = bench.table("assembly_summary")
    assert summary["assembler"] == "Oatk 1.0"
    assert summary["assembler_programs"] == "HMMER 3.4"
    assert summary["seed_database"] == "OatkDB v20230921"
    assert summary["status"] == "complete"
    assert summary["read_support"] == "supported"
    assert summary["graph_call"] == "resolved"
    assert summary["word_limit"] == "not applicable"
    assert int(summary["reads_used"]) == 590
    assert int(summary["ir_bp"]) == _IR
    [artifact] = bench.manifest()["artifacts"]
    text = (bench.run_dir() / artifact["path"]).read_text(encoding="utf-8")
    assert "".join(text.splitlines()[1:]) == SEQUENCE


def test_the_second_isomer_comes_from_the_graph(bench: Bench) -> None:
    assert bench.run() == 0
    paths = {row["path"]: row for row in bench.table("assembly_paths")}
    assert set(paths) == {"ctg000001c", "ctg000001c.isomer"}
    [graph] = bench.table("assembly_graph")
    assert (graph["structures"], graph["paths_written"], graph["isomers"]) == (2, 2, "yes")
    # The repeat comes out of the walk, carried past the segment's ends by the bases.
    assert graph["repeat"] == f"{_IR} bp, inverted"
    text = (bench.run_dir() / "01_assembly" / "assembly.paths.fasta").read_text(encoding="utf-8")
    assert text.count(">") == 2


def test_the_command_lines_are_the_measured_ones(bench: Bench) -> None:
    assert bench.run("--set", "threads=2") == 0
    module = bench.manifest()["modules"][0]
    oatk, mapper = module["tools"][0]["argv"], module["tools"][1]["argv"]
    assert oatk[oatk.index("-k") + 1] == "1001"
    assert oatk[oatk.index("-c") + 1] == "30"
    assert oatk[oatk.index("-D") + 1] == "1500000000"
    assert oatk[oatk.index("-t") + 1] == "2"
    assert oatk[oatk.index("-m") + 1] == str(bench.database / "embryophyta_mito.fam")
    assert oatk[oatk.index("-p") + 1] == str(bench.database / "embryophyta_pltd.fam")
    assert mapper[mapper.index("-ax") + 1] == "map-hifi"
    assert [tool["name"] for tool in module["tools"]] == ["oatk", "minimap2"]


def test_the_profiles_are_recorded_by_digest(bench: Bench) -> None:
    assert bench.run() == 0
    [record] = bench.manifest()["reference_sets"]
    assert record["format"] == "oatk_database"
    assert {entry["database"] for entry in record["databases"]} == {
        "embryophyta_pltd.fam",
        "embryophyta_mito.fam",
    }
    assert all(entry["version"] == "v20230921" for entry in record["databases"])


def test_the_reason_oatk_stopped_reaches_the_person(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # It stands at the end of stderr, past the excerpt a tool error carries.
    padding = "x" * 400
    body = (
        f'echo "[M::syncasm] {padding}" >&2\n'
        'echo "[E::syncasm] empty syncmer graph" >&2\n'
        'echo "[E::main] syncasm assembly program failed" >&2\n'
        "exit 1\n"
    )
    bench = Bench(tmp_path, oatk_body=body)
    assert bench.run() == 3
    err = capsys.readouterr().err
    assert "empty syncmer graph" in err
    assert "plastid assemble prepare --hifi" in err


def test_a_missing_tool_stops_before_anything_runs(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    (bench.tools / "oatk").unlink()
    assert bench.run() == 3
    err = capsys.readouterr().err
    assert "oatk not found" in err
    assert "docs/tools/oatk.md" in err
    assert "plastid assemble prepare --hifi" in err
    assert not bench.out.exists()


def test_another_version_of_oatk_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _script(bench.tools / "oatk", "echo 0.9\n")
    assert bench.run() == 3
    assert "expected Oatk 1.0" in capsys.readouterr().err


def test_another_version_of_hmmer_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    _script(bench.tools / "nhmmscan", 'echo "# HMMER 3.3 (Nov 2019)"\n')
    assert bench.run() == 3
    assert "expected HMMER 3.4" in capsys.readouterr().err


def test_profiles_that_were_not_pressed_are_named(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    (bench.database / "embryophyta_pltd.fam.h3m").unlink()
    assert bench.run() == 3
    assert "embryophyta_pltd.fam.h3m" in capsys.readouterr().err


def test_no_database_is_a_named_stop(
    bench: Bench, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("OATK_DB", raising=False)
    code = main(
        ["assemble", "run", *bench.given, "-o", str(bench.out), "--set", f"tools={bench.tools}"]
    )
    assert code == 3
    assert "no Oatk database" in capsys.readouterr().err


def test_reads_named_from_where_the_command_was_typed_are_found(
    bench: Bench, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(bench.root)
    reads = ["--hifi", "hifi.fastq.gz", "--sample", "s"]
    sets = ["--set", "tools=tools", "--set", "oatk_db=database"]
    assert main(["assemble", "run", *reads, "-o", "relative", *sets]) == 0


def test_hifi_replaces_the_short_reads(capsys: pytest.CaptureFixture[str], bench: Bench) -> None:
    code = main(["assemble", "run", *bench.given, "-1", "r1", "-2", "r2", *bench.sets])
    assert code == 2
    assert "--hifi replaces -1 and -2" in capsys.readouterr().err


def test_hifi_reads_need_a_sample_name(capsys: pytest.CaptureFixture[str], bench: Bench) -> None:
    assert main(["assemble", "run", "--hifi", str(bench.reads), *bench.sets]) == 2
    assert "--sample required" in capsys.readouterr().err


def test_a_run_that_does_not_close_is_no_artifact(bench: Bench) -> None:
    pieces = _contig("u1+", SEQUENCE[:20_000], circular=False, name="ctg000001l")
    pieces += _contig("u3+", SEQUENCE[20_000:30_000], circular=False, name="ctg000002l")
    bench.canned_files(pieces, _gfa().split("L\t")[0])
    assert bench.run() == 1
    [summary] = bench.table("assembly_summary")
    assert summary["status"] == "fragmented"
    assert summary["assembler_says"] == "2 contig(s)"
    assert summary["graph_call"] == "unresolved"
    assert bench.manifest()["artifacts"] == []
    assert (bench.run_dir() / "01_assembly" / "scaffolds.fasta").is_file()


def test_nothing_plastid_in_a_finished_run_is_no_result(bench: Bench) -> None:
    (bench.canned / "oatk.pltd.ctg.fasta").write_text("", encoding="utf-8")
    assert bench.run() == 1
    [summary] = bench.table("assembly_summary")
    assert summary["status"] == "no_result"
    assert summary["graph_call"] == "not measured"


def test_a_contig_its_graph_cannot_spell_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    changed = ("A" if SEQUENCE[5_000] != "A" else "C").join((SEQUENCE[:5_000], SEQUENCE[5_001:]))
    bench.canned_files(_contig("u1+,u2+,u3+,u2-", changed), _gfa())
    assert bench.run() == 1
    assert "not the sequence its walk spells" in capsys.readouterr().err


def _prepare(bench: Bench) -> Path:
    command = ["assemble", "prepare", *bench.given, "-o", str(bench.root / "p"), *bench.sets]
    assert main(command) == 0
    [prepared] = list((bench.root / "p").glob("*/*/01_assembly/prepared.json"))
    return prepared


def _result(bench: Bench) -> Path:
    result = bench.root / "result"
    result.mkdir()
    for path in bench.canned.iterdir():
        (result / path.name).write_bytes(path.read_bytes())
    (result / "oatk.log").write_text("[M::main] Version: 1.0\n", encoding="utf-8")
    return result


def _import(bench: Bench, prepared: Path, result: Path) -> int:
    return main(
        [
            "assemble",
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
        ]
    )


def test_prepare_runs_no_tool_and_prints_the_command(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    prepared = _prepare(bench)
    payload = json.loads(prepared.read_text(encoding="utf-8"))
    assert payload["tool"]["name"] == "oatk"
    assert payload["read_qc"] is None
    assert payload["raw_reads"][0]["path"] == str(bench.reads.resolve())
    assert "oatk -k 1001 -c 30 -D 1500000000" in capsys.readouterr().out


def test_an_imported_run_gives_the_local_result(bench: Bench) -> None:
    prepared = _prepare(bench)
    assert _import(bench, prepared, _result(bench)) == 0
    [summary] = bench.table("assembly_summary")
    assert summary["assembler"] == "Oatk 1.0"
    assert summary["read_support"] == "supported"
    assert summary["graph_call"] == "resolved"
    assert summary["backend"] == "import"
    assert summary["assembler_programs"] == "not stated"
    [artifact] = bench.manifest()["artifacts"]
    assert artifact["format"] == "fasta"


def test_an_import_with_no_log_is_refused(bench: Bench, capsys: pytest.CaptureFixture[str]) -> None:
    prepared = _prepare(bench)
    result = _result(bench)
    (result / "oatk.log").unlink()
    assert _import(bench, prepared, result) == 2
    assert "keep Oatk's stderr" in capsys.readouterr().err


def test_an_import_of_changed_reads_is_refused(
    bench: Bench, capsys: pytest.CaptureFixture[str]
) -> None:
    prepared = _prepare(bench)
    bench.reads.write_bytes(gzip.compress(b"@r\nC\n+\nI\n"))
    assert _import(bench, prepared, _result(bench)) == 2
    assert "not the input the imported result was asked for" in capsys.readouterr().err
