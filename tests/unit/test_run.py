"""Run layer: layout · manifest · run_tool · cache · runner (P1-05)."""

from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Genome
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout, slugify
from plastid.run.manifest import (
    MANIFEST_SCHEMA_VERSION,
    ManifestError,
    build_manifest,
    validate,
    write_manifest,
)
from plastid.run.runner import _input_hash, _module_source_hash, _params_hash, execute_step
from plastid.run.tools import ToolError, run_tool, tool_version, which

WHEN = datetime(2026, 10, 2, 15, 30, 0)


def test_slugify_strips_diacritics() -> None:
    assert slugify("Iris domestica — NC_050833.1") == "iris-domestica-nc_050833.1"
    assert slugify("chạy Trí") == "chay-tri"  # PROJECT.md §9.30: Vietnamese paths caused bugs
    assert slugify("iris_domestica") == "iris_domestica"


def test_layout_builds_run_scoped_paths(tmp_path: Path) -> None:
    layout = RunLayout.create(tmp_path, "Iris domestica", fingerprint="abcdef1234567890", now=WHEN)
    assert layout.dataset == "iris-domestica"
    assert layout.run_id == "20261002-153000-abcdef12"
    assert layout.run_dir == tmp_path / "iris-domestica" / "20261002-153000-abcdef12"
    layout.make_dirs()
    assert layout.run_info_dir.is_dir()
    assert layout.logs_dir.is_dir()
    assert layout.module_dir(1, "genome_stats").name == "01_genome_stats"


def test_two_runs_never_share_a_directory_however_the_clock_behaves(tmp_path: Path) -> None:
    """One run = one directory, even in the same second — and even when the clock steps back.

    Measured 08/10/2026 (PROJECT.md §9.C-37): this machine's clock moved ~1 s *backwards* between
    two runs, so a later run was stamped earlier than the one before it and the run directories
    no longer sorted into the order they were made. The guarantee that survives a clock that is
    not monotonic is not "the names sort" but "no run overwrites another".
    """
    first = RunLayout.create(tmp_path, "iris", fingerprint="abcdef1234567890", now=WHEN)
    first.make_dirs()
    same_second = RunLayout.create(tmp_path, "iris", fingerprint="abcdef1234567890", now=WHEN)
    stepped_back = RunLayout.create(
        tmp_path, "iris", fingerprint="abcdef1234567890", now=datetime(2026, 10, 2, 15, 29, 59)
    )
    assert len({first.run_id, same_second.run_id, stepped_back.run_id}) == 3
    assert stepped_back.run_dir.parent == first.run_dir.parent  # same dataset, different run


def test_manifest_validates_and_round_trips(tmp_path: Path) -> None:
    manifest = build_manifest(
        run_id="r1",
        command=["plastid", "run"],
        package_version="0.1.0",
        started_at="2026-10-02T15:30:00",
        params={},
        inputs=[],
        coordinate_system={"name": "genbank"},
    )
    validate(manifest)
    assert manifest["schema_version"] == MANIFEST_SCHEMA_VERSION

    broken = dict(manifest)
    del broken["inputs"]
    with pytest.raises(ManifestError, match="inputs"):
        validate(broken)

    path = tmp_path / "manifest.json"
    write_manifest(path, manifest)
    assert json.loads(path.read_text(encoding="utf-8"))["run_id"] == "r1"


def test_run_tool_records_version_and_logs(tmp_path: Path) -> None:
    record = run_tool(
        "python",
        [sys.executable, "-c", "print('hello from tool')"],
        kind="binary",
        log_dir=tmp_path,
    )
    assert record.rc == 0
    assert record.version  # e.g. "Python 3.12.3"
    assert Path(record.stdout_file).read_text(encoding="utf-8").strip() == "hello from tool"
    assert record.duration_s >= 0
    assert record.stdout_tail.strip() == "hello from tool"


def test_run_tool_raises_with_tail_on_failure(tmp_path: Path) -> None:
    with pytest.raises(ToolError, match="boom"):
        run_tool(
            "failing-tool",
            [sys.executable, "-c", "import sys; print('boom'); sys.exit(3)"],
            kind="binary",
            log_dir=tmp_path,
        )


def test_which_and_pinned_version() -> None:
    assert which("definitely-not-a-real-tool-xyz") is None
    assert tool_version("misa", "vendored", path=None, pinned="2.1") == ("2.1", "pinned")


def test_cache_round_trip(tmp_path: Path) -> None:
    cache = Cache(tmp_path)
    key = Cache.key(module_id="m", module_version="1", params_hash="p", input_sha256s=("a",))
    assert cache.get(key) is None
    cache.put(key, {"tables/x.csv": "a,b\n"})
    assert cache.get(key) == {"tables/x.csv": "a,b\n"}


def test_cache_round_trips_bytes_without_touching_them(tmp_path: Path) -> None:
    """A .docx is a zip: JSON has no bytes, and a cache hit has to restore the exact file.

    The payload is stored as base64 — any lossy step here would show up as a document whose
    bytes differ between a fresh run and a cached one, which is exactly what P2-09 forbids.
    """
    cache = Cache(tmp_path)
    blob = bytes(range(256)) * 4
    key = Cache.key(module_id="m", module_version="1", params_hash="p", input_sha256s=("a",))
    cache.put(key, {"tables/x.csv": "a,b\n", "report.docx": blob})
    restored = cache.get(key)
    assert restored == {"tables/x.csv": "a,b\n", "report.docx": blob}
    assert isinstance(restored["report.docx"], bytes)
    # the cache file itself stays valid JSON — no bytes leaked into it
    assert "report.docx" in (tmp_path / f"{key}.json").read_text(encoding="utf-8")


def _stub_result() -> Result:
    return Result(
        module="stub",
        module_version="0.1.0",
        tables={
            "stub": Table(
                id="stub",
                title="Stub",
                columns=(Column("n", "N", dtype="int"),),
                rows=({"n": 1},),
                scope=Scope(label="synthetic"),
            )
        },
    )


def _stub_genome() -> Genome:
    return Genome(name="synthetic", sequence="ACGT" * 25)


def test_runner_refuses_the_wrong_number_of_genomes(tmp_path: Path) -> None:
    """API v0.2: a module's declared arity is enforced before anything runs."""
    layout = RunLayout.create(tmp_path, "dataset", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    with pytest.raises(ValueError, match="at least 2 genomes"):
        execute_step(
            module_id="stub",
            module_version="0.1.0",
            run_fn=lambda genomes, params, ctx: _stub_result(),
            genomes=(_stub_genome(),),
            arity=2,
            params={},
            ctx=RunContext(layout=layout, params={}, module_id="stub", module_order=1),
            cache=Cache(layout.cache_dir),
        )


def test_runner_caches_and_restores_files_byte_identically(tmp_path: Path) -> None:
    calls: list[int] = []

    def run_fn(genome: Genome, params: dict, ctx: RunContext) -> Result:
        calls.append(1)
        return _stub_result()

    layout = RunLayout.create(tmp_path, "dataset", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    cache = Cache(layout.cache_dir)

    def context() -> RunContext:
        return RunContext(layout=layout, params={}, module_id="stub", module_order=1)

    first = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
    )
    assert first.status == "ok"
    assert len(calls) == 1
    csv_path = layout.module_dir(1, "stub") / "tables" / "stub.csv"
    assert csv_path.is_file()
    original = csv_path.read_bytes()

    second = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
    )
    assert second.status == "cached"
    assert len(calls) == 1  # the module never ran again
    assert csv_path.read_bytes() == original

    csv_path.unlink()
    third = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
    )
    assert third.status == "cached"
    assert csv_path.read_bytes() == original  # restored byte-identically


def test_require_ir_fails_loudly() -> None:
    from plastid.core.ir import IRGateError, IRResult
    from plastid.core.model import Interval

    without_ir = RunContext(layout=None, params={}, module_id="m", module_order=1)  # type: ignore[arg-type]
    with pytest.raises(IRGateError, match="never determined"):
        without_ir.require_ir()

    failing_ir = IRResult(
        ira=Interval(0, 10),
        irb=Interval(20, 30),
        identity=0.92,
        mismatches=8,
        compared_length=10,
        exact_core=1,
        mirror_constant=29,
        source="manual",
    )
    gated = RunContext(layout=None, params={}, module_id="m", module_order=1, ir=failing_ir)  # type: ignore[arg-type]
    with pytest.raises(IRGateError, match="IR gate failed"):
        gated.require_ir()


def test_the_version_probe_knows_the_single_dash_spelling(tmp_path: Path) -> None:
    """Measured 03/10/2026 on BLAST+ 2.17.0: ``blastx --version`` exits 1 with a usage message and
    ``blastx -version`` prints the version. A probe that only tried the double dash reported the
    tool that now annotates mitogenomes as "unknown" — a silent hole in the tool report."""
    script = tmp_path / "blastx"
    script.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "-version" ]; then echo "blastx: 2.17.0+"; exit 0; fi\n'
        "echo USAGE; exit 1\n",
        encoding="utf-8",
    )
    script.chmod(0o755)

    version, source = tool_version("blastx", "binary", str(script))

    assert (version, source) == ("blastx: 2.17.0+", "flag")


def test_a_tool_that_answers_a_version_flag_with_an_error_is_unknown(tmp_path: Path) -> None:
    """Measured 03/10/2026 on ARAGORN: ``aragorn -V`` exits **0** and prints "No sequence file
    specified, type aragorn -h for help". Recorded verbatim, that reads like a version in a
    manifest — a provenance entry nobody can notice is wrong."""
    script = tmp_path / "aragorn"
    script.write_text(
        "#!/bin/sh\necho 'No sequence file specified, type aragorn -h for help' >&2\nexit 0\n",
        encoding="utf-8",
    )
    script.chmod(0o755)

    version, source = tool_version("aragorn", "binary", str(script))

    assert (version, source) == (None, "unknown")


def test_the_version_probe_never_runs_a_perl_tool(tmp_path: Path, monkeypatch) -> None:
    """Measured 2026-10-03: probing PGA for a version with an unknown flag made it run itself with
    its *default* arguments, and its first act is `rm -rf gb` — deleting a directory named ``gb``
    in whatever directory the toolkit ran from. A Perl tool is now asked for its interpreter's
    version instead, and the probe never executes the script."""
    from plastid.run.tools import tool_version

    victim = tmp_path / "gb"
    victim.mkdir()
    (victim / "important.txt").write_text("do not delete\n", encoding="utf-8")

    script = tmp_path / "tool.pl"
    script.write_text(
        "#!/usr/bin/env perl\n"
        "use strict;\n"
        'system("rm -rf gb");\n'  # what PGA does on its first run
        'print "v9.9\\n";\n',
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.chdir(tmp_path)

    version, source = tool_version("tool", "perl", str(script))

    assert (victim / "important.txt").is_file(), "the probe executed the tool"
    assert source == "stdout-header"  # the interpreter's version, which is the honest answer
    assert version and version.startswith("v5")


def test_a_data_label_reaches_every_exported_table(tmp_path: Path) -> None:
    """A stand-in dataset never travels unlabelled (CLAUDE.md) — and the label is applied once,
    in the runner, so no module has to remember it. Every table's scope line carries it, which is
    the line the CSV header, the JSON export and the report all read."""
    from plastid.run.runner import DATA_LABEL_PARAM, label_result

    result = _stub_result()
    labelled = label_result(result, "stand-in (public records)")
    for table in labelled.tables.values():
        assert table.scope.label.startswith("stand-in (public records) · ")
        assert table.scope.coord_system == result.tables[table.id].scope.coord_system
    assert labelled.tables["stub"].rows == result.tables["stub"].rows
    # no label, no change — the parameter is what decides, not a default
    assert label_result(result, None).tables["stub"].scope.label == "synthetic"
    assert DATA_LABEL_PARAM == "data_label"


def test_the_data_label_reaches_the_rendered_report(tmp_path: Path) -> None:
    """The tables carry the label in their scope lines; the prose carries it here, so a report
    read on its own says where its numbers came from."""
    from plastid.cli import _with_data_label

    def renderer(result):
        return "## Genome statistics — X\n\nNumbers.\n\n*Draft — needs human review.*\n"

    with_label = _with_data_label(renderer, {"data_label": "stand-in (public records)"})(
        _stub_result()
    )
    assert "> **Data:** stand-in (public records)" in with_label
    assert with_label.startswith("## Genome statistics — X\n")
    # no label, no banner
    assert _with_data_label(renderer, {})(_stub_result()) == renderer(_stub_result())


def test_every_tool_call_keeps_its_own_log(tmp_path: Path) -> None:
    """Measured 2026-10-04: the log index was counted with a glob that never matched the names
    being written (``*__name__*.stdout.log`` vs ``001__name__stdout.log``), so **every call
    overwrote the previous one's logs** — each run kept the last tool call's output and nothing
    else, and "the raw output travels with the run" was false for every tool but the last."""
    layout = RunLayout.create(tmp_path, "logs", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="stub", module_order=1)
    for index in range(3):
        ctx.tool("echo", ["echo", f"call {index}"], cwd=tmp_path)
    names = sorted(path.name for path in layout.logs_dir.iterdir())
    assert names == [
        "001__echo__stderr.log",
        "001__echo__stdout.log",
        "002__echo__stderr.log",
        "002__echo__stdout.log",
        "003__echo__stderr.log",
        "003__echo__stdout.log",
    ]
    # the first call's output is still there — that is the point
    assert (layout.logs_dir / "001__echo__stdout.log").read_text(
        encoding="utf-8"
    ).strip() == "call 0"
    assert [record.stdout_file for record in ctx.tools] == [
        str(layout.logs_dir / f"{index:03d}__echo__stdout.log") for index in (1, 2, 3)
    ]


def test_two_annotations_of_one_sequence_are_two_cache_keys() -> None:
    """H-15, measured 05/10/2026: the corpus's own pair shares a sequence and not a gene set.

    The published and PGA records of *Iris domestica* carry the same 153,724 bp; they differ in
    their annotation (86 vs 85 CDS). Hashing the sequence alone gave them one cache key, and six
    modules take ``genomes[0]`` as their reference — so the second run was served the first's
    tables. The premise is asserted too: if the sequences ever differed, this test would stop
    measuring what it says it measures.
    """
    corpus = Path(__file__).resolve().parents[1] / "data" / "corpus"
    published = read_genbank(corpus / "Iris_domestica_NC_050833.1.pga.gb")
    golden = read_genbank(Path(__file__).resolve().parents[1] / "data" / "NC_050833.1.gb")
    assert published.sequence == golden.sequence, "the pair must carry one sequence"
    assert len(published.genes) != len(golden.genes), "and two different annotations"
    assert _input_hash(published) != _input_hash(golden)


def test_the_cache_key_keeps_the_input_order() -> None:
    """Input order is meaning: six modules read ``genomes[0]`` as their reference."""
    first = Cache.key(
        module_id="kaks",
        module_version="0.1.0",
        params_hash="p",
        input_sha256s=("aa", "bb"),
    )
    swapped = Cache.key(
        module_id="kaks",
        module_version="0.1.0",
        params_hash="p",
        input_sha256s=("bb", "aa"),
    )
    assert first != swapped


def test_a_parameter_that_names_a_file_is_content_addressed(tmp_path: Path) -> None:
    """D5-09 §7.2: a path is not an identity — the ★ record's sha256 has to reach the key."""
    reference = tmp_path / "naming_reference.gb"
    reference.write_text("LOCUS       NC_050833\n", encoding="utf-8")
    before = _params_hash({"naming_reference": str(reference)})
    assert before == _params_hash({"naming_reference": str(reference)})
    reference.write_text("LOCUS       NC_024936\n", encoding="utf-8")
    assert _params_hash({"naming_reference": str(reference)}) != before
    # a parameter that is not a path is hashed as itself, and stays stable
    assert _params_hash({"mode": "normalize"}) == _params_hash({"mode": "normalize"})
    assert _params_hash({"mode": "normalize"}) != _params_hash({"mode": "verify"})


def test_a_changed_module_source_changes_the_key(tmp_path: Path) -> None:
    """H-08: the key used to carry the version and nothing else about the code.

    The module file is loaded for real rather than stubbed: what the guard hashes is *the file
    the function came from*, so a stub defined inside this test would hash this test.
    """

    def load(text: str):
        path = tmp_path / "stub_module.py"
        path.write_text(text, encoding="utf-8")
        spec = importlib.util.spec_from_file_location(f"stub_module_{len(text)}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.run

    first = _module_source_hash(load("def run():\n    return 1\n"))
    same = _module_source_hash(load("def run():\n    return 1\n"))
    changed = _module_source_hash(load("def run():\n    return 2\n"))
    assert first == same
    assert changed != first


def test_a_second_run_asking_for_a_document_is_not_served_the_cached_payload(
    tmp_path: Path,
) -> None:
    """H-14, measured 05/10/2026: a `--docx` run came back `cached` with 0 .docx files."""
    layout = RunLayout.create(tmp_path, "dataset", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    cache = Cache(layout.cache_dir)
    calls: list[int] = []

    def run_fn(genome: Genome, params: dict, ctx: RunContext) -> Result:
        calls.append(1)
        return _stub_result()

    def context() -> RunContext:
        return RunContext(layout=layout, params={}, module_id="stub", module_order=1)

    first = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
        render_document=lambda result, report: b"a document",
    )
    assert first.status == "ok" and len(calls) == 1

    # the same run, now also asked for the Word document: a new key, and the document appears
    second = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
        render_document=lambda result, report: b"a document",
    )
    assert second.status == "cached"

    third = execute_step(
        module_id="stub",
        module_version="0.1.0",
        run_fn=run_fn,
        genomes=(_stub_genome(),),
        arity=1,
        params={},
        ctx=context(),
        cache=cache,
        render_document=None,
    )
    assert third.status == "ok", "a payload without the document must not answer a --docx run"
