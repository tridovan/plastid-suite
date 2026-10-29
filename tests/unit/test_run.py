"""Run layer: layout · manifest · run_tool · cache · runner (P1-05)."""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import pytest

from plastid.core.genbank import read_genbank
from plastid.core.model import Gene, Genome, Interval
from plastid.result.result import Column, Result, Scope, Table
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.contract import (
    ARTIFACT_SCHEMA_VERSION,
    ContractError,
    describe_artifact,
    describe_reference_set,
    validate_artifact,
)
from plastid.run.layout import RunLayout, slugify
from plastid.run.manifest import (
    MANIFEST_SCHEMA_VERSION,
    ManifestError,
    build_manifest,
    validate,
    write_manifest,
)
from plastid.run.runner import (
    _input_hash,
    _module_source_hash,
    _package_source_hash,
    _params_hash,
    execute_step,
)
from plastid.run.tools import ToolError, run_tool, tool_version, which

WHEN = datetime(2026, 10, 2, 15, 30, 0)


def test_slugify_strips_diacritics() -> None:
    assert slugify("Iris domestica — NC_050833.1") == "iris-domestica-nc_050833.1"
    assert slugify("chạy Trí") == "chay-tri"  # measured: Vietnamese paths caused bugs
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

    Measured 08/10/2026: this machine's clock moved ~1 s *backwards* between
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


def test_run_info_names_the_original_computation_of_a_cached_step() -> None:
    """A hit is not a computation that ran now: say when it ran and which tools it used."""
    from plastid.run.runinfo import render_run_info

    manifest = build_manifest(
        run_id="r1",
        command=["plastid", "run"],
        package_version="1.0.0",
        started_at="2026-10-26T12:00:00",
        params={},
        inputs=[],
        coordinate_system={"name": "genbank"},
    )
    manifest["modules"] = [
        {
            "id": "ssr",
            "version": "0.1.0",
            "status": "cached",
            "files": [],
            "provenance": {
                "computed_at": "2026-10-25T09:00:00",
                "tools": [{"name": "misa", "version": "2.1"}],
            },
        }
    ]
    text = render_run_info(manifest, ())
    assert "payload computed 2026-10-25T09:00:00 by misa 2.1" in text


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
    provenance = {"module": "m", "tools": [], "computed_at": "2026-10-26T00:00:00"}
    cache.put(key, {"tables/x.csv": "a,b\n"}, provenance=provenance)
    entry = cache.get(key)
    assert entry is not None
    assert entry.files == {"tables/x.csv": "a,b\n"}
    assert entry.provenance == provenance


def test_an_entry_without_provenance_is_not_reused(tmp_path: Path) -> None:
    """A payload whose origin cannot be stated is not one this runner serves (M11.12).

    Old entries are the concrete case: the key changes too, but the schema guard is what makes
    "old payload not reused" an explicit rule rather than a property of one hash function.
    """
    cache = Cache(tmp_path)
    key = Cache.key(module_id="m", module_version="1", params_hash="p", input_sha256s=("a",))
    path = tmp_path / f"{key}.json"
    path.write_text('{"files": {"tables/x.csv": "a,b\\n"}}\n', encoding="utf-8")
    assert cache.get(key) is None
    path.write_text(
        '{"schema": 1, "provenance": {}, "files": {"tables/x.csv": "a,b\\n"}}\n',
        encoding="utf-8",
    )
    assert cache.get(key) is None


def test_cache_round_trips_bytes_without_touching_them(tmp_path: Path) -> None:
    """A .docx is a zip: JSON has no bytes, and a cache hit has to restore the exact file.

    The payload is stored as base64 — any lossy step here would show up as a document whose
    bytes differ between a fresh run and a cached one, which is exactly what P2-09 forbids.
    """
    cache = Cache(tmp_path)
    blob = bytes(range(256)) * 4
    key = Cache.key(module_id="m", module_version="1", params_hash="p", input_sha256s=("a",))
    cache.put(
        key, {"tables/x.csv": "a,b\n", "report.docx": blob}, provenance={"module": "m", "tools": []}
    )
    entry = cache.get(key)
    assert entry is not None
    restored = entry.files
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


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("transl_table", 11),
        ("codon_start", 2),
        ("exception", "RNA editing"),
        ("trans_splicing", True),
    ],
)
def test_each_annotation_qualifier_moves_the_input_hash(field: str, value: object) -> None:
    """K-02, the unit half: every field an analysis can read is part of the input's identity."""
    genome = Genome(
        "probe",
        "GTGAAATAA",
        genes=(Gene("x", "CDS", (Interval(0, 9),), transl_table=1),),
    )
    changed = replace(genome, genes=(replace(genome.genes[0], **{field: value}),))
    assert _input_hash(genome) != _input_hash(changed)


def test_a_changed_translation_table_is_not_served_the_cached_orf_result(
    tmp_path: Path,
) -> None:
    """K-02 through the runner, with the module really running and the files really read.

    The same sequence and CDS at table 1 and table 11 must produce different status tables. A
    hash equality check alone would have missed the defect: this asserts the exported cells,
    then asserts that the *same* input is still a cache hit with the original provenance.
    """
    from plastid.modules.registry import get as get_module

    spec = get_module("orf")
    first = Genome(
        "probe",
        "GTGAAATAA",
        genes=(Gene("x", "CDS", (Interval(0, 9),), transl_table=1),),
    )
    second = replace(first, genes=(replace(first.genes[0], transl_table=11),))
    layout = RunLayout.create(tmp_path, "dataset", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    cache = Cache(layout.cache_dir)

    def run(genome: Genome):
        ctx = RunContext(layout=layout, params={}, module_id="orf", module_order=1)
        return execute_step(
            module_id=spec.id,
            module_version=spec.version,
            run_fn=spec.run,
            genomes=(genome,),
            arity=1,
            params={},
            ctx=ctx,
            cache=cache,
        )

    first_outcome = run(first)
    assert first_outcome.status == "ok"
    status_csv = layout.module_dir(1, "orf") / "tables" / "orf_status.csv"
    assert _status_count(status_csv, "no_start") == 1
    assert _status_count(status_csv, "alternative_start") == 0

    second_outcome = run(second)
    assert second_outcome.status == "ok", "the table-1 payload answered a table-11 run"
    assert _status_count(status_csv, "no_start") == 0
    assert _status_count(status_csv, "alternative_start") == 1

    hit = run(second)
    assert hit.status == "cached"
    assert hit.provenance == second_outcome.provenance
    assert hit.provenance is not None and hit.provenance["tools"] == []


def _status_count(path: Path, status: str) -> int:
    for line in path.read_text(encoding="utf-8").splitlines():
        cells = line.split(",")
        if cells and cells[0] == status:
            return int(cells[-1])
    raise AssertionError(f"{status} not found in {path}")


def test_shared_package_code_and_resources_move_the_key(tmp_path: Path) -> None:
    """The module file is not the whole code (M11.12).

    ``core/`` helpers and packaged resources (``defaults.yaml``, templates, name lists) can change
    a module's output without the module's own file changing. The shared digest covers both, and
    a new file is a change too — not just editing an existing one.
    """
    package = tmp_path / "plastid"
    (package / "core").mkdir(parents=True)
    (package / "data").mkdir()
    (package / "core" / "helper.py").write_text("VALUE = 1\n", encoding="utf-8")
    (package / "data" / "defaults.yaml").write_text("orf: 1\n", encoding="utf-8")
    baseline = _package_source_hash(package)

    (package / "data" / "defaults.yaml").write_text("orf: 11\n", encoding="utf-8")
    edited = _package_source_hash(package)
    assert edited != baseline, "editing a packaged resource did not change the shared digest"

    (package / "core" / "new.py").write_text("X = 2\n", encoding="utf-8")
    added = _package_source_hash(package)
    assert added not in {baseline, edited}, "adding a shared source file did not change it"


def test_a_changed_tool_version_invalidates_a_cached_payload(tmp_path: Path) -> None:
    """A hit restores a payload, not the tool that made it (M11.12).

    The recorded provenance is revalidated on every hit. When the tool now reports another
    version, the stored payload is not served as this run's result; the current run's tool call
    is recorded instead, while the hit's own provenance keeps saying what originally ran.
    """
    calls: list[int] = []
    script = tmp_path / "faketool"
    script.write_text(_fake_tool(version="1.0"), encoding="utf-8")
    script.chmod(0o755)

    def run_fn(genome: Genome, params: dict, ctx: RunContext) -> Result:
        calls.append(1)
        ctx.tool("faketool", [str(script), "work"])
        return _stub_result()

    layout = RunLayout.create(tmp_path, "dataset", fingerprint="0" * 16, now=WHEN)
    layout.make_dirs()
    cache = Cache(layout.cache_dir)

    def run():
        ctx = RunContext(layout=layout, params={}, module_id="stub", module_order=1)
        return execute_step(
            module_id="stub",
            module_version="0.1.0",
            run_fn=run_fn,
            genomes=(_stub_genome(),),
            arity=1,
            params={},
            ctx=ctx,
            cache=cache,
        )

    first = run()
    assert first.status == "ok"
    assert first.provenance is not None
    assert first.provenance["tools"][0]["version"] == "1.0"

    same = run()
    assert same.status == "cached"
    assert len(calls) == 1
    assert same.provenance == first.provenance

    script.write_text(_fake_tool(version="2.0"), encoding="utf-8")
    script.chmod(0o755)
    changed = run()
    assert changed.status == "ok", "a payload made by tool 1.0 was served after the upgrade"
    assert len(calls) == 2
    assert changed.provenance is not None
    assert changed.provenance["tools"][0]["version"] == "2.0"


def _fake_tool(*, version: str) -> str:
    # `tool_version` probes with each supported flag; the probe runs in a throwaway directory, so
    # the script must not rely on a relative path.
    return (
        "#!/bin/sh\n"
        'case "$1" in\n'
        f"  --version|-version|-v|-V) echo '{version}'; exit 0;;\n"
        "esac\n"
        "echo work\n"
    )


def test_a_parameter_that_names_a_file_is_content_addressed(tmp_path: Path) -> None:
    """D5-09: a path is not an identity — the ★ record's sha256 has to reach the key."""
    reference = tmp_path / "naming_reference.gb"
    reference.write_text("LOCUS       NC_050833\n", encoding="utf-8")
    before = _params_hash({"naming_reference": str(reference)})
    assert before == _params_hash({"naming_reference": str(reference)})
    reference.write_text("LOCUS       NC_024936\n", encoding="utf-8")
    assert _params_hash({"naming_reference": str(reference)}) != before
    # a parameter that is not a path is hashed as itself, and stays stable
    assert _params_hash({"mode": "normalize"}) == _params_hash({"mode": "normalize"})
    assert _params_hash({"mode": "normalize"}) != _params_hash({"mode": "verify"})


def test_a_parameter_that_names_a_directory_is_content_addressed(tmp_path: Path) -> None:
    """The same rule for a *directory* — the annotate step's reference set (P3B-09).

    Measured 09/10/2026: ``--set reference=<directory>`` was hashed as that string, so editing or
    re-pinning a record served the previous run's tables out of the cache with nothing to say the
    reference had moved — H-15's failure mode, one level up. The digest is over the file names
    *and* their contents, so all three ways a set changes (add, remove, edit) move the key.
    """
    reference = tmp_path / "reference"
    reference.mkdir()
    first = reference / "NC_050833.1.gb"
    first.write_text("LOCUS       NC_050833\n", encoding="utf-8")
    before = _params_hash({"reference": str(reference)})
    assert before == _params_hash({"reference": str(reference)})

    first.write_text("LOCUS       NC_050833  edited\n", encoding="utf-8")
    edited = _params_hash({"reference": str(reference)})
    assert edited != before, "editing a record did not move the key"

    (reference / "NC_024936.1.gb").write_text("LOCUS       NC_024936\n", encoding="utf-8")
    added = _params_hash({"reference": str(reference)})
    assert added not in {before, edited}, "adding a record did not move the key"

    (reference / "NC_024936.1.gb").unlink()
    assert _params_hash({"reference": str(reference)}) == edited, "removing restored the state"
    assert _params_hash({"reference": str(reference)}) != added


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


# ------------------------------------------------- artifact and reference-set contract (M11.16)


def _artifact_file(tmp_path: Path) -> Path:
    path = tmp_path / "plastid_annotation.gb"
    path.write_text("LOCUS       TEST                    36 bp    DNA\n//\n", encoding="utf-8")
    return path


def test_the_artifact_contract_names_what_a_next_step_reads(tmp_path: Path) -> None:
    path = _artifact_file(tmp_path)
    entry = describe_artifact(
        path,
        relative_to=tmp_path,
        sample_id="SYNTH.1",
        topology="circular",
        status="ok",
        reference_sha256="a" * 64,
    )
    assert entry["path"] == "plastid_annotation.gb"
    assert entry["format"] == "genbank"
    assert entry["artifact_schema_version"] == ARTIFACT_SCHEMA_VERSION
    assert (entry["sample_id"], entry["topology"], entry["frame"], entry["status"]) == (
        "SYNTH.1",
        "circular",
        "genbank",
        "ok",
    )
    assert entry["reference_sha256"] == "a" * 64
    validate_artifact(entry)

    for broken in (
        {**entry, "status": "failed"},  # a failed step's output may not be consumed
        {**entry, "format": "fasta"},  # a format this build cannot read as an artifact
        {key: value for key, value in entry.items() if key != "sha256"},
    ):
        with pytest.raises(ContractError):
            validate_artifact(broken)


def test_a_reference_set_digest_names_the_records_a_run_searched(tmp_path: Path) -> None:
    directory = tmp_path / "references"
    directory.mkdir()
    first = directory / "ONE.gb"
    second = directory / "TWO.gb"
    first.write_text("LOCUS       ONE\n//\n", encoding="utf-8")
    second.write_text("LOCUS       TWO\n//\n", encoding="utf-8")

    described = describe_reference_set(directory)
    assert described["role"] == "annotation_db"
    assert described["records"] == 2
    assert described["format"] == "genbank"

    second.write_text("LOCUS       TWO changed\n//\n", encoding="utf-8")
    assert describe_reference_set(directory)["sha256"] != described["sha256"], (
        "a changed record must change the reference-set digest"
    )

    with pytest.raises(ContractError, match="no GenBank records"):
        describe_reference_set(tmp_path)


def test_the_manifest_always_carries_the_contract_lists() -> None:
    manifest = build_manifest(
        run_id="r1",
        command=["plastid", "run"],
        package_version="1.0.0",
        started_at="2026-10-27T12:00:00",
        params={},
        inputs=[],
        coordinate_system={"name": "genbank"},
    )
    assert manifest["artifacts"] == []
    assert manifest["reference_sets"] == []
    validate(manifest)

    without = dict(manifest)
    del without["artifacts"]
    with pytest.raises(ManifestError, match="artifacts"):
        validate(without)


def test_the_ir_manifest_records_the_policy_that_measured_identity() -> None:
    """M11.16: an IR identity is a measurement, and the artifact must say how it was measured."""
    from plastid.cli import _ir_manifest_fields
    from plastid.core.ir import IR_SCORING_POLICY, IR_SCORING_POLICY_VERSION, IRResult
    from plastid.core.model import Interval

    ir = IRResult(
        ira=Interval(0, 10),
        irb=Interval(20, 30),
        identity=0.999,
        mismatches=1,
        compared_length=10,
        exact_core=1,
        mirror_constant=29,
        source="manual",
        matched=7,
        ambiguous_columns=1,
        internal_gaps=1,
        terminal_gaps=0,
        known_columns=8,
    )
    fields = _ir_manifest_fields(ir)
    assert fields["scoring_policy"] == IR_SCORING_POLICY
    assert fields["scoring_policy_version"] == IR_SCORING_POLICY_VERSION
    # The four column classes partition the compared columns, and each is named on its own.
    assert fields["matched_columns"] == 7
    assert fields["mismatched_columns"] == 1
    assert fields["ambiguous_columns"] == 1
    assert fields["gap_columns"] == 1
    assert (
        fields["matched_columns"]
        + fields["mismatched_columns"]
        + fields["ambiguous_columns"]
        + fields["gap_columns"]
        == fields["compared_columns"]
    )
    assert fields["requested_source"] == "annotation"


def test_the_ir_manifest_names_the_sequence_path_and_a_refused_annotation_differently() -> None:
    """M11.16: four IR outcomes, four states — never one "IR ok".

    The de-novo scan runs both when the operator asked for it and when the annotation was absent
    or refused, so the requested source and the presence of a refusal note are what separate the
    states. A gate failure stays visible through `gate_passed`, not by renaming the path.
    """
    from plastid.cli import _ir_status
    from plastid.core.ir import IRResult
    from plastid.core.model import Interval

    base = {
        "ira": Interval(0, 10),
        "irb": Interval(20, 30),
        "identity": 0.999,
        "mismatches": 0,
        "compared_length": 10,
        "exact_core": 1,
        "mirror_constant": 29,
    }
    assert _ir_status(IRResult(source="annotation", **base), "annotation") == "annotation"
    assert _ir_status(IRResult(source="manual", **base), "annotation") == "manual"
    assert (
        _ir_status(IRResult(source="denovo_diag", annotated_lengths=(10, 10), **base), "sequence")
        == "sequence_annotation_present"
    )
    assert (
        _ir_status(
            IRResult(
                source="denovo_diag",
                annotation_note="expected exactly two repeat_region copies",
                **base,
            ),
            "sequence",
        )
        == "sequence_annotation_refused"
    )
    assert (
        _ir_status(IRResult(source="denovo_diag", **base), "annotation") == "sequence_no_annotation"
    )
