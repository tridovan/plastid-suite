"""Run layer: layout · manifest · run_tool · cache · runner (P1-05)."""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import pytest

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
from plastid.run.runner import execute_step
from plastid.run.tools import ToolError, run_tool, tool_version, which

WHEN = datetime(2026, 10, 2, 15, 30, 0)


def test_slugify_strips_diacritics() -> None:
    assert slugify("Iris domestica — NC_050833.1") == "iris-domestica-nc_050833.1"
    assert slugify("chạy Trí") == "chay-tri"  # PROJECT.md §9.16: Vietnamese paths caused bugs
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
