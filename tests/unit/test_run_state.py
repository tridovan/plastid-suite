"""What a run leaves on disk while it runs and when it is stopped.

A viewer that opens a run directory has only the directory. Before these rules a run wrote its
manifest as its last act, so a run in progress, a run that was killed and a directory somebody
made by hand all looked the same, and `subprocess.run` killed only the direct child, leaving the
program a wrapper had started still computing.

These tests run the real command in its own process, because a signal is what is under test, and
stand in for the tool with a shell script, because the fast suite has no external tools. The script
is shaped like `mafft`: a wrapper waiting on a child that does the work.
"""

from __future__ import annotations

import json
import os
import signal
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest

from plastid.cli import main
from plastid.run.cache import Cache
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.state import EXIT_STATES, exit_state, list_runs, read_state
from plastid.run.tools import RunCancelled, run_tool

DATA = Path(__file__).parents[1] / "data"
GENOME = str(DATA / "NC_050833.1.gb")
ENTRY = "import sys; from plastid.cli import main; sys.exit(main(sys.argv[1:]))"

_WRAPPER = """#!/bin/sh
sleep 120 &
echo $! > "$FAKE_CHILD"
wait
"""


def _tool_dir(path: Path) -> Path:
    """`repeat-match` as a wrapper that never finishes, and the sibling its version comes from."""
    path.mkdir()
    (path / "repeat-match").write_text(_WRAPPER, encoding="utf-8")
    (path / "mummer").write_text("#!/bin/sh\necho 4.0.1\n", encoding="utf-8")
    for name in ("repeat-match", "mummer"):
        (path / name).chmod((path / name).stat().st_mode | stat.S_IEXEC)
    return path


def _alive(pid: int) -> bool:
    """A zombie has finished: it only waits for a parent nobody here controls."""
    try:
        state = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0]
    except OSError:
        return False
    return state != "Z"


def _wait_for(check, *, what: str, seconds: float = 30.0):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        found = check()
        if found:
            return found
        time.sleep(0.02)
    raise AssertionError(f"timed out waiting for {what}")


def _start(tmp_path: Path, argv: list[str]) -> tuple[subprocess.Popen, Path]:
    child_file = tmp_path / "child.pid"
    tools = _tool_dir(tmp_path / "tools")
    environment = {
        **os.environ,
        "PATH": f"{tools}{os.pathsep}{os.environ['PATH']}",
        "FAKE_CHILD": str(child_file),
    }
    process = subprocess.Popen(
        [sys.executable, "-c", ENTRY, *argv],
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return process, child_file


def _three_steps(tmp_path: Path) -> list[str]:
    return [
        "run",
        "-m",
        "genome_stats",
        "-m",
        "repeats",
        "-m",
        "structure",
        "-i",
        GENOME,
        "-o",
        str(tmp_path / "runs"),
        "--dataset",
        "stopped",
    ]


def _manifest_path(tmp_path: Path) -> Path:
    return next((tmp_path / "runs").glob("stopped/*/manifest.json"))


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _running_child(tmp_path: Path, process: subprocess.Popen, child_file: Path) -> int:
    def ready() -> int | None:
        if process.poll() is not None:
            raise AssertionError(f"the run ended early: {process.stderr.read()}")
        text = child_file.read_text().strip() if child_file.exists() else ""
        return int(text) if text else None

    return _wait_for(ready, what="the tool to start")


def test_a_run_in_progress_says_so_on_disk(tmp_path: Path) -> None:
    process, child_file = _start(tmp_path, _three_steps(tmp_path))
    try:
        _running_child(tmp_path, process, child_file)
        manifest = _read(_manifest_path(tmp_path))

        assert manifest["status"] == "running"
        assert manifest["pid"] == process.pid
        assert manifest["finished_at"] is None
        assert [entry["id"] for entry in manifest["modules"]] == ["genome_stats"]
        assert manifest["current"]["module"] == "repeats"
        assert manifest["current"]["order"] == 2
        assert manifest["current"]["of"] == 3
        assert read_state(_manifest_path(tmp_path).parent).status == "running"
    finally:
        process.send_signal(signal.SIGTERM)
        process.communicate(timeout=20)


@pytest.mark.parametrize("sent", [signal.SIGTERM, signal.SIGINT])
def test_a_cancelled_run_stops_its_tool_and_says_what_happened(
    tmp_path: Path, sent: signal.Signals
) -> None:
    process, child_file = _start(tmp_path, _three_steps(tmp_path))
    child = _running_child(tmp_path, process, child_file)
    process.send_signal(sent)
    _, stderr = process.communicate(timeout=20)

    assert process.returncode == 130, stderr
    assert "Traceback" not in stderr
    assert "cancelled" in stderr
    _wait_for(lambda: not _alive(child), what="the tool's own child to stop", seconds=10)

    run_dir = _manifest_path(tmp_path).parent
    manifest = _read(_manifest_path(tmp_path))
    assert manifest["status"] == "cancelled"
    assert manifest["exit_code"] == 130
    assert manifest["finished_at"] is not None
    assert "current" not in manifest
    assert [(entry["id"], entry["status"]) for entry in manifest["modules"]] == [
        ("genome_stats", "ok"),
        ("repeats", "cancelled"),
        ("structure", "skipped"),
    ]
    assert (run_dir / "00_run_info" / "run_info.md").is_file()
    assert read_state(run_dir).status == "cancelled"


def test_a_signal_already_ignored_stays_ignored(tmp_path: Path) -> None:
    """`nohup` keeps a long run alive past its terminal by ignoring SIGHUP before it starts."""
    child_file = tmp_path / "child.pid"
    tools = _tool_dir(tmp_path / "tools")
    process = subprocess.Popen(
        [sys.executable, "-c", ENTRY, *_three_steps(tmp_path)],
        env={
            **os.environ,
            "PATH": f"{tools}{os.pathsep}{os.environ['PATH']}",
            "FAKE_CHILD": str(child_file),
        },
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        preexec_fn=lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN),
    )
    try:
        _running_child(tmp_path, process, child_file)
        process.send_signal(signal.SIGHUP)
        time.sleep(0.3)
        assert process.poll() is None
        assert _read(_manifest_path(tmp_path))["status"] == "running"
    finally:
        process.send_signal(signal.SIGTERM)
        process.communicate(timeout=20)


def test_a_killed_run_reads_as_interrupted_not_as_running(tmp_path: Path) -> None:
    """SIGKILL gives the run no chance to write anything, so the reader has to work it out."""
    process, child_file = _start(tmp_path, _three_steps(tmp_path))
    child = _running_child(tmp_path, process, child_file)
    process.kill()
    process.communicate(timeout=20)
    try:
        run_dir = _manifest_path(tmp_path).parent
        assert _read(run_dir / "manifest.json")["status"] == "running"
        state = read_state(run_dir)
        assert state.status == "interrupted"
        assert state.reason
    finally:
        os.kill(child, signal.SIGKILL)


def test_a_timeout_takes_the_tools_own_children_with_it(tmp_path: Path) -> None:
    """The same defect from the other side: `subprocess.run` killed the wrapper and nothing else."""
    tool = tmp_path / "wrapper.sh"
    child_file = tmp_path / "child.pid"
    tool.write_text(_WRAPPER, encoding="utf-8")
    tool.chmod(tool.stat().st_mode | stat.S_IEXEC)
    os.environ["FAKE_CHILD"] = str(child_file)
    try:
        with pytest.raises(subprocess.TimeoutExpired):
            run_tool(
                "wrapper", [str(tool)], log_dir=tmp_path / "logs", timeout=1, pinned_version="1"
            )
    finally:
        del os.environ["FAKE_CHILD"]
    child = int(child_file.read_text().strip())
    _wait_for(lambda: not _alive(child), what="the child to stop", seconds=10)


def test_a_cancelled_batch_says_which_samples_never_ran(tmp_path: Path) -> None:
    process, child_file = _start(
        tmp_path,
        [
            "batch",
            "-m",
            "repeats",
            "-i",
            GENOME,
            "-i",
            GENOME,
            "-o",
            str(tmp_path / "runs"),
            "--dataset",
            "set",
        ],
    )
    _running_child(tmp_path, process, child_file)
    path = next((tmp_path / "runs").glob("set/*/batch.json"))
    during = _read(path)
    assert during["status"] == "running"
    assert during["pid"] == process.pid
    assert during["samples_total"] == 2
    assert read_state(path.parent).status == "running"

    process.send_signal(signal.SIGTERM)
    _, stderr = process.communicate(timeout=20)

    assert process.returncode == 130, stderr
    after = _read(path)
    assert after["status"] == "cancelled"
    assert after["exit_code"] == 130
    assert [sample["status"] for sample in after["samples"]] == ["cancelled"]
    assert after["not_run"] == [GENOME]
    assert read_state(path.parent).status == "cancelled"


def _manifest() -> dict:
    return build_manifest(
        run_id="20260101-000000-aaaaaaaa",
        command=["plastid", "run"],
        package_version="0",
        started_at="2026-01-01T00:00:00",
        params={},
        inputs=[],
        coordinate_system={},
    )


def test_the_manifest_is_replaced_whole(tmp_path: Path) -> None:
    """A viewer reads the file while the run rewrites it, so a half-written file must not exist."""
    path = tmp_path / "manifest.json"
    write_manifest(path, _manifest())
    write_manifest(path, _manifest())
    assert [item.name for item in tmp_path.iterdir()] == ["manifest.json"]


def test_a_cache_entry_is_written_whole_or_not_at_all(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A stop while an entry was being written left cut JSON that every later run then read."""
    cache = Cache(tmp_path)

    def stopped(*_args: object) -> None:
        raise RunCancelled("SIGTERM")

    monkeypatch.setattr("plastid.run.cache.os.replace", stopped)
    with pytest.raises(RunCancelled):
        cache.put("key", {"tables/a.csv": "text"}, provenance={})
    assert cache.get("key") is None


def test_a_manifest_names_the_process_and_machine_it_ran_on() -> None:
    """A process number means nothing on another machine, or after a container restarts."""
    manifest = _manifest()
    assert manifest["pid"] == os.getpid()
    assert manifest["host"]


def test_a_directory_without_a_manifest_is_named_as_unknown(tmp_path: Path) -> None:
    """Runs killed before this rule existed left exactly this."""
    run_dir = tmp_path / "20260101-000000-aaaaaaaa"
    (run_dir / "00_run_info").mkdir(parents=True)
    state = read_state(run_dir)
    assert state.status == "unknown"
    assert "no manifest.json" in state.reason


def test_runs_are_listed_newest_first(tmp_path: Path) -> None:
    for name, status in (
        ("20260101-000000-aaaaaaaa", "ok"),
        ("20260102-000000-bbbbbbbb", "failed"),
    ):
        run_dir = tmp_path / "set" / name
        run_dir.mkdir(parents=True)
        manifest = _manifest()
        manifest["run_id"] = name
        manifest["status"] = status
        write_manifest(run_dir / "manifest.json", manifest)
    (tmp_path / ".cache").mkdir()

    assert [(state.run_id, state.status) for state in list_runs(tmp_path)] == [
        ("20260102-000000-bbbbbbbb", "failed"),
        ("20260101-000000-aaaaaaaa", "ok"),
    ]


def test_every_exit_code_has_one_name() -> None:
    codes = [state.code for state in EXIT_STATES]
    names = [state.status for state in EXIT_STATES]
    assert len(set(codes)) == len(codes)
    assert len(set(names)) == len(names)
    assert exit_state(0).status == "ok"
    assert exit_state(2).status == "refused"
    assert exit_state(130).status == "cancelled"


def test_a_refused_run_is_recorded_as_refused(tmp_path: Path) -> None:
    """Exit 2 means the run never started on the work, which `failed` does not say."""
    gated = str(DATA / "NC_024936.1.gb")
    out = tmp_path / "runs"
    assert main(["run", "-m", "structure", "-i", gated, "-o", str(out), "--dataset", "gate"]) == 2
    manifest = _read(next(out.glob("gate/*/manifest.json")))
    assert manifest["status"] == "refused"
    assert manifest["exit_code"] == 2
    assert "IR gate" in manifest["modules"][0]["error"]


def test_the_manifest_records_the_command_that_was_typed(tmp_path: Path) -> None:
    """The recorded command must run the same analysis again, so nothing may fall out of it."""
    out = tmp_path / "runs"
    argv = [
        "run",
        "-m",
        "structure",
        "-i",
        GENOME,
        "-o",
        str(out),
        "--dataset",
        "typed",
        "--label",
        "stand-in",
        "--set",
        "window=2000",
        "--docx",
    ]
    assert main(argv) == 0
    manifest = _read(next(out.glob("typed/*/manifest.json")))
    assert manifest["command"] == ["plastid", *argv]


def test_the_output_commands_read_a_cancelled_run(tmp_path: Path) -> None:
    """A stopped run keeps what it finished, so its record and tables have to stay readable."""
    process, child_file = _start(tmp_path, _three_steps(tmp_path))
    _running_child(tmp_path, process, child_file)
    process.send_signal(signal.SIGTERM)
    process.communicate(timeout=20)
    run_dir = _manifest_path(tmp_path).parent

    assert main(["qc", str(run_dir)]) == 0
    assert main(["figures", str(run_dir), "--list"]) == 0


@pytest.mark.parametrize("command", ["qc", "package", "figures"])
def test_a_run_in_progress_is_not_read_as_a_result(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], command: str
) -> None:
    """Its tables are not finished, and these commands would rewrite the record it is writing."""
    process, child_file = _start(tmp_path, _three_steps(tmp_path))
    try:
        _running_child(tmp_path, process, child_file)
        run_dir = _manifest_path(tmp_path).parent
        assert main([command, str(run_dir)]) == 2
        assert "still running" in capsys.readouterr().err
    finally:
        process.send_signal(signal.SIGTERM)
        process.communicate(timeout=20)


def test_a_cancelled_assembly_step_says_so(tmp_path: Path) -> None:
    """A command other than `plastid run` holds its own manifest, and closes it the same way."""
    tools = tmp_path / "tools"
    tools.mkdir()
    fastp = tools / "fastp"
    fastp.write_text(
        '#!/bin/sh\nif [ "$1" = "--version" ]; then echo "fastp 1.4.0"; exit 0; fi\n'
        + _WRAPPER.split("\n", 1)[1],
        encoding="utf-8",
    )
    fastp.chmod(fastp.stat().st_mode | stat.S_IEXEC)
    for name in ("r_1.fastq.gz", "r_2.fastq.gz"):
        (tmp_path / name).write_bytes(b"reads")
    child_file = tmp_path / "child.pid"
    process = subprocess.Popen(
        [
            sys.executable,
            "-c",
            ENTRY,
            "assemble",
            "prepare",
            "-1",
            str(tmp_path / "r_1.fastq.gz"),
            "-2",
            str(tmp_path / "r_2.fastq.gz"),
            "--sample",
            "stopped",
            "-o",
            str(tmp_path / "runs"),
        ],
        env={
            **os.environ,
            "PATH": f"{tools}{os.pathsep}{os.environ['PATH']}",
            "FAKE_CHILD": str(child_file),
        },
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    child = _running_child(tmp_path, process, child_file)
    assert _read(_manifest_path(tmp_path))["status"] == "running"
    process.send_signal(signal.SIGTERM)
    _, stderr = process.communicate(timeout=20)

    assert process.returncode == 130, stderr
    _wait_for(lambda: not _alive(child), what="the tool's own child to stop", seconds=10)
    manifest = _read(_manifest_path(tmp_path))
    assert manifest["status"] == "cancelled"
    assert [entry["status"] for entry in manifest["modules"]] == ["cancelled"]


def test_a_run_that_finished_no_step_is_read_as_what_it_is(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A refused or stopped run is a run: it has a report, and nothing to draw."""
    out = tmp_path / "runs"
    gated = str(DATA / "NC_024936.1.gb")
    assert main(["run", "-m", "structure", "-i", gated, "-o", str(out), "--dataset", "gate"]) == 2
    run_dir = next(out.glob("gate/*"))
    capsys.readouterr()

    assert main(["qc", str(run_dir)]) == 0
    assert main(["figures", str(run_dir), "--list"]) == 0
    capsys.readouterr()
    for command in ("figures", "package"):
        assert main([command, str(run_dir)]) == 2
        assert "finished no step (status refused)" in capsys.readouterr().err
    state = read_state(run_dir)
    assert (state.done, state.total) == (0, 1)


def test_runs_with_no_runs_yet_says_so(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """The first command of a new user, in an empty directory."""
    monkeypatch.chdir(tmp_path)
    assert main(["runs"]) == 0
    assert "no runs" in capsys.readouterr().out
