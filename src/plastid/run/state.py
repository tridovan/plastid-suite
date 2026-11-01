"""A run's state: the vocabulary shared by exit codes, manifests and whoever reads a run directory.

A command ends in an exit code and a run directory holds a status, and the two are one word: a
script reads the number, a person or a viewer reads the name, and neither has to translate. The
states a directory can only be *read* as, because no process was left to write them, are worked out
here too, so every reader draws the same conclusion from the same files.
"""

from __future__ import annotations

import json
import os
import socket
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "EXIT_STATES",
    "ExitState",
    "RunState",
    "exit_state",
    "host_name",
    "list_runs",
    "read_state",
]


@dataclass(frozen=True, slots=True)
class ExitState:
    code: int
    status: str
    meaning: str


EXIT_STATES: tuple[ExitState, ...] = (
    ExitState(0, "ok", "finished"),
    ExitState(1, "failed", "the work started and did not finish, or a check found a problem"),
    ExitState(
        2,
        "refused",
        "nothing was computed: a wrong command line, an input that cannot be read, the IR gate or "
        "the reference-set audit",
    ),
    ExitState(3, "tool_error", "a tool or an optional package is missing or failed, by name"),
    ExitState(4, "conflict", "the data disagrees with the chosen naming reference"),
    ExitState(130, "cancelled", "stopped on request; what had finished stays in the run directory"),
)
_BY_CODE = {state.code: state for state in EXIT_STATES}
# The statuses of a step that left its results. A failed step was reached, and is not done.
_FINISHED_STEPS = ("ok", "cached", "prepared")


def exit_state(code: int) -> ExitState:
    return _BY_CODE[code]


def host_name() -> str:
    return socket.gethostname()


@dataclass(frozen=True, slots=True)
class RunState:
    """One run directory, as a reader finds it."""

    run_dir: Path
    run_id: str
    dataset: str
    kind: str  # "run" | "batch" | "unknown"
    # An exit-state name, or `running`, `incomplete` (an assembly that left no usable sequence),
    # `partial` (a batch with some samples ok), `interrupted` (the process died without saying so)
    # or `unknown` (the directory holds no record).
    status: str
    started_at: str | None = None
    finished_at: str | None = None
    exit_code: int | None = None
    done: int = 0
    total: int | None = None
    current: str | None = None
    command: tuple[str, ...] = ()
    reason: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {**asdict(self), "run_dir": str(self.run_dir), "command": list(self.command)}


def _process_is_this_run(payload: dict[str, Any]) -> bool:
    """Whether the process a record names is still alive, and is still this toolkit.

    A process number is only meaningful on the machine that issued it, and a container restart
    reissues them from the start, so the host has to match. A number can also be reused, which is
    what the command-line check is for where the system can answer it.
    """
    pid = payload.get("pid")
    if not isinstance(pid, int) or payload.get("host") != host_name():
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    try:
        command_line = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return True
    return b"plastid" in command_line


def _read_record(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def read_state(run_dir: Path) -> RunState:
    """What a directory says about its run. ``running`` is believed only while the process lives."""
    manifest_path = run_dir / "manifest.json"
    batch_path = run_dir / "batch.json"
    kind = "batch" if batch_path.is_file() else "run"
    path = batch_path if kind == "batch" else manifest_path
    dataset = run_dir.parent.name
    if not path.is_file():
        return RunState(
            run_dir,
            run_dir.name,
            dataset,
            "unknown",
            "unknown",
            reason="no manifest.json or batch.json in this directory",
        )
    payload = _read_record(path)
    if payload is None:
        return RunState(
            run_dir,
            run_dir.name,
            dataset,
            "unknown",
            "unknown",
            reason=f"{path.name} cannot be read",
        )

    status = str(payload.get("status", "unknown"))
    reason = ""
    if status == "running" and not _process_is_this_run(payload):
        status = "interrupted"
        reason = "the record says running and its process is gone: it stopped without finishing"
    current = payload.get("current")
    if kind == "batch":
        samples = payload.get("samples", [])
        done = sum(1 for sample in samples if sample.get("status") == "ok")
        total = payload.get("samples_total", len(samples))
    else:
        entries = payload.get("modules", [])
        done = sum(1 for entry in entries if entry.get("status") in _FINISHED_STEPS)
        total = current.get("of") if isinstance(current, dict) else len(entries)
    exit_code = payload.get("exit_code")
    return RunState(
        run_dir=run_dir,
        run_id=str(payload.get("run_id", run_dir.name)),
        dataset=str(payload.get("dataset") or dataset),
        kind=kind,
        status=status,
        started_at=payload.get("started_at"),
        finished_at=payload.get("finished_at"),
        exit_code=exit_code if isinstance(exit_code, int) else None,
        done=done,
        total=total if isinstance(total, int) else None,
        current=current.get("module") if isinstance(current, dict) else None,
        command=tuple(str(part) for part in payload.get("command", ())),
        reason=reason,
    )


def _is_run_directory(path: Path) -> bool:
    return (
        path.is_dir()
        and not path.name.startswith(".")
        and any((path / name).exists() for name in ("manifest.json", "batch.json", "00_run_info"))
    )


def list_runs(root: Path) -> list[RunState]:
    """Every run under an output root, newest first. A run id starts with its date and time."""
    found = (
        read_state(run_dir)
        for dataset in sorted(root.iterdir())
        if dataset.is_dir() and not dataset.name.startswith(".")
        for run_dir in sorted(dataset.iterdir())
        if _is_run_directory(run_dir)
    )
    return sorted(found, key=lambda state: state.run_id, reverse=True)
