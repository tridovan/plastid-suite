"""The single door for external tools (seam 3).

``run_tool()`` is the only place in the package allowed to call ``subprocess`` and the only
place that records provenance: name · version · path · argv · time · exit code · logs.
CGAS's failure mode was six different call styles and tool versions thrown away.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

__all__ = ["ToolError", "ToolRecord", "run_tool", "tool_version", "which"]

# Flags tried in order. ``-version`` was added 03/10/2026 after measuring BLAST+: it accepts only
# the single-dash spelling and exits 1 for ``--version``, ``-v`` and ``-V`` — so the version of the
# tool that now annotates mitogenomes was recorded as "unknown" while every probe "succeeded" in
# failing. The probe runs in a throwaway directory for exactly this class of surprise.
_VERSION_FLAGS = ("--version", "-version", "-v", "-V")
_TAIL_LINES = 20


class ToolError(Exception):
    """An external tool exited outside its expected codes."""


@dataclass(frozen=True, slots=True)
class ToolRecord:
    """Provenance for one tool invocation — this is what Materials & Methods needs."""

    name: str
    kind: str  # "binary" | "perl" | "vendored"
    version: str | None
    version_source: str  # "flag" | "stdout-header" | "pinned" | "unknown"
    path: str | None
    argv: tuple[str, ...]
    cwd: str
    rc: int
    started_at: str
    duration_s: float
    stdout_file: str
    stderr_file: str
    stdout_tail: str
    interpreter: str | None = None  # e.g. "perl 5.38.2" for script-based tools


def which(name: str, extra: Sequence[Path] = ()) -> str | None:
    """Locate a tool on PATH or in explicit extra directories.

    PATH alone is not enough (``misa.pl`` lives outside it).
    """
    found = shutil.which(name)
    if found:
        return found
    for directory in extra:
        candidate = Path(directory) / name
        if candidate.is_file():
            return str(candidate)
    return None


def tool_version(
    name: str, kind: str, path: str | None, pinned: str | None = None
) -> tuple[str | None, str]:
    """Best-effort version: a version flag for binaries, ``$^V`` for perl, pinned for vendored.

    The probe *runs the tool*, so two rules apply, both measured 2026-10-03 on PGA:

    * a Perl script is asked for the **interpreter's** version, never probed with a flag. Passing
      an unknown flag makes the script run itself with its default arguments: PGA's next act is
      ``rm -rf gb`` — deleting a directory named ``gb`` in whatever directory the toolkit happened
      to be running from;
    * a binary is probed in a throwaway directory, for the same reason: ``--version`` is only a
      convention, and a tool that ignores it may do real work where it stands.
    """
    if kind == "vendored" or pinned is not None:
        return pinned, "pinned" if pinned is not None else "unknown"
    if kind == "perl":
        binary = which("perl") or "perl"
        try:
            result = subprocess.run(
                [binary, "-e", "print $^V"],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
        except OSError:
            return None, "unknown"
        output = result.stdout.strip()
        return (output, "stdout-header") if output else (None, "unknown")
    binary = path or which(name)
    if binary is None:
        return None, "unknown"
    with tempfile.TemporaryDirectory(prefix="version_probe_") as probe:
        for flag in _VERSION_FLAGS:
            try:
                result = subprocess.run(
                    [binary, flag],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False,
                    cwd=probe,
                )
            except OSError:
                return None, "unknown"
            if result.returncode == 0:
                # A line is a version only if it contains a digit. Measured 03/10/2026 on
                # ARAGORN: its ``-V`` is not a version flag, it exits **0** and prints
                # "No sequence file specified, type aragorn -h for help" — which was recorded as
                # the version of the tool that calls tRNA genes. A plausible-looking string in a
                # provenance record is worse than "unknown": it cannot be noticed.
                for line in (result.stdout or result.stderr).splitlines():
                    if line.strip() and any(character.isdigit() for character in line):
                        return line.strip(), "flag"
    return None, "unknown"


def run_tool(
    name: str,
    argv: Sequence[str],
    *,
    kind: str = "binary",
    log_dir: str | Path,
    pinned_version: str | None = None,
    timeout: int | None = None,
    expect: tuple[int, ...] = (0,),
    cwd: str | Path | None = None,
    tool_path: str | None = None,
    interpreter: str | None = None,
) -> ToolRecord:
    """Run one external tool, capture everything, raise :class:`ToolError` on a bad exit code.

    ``cwd`` defaults to the run directory; ``tool_path`` records the tool's own path when the
    executed binary is an interpreter (e.g. ``perl misa.pl``).
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    workdir = str(cwd) if cwd is not None else str(log_dir.parent.parent)
    resolved = tool_path or (which(argv[0]) if not Path(argv[0]).is_absolute() else argv[0])
    version, version_source = tool_version(name, kind, resolved, pinned_version)

    started = datetime.now()
    clock = time.monotonic()
    result = subprocess.run(
        list(argv),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
        cwd=workdir,
    )
    duration = time.monotonic() - clock

    safe_name = "".join(
        character if character.isalnum() or character in "-_" else "_" for character in name
    )
    # The index counts the calls already logged for this tool. The pattern must match the names
    # written below exactly — measured 2026-10-04: the first version asked for
    # ``*__{name}__*.stdout.log``, which does not match ``001__mafft__stdout.log`` (there is no
    # dot before ``stdout``), so the count stayed 0 and **every call overwrote the previous
    # one's logs**. Every run kept the last tool call's output and nothing else — which quietly
    # made "the raw output travels with the run" false for every tool but the last.
    index = len(list(log_dir.glob(f"*__{safe_name}__stdout.log"))) + 1
    stdout_file = log_dir / f"{index:03d}__{safe_name}__stdout.log"
    stderr_file = log_dir / f"{index:03d}__{safe_name}__stderr.log"
    stdout_file.write_text(result.stdout or "", encoding="utf-8")
    stderr_file.write_text(result.stderr or "", encoding="utf-8")

    tail = "\n".join((result.stdout or "").splitlines()[-_TAIL_LINES:])
    record = ToolRecord(
        name=name,
        kind=kind,
        version=version,
        version_source=version_source,
        path=resolved,
        argv=tuple(str(argument) for argument in argv),
        cwd=workdir,
        interpreter=interpreter,
        rc=result.returncode,
        started_at=started.isoformat(timespec="seconds"),
        duration_s=round(duration, 3),
        stdout_file=str(stdout_file),
        stderr_file=str(stderr_file),
        stdout_tail=tail,
    )
    if result.returncode not in expect:
        detail = tail or (result.stderr or "").strip()
        raise ToolError(
            f"{name} exited {result.returncode} (expected {list(expect)}): {detail[:300]}"
        )
    return record


def python_self() -> str:
    """The interpreter running the toolkit — recorded in run info."""
    return sys.executable
