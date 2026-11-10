"""What the server does for the page, with no HTTP in it.

The workspace is the one folder the interface can see, and every path that crosses this layer is
checked against it. A run is started as the `plastid` command line in its own process, so closing
the page stops nothing, and the run directory on disk stays the only record of what happened.
"""

from __future__ import annotations

import json
import os
import shutil
import signal
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock
from typing import Any
from uuid import uuid4

from plastid.core import plastid_reference
from plastid.run.state import EXIT_STATES, exit_state, list_runs, read_state
from plastid.run.tools import ToolError, plastid_output, start_plastid
from plastid.ui import commands, reader
from plastid.ui.catalogue import PAGES
from plastid.ui.records import RECORD_SUFFIXES, RecordBook, brief

__all__ = ["App", "NotFound", "Refused"]

_PLACES = ("records", "reads", "references", "runs", "tools")
_STATE_DIR = ".plastid-ui"
_TEXT_SUFFIXES = (".md", ".txt", ".log", ".tsv", ".cfg", ".csv", ".json")
_TEXT_LIMIT = 200_000  # characters sent to the page; the file itself holds them all
_LOG_TAIL = 4_000  # characters of a command's output kept on the run page
_PREVIEW_WIDTH = 1600  # px; an R figure's PNG can be 20,000 px wide
_SEEN_LIMIT_S = 5  # a run started this close to a command is that command's run
_ON_A_RUN = ("figures", "package", "qc")  # commands that write into a run, not a new directory
_WRITES_A_RUN = (
    "run",
    "annotate",
    "assemble run",
    "assemble prepare",
    "assemble import",
    "reference audit",
    "editing-reference build",
    "ir-junction",
)


class NotFound(LookupError):
    """Nothing at that address."""


class Refused(ValueError):
    """A request the interface will not act on, in words the page can show."""


@dataclass(slots=True)
class _Job:
    id: str
    argv: list[str]
    process: Any
    log: Path
    started: float = field(default_factory=time.time)


class App:
    def __init__(self, workspace: Path, *, version: str = "", shown_as: str | None = None) -> None:
        self.workspace = workspace.resolve()
        self.version = version
        self.shown_as = shown_as or str(self.workspace)  # the name of the folder outside an image
        self.runs_dir = self.workspace / "runs"
        for place in _PLACES:
            (self.workspace / place).mkdir(parents=True, exist_ok=True)
        self.state_dir = self.workspace / _STATE_DIR
        self.state_dir.mkdir(exist_ok=True)
        self._records = RecordBook()
        self._cached: dict[str, Any] = {}
        self._jobs: dict[str, _Job] = {}
        self._lock = Lock()
        self._preview_lock = Lock()

    # Paths -----------------------------------------------------------------------------------

    def _inside(self, relative: str) -> Path:
        """A path the page named, only if it lies in the workspace."""
        path = (self.workspace / relative).resolve()
        if not path.is_relative_to(self.workspace):
            raise Refused("that path is outside the workspace")
        return path

    def relative(self, path: Path) -> str:
        return path.resolve().relative_to(self.workspace).as_posix()

    def readable(self, relative: str) -> Path:
        path = self._inside(relative)
        if not path.is_file():
            raise NotFound(f"no file at {relative}")
        return path

    def upload_target(self, name: str, to: str) -> Path:
        """Where a dropped file lands. Its name is checked, and nothing is overwritten."""
        bare = Path(name).name
        if not bare or bare != name or bare.startswith("."):
            raise Refused("that file name cannot be used")
        folder = self._inside(to or "records")
        if not folder.is_dir():
            raise Refused(f"no folder {to}")
        target = folder / bare
        if target.exists():
            raise Refused(f"{bare} is already there; files here are not overwritten")
        return target

    def _tool_dirs(self) -> list[Path]:
        """Where the person keeps a tool of their own: `tools`, and a program's `bin` inside it.

        The image cannot carry a tool with no licence, and someone on Windows has no PATH to set.
        """
        tools = self.workspace / "tools"
        return [tools, *sorted(tools.glob("*/bin"))]

    # The engine's own lists ------------------------------------------------------------------

    def about(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "workspace": self.shown_as,
            "runs": self.relative(self.runs_dir),
            "pages": PAGES,
            "exit_states": [
                {"code": state.code, "status": state.status, "meaning": state.meaning}
                for state in EXIT_STATES
            ],
        }

    def _listing(self, key: str, args: tuple[str, ...]) -> Any:
        with self._lock:
            if key in self._cached:
                return self._cached[key]
        value = json.loads(plastid_output(args, cwd=self.workspace, tool_dirs=self._tool_dirs()))
        with self._lock:
            self._cached[key] = value
        return value

    def modules(self) -> Any:
        return self._listing("modules", ("modules", "--json"))

    def tools(self) -> Any:
        return self._listing("tools", ("tools", "--json"))

    def refresh_tools(self, form: Mapping[str, Any]) -> Any:
        with self._lock:
            self._cached.pop("tools", None)
        return self.tools()

    # Records ---------------------------------------------------------------------------------

    def records(self, folder: str = "records") -> dict[str, Any]:
        """The records under one folder of the workspace, for its list or a reference set."""
        base = self._inside(folder)
        if not base.is_dir():
            raise NotFound(f"no folder {folder}")
        return {
            "records": [
                {"path": self.relative(path), "name": path.name, "bytes": path.stat().st_size}
                for path in sorted(base.rglob("*"))
                if path.is_file() and path.suffix.lower() in RECORD_SUFFIXES
            ]
        }

    def reference_set(self, folder: str = "references") -> dict[str, Any]:
        """The records a reference set holds, as `annotate` will search them, for Annotate to tick.

        Only the files a run searches are listed (``plastid_reference.record_files``), so every box
        the page draws is a record the run reads. A row's ``key`` is the name ``--exclude`` takes:
        the file's stem, which is not the accession when a file is named otherwise.
        """
        base = self._inside(folder)
        if base.is_file():
            raise Refused(f"{folder} is a file; choose the folder that holds it")
        if not base.is_dir():
            raise NotFound(f"no folder {folder}")
        rows: list[dict[str, Any]] = []
        for path in plastid_reference.record_files(base):
            row: dict[str, Any] = {"path": self.relative(path), "name": path.name, "key": path.stem}
            try:
                row |= brief(path)
            except (OSError, ValueError) as error:
                row["unreadable"] = str(error)
            rows.append(row)
        return {"records": rows}

    def record(self, path: str) -> dict[str, Any]:
        file = self.readable(path)
        if file.suffix.lower() not in RECORD_SUFFIXES:
            raise Refused(f"{file.name} is not a GenBank or FASTA record")
        try:
            return {"path": path, **self._records.read(file)}
        except ValueError as error:
            return {"path": path, "unreadable": str(error)}

    def adopt(self, form: Mapping[str, Any]) -> dict[str, str]:
        """Copy a record a run wrote into the records, so the route can go on from it."""
        source = self.readable(str(form.get("path", "")))
        name = str(form.get("name", ""))
        if Path(name).suffix.lower() not in RECORD_SUFFIXES:
            raise Refused(f"{name} is not a GenBank or FASTA record")
        target = self.upload_target(name, "records")
        shutil.copyfile(source, target)
        return {"path": self.relative(target)}

    def browse(self, path: str = "") -> dict[str, Any]:
        """One folder of the workspace, for the file picker: the server reads the disk."""
        folder = self._inside(path)
        if not folder.is_dir():
            raise NotFound(f"no folder {path}")
        entries = []
        for item in sorted(folder.iterdir(), key=lambda p: (p.is_file(), p.name.lower())):
            if item.name.startswith("."):
                continue
            entries.append(
                {
                    "name": item.name,
                    "path": self.relative(item),
                    "folder": item.is_dir(),
                    "bytes": item.stat().st_size if item.is_file() else None,
                }
            )
        return {
            "path": self.relative(folder) if folder != self.workspace else "",
            "entries": entries,
        }

    # Runs ------------------------------------------------------------------------------------

    def runs(self) -> dict[str, Any]:
        if not self.runs_dir.is_dir():
            return {"runs": []}
        return {"runs": [reader.state_row(state) for state in list_runs(self.runs_dir)]}

    def _run_dir(self, dataset: str, run_id: str) -> Path:
        folder = (self.runs_dir / dataset / run_id).resolve()
        if not folder.is_relative_to(self.runs_dir) or not folder.is_dir():
            raise NotFound(f"no run {dataset}/{run_id}")
        return folder

    def run(self, dataset: str, run_id: str) -> dict[str, Any]:
        return reader.run_detail(self._run_dir(dataset, run_id), self.workspace)

    def cancel_run(self, form: Mapping[str, Any], dataset: str, run_id: str) -> dict[str, Any]:
        """Ask a run to stop, wherever it was started. It closes its own record, as for Ctrl-C."""
        folder = self._run_dir(dataset, run_id)
        state = read_state(folder)
        if state.status != "running":
            raise Refused(f"that run is {state.status}, not running")
        pid = reader.record_of(folder).get("pid")
        if not isinstance(pid, int):
            raise Refused("that run does not name its process")
        os.kill(pid, signal.SIGTERM)
        return {"asked": True}

    def figures(self, run: str) -> Any:
        """What a run can draw, why not, and the drawing options: the engine's own list."""
        folder = self._inside(run)
        if not folder.is_dir():
            raise NotFound(f"no run {run}")
        args = ("figures", self.relative(folder), "--list", "--json")
        try:
            return json.loads(plastid_output(args, cwd=self.workspace, tool_dirs=self._tool_dirs()))
        except ToolError as error:
            raise Refused(str(error)) from None

    def table(self, path: str) -> dict[str, Any]:
        return reader.table(self.readable(path))

    def text(self, path: str) -> dict[str, Any]:
        file = self.readable(path)
        if file.suffix.lower() not in _TEXT_SUFFIXES:
            raise Refused(f"{file.name} is not a text file")
        return {"text": file.read_text(encoding="utf-8", errors="replace")[:_TEXT_LIMIT]}

    def preview(self, relative: str) -> Path:
        """A figure small enough for a browser. Without Pillow the file itself is sent."""
        path = self.readable(relative)
        if path.suffix.lower() != ".png":
            return path
        try:
            from PIL import Image
        except ImportError:
            return path
        Image.MAX_IMAGE_PIXELS = None  # the toolkit's own R figures reach 235 million pixels
        target = self.state_dir / "previews" / f"{path.stat().st_mtime_ns}-{path.name}"
        # One at a time: that many pixels take most of a gigabyte each while they are decoded.
        with self._preview_lock:
            if target.is_file():
                return target
            target.parent.mkdir(exist_ok=True)
            with Image.open(path) as image:
                if image.width <= _PREVIEW_WIDTH:
                    return path
                height = max(1, image.height * _PREVIEW_WIDTH // image.width)
                image.convert("RGB").resize((_PREVIEW_WIDTH, height)).save(target)
        return target

    # Starting work ---------------------------------------------------------------------------

    def _with_defaults(self, command: str, form: Mapping[str, Any]) -> dict[str, Any]:
        filled = dict(form)
        if command in _WRITES_A_RUN:
            filled.setdefault("out", self.relative(self.runs_dir))
        return filled

    def command(self, form: Mapping[str, Any]) -> dict[str, Any]:
        """The line a form would run, shown beside it as it is filled in."""
        command = str(form.get("command", ""))
        try:
            args = commands.build(command, self._with_defaults(command, form.get("form") or {}))
        except commands.FormError as error:
            return {"incomplete": str(error)}
        return {"line": commands.shown(args), "args": args}

    def start_job(self, form: Mapping[str, Any]) -> dict[str, Any]:
        command = str(form.get("command", ""))
        try:
            args = commands.build(command, self._with_defaults(command, form.get("form") or {}))
        except commands.FormError as error:
            raise Refused(str(error)) from None
        job_id = uuid4().hex[:12]
        log = self.state_dir / "jobs" / f"{job_id}.log"
        log.parent.mkdir(exist_ok=True)
        process = start_plastid(args, cwd=self.workspace, log=log, tool_dirs=self._tool_dirs())
        job = _Job(job_id, args, process, log)
        with self._lock:
            self._jobs[job_id] = job
        return self._job_view(job)

    def _job_run(self, job: _Job) -> dict[str, str] | None:
        """The run this command made or worked on.

        A new run is found by the process number it recorded. A batch comes first, because each
        of its samples is a run under the same process.
        """
        if job.argv[0] in _ON_A_RUN:
            parts = Path(job.argv[1]).parts[-2:]
            return {"dataset": parts[0], "run_id": parts[1]} if len(parts) == 2 else None
        if not self.runs_dir.is_dir():
            return None
        for name in ("batch.json", "manifest.json"):
            for path in self.runs_dir.glob(f"*/*/{name}"):
                if path.stat().st_mtime < job.started - _SEEN_LIMIT_S:
                    continue
                try:
                    record = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    continue
                if record.get("pid") == job.process.pid:
                    return {"dataset": path.parent.parent.name, "run_id": path.parent.name}
        return None

    def _job_view(self, job: _Job) -> dict[str, Any]:
        code = job.process.poll()
        text = job.log.read_text(encoding="utf-8", errors="replace") if job.log.is_file() else ""
        return {
            "id": job.id,
            "line": commands.shown(job.argv),
            "running": code is None,
            "exit_code": code,
            "status": "running" if code is None else exit_state(code).status,
            "output": text[-_LOG_TAIL:],
            "run": self._job_run(job),
        }

    def jobs(self) -> dict[str, Any]:
        with self._lock:
            jobs = list(self._jobs.values())
        return {"jobs": [self._job_view(job) for job in jobs]}

    def job(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs.get(job_id)
        if job is None:
            raise NotFound(f"no job {job_id}")
        return self._job_view(job)
