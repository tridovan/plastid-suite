"""The per-module run context: output paths, the IR gate, and the tool door.

A module receives exactly this and nothing else — no cwd, no globals (seam 1/5).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from plastid.core.ir import IRGateError, IRResult, identity_percent
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolRecord, run_tool

__all__ = ["ComputedOutput", "RunContext"]


@dataclass(frozen=True, slots=True)
class ComputedOutput:
    """A tool's output this run already holds, with the call that produced it and the step's id."""

    stdout: str
    record: ToolRecord
    module_id: str


@dataclass
class RunContext:
    layout: RunLayout
    params: dict[str, Any]
    module_id: str
    module_order: int
    ir: IRResult | None = None
    tools: list[ToolRecord] = field(default_factory=list)
    # One tool of a long step can take hours, and a run that says nothing looks stopped.
    announce_tools: bool = False
    # Outputs this run's steps have already computed, by a digest of the call (tool, binary,
    # arguments, input). The runner hands every step of one run the same dict: three comparison
    # modules ask MAFFT for one alignment of the same rows, and one call serves them all. It lives
    # for one run only, so an output never outlives the binary that produced it.
    computed: dict[str, ComputedOutput] = field(default_factory=dict)
    # This step's uses of an output an earlier step computed — recorded beside ``tools``, which
    # keeps meaning "the calls this step made".
    reused: list[dict[str, Any]] = field(default_factory=list)

    def out(self, logical: str, ext: str) -> Path:
        """Where a logical output goes — the layout decides, the module never names paths."""
        path = self.layout.module_dir(self.module_order, self.module_id) / f"{logical}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def require_ir(self) -> IRResult:
        """Every position-based module calls this first; a failed gate stops the run."""
        if self.ir is None:
            raise IRGateError(
                f"{self.module_id}: needs an IR result — the run never determined one"
            )
        if not self.ir.gate_passed:
            raise IRGateError(
                f"{self.module_id}: IR gate failed "
                f"({identity_percent(self.ir.identity):.2f} % identity) — "
                f"position-based analysis must not run"
            )
        return self.ir

    def tool(self, name: str, argv: Sequence[str], **kwargs: Any) -> ToolRecord:
        """Run an external tool through the single door and keep its provenance."""
        if self.announce_tools and name not in {record.name for record in self.tools}:
            print(f"[{name}]", flush=True)
        record = run_tool(name, argv, log_dir=self.layout.logs_dir, **kwargs)
        self.tools.append(record)
        return record

    def computed_output(self, key: str) -> str | None:
        """An output an earlier step of this run computed for the same call, recorded as reused."""
        found = self.computed.get(key)
        if found is None:
            return None
        self.reused.append(
            {
                "name": found.record.name,
                "argv": list(found.record.argv),
                "computed_by": found.module_id,
                "stdout_file": found.record.stdout_file,
            }
        )
        return found.stdout

    def keep_output(self, key: str, stdout: str, record: ToolRecord) -> None:
        """Hold one call's output for the later steps of this run."""
        self.computed[key] = ComputedOutput(stdout=stdout, record=record, module_id=self.module_id)
