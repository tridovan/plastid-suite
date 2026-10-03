"""The per-module run context: output paths, the IR gate, and the tool door.

A module receives exactly this and nothing else — no cwd, no globals (seam 1/5).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from plastid.core.ir import IRGateError, IRResult
from plastid.run.layout import RunLayout
from plastid.run.tools import ToolRecord, run_tool

__all__ = ["RunContext"]


@dataclass
class RunContext:
    layout: RunLayout
    params: dict[str, Any]
    module_id: str
    module_order: int
    ir: IRResult | None = None
    tools: list[ToolRecord] = field(default_factory=list)

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
                f"{self.module_id}: IR gate failed ({self.ir.identity * 100:.2f} % identity) — "
                f"position-based analysis must not run"
            )
        return self.ir

    def tool(self, name: str, argv: Sequence[str], **kwargs: Any) -> ToolRecord:
        """Run an external tool through the single door and keep its provenance."""
        record = run_tool(name, argv, log_dir=self.layout.logs_dir, **kwargs)
        self.tools.append(record)
        return record
