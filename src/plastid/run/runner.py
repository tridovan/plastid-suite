"""Step execution: cache lookup → run the module → export → store.

The runner is the only place that turns a :class:`Result` into files; modules stay pure.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from plastid.core.ir import IRResult
from plastid.core.model import Genome
from plastid.result.export import export_csv, export_json
from plastid.result.result import Result
from plastid.run.cache import Cache
from plastid.run.context import RunContext

__all__ = ["StepOutcome", "execute_step"]

ReportRenderer = Callable[[Result], str]


@dataclass(frozen=True, slots=True)
class StepOutcome:
    status: str  # "ok" | "cached"
    key: str
    files: tuple[str, ...]  # absolute paths, deterministic order


def _params_hash(params: dict) -> str:
    payload = json.dumps(params, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _input_hash(genome: Genome) -> str:
    return hashlib.sha256(genome.sequence.encode("ascii", "replace")).hexdigest()


def _ir_fingerprint(ir: IRResult | None) -> str:
    if ir is None:
        return ""
    return f"{ir.length}:{ir.identity:.6f}:{ir.source}"


def export_files(result: Result, render_report: ReportRenderer | None = None) -> dict[str, str]:
    """Logical name → text. Deterministic: table order, then result.json, then report."""
    files: dict[str, str] = {}
    for table in result.tables.values():
        files[f"tables/{table.id}.csv"] = export_csv(table)
    files["result.json"] = export_json(result)
    if render_report is not None:
        files["report.md"] = render_report(result)
    return files


def _write(ctx: RunContext, files: dict[str, str]) -> tuple[str, ...]:
    written: list[Path] = []
    for logical_with_ext, text in files.items():
        logical, _, ext = logical_with_ext.rpartition(".")
        path = ctx.out(logical, ext)
        path.write_text(text, encoding="utf-8")
        written.append(path)
    return tuple(str(path) for path in written)


def execute_step(
    *,
    module_id: str,
    module_version: str,
    run_fn: Callable[[Genome, dict, RunContext], Result],
    genome: Genome,
    params: dict,
    ctx: RunContext,
    cache: Cache,
    force: bool = False,
    render_report: ReportRenderer | None = None,
) -> StepOutcome:
    """Run one module step, or restore it from the cache when nothing relevant changed."""
    key = Cache.key(
        module_id=module_id,
        module_version=module_version,
        params_hash=_params_hash(params),
        input_sha256s=(_input_hash(genome),),
        ir_fingerprint=_ir_fingerprint(ctx.ir),
    )
    if not force:
        cached = cache.get(key)
        if cached is not None:
            return StepOutcome(status="cached", key=key, files=_write(ctx, cached))

    result = run_fn(genome, params, ctx)
    if result.module != module_id:
        raise ValueError(f"module {module_id!r} returned a Result for {result.module!r}")
    payload = export_files(result, render_report)
    files = _write(ctx, payload)
    cache.put(key, payload)
    return StepOutcome(status="ok", key=key, files=files)
