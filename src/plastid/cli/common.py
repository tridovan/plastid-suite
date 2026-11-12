"""What every command shares: exit codes, run-directory flags, `--set` and manifest bookkeeping."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any

from plastid.core.digest import sha256_file
from plastid.modules import registry
from plastid.result.result import Result
from plastid.run.context import RunContext
from plastid.run.figures import (
    read_run_results,
)
from plastid.run.layout import RunLayout
from plastid.run.manifest import write_manifest
from plastid.run.runinfo import render_run_info
from plastid.run.runner import (
    DATA_LABEL_PARAM,
)
from plastid.run.state import exit_state, read_state

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_TOOL = 3
# 4 — a naming conflict stopped the step. Its own number, because it is not a usage error and not a
# missing tool: the run was well-formed and the data disagrees with the reference.
EXIT_CONFLICT = 4
# 130 — what a shell reports for a job stopped by Ctrl-C, here for every way of stopping a run.
EXIT_CANCELLED = 130

_OUT_HELP = "directory that holds every run (default: .runs)"
_LABEL_HELP = "where the data came from, printed with every table, e.g. 'NCBI RefSeq, 2026-10'"
_NO_CACHE_HELP = "compute every step again instead of reusing a stored result"
_SET_HELP = "a module parameter, repeatable; `plastid modules <module>` lists them"
_JSON_HELP = "print JSON instead of text"
_RUN_DIR_HELP = "a run directory, as `plastid runs` lists them"
_REFERENCE_SET_HELP = "directory of pinned GenBank records (default: .panel/plastid-reference)"


def _add_run_directory(
    parser: argparse.ArgumentParser, dataset_default: str, *, fixed: str | None = None
) -> None:
    """The three flags of every command that makes a run: where it goes and what it is."""
    parser.add_argument("-o", "--out", default=".runs", help=_OUT_HELP)
    parser.add_argument(
        "--dataset",
        default=fixed,
        metavar="NAME",
        help=f"name of this run's directory inside --out (default: {dataset_default})",
    )
    parser.add_argument("--label", metavar="TEXT", help=_LABEL_HELP)


def _parse_sets(pairs: list[str]) -> dict:
    params: dict = {}
    for pair in pairs:
        key, separator, raw = pair.partition("=")
        if not separator or not key:
            raise ValueError(f"--set expects KEY=VALUE, got {pair!r}")
        params[key] = _coerce(raw)
    return params


def _checked_sets(pairs: Sequence[str], module_ids: Sequence[str]) -> dict | None:
    """The ``--set`` values of a run, or None after saying what is wrong with every one of them.

    Measured before this check: ``--set this_key_does_not_exist=42`` ran and exited 0, so a
    misspelt ``min_length_bp`` looked like a successful override. A key a module accepts but
    ignores is the same failure from the other side, which is why each module declares what it
    reads and a test compares that declaration with its source. A value is checked here too, so
    a wrong one is a refused command line and not a Python error from inside a module.
    """
    try:
        params = _parse_sets(list(pairs))
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return None
    problems = [
        *registry.unknown_params(module_ids, params),
        *registry.refused_values(module_ids, params),
    ]
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return None
    return registry.anchored_paths(registry.declared_params(module_ids), params)


def _refused_for_own_command(specs: Sequence[registry.ModuleSpec]) -> bool:
    """Whether any of these modules has to run through a command of its own, saying which."""
    for spec in specs:
        if spec.command:
            print(
                f"{spec.id} runs through `{spec.command}`, which checks the input and records "
                f"the reference set — `plastid modules {spec.id}` shows how",
                file=sys.stderr,
            )
    return any(spec.command for spec in specs)


def _coerce(raw: str) -> object:
    for cast in (int, float):
        try:
            return cast(raw)
        except ValueError:
            continue
    return raw


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in ("1", "true", "yes", "on")


def _finish(manifest: dict, layout: RunLayout, ctx: RunContext, code: int) -> int:
    """Write the manifest and run info whatever happened — a failed run stays traceable."""
    manifest["finished_at"] = datetime.now().isoformat(timespec="seconds")
    manifest["exit_code"] = code
    manifest.pop("current", None)
    write_manifest(layout.run_dir / "manifest.json", manifest)
    (layout.run_info_dir / "run_info.md").write_text(
        render_run_info(manifest, ctx.tools), encoding="utf-8"
    )
    if code == EXIT_OK:
        steps = ", ".join(f"{entry['id']}: {entry['status']}" for entry in manifest["modules"])
        print(f"{steps} — {layout.run_dir}")
    return code


def _downstream(specs: Sequence[registry.ModuleSpec], order: int) -> list[tuple[str, str]]:
    """The modules still to run after the one at 1-based ``order`` — recorded as skipped."""
    return [(spec.id, spec.version) for spec in specs[order:]]


def _fail(
    manifest: dict,
    layout: RunLayout,
    ctx: RunContext,
    version: str,
    message: str,
    code: int,
    downstream: Sequence[tuple[str, str]] = (),
) -> int:
    """Record the failed step, then every step that therefore did not run.

    A failure stops the remaining modules; leaving them out of the manifest made "the run never
    reached them" look like "they were never asked for". ``downstream`` carries
    ``(module_id, version)`` for each, recorded as ``skipped`` with the reason.
    """
    print(message, file=sys.stderr)
    manifest["modules"].append(
        {
            "id": ctx.module_id,
            "version": version,
            "status": "failed",
            "error": message,
            "tools": [asdict(record) for record in ctx.tools],
            "files": [],
        }
    )
    for module_id, module_version in downstream:
        manifest["modules"].append(
            {
                "id": module_id,
                "version": module_version,
                "status": "skipped",
                "reason": f"skipped: {ctx.module_id} failed, so this step did not run",
                "tools": [],
                "files": [],
            }
        )
    manifest["status"] = exit_state(code).status
    return _finish(manifest, layout, ctx, code)


def _figures_stopped(
    manifest: dict, layout: RunLayout, ctx: RunContext, message: str, code: int
) -> int:
    """Figures that could not start. Every analysis finished, so none of them is recorded failed."""
    print(message, file=sys.stderr)
    print(f"the analyses finished — {layout.run_dir}", file=sys.stderr)
    manifest["figures"] = {"status": "failed", "error": message}
    manifest["status"] = exit_state(code).status
    return _finish(manifest, layout, ctx, code)


def _cancelled(
    manifest: dict,
    layout: RunLayout,
    ctx: RunContext,
    steps: Sequence[tuple[str, str]],
) -> int:
    """Close a run that was asked to stop. ``steps`` is every ``(id, version)`` the run planned.

    The first step with no entry is the one the signal reached. Its tool calls so far stay on the
    record, and what had finished stays in the directory, because a long run stopped near its end
    still holds work somebody may want to read.
    """
    print(f"cancelled — what had finished is in {layout.run_dir}", file=sys.stderr)
    done = {entry["id"] for entry in manifest["modules"]}
    waiting = [step for step in steps if step[0] not in done]
    for position, (step_id, version) in enumerate(waiting):
        manifest["modules"].append(
            {
                "id": step_id,
                "version": version,
                "status": "skipped" if position else "cancelled",
                "reason": (
                    "skipped: the run was cancelled before this step"
                    if position
                    else "cancelled while this step ran"
                ),
                "tools": [] if position else [asdict(record) for record in ctx.tools],
                "files": [],
            }
        )
    manifest["status"] = exit_state(EXIT_CANCELLED).status
    return _finish(manifest, layout, ctx, EXIT_CANCELLED)


def _file_entry(path: Path, *, relative_to: Path | None = None, module: str | None = None) -> dict:
    """One file as a manifest lists it: path (relative to a base when given), step, digest, size."""
    entry: dict[str, object] = {
        "path": path.relative_to(relative_to).as_posix() if relative_to is not None else str(path)
    }
    if module is not None:
        entry["module"] = module
    entry["sha256"] = sha256_file(path)
    entry["bytes"] = path.stat().st_size
    return entry


def _stored_results(run_dir: Path) -> dict[str, Result]:
    """What a run stored. A run stopped or refused before its first step finished holds nothing."""
    return read_run_results(run_dir) if any(run_dir.glob("*/result.json")) else {}


def _nothing_finished(run_dir: Path) -> bool:
    """Whether this run finished no step, saying so: it is a state of the run, not a wrong path."""
    if any(run_dir.glob("*/result.json")):
        return False
    status = read_state(run_dir).status
    print(f"{run_dir} finished no step (status {status}), so it holds no tables", file=sys.stderr)
    return True


def _run_dir_or_error(path: str, *, allow_batch: bool = False) -> Path | None:
    """A run directory (or, where the caller can use one, a batch directory)."""
    run_dir = Path(path)
    is_run = (run_dir / "manifest.json").is_file()
    is_batch = allow_batch and (run_dir / "batch.json").is_file()
    if (is_run or is_batch) and read_state(run_dir).status == "running":
        # The run is rewriting its record, and what its tables hold so far is not a result.
        print(f"{run_dir} is still running — wait for it, or stop it, first", file=sys.stderr)
        return None
    if is_run or is_batch:
        return run_dir
    print(
        f"{run_dir} is not a run directory — no manifest.json"
        + (" or batch.json" if allow_batch else "")
        + " in it (run `plastid run` first)",
        file=sys.stderr,
    )
    return None


def _manifest_provenance(manifest: Mapping[str, Any]) -> tuple[str, ...]:
    """The run-level lines every figure and package must carry.

    A figure or a supplement made from third-party records must say whose records, with which
    data label, and on which IR path — not leave the reader to open ``manifest.json``. Only what
    the run recorded is printed; nothing is inferred.
    """
    lines: list[str] = []
    label = (manifest.get("params") or {}).get(DATA_LABEL_PARAM)
    if label:
        lines.append(f"data label: {label}")
    for entry in manifest.get("inputs", []):
        digest = str(entry.get("sha256", ""))
        lines.append(f"input: {entry.get('path')} (sha256 {digest[:16]}…)")
    ir = manifest.get("ir") or {}
    if ir:
        gate = "passed" if ir.get("gate_passed") else "failed"
        lines.append(f"IR: source {ir.get('source')} · identity {ir.get('identity')} · gate {gate}")
    return tuple(lines)


def _batch_provenance(payload: Mapping[str, Any]) -> tuple[str, ...]:
    """The batch-level lines: the data label, the accessions, and the IR statuses.

    A multi-sample figure cannot name one accession, so it names the set. The per-sample mapping
    (accession · IR source · gate) is ``samples.csv`` beside ``batch.json``; the summary here is
    what fits in a caption without pretending the set is one record.
    """
    lines: list[str] = []
    label = payload.get("data_label") or (payload.get("params") or {}).get(DATA_LABEL_PARAM)
    if label:
        lines.append(f"data label: {label}")
    samples = list(payload.get("samples", []))
    lines.append(f"batch: {payload.get('dataset')} — {len(samples)} sample(s)")
    accessions = [str(sample.get("accession")) for sample in samples if sample.get("accession")]
    if accessions:
        lines.append("accessions: " + ", ".join(accessions))
    for key, label_text in (("ir_source", "IR source"), ("ir_gate", "IR gate")):
        counts: dict[str, int] = {}
        for sample in samples:
            value = str(sample.get(key) or "not recorded")
            counts[value] = counts.get(value, 0) + 1
        if counts:
            lines.append(
                f"{label_text}: "
                + ", ".join(f"{value} ×{count}" for value, count in counts.items())
            )
    lines.append("per-sample accession/IR detail: samples.csv")
    return tuple(lines)


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _print_json(payload: Any) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def _print_table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    widths = [max(len(row[i]) for row in (header, *rows)) for i in range(len(header) - 1)]
    for row in (header, *rows):
        cells = [cell.ljust(width) for cell, width in zip(row, widths, strict=False)]
        print("  ".join((*cells, row[-1])).rstrip())
