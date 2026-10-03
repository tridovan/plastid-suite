"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by returning an
exit code to the console-script wrapper). Library code raises exceptions instead.

Exit codes:
    0 — success
    1 — general error (bad input, unknown module, internal failure)
    2 — usage error, or a failed gate (e.g. the IR gate)
    3 — missing or failed external tool
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

from plastid import __version__
from plastid.core.genbank import read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, detect_ir, enforce_ir_gate
from plastid.core.normalize import normalize
from plastid.core.panel import (
    IndexSnapshot,
    entry_from_summary,
    load_index,
    pin,
    propose,
    save_index,
)
from plastid.modules import registry
from plastid.report.citations import load_citations
from plastid.report.engine import load_template, render
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.runinfo import render_run_info
from plastid.run.runner import execute_step
from plastid.run.tools import ToolError

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_TOOL = 3


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plastid",
        description="PlastidToolSuite — plastid and mitochondrial genome analysis.",
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    run_parser = subparsers.add_parser("run", help="run one analysis module on one genome")
    run_parser.add_argument("-m", "--module", required=True, help="module id (e.g. genome_stats)")
    run_parser.add_argument("-i", "--input", required=True, help="GenBank file (one record)")
    run_parser.add_argument("-o", "--out", default=".runs", help="output root (default: .runs)")
    run_parser.add_argument("--dataset", default=None, help="dataset label (default: accession)")
    run_parser.add_argument("--force", action="store_true", help="ignore the cache and re-run")
    run_parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="module parameter (repeatable)",
    )

    panel = subparsers.add_parser("panel", help="Reference Panel: propose and pin genomes")
    panel_actions = panel.add_subparsers(dest="panel_action", metavar="ACTION")
    refresh = panel_actions.add_parser("refresh", help="rebuild the RefSeq index from NCBI")
    refresh.add_argument("--taxon", required=True, help='e.g. "Iridaceae"')
    refresh.add_argument("--out", default=".panel", help="index directory (default: .panel)")
    refresh.add_argument("--retmax", type=int, default=200)
    propose_parser = panel_actions.add_parser("propose", help="list candidates (no network)")
    propose_parser.add_argument("--index", default=".panel/index.json")
    propose_parser.add_argument("--target", required=True, help="accession of the target genome")
    propose_parser.add_argument("--mode", choices=["A", "B"], default="B")
    propose_parser.add_argument("--level", default="family")
    propose_parser.add_argument("--per-tier", type=int, default=3)
    pin_parser = panel_actions.add_parser("pin", help="download and checksum the chosen genomes")
    pin_parser.add_argument("--index", default=".panel/index.json")
    pin_parser.add_argument("--accessions", required=True, help="comma-separated accessions")
    pin_parser.add_argument("--out", default=".panel/genomes")
    return parser


def _parse_sets(pairs: list[str]) -> dict:
    params: dict = {}
    for pair in pairs:
        key, separator, raw = pair.partition("=")
        if not separator or not key:
            raise ValueError(f"--set expects KEY=VALUE, got {pair!r}")
        params[key] = _coerce(raw)
    return params


def _coerce(raw: str) -> object:
    for cast in (int, float):
        try:
            return cast(raw)
        except ValueError:
            continue
    return raw


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _finish(manifest: dict, layout: RunLayout, ctx: RunContext, version: str, code: int) -> int:
    """Write the manifest and run info whatever happened — a failed run stays traceable."""
    manifest["finished_at"] = datetime.now().isoformat(timespec="seconds")
    write_manifest(layout.run_dir / "manifest.json", manifest)
    (layout.run_info_dir / "run_info.md").write_text(
        render_run_info(manifest, ctx.tools), encoding="utf-8"
    )
    if code == EXIT_OK:
        print(f"{ctx.module_id}: {manifest['modules'][-1]['status']} — {layout.run_dir}")
    return code


def _fail(
    manifest: dict,
    layout: RunLayout,
    ctx: RunContext,
    version: str,
    message: str,
    code: int,
) -> int:
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
    manifest["status"] = "failed"
    return _finish(manifest, layout, ctx, version, code)


def _cmd_run(args: argparse.Namespace) -> int:
    try:
        spec = registry.get(args.module)
    except registry.UnknownModuleError as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR

    try:
        genome = read_genbank(args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read {args.input}: {error}", file=sys.stderr)
        return EXIT_ERROR
    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE

    overlay = json.dumps(params, sort_keys=True)
    fingerprint = hashlib.sha256((genome.sequence + overlay).encode("ascii", "replace")).hexdigest()
    dataset = args.dataset or genome.accession or genome.name
    layout = RunLayout.create(args.out, dataset, fingerprint=fingerprint)
    layout.make_dirs()

    ctx = RunContext(layout=layout, params=params, module_id=spec.id, module_order=1)
    manifest = build_manifest(
        run_id=layout.run_id,
        command=["plastid", "run", "--module", spec.id, "--input", str(args.input)],
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params=params,
        inputs=[
            {
                "path": str(args.input),
                "sha256": _sha256_file(Path(args.input)),
                "bytes": Path(args.input).stat().st_size,
            }
        ],
        coordinate_system={
            "name": "genbank",
            "rotation_offset": 0,
            "note": (
                "per-genome tables use GenBank coordinates; cross-genome comparisons use the "
                "normalized frame (LSC -> IRa -> SSC -> IRb)"
            ),
        },
    )

    if spec.requires_ir:
        try:
            ir = enforce_ir_gate(detect_ir(genome))
        except (IRGateError, IRDetectionError) as error:
            manifest["ir"] = {"gate_passed": False, "error": str(error)}
            return _fail(manifest, layout, ctx, spec.version, f"IR gate: {error}", EXIT_USAGE)
        ctx.ir = ir
        manifest["ir"] = {
            "source": ir.source,
            "length": ir.length,
            "identity": ir.identity,
            "gate_passed": ir.gate_passed,
            "mirror_constant": ir.mirror_constant,
            "denovo_length": ir.denovo_length,
            "denovo_delta": ir.denovo_delta,
        }
        manifest["coordinate_system"]["rotation_offset"] = normalize(genome, ir).rotation_offset

    try:
        template = load_template(f"{spec.id}.en.md")
    except FileNotFoundError:
        template = None
    citations = load_citations()
    render_report = (
        (lambda result: render(template, result, citations=citations)) if template else None
    )

    try:
        outcome = execute_step(
            module_id=spec.id,
            module_version=spec.version,
            run_fn=spec.run,
            genome=genome,
            params=params,
            ctx=ctx,
            cache=Cache(layout.cache_dir),
            force=args.force,
            render_report=render_report,
        )
    except IRGateError as error:
        return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_USAGE)
    except ToolError as error:
        return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_TOOL)
    except Exception as error:  # noqa: BLE001 — report any module failure, then exit 1
        message = f"{type(error).__name__}: {error}"
        return _fail(manifest, layout, ctx, spec.version, message, EXIT_ERROR)

    files = []
    for path_text in outcome.files:
        path = Path(path_text)
        files.append(
            {
                "path": str(path.relative_to(layout.run_dir)),
                "module": spec.id,
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    manifest["modules"].append(
        {
            "id": spec.id,
            "version": spec.version,
            "status": outcome.status,
            "key": outcome.key,
            "tools": [asdict(record) for record in ctx.tools],
            "files": [entry["path"] for entry in files],
        }
    )
    manifest["files"] = files
    manifest["status"] = "ok"
    return _finish(manifest, layout, ctx, spec.version, EXIT_OK)


def _cmd_panel(args: argparse.Namespace) -> int:
    action = getattr(args, "panel_action", None)
    if action == "refresh":
        return _panel_refresh(args)
    if action == "propose":
        return _panel_propose(args)
    if action == "pin":
        return _panel_pin(args)
    print("panel: choose an action — refresh | propose | pin", file=sys.stderr)
    return EXIT_USAGE


def _panel_refresh(args: argparse.Namespace) -> int:
    client = NcbiClient()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        ids = client.search_refseq_plastomes(args.taxon, retmax=args.retmax)
        summaries = client.summaries(ids)
        lineages = client.lineages([str(record.get("taxid", "")) for record in summaries])
    except (NcbiError, OSError) as error:
        print(f"NCBI: {error}", file=sys.stderr)
        return EXIT_ERROR

    entries = []
    for summary in summaries:
        lineage = lineages.get(str(summary.get("taxid", "")), "")
        try:
            entries.append(entry_from_summary(summary, lineage))
        except ValueError as error:
            print(f"skipped: {error}", file=sys.stderr)
    snapshot = IndexSnapshot(
        retrieved_at=datetime.now().isoformat(timespec="seconds"), entries=tuple(entries)
    )
    save_index(out_dir / "index.json", snapshot)
    print(f"panel: {len(entries)} RefSeq genome(s) for {args.taxon} → {out_dir / 'index.json'}")
    return EXIT_OK


def _panel_propose(args: argparse.Namespace) -> int:
    try:
        snapshot = load_index(args.index)
    except (OSError, KeyError) as error:
        print(
            f"cannot read {args.index}: {error} — run 'plastid panel refresh' first",
            file=sys.stderr,
        )
        return EXIT_ERROR
    try:
        proposal = propose(
            snapshot, args.target, mode=args.mode, level=args.level, per_tier=args.per_tier
        )
    except (ValueError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR

    print(
        f"target: {proposal.target.accession} ({proposal.target.organism}) — mode {proposal.mode}"
    )
    for candidate in proposal.candidates:
        entry = candidate.entry
        print(
            f"  tier {candidate.tier}  {candidate.tier_label:30}  {entry.accession:14}  "
            f"{entry.organism:35}  {entry.length} bp"
        )
    for note in proposal.notes:
        print(f"  note: {note}")
    print("the engine only proposes — pick the final list, then run 'plastid panel pin'")
    return EXIT_OK


def _panel_pin(args: argparse.Namespace) -> int:
    try:
        snapshot = load_index(args.index)
    except (OSError, KeyError) as error:
        print(
            f"cannot read {args.index}: {error} — run 'plastid panel refresh' first",
            file=sys.stderr,
        )
        return EXIT_ERROR
    wanted = [accession.strip() for accession in args.accessions.split(",") if accession.strip()]
    by_accession = {entry.accession: entry for entry in snapshot.entries}
    missing = [accession for accession in wanted if accession not in by_accession]
    if missing:
        print(f"not in the index: {missing} — refresh the panel first", file=sys.stderr)
        return EXIT_ERROR

    client = NcbiClient()
    try:
        pinned = pin(
            (by_accession[accession] for accession in wanted),
            client.fetch_genbank,
            args.out,
            retrieved_at=snapshot.retrieved_at,
        )
    except (NcbiError, OSError) as error:
        print(f"NCBI: {error}", file=sys.stderr)
        return EXIT_ERROR
    for record in pinned:
        detail = (
            f"{record.length} bp, {record.cds} CDS"
            if record.length is not None
            else f"UNVERIFIED: {record.parse_error}"
        )
        digest = record.sha256[:12]
        print(f"pinned {record.accession}: {record.bytes} bytes · sha256 {digest}… ({detail})")
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(f"plastid {__version__}")
        return EXIT_OK

    if args.command is None:
        parser.print_help()
        return EXIT_OK

    if args.command == "run":
        return _cmd_run(args)

    if args.command == "panel":
        return _cmd_panel(args)

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
