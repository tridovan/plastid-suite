"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by returning an
exit code to the console-script wrapper). Library code raises exceptions instead.

Exit codes:
    0 — success
    1 — general error (bad input, unknown module, internal failure)
    2 — usage error, or a failed gate (e.g. the IR gate)
    3 — missing or failed external tool, or a missing optional dependency
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
from plastid.core.ir import IRDetectionError, IRGateError, IRResult, detect_ir, enforce_ir_gate
from plastid.core.normalize import normalize
from plastid.core.panel import (
    IndexSnapshot,
    PanelEntry,
    entries_from_summaries,
    load_index,
    pin,
    propose,
    save_index,
)
from plastid.core.samples import group_of, load_groups
from plastid.modules import registry
from plastid.report.citations import Citation, format_citation, load_citations
from plastid.report.engine import CitationError, load_template, render
from plastid.result.docx_export import DocxMeta, ExportDependencyError, export_docx
from plastid.result.result import Result
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.runinfo import render_run_info
from plastid.run.runner import DocumentRenderer, execute_step
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
    run_parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        help="GenBank file (one record); repeat the flag for a multi-genome module",
    )
    run_parser.add_argument("-o", "--out", default=".runs", help="output root (default: .runs)")
    run_parser.add_argument("--dataset", default=None, help="dataset label (default: accession)")
    run_parser.add_argument(
        "--groups",
        default=None,
        metavar="FILE",
        help="sample-group YAML ({group: [accessions]}) — labels every input (P2-11)",
    )
    run_parser.add_argument(
        "--docx",
        action="store_true",
        help="also write report.docx — a publication document (needs the 'docx' extra)",
    )
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
    pin_parser.add_argument(
        "--index",
        default=None,
        help="panel index to take entries from; without it every accession is resolved from NCBI",
    )
    wanted = pin_parser.add_mutually_exclusive_group(required=True)
    wanted.add_argument("--accessions", help="comma-separated accessions")
    wanted.add_argument(
        "--groups",
        metavar="FILE",
        help="sample-group YAML ({group: [accessions]}) — pins every accession in it",
    )
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


def _document_renderer(
    *,
    spec: registry.ModuleSpec,
    layout: RunLayout,
    manifest: dict,
    citations: dict[str, Citation],
    accession: str,
) -> DocumentRenderer:
    """Build the ``--docx`` renderer: provenance from the run, references from the registry.

    Everything that reaches the page is either a table cell or one of these labelled strings, so
    the document cannot carry a number that no run produced.
    """
    inputs = [
        ("Input", f"{entry['path']} · sha256 {entry['sha256'][:16]}…")
        for entry in manifest["inputs"]
    ]
    coordinates = manifest["coordinate_system"]
    provenance = (
        ("Dataset", layout.dataset),
        ("Run", layout.run_id),
        *inputs,
        (
            "Coordinates",
            f"{coordinates['name']} · rotation offset {coordinates['rotation_offset']}",
        ),
        ("Module", f"{spec.id} {spec.version}"),
        ("Generated by", f"PlastidToolSuite {__version__}"),
    )

    def render_document(result: Result, prose: str | None) -> bytes:
        references = []
        for key in result.citations:
            if key not in citations:
                raise CitationError(f"citation key {key!r} is not in the registry")
            references.append(format_citation(citations[key]))
        return export_docx(
            result,
            DocxMeta(
                title=spec.title,
                label=f"{accession} · {spec.id}",
                prose=prose or "",
                provenance=provenance,
                references=tuple(references),
            ),
        )

    return render_document


def _cmd_run(args: argparse.Namespace) -> int:
    try:
        spec = registry.get(args.module)
    except registry.UnknownModuleError as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR

    try:
        genomes = tuple(read_genbank(text) for text in args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read the input: {error}", file=sys.stderr)
        return EXIT_ERROR
    if len(genomes) != spec.arity:
        print(
            f"module {spec.id!r} takes {spec.arity} genome(s), got {len(genomes)} — "
            f"pass --input once per genome",
            file=sys.stderr,
        )
        return EXIT_USAGE
    genome = genomes[0]
    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE

    if args.groups:
        # Shared group labels for every comparison output (P2-11). Loaded here because a module
        # never touches the filesystem; the module receives plain data in params.
        try:
            groups = load_groups(args.groups)
        except (OSError, ValueError) as error:
            print(f"cannot read {args.groups}: {error}", file=sys.stderr)
            return EXIT_ERROR
        labels: dict[str, str] = {}
        unlabelled = []
        for item in genomes:
            key = item.accession or item.name
            label = group_of(groups, key)
            if label is None:
                unlabelled.append(key)
            else:
                labels[key] = label
        if unlabelled:
            print(
                f"{args.groups}: no group for {', '.join(unlabelled)} — every input needs one, "
                f"or the comparison output would carry an unlabelled row",
                file=sys.stderr,
            )
            return EXIT_ERROR
        params["groups"] = labels

    overlay = json.dumps(params, sort_keys=True)
    sequences = "".join(item.sequence + "\x00" for item in genomes)  # ordered: order is meaning
    fingerprint = hashlib.sha256((sequences + overlay).encode("ascii", "replace")).hexdigest()
    dataset = args.dataset or genome.accession or genome.name
    layout = RunLayout.create(args.out, dataset, fingerprint=fingerprint)
    layout.make_dirs()

    ctx = RunContext(layout=layout, params=params, module_id=spec.id, module_order=1)
    manifest = build_manifest(
        run_id=layout.run_id,
        command=[
            "plastid",
            "run",
            "--module",
            spec.id,
            *[part for text in args.input for part in ("--input", text)],
        ],
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params=params,
        inputs=[
            {
                "path": str(text),
                "sha256": _sha256_file(Path(text)),
                "bytes": Path(text).stat().st_size,
            }
            for text in args.input
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
        gated: list[IRResult] = []
        for candidate in genomes:
            try:
                gated.append(enforce_ir_gate(detect_ir(candidate)))
            except (IRGateError, IRDetectionError) as error:
                # every genome in a comparison must pass: one diverged IR corrupts them all
                manifest["ir"] = {"gate_passed": False, "error": str(error)}
                return _fail(manifest, layout, ctx, spec.version, f"IR gate: {error}", EXIT_USAGE)
        ctx.ir = gated[0]
        manifest["ir"] = {
            "source": gated[0].source,
            "length": gated[0].length,
            "identity": gated[0].identity,
            "gate_passed": gated[0].gate_passed,
            "mirror_constant": gated[0].mirror_constant,
            "denovo_length": gated[0].denovo_length,
            "denovo_delta": gated[0].denovo_delta,
        }
        manifest["coordinate_system"]["rotation_offset"] = normalize(
            genome, gated[0]
        ).rotation_offset
        if len(gated) > 1:
            manifest["ir"]["genomes"] = [
                {
                    "accession": candidate.accession or candidate.name,
                    "length": result.length,
                    "identity": result.identity,
                    "source": result.source,
                }
                for candidate, result in zip(genomes, gated, strict=True)
            ]

    try:
        template = load_template(f"{spec.id}.en.md")
    except FileNotFoundError:
        template = None
    citations = load_citations()
    render_report = (
        (lambda result: render(template, result, citations=citations)) if template else None
    )
    render_document = (
        _document_renderer(
            spec=spec,
            layout=layout,
            manifest=manifest,
            citations=citations,
            accession=genome.accession or genome.name,
        )
        if args.docx
        else None
    )

    try:
        outcome = execute_step(
            module_id=spec.id,
            module_version=spec.version,
            run_fn=spec.run,
            genomes=genomes,
            arity=spec.arity,
            params=params,
            ctx=ctx,
            cache=Cache(layout.cache_dir),
            force=args.force,
            render_report=render_report,
            render_document=render_document,
        )
    except IRGateError as error:
        return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_USAGE)
    except ExportDependencyError as error:
        return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_TOOL)
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

    entries, problems = entries_from_summaries(
        (summary, lineages.get(str(summary.get("taxid", "")), "")) for summary in summaries
    )
    for problem in problems:
        print(f"skipped: {problem}", file=sys.stderr)
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
    """Pin chosen genomes; accessions a panel index cannot hold are resolved live from NCBI.

    A RefSeq-only index can never contain a GenBank submission with no RefSeq twin, nor a
    mitogenome — yet those are exactly the genomes a manuscript may need (measured 2026-10-03:
    neither *Oldenlandia diffusa* nor *Coleus amboinicus* has a RefSeq plastome).
    """
    if args.groups:
        try:
            wanted = [a for group in load_groups(args.groups) for a in group.accessions]
        except (OSError, ValueError) as error:
            print(f"cannot read {args.groups}: {error}", file=sys.stderr)
            return EXIT_ERROR
    else:
        wanted = [a.strip() for a in args.accessions.split(",") if a.strip()]
    duplicates = sorted({a for a in wanted if wanted.count(a) > 1})
    if duplicates:
        print(f"panel pin: duplicated accession(s): {', '.join(duplicates)}", file=sys.stderr)
        return EXIT_USAGE
    if not wanted:
        print("panel pin: no accession to pin", file=sys.stderr)
        return EXIT_USAGE

    by_accession: dict[str, PanelEntry] = {}
    retrieved_at: str | None = None
    if args.index:
        try:
            snapshot = load_index(args.index)
        except (OSError, KeyError) as error:
            print(f"cannot read {args.index}: {error}", file=sys.stderr)
            return EXIT_ERROR
        by_accession = {entry.accession: entry for entry in snapshot.entries}
        retrieved_at = snapshot.retrieved_at

    client = NcbiClient()
    entries: list[PanelEntry] = []
    failed: list[str] = []
    for accession in wanted:
        entry = by_accession.get(accession)
        if entry is None:
            try:
                summary, lineage = client.resolve(accession)
            except (NcbiError, OSError) as error:
                print(f"{accession}: {error}", file=sys.stderr)
                failed.append(accession)
                continue
            resolved, problems = entries_from_summaries([(summary, lineage)])
            for problem in problems:
                print(f"{accession}: {problem}", file=sys.stderr)
            if not resolved:
                failed.append(accession)
                continue
            entry = resolved[0]
            print(f"{accession}: not in the index — resolved live ({entry.organism})")
        entries.append(entry)

    if not entries:
        print("panel pin: nothing resolved, nothing pinned", file=sys.stderr)
        return EXIT_ERROR
    try:
        pinned = pin(entries, client.fetch_genbank, args.out, retrieved_at=retrieved_at)
    except (NcbiError, OSError) as error:
        print(f"NCBI: {error}", file=sys.stderr)
        return EXIT_ERROR
    problems = list(failed)
    for record in pinned:
        digest = record.sha256[:12]
        if record.parse_error is not None:
            print(
                f"pinned {record.accession}: {record.bytes} bytes · sha256 {digest}… "
                f"(UNVERIFIED: {record.parse_error})"
            )
            problems.append(record.accession)
            continue
        genes = f"{record.cds} CDS + {record.trna} tRNA + {record.rrna} rRNA"
        ir = (
            f"IR {record.ir_length} bp {record.ir_source} @ {record.ir_identity * 100:.2f} %"
            if record.ir_error is None
            else f"IR: {record.ir_error}"
        )
        print(
            f"pinned {record.accession}: {record.bytes} bytes · sha256 {digest}… "
            f"({record.length} bp · {genes} · {record.topology} · {ir})"
        )
    if problems:
        print(f"panel pin: {len(problems)} problem(s): {', '.join(problems)}", file=sys.stderr)
        return EXIT_ERROR
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
