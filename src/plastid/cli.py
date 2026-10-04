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
from collections.abc import Sequence
from dataclasses import asdict
from datetime import datetime
from functools import partial
from pathlib import Path

from plastid import __version__
from plastid.core.genbank import read_fasta, read_genbank
from plastid.core.ir import IRDetectionError, IRGateError, IRResult, detect_ir, enforce_ir_gate
from plastid.core.model import Genome
from plastid.core.normalize import normalize
from plastid.core.panel import (
    IndexSnapshot,
    PanelEntry,
    entries_from_summaries,
    load_index,
    pin,
    pin_files,
    propose,
    save_index,
)
from plastid.core.samples import group_of, load_groups
from plastid.modules import registry
from plastid.report.citations import Citation, format_citation, load_citations
from plastid.report.engine import CitationError, load_template, render
from plastid.result.docx_export import (
    DocxMeta,
    ExportDependencyError,
    export_docx,
    merge_documents,
)
from plastid.result.result import Result
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.runinfo import render_run_info
from plastid.run.runner import DATA_LABEL_PARAM, DocumentRenderer, execute_step
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
    run_parser.add_argument(
        "-m",
        "--module",
        required=True,
        action="append",
        help="module id (e.g. genome_stats); repeat the flag to run several in one run",
    )
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
        help=(
            "also write report.docx in each module directory, and a combined one in the run's "
            "report directory (needs the 'docx' extra)"
        ),
    )
    run_parser.add_argument(
        "--label",
        default=None,
        metavar="TEXT",
        help=(
            "data-provenance label written into every exported table's scope — a stand-in "
            "dataset never travels unlabelled (e.g. --label stand-in)"
        ),
    )
    run_parser.add_argument("--force", action="store_true", help="ignore the cache and re-run")
    run_parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="module parameter (repeatable)",
    )

    mito = subparsers.add_parser("mito", help="mitochondrial annotation (Phase 3A prototype)")
    mito_actions = mito.add_subparsers(dest="mito_action", metavar="ACTION")
    annotate = mito_actions.add_parser(
        "annotate",
        help="annotate a mitogenome by transferring a gene model from a reference set",
    )
    annotate.add_argument(
        "-i",
        "--input",
        required=True,
        help="FASTA file, one record — a bare assembled sequence is the usual input",
    )
    annotate.add_argument(
        "--reference",
        default=None,
        help=(
            "directory of GenBank references to transfer from (default: the packaged mito "
            "reference set built in P3A-02)"
        ),
    )
    annotate.add_argument("-o", "--out", default=".runs", help="output root (default: .runs)")
    annotate.add_argument("--dataset", default=None, help="dataset label (default: accession)")
    annotate.add_argument("--force", action="store_true", help="ignore the cache and re-run")
    annotate.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help="module parameter"
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
    wanted.add_argument(
        "--file",
        action="append",
        metavar="PATH",
        help="pin a record produced elsewhere (a web service, a cluster) — repeat for several; "
        "the answer of a heavy tool is a file, and a file that is not pinned cannot be reproduced",
    )
    pin_parser.add_argument(
        "--note",
        default=None,
        help="with --file: which service produced it and when, e.g. 'PMGA web, 2026-10-03'",
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


def _with_data_label(renderer, params: dict):
    """Put the run's data label at the top of every report — a stand-in never travels unlabelled.

    The tables carry it in their scope lines; the prose gets it here, so a report read on its own
    says where its numbers come from without the reader having to open the manifest.
    """
    label = params.get(DATA_LABEL_PARAM)

    def render(result: Result) -> str:
        text = renderer(result)
        if not label:
            return text
        title, separator, body = text.partition("\n")
        banner = f"> **Data:** {label}"
        return f"{title}\n\n{banner}\n{body}" if separator else f"{title}\n\n{banner}\n"

    return render


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
    specs: list[registry.ModuleSpec] = []
    for module_id in args.module:
        try:
            specs.append(registry.get(module_id))
        except registry.UnknownModuleError as error:
            print(str(error), file=sys.stderr)
            return EXIT_ERROR
    repeated = sorted({spec.id for spec in specs if [s.id for s in specs].count(spec.id) > 1})
    if repeated:
        print(
            f"run: module(s) named twice: {', '.join(repeated)} — one run holds each module once",
            file=sys.stderr,
        )
        return EXIT_USAGE

    try:
        genomes = tuple(read_genbank(text) for text in args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read the input: {error}", file=sys.stderr)
        return EXIT_ERROR
    for spec in specs:
        # API v0.2: arity ≥ 2 means *at least* that many genomes — an estimate over N accessions
        # wants all of them. arity 1 still takes exactly one. The runner enforces "at least"
        # too; this check is what turns a wrong count into a usage error before anything runs.
        if len(genomes) < spec.arity or (spec.arity == 1 and len(genomes) != 1):
            wanted = "exactly one genome" if spec.arity == 1 else f"at least {spec.arity} genomes"
            print(
                f"module {spec.id!r} takes {wanted}, got {len(genomes)} — "
                f"pass --input once per genome",
                file=sys.stderr,
            )
            return EXIT_USAGE
    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
    if args.label:
        params[DATA_LABEL_PARAM] = args.label

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

    return _execute_run(
        genomes=genomes,
        specs=specs,
        params=params,
        out=args.out,
        dataset=args.dataset,
        force=args.force,
        docx=args.docx,
        command=[
            "plastid",
            "run",
            *[part for spec in specs for part in ("--module", spec.id)],
            *[part for text in args.input for part in ("--input", text)],
            *(["--label", args.label] if args.label else []),
        ],
        inputs=list(args.input),
    )


def _execute_run(
    *,
    genomes: Sequence[Genome],
    specs: Sequence[registry.ModuleSpec],
    params: dict,
    out: str,
    dataset: str | None,
    force: bool,
    docx: bool,
    command: Sequence[str],
    inputs: Sequence[str],
) -> int:
    """One run's shared body: both ``plastid run`` and ``plastid mito annotate`` end up here.

    The point of sharing it is seam 3 — a second command must not become a second code path.
    Every output, every cache entry and every provenance record is written by the same lines,
    so an annotation produced by the mito command is traceable exactly like any other run.
    """

    overlay = json.dumps(params, sort_keys=True)
    sequences = "".join(item.sequence + "\x00" for item in genomes)  # ordered: order is meaning
    fingerprint = hashlib.sha256((sequences + overlay).encode("ascii", "replace")).hexdigest()
    label = dataset or genomes[0].accession or genomes[0].name
    layout = RunLayout.create(out, label, fingerprint=fingerprint)
    layout.make_dirs()

    ctx = RunContext(layout=layout, params=params, module_id=specs[0].id, module_order=1)
    manifest = build_manifest(
        run_id=layout.run_id,
        command=command,
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params=params,
        inputs=[
            {
                "path": str(text),
                "sha256": _sha256_file(Path(text)),
                "bytes": Path(text).stat().st_size,
            }
            for text in inputs
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

    if any(spec.requires_ir for spec in specs):
        gated: list[IRResult] = []
        for candidate in genomes:
            try:
                gated.append(enforce_ir_gate(detect_ir(candidate)))
            except (IRGateError, IRDetectionError) as error:
                # every genome in a comparison must pass: one diverged IR corrupts them all
                manifest["ir"] = {"gate_passed": False, "error": str(error)}
                return _fail(
                    manifest, layout, ctx, specs[0].version, f"IR gate: {error}", EXIT_USAGE
                )
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
            genomes[0], gated[0]
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

    citations = load_citations()
    gated_ir = gated[0] if any(spec.requires_ir for spec in specs) else None
    documents: list[Path] = []
    files: list[dict] = []
    for order, spec in enumerate(specs, start=1):
        ctx = RunContext(layout=layout, params=params, module_id=spec.id, module_order=order)
        ctx.ir = gated_ir
        try:
            template = load_template(f"{spec.id}.en.md")
        except FileNotFoundError:
            template = None
        render_report = (
            _with_data_label(partial(render, template, citations=citations), params)
            if template
            else None
        )
        render_document = (
            _document_renderer(
                spec=spec,
                layout=layout,
                manifest=manifest,
                citations=citations,
                accession=genomes[0].accession or genomes[0].name,
            )
            if docx
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
                force=force,
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

        module_files = []
        for path_text in outcome.files:
            path = Path(path_text)
            entry = {
                "path": str(path.relative_to(layout.run_dir)),
                "module": spec.id,
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
            module_files.append(entry)
            files.append(entry)
            if path.suffix == ".docx":
                documents.append(path)
        manifest["modules"].append(
            {
                "id": spec.id,
                "version": spec.version,
                "status": outcome.status,
                "key": outcome.key,
                "tools": [asdict(record) for record in ctx.tools],
                "files": [entry["path"] for entry in module_files],
            }
        )

    manifest["files"] = files
    if len(documents) > 1:
        # One run, one report: the modules' own documents are joined rather than re-rendered,
        # because a cache hit restores files and not results (see result/docx_export).
        try:
            combined = _combine_documents(layout, manifest, documents, label)
        except (OSError, ValueError) as error:
            return _fail(
                manifest,
                layout,
                ctx,
                specs[-1].version,
                f"cannot combine documents: {error}",
                EXIT_ERROR,
            )
        manifest["files"].append(
            {
                "path": str(combined.relative_to(layout.run_dir)),
                "module": "run",
                "sha256": _sha256_file(combined),
                "bytes": combined.stat().st_size,
            }
        )
    manifest["status"] = "ok"
    return _finish(manifest, layout, ctx, specs[-1].version, EXIT_OK)


def _combine_documents(
    layout: RunLayout, manifest: dict, documents: list[Path], dataset: str
) -> Path:
    """Join the per-module documents into the run's own report — one file per run."""
    blob = merge_documents(
        [path.read_bytes() for path in documents],
        title=f"{dataset} — analysis report",
        label=f"{dataset} · {len(documents)} modules",
    )
    path = layout.report_dir / "report.docx"
    path.write_bytes(blob)
    return path


MITO_MODULE = "mito_annotate"


def _cmd_mito(args: argparse.Namespace) -> int:
    action = getattr(args, "mito_action", None)
    if action == "annotate":
        return _cmd_mito_annotate(args)
    print("mito: choose an action — annotate", file=sys.stderr)
    return EXIT_USAGE


def _cmd_mito_annotate(args: argparse.Namespace) -> int:
    """Annotate an assembled mitogenome: FASTA in, GenBank out, the standard run report.

    This command exists because the phase's input is a **bare sequence** — an assembler's output
    has no features to read, and ``plastid run`` takes GenBank. Everything after the input is the
    same body ``plastid run`` uses (:func:`_execute_run`), so the annotation lands in a run
    directory with a manifest, a tool record and a report like every other output here: a second
    command must not become a second code path (seam 3).

    The GenBank it produces is written by :func:`plastid.core.genbank.to_genbank` and re-read by
    the module before it is returned, so the file can be handed straight to
    ``plastid run -m orf -m mito_genes``.
    """
    try:
        spec = registry.get(MITO_MODULE)
    except registry.UnknownModuleError as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR

    try:
        genome = read_fasta(args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read the input: {error}", file=sys.stderr)
        return EXIT_ERROR

    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
    if args.reference:
        params["reference"] = args.reference

    return _execute_run(
        genomes=(genome,),
        specs=(spec,),
        params=params,
        out=args.out,
        dataset=args.dataset,
        force=args.force,
        docx=False,
        command=[
            "plastid",
            "mito",
            "annotate",
            "--input",
            args.input,
            *(["--reference", args.reference] if args.reference else []),
        ],
        inputs=[args.input],
    )


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


def _panel_pin_files(args: argparse.Namespace) -> int:
    """Pin records produced outside this machine — the other half of the backend rule (P1-12)."""
    try:
        pinned = pin_files(args.file, args.out, note=args.note)
    except OSError as error:
        print(f"panel pin: {error}", file=sys.stderr)
        return EXIT_ERROR
    problems = []
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
        print(
            f"pinned {record.accession}: {record.bytes} bytes · sha256 {digest}… "
            f"({record.length} bp · {genes} · {record.topology})"
            + (f" · {record.note}" if record.note else "")
        )
    if problems:
        print(f"panel pin: {len(problems)} problem(s): {', '.join(problems)}", file=sys.stderr)
        return EXIT_ERROR
    return EXIT_OK


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

    A heavy tool runs where it fits (PROJECT.md §3, P1-12) — this machine when it can, a web
    service or a cluster when it cannot — so ``--file`` pins what such a tool produced: the same
    checksums, the same read-back, and the note saying where it came from.

    A RefSeq-only index can never contain a GenBank submission with no RefSeq twin, nor a
    mitogenome — yet those are exactly the genomes a manuscript may need (measured 2026-10-03:
    neither *Oldenlandia diffusa* nor *Coleus amboinicus* has a RefSeq plastome).
    """
    if args.file:
        return _panel_pin_files(args)

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

    if args.command == "mito":
        return _cmd_mito(args)

    if args.command == "panel":
        return _cmd_panel(args)

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
