"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by returning an
exit code to the console-script wrapper). Library code raises exceptions instead.

Exit codes:
    0 — success
    1 — general error (bad input, unknown module, internal failure)
    2 — usage error, or a failed gate (e.g. the IR gate, the reference audit)
    3 — missing or failed external tool, or a missing optional dependency
    4 — a naming conflict stopped the step (gene_normalize, D5-09 §5)
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
from plastid.core import organelle, plastid_reference, reference_audit
from plastid.core.deps import ExportDependencyError
from plastid.core.genbank import read_fasta, read_genbank
from plastid.core.gene_names import NamingConflictError, load_recorded_reference, repo_root
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
from plastid.report.package import build_package
from plastid.report.qc import render_qc_for_run
from plastid.result.docx_export import DocxMeta, export_docx, merge_documents
from plastid.result.result import Result
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.figures import read_run_results, render_run_figures
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.runinfo import render_run_info
from plastid.run.runner import (
    DATA_LABEL_PARAM,
    DocumentRenderer,
    execute_step,
    export_files,
    write_payload,
)
from plastid.run.tools import ToolError

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_TOOL = 3
# 4 — a naming conflict stopped the step (D5-09 §5). Its own number, because it is not a usage
# error and not a missing tool: the run was well-formed and the data disagrees with the reference.
EXIT_CONFLICT = 4


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
        help="GenBank or FASTA file (one record); repeat the flag for a multi-genome module. "
        "A bare FASTA is a genome with no features yet — what `plastid annotate` takes",
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
        "--figures",
        action="store_true",
        help=(
            "also render the figures whose source tables this run produces, into the run's "
            "figures directory (needs the 'figures' extra)"
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

    package_parser = subparsers.add_parser(
        "package",
        help="regenerate the Supplementary package from a run directory (tables + figures + index)",
    )
    package_parser.add_argument("run", metavar="RUN_DIR", help="a run directory")

    qc_parser = subparsers.add_parser(
        "qc", help="write the technical QC report for a run (what ran, versions, inputs, audits)"
    )
    qc_parser.add_argument("run", metavar="RUN_DIR", help="a run directory")

    figures_parser = subparsers.add_parser(
        "figures",
        help="re-render a run's figures from the tables it stored (no analysis is re-run)",
    )
    figures_parser.add_argument("run", metavar="RUN_DIR", help="a run directory")
    figures_parser.add_argument(
        "--figure",
        action="append",
        default=None,
        metavar="ID",
        help="render only this figure (repeatable; default: every figure the run's tables feed)",
    )

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="annotate plastomes by transfer from the pinned plastid reference set",
    )
    annotate_parser.add_argument(
        "-i",
        "--input",
        action="append",
        default=[],
        help="FASTA or GenBank sample, one record — repeat the flag to annotate a batch, "
        "one run directory and one artifact per sample (not needed with --audit-only)",
    )
    annotate_parser.add_argument(
        "--reference",
        default=None,
        help="directory of GenBank records to transfer from (default: the pinned plastid set)",
    )
    annotate_parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="ACCESSION",
        help=(
            "leave this reference record out of the search (repeatable) — the mechanical form of "
            "D5-08's 'fix or drop the bad file': the run names what it did not search"
        ),
    )
    annotate_parser.add_argument(
        "--audit-only",
        action="store_true",
        help="audit the reference set and stop — nothing is annotated (P3B-03's flag)",
    )
    annotate_parser.add_argument(
        "--no-confirm",
        action="store_true",
        help="with --audit-only: skip the second reading; the report then says not measured",
    )
    annotate_parser.add_argument(
        "-o", "--out", default=".runs", help="output root (default: .runs)"
    )
    annotate_parser.add_argument(
        "--dataset", default=None, help="dataset label (default: each sample's accession)"
    )
    annotate_parser.add_argument(
        "--docx",
        action="store_true",
        help="also write report.docx in each run's module directory (needs the 'docx' extra)",
    )
    annotate_parser.add_argument("--force", action="store_true", help="ignore the cache and re-run")
    annotate_parser.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help="module parameter"
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

    reference = subparsers.add_parser(
        "reference",
        help="the plastid reference set: re-verify its pin, write its leave-one-out design",
    )
    reference_actions = reference.add_subparsers(dest="reference_action", metavar="ACTION")
    ref_verify = reference_actions.add_parser(
        "verify", help="re-hash every pinned record and name every file that changed"
    )
    ref_verify.add_argument(
        "--panel",
        default=".panel/plastid-reference",
        help="the pinned reference set (default: .panel/plastid-reference)",
    )
    ref_loo = reference_actions.add_parser(
        "loo", help="write the leave-one-out manifest: which records may annotate which target"
    )
    ref_loo.add_argument("--panel", default=".panel/plastid-reference")
    ref_loo.add_argument(
        "--family",
        required=True,
        help='the family the taxonomy rule was applied at, e.g. "Iridaceae"',
    )
    ref_loo.add_argument(
        "--rule", default="same family first (D-b)", help="the sentence the set was chosen by"
    )
    ref_audit = reference_actions.add_parser(
        "audit",
        help="audit the reference set: block what cannot be used, report what a person must know",
    )
    ref_audit.add_argument("--panel", default=".panel/plastid-reference")
    ref_audit.add_argument("-o", "--out", default=".runs", help="output root (default: .runs)")
    ref_audit.add_argument(
        "--exceptions",
        default=None,
        help=(
            "the pinned exceptions (default: <panel>/AUDIT_EXCEPTIONS.json); an exception is "
            "bound to the file's digest, so editing the file voids it"
        ),
    )
    ref_audit.add_argument(
        "--no-confirm",
        action="store_true",
        help="skip the second reading (the aligner call); the report then says not measured",
    )
    ref_audit.add_argument(
        "--confirm-exception",
        default=None,
        metavar="FILE:CHECK",
        help=(
            "pin a confirmed exception for one reported check and exit — the digest is written "
            "from the file as it is now, never typed"
        ),
    )
    ref_audit.add_argument("--confirmed-by", default=None, help="with --confirm-exception: who")
    ref_audit.add_argument(
        "--reason", default=None, help="with --confirm-exception: why the finding is accepted"
    )
    ref_loo.add_argument(
        "--out",
        default=None,
        help="where to write it (default: <panel>/REFERENCE_SET.json; point it at "
        "src/plastid/data/ to commit the offline copy)",
    )
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
        genomes = tuple(_read_input(path) for path in args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read the input: {error}", file=sys.stderr)
        return EXIT_ERROR
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

    try:
        reference, genomes, inputs = _resolve_reference(specs, params, genomes, list(args.input))
    except (OSError, ValueError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR
    for spec in specs:
        # API v0.2: arity ≥ 2 means *at least* that many genomes — an estimate over N accessions
        # wants all of them. arity 1 still takes exactly one. The runner enforces "at least"
        # too; this check is what turns a wrong count into a usage error before anything runs.
        # The check runs *after* the reference is resolved, because a module that declares one
        # counts it among its genomes: the naming reference is genomes[0], not a parameter the
        # module never sees.
        if len(genomes) < spec.arity or (spec.arity == 1 and len(genomes) != 1):
            wanted = "exactly one genome" if spec.arity == 1 else f"at least {spec.arity} genomes"
            extra = (
                " (plus the naming reference, which the run adds itself)"
                if spec.reference_param
                else ""
            )
            print(
                f"module {spec.id!r} takes {wanted}, got {len(genomes)} — "
                f"pass --input once per genome{extra}",
                file=sys.stderr,
            )
            return EXIT_USAGE

    return _execute_run(
        genomes=genomes,
        specs=specs,
        params=params,
        out=args.out,
        dataset=args.dataset,
        force=args.force,
        docx=args.docx,
        figures=args.figures,
        command=[
            "plastid",
            "run",
            *[part for spec in specs for part in ("--module", spec.id)],
            *[part for text in args.input for part in ("--input", text)],
            *(["--label", args.label] if args.label else []),
            *(["--figures"] if args.figures else []),
        ],
        inputs=inputs,
        references=reference,
    )


def _resolve_reference(
    specs: Sequence[registry.ModuleSpec],
    params: dict,
    genomes: tuple[Genome, ...],
    inputs: list[str],
) -> tuple[dict | None, tuple[Genome, ...], list[str]]:
    """Put a module's declared reference in front of the run's inputs — read, audited, recorded.

    A module that names a ``reference_param`` is working against a *privileged record*, not
    against one of the run's equal inputs (D5-09 §7: the naming reference, chosen by the end user
    and fixed for the project). Three things follow, and all three happen here rather than inside
    the module — modules never touch the filesystem:

    * the path comes from the parameter, or from the project's own recorded standard
      (``data/naming_reference.json``) when the run does not override it;
    * the record is **audited before it is used**: a reference that cannot be transferred from
      would make every name in the run wrong in the same way, which is D5-08's gate applied at
      the one place it is load-bearing rather than decorative;
    * it becomes ``genomes[0]`` — the convention six modules already follow — and its digest goes
      into the parameter, so the cache key carries the bytes (P3B-05) and a run that used a
      different standard says so.
    """
    declaring = [spec for spec in specs if spec.reference_param]
    if not declaring:
        return None, genomes, inputs
    if len(declaring) > 1:
        raise ValueError(
            f"one run cannot serve two privileged references: "
            f"{', '.join(sorted(spec.id for spec in declaring))} each declare one"
        )
    spec = declaring[0]
    key = spec.reference_param or ""
    override = params.get(key)
    if override:
        path = Path(str(override))
        recorded = None
        if not path.is_file():
            raise ValueError(f"{key}={override} is not a file — the naming reference has to be one")
    else:
        recorded = load_recorded_reference()
        path = recorded.resolve(root=repo_root())

    digest = _sha256_file(path)
    audit = reference_audit.audit_files([path], confirm=None)
    if not audit.ok:
        raise ValueError(
            f"{path} cannot be the naming reference — every name in the run would inherit the "
            f"fault:\n" + "\n".join(finding.line for finding in audit.blocked)
        )
    if recorded is not None and digest != recorded.sha256:
        raise ValueError(
            f"{path} is not the bytes the recorded naming reference was recorded at "
            f"({digest[:12]}… against {recorded.sha256[:12]}…) — re-record it, or pass "
            f"--set {key}=<the file you mean>"
        )
    reference = {
        "path": str(path),
        "accession": path.stem if recorded is None else recorded.accession,
        "sha256": digest,
        "provisional": bool(recorded.provisional) if recorded is not None else True,
    }
    params[key] = str(path)
    params[f"{key}_sha256"] = digest
    return reference, (_read_input(str(path)), *genomes), [str(path), *inputs]


def _read_input(path: str) -> Genome:
    """One input file as a genome — GenBank or FASTA, decided by the file's own first line.

    A sample that has just come off an assembler is a bare sequence, and annotating it is exactly
    what ``plastid annotate`` is for, so ``plastid run`` has to accept one (P4a-04). The
    choice is made on content rather than on the suffix: ``.fa``, ``.fasta``, ``.fna`` and ``.fas``
    all occur in the wild, and a sequencing facility's file naming is not ours to require.

    An empty file is refused here rather than by whichever reader happens to be picked — the
    GenBank reader's own message ("no GenBank records found") would name the wrong format.
    """
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                return read_fasta(path) if line.lstrip().startswith(">") else read_genbank(path)
    raise ValueError(f"{path}: empty file")


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
    figures: bool = False,
    references: dict | None = None,
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
        dataset=label,
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
        except reference_audit.ReferenceAuditError as error:
            # 2, the gate's own number, next to the IR gate: the run was well-formed and it is the
            # reference set that cannot be used. The message is the audit's own `file — reason`
            # lines; the way out (drop the file, or re-check the set) is said here, once.
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                f"{error}\n`plastid annotate --audit-only` re-checks the set without annotating "
                f"anything; `--exclude <accession>` drops a record from the search, and the run "
                f"says which records it did not search",
                EXIT_USAGE,
            )
        except NamingConflictError as error:
            # 4, not 2: the run was well-formed and the data disagrees with the reference — a
            # number of its own keeps it apart from the IR gate and the audit (D5-09 §5, §9.7).
            return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_CONFLICT)
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
        module_entry = {
            "id": spec.id,
            "version": spec.version,
            "status": outcome.status,
            "key": outcome.key,
            "tools": [asdict(record) for record in ctx.tools],
            "files": [entry["path"] for entry in module_files],
        }
        if spec.reference_param and references is not None:
            # The identity of the privileged record this step worked against — path *and*
            # digest, because a path alone cannot say whether the reference changed (D5-09 §7.2).
            module_entry["reference"] = references
        manifest["modules"].append(module_entry)

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
    if figures:
        # Rendering reads the results the run has just written, not the ones in memory: a module
        # served from the cache has no result object left, and a figure must not depend on whether
        # the analysis ran or was restored (P5-01).
        try:
            outcome = render_run_figures(layout.run_dir)
        except ExportDependencyError as error:
            return _fail(manifest, layout, ctx, specs[-1].version, str(error), EXIT_TOOL)
        except (OSError, ValueError) as error:
            return _fail(manifest, layout, ctx, specs[-1].version, str(error), EXIT_ERROR)
        for path in outcome.files:
            manifest["files"].append(
                {
                    "path": str(path.relative_to(layout.run_dir)),
                    "module": "figures",
                    "sha256": _sha256_file(path),
                    "bytes": path.stat().st_size,
                }
            )
        manifest["figures"] = {
            "rendered": [path.name for path in outcome.files],
            "skipped": [
                {"figure": figure_id, "reason": reason} for figure_id, reason in outcome.skipped
            ],
        }
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


def _run_dir_or_error(path: str) -> Path | None:
    """A run directory, or None after saying what is missing."""
    run_dir = Path(path)
    if not (run_dir / "manifest.json").is_file():
        print(
            f"{run_dir} is not a run directory — no manifest.json in it (run `plastid run` first)",
            file=sys.stderr,
        )
        return None
    return run_dir


def _cmd_package(args: argparse.Namespace) -> int:
    """The Supplementary package, rebuilt from the run's own stored results."""
    run_dir = _run_dir_or_error(args.run)
    if run_dir is None:
        return EXIT_ERROR
    try:
        results = read_run_results(run_dir)
        figures = render_run_figures(run_dir)
        package = build_package(
            run_dir, results=results, figures=figures.files, skipped=figures.skipped
        )
    except ExportDependencyError as error:
        print(str(error), file=sys.stderr)
        return EXIT_TOOL
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot build the package: {error}", file=sys.stderr)
        return EXIT_ERROR
    for path in package.files:
        print(f"wrote {path.relative_to(run_dir)}")
    for figure_id, reason in figures.skipped:
        print(f"skipped figure {figure_id}: {reason}", file=sys.stderr)
    return EXIT_OK


def _cmd_qc(args: argparse.Namespace) -> int:
    """The technical QC report — manifest plus the queryable-number audit of every report."""
    run_dir = _run_dir_or_error(args.run)
    if run_dir is None:
        return EXIT_ERROR
    try:
        figures = render_run_figures(run_dir, formats=("pdf",))
    except ExportDependencyError:
        figures = None  # the QC report is text; it must not need a plotting library
    except (OSError, ValueError, KeyError):
        figures = None
    rendered = tuple(path.name for path in figures.files) if figures else ()
    skipped = figures.skipped if figures else ()
    try:
        _, target = render_qc_for_run(run_dir, figures=rendered, skipped=skipped)
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot write the QC report: {error}", file=sys.stderr)
        return EXIT_ERROR
    print(f"wrote {target.relative_to(run_dir)}")
    return EXIT_OK


def _cmd_figures(args: argparse.Namespace) -> int:
    """Re-render a run's figures from the tables it stored — no module runs again."""
    run_dir = _run_dir_or_error(args.run)
    if run_dir is None:
        return EXIT_ERROR
    manifest_path = run_dir / "manifest.json"
    try:
        outcome = render_run_figures(run_dir, figures=args.figure)
    except ExportDependencyError as error:
        print(str(error), file=sys.stderr)
        return EXIT_TOOL
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot render figures: {error}", file=sys.stderr)
        return EXIT_ERROR

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    known = {entry["path"] for entry in manifest.get("files", [])}
    for path in outcome.files:
        relative = path.relative_to(run_dir).as_posix()
        if relative not in known:  # a re-render replaces the file; it is not listed twice
            manifest.setdefault("files", []).append(
                {
                    "path": relative,
                    "module": "figures",
                    "sha256": _sha256_file(path),
                    "bytes": path.stat().st_size,
                }
            )
        print(f"wrote {relative}")
    for figure_id, reason in outcome.skipped:
        print(f"skipped {figure_id}: {reason}", file=sys.stderr)
    if not outcome.files:
        print("no figure could be rendered from this run's tables", file=sys.stderr)
        return EXIT_ERROR
    write_manifest(manifest_path, manifest)
    return EXIT_OK


MITO_MODULE = "mito_annotate"
PLASTID_MODULE = "plastid_annotate"


def _cmd_annotate(args: argparse.Namespace) -> int:
    """Annotate 1..N plastomes from the reference set — one run per sample (P3B-09, M3B.4).

    **The batch is a batch of runs, not one run with N samples.** Each sample gets the run
    directory, manifest, tool record and report that every other output in this toolkit has
    (PROJECT.md §2: output is split by run/dataset and never overwritten in place), so a sample's
    annotation is traceable on its own — and the command's own numbers are the reference set it
    searched, recorded in each manifest, plus its exit code.

    **A sample that fails does not kill the batch.** It is named on stderr, its own run directory
    records the failure and the reason, and the command exits non-zero so a script cannot mistake
    a partial batch for a whole one (D5-08: "a failed sample is marked and named without killing
    the run, but the exit code reflects that something failed"). The first failure's code is what
    the command returns: the number is about *why* the run was refused, and the samples that
    failed for the most specific reason are the ones that ran furthest.

    **One organelle, declared at input.** A plastid sample that declares ``/organelle`` otherwise
    is refused by name, as is a reference set of mito records, and so is a batch that mixes the
    two — the path declares for a bare sequence (D5-06), and where a record declares, the record
    is what is read. ``plastid mito annotate`` is the path that serves a mitogenome.
    """
    try:
        spec = registry.get(PLASTID_MODULE)
    except registry.UnknownModuleError as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR

    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
    if args.reference:
        params["reference"] = args.reference
    planted = params.get("exclude")
    exclusions = [planted] if isinstance(planted, str) else list(planted or ())
    if args.exclude:
        exclusions.extend(args.exclude)
    if exclusions:
        params["exclude"] = exclusions

    directory = plastid_reference.resolve_directory(params.get("reference"))
    if args.audit_only:
        return _audit_reference_directory(
            panel=directory,
            out=args.out,
            no_confirm=args.no_confirm,
            exceptions_path=directory / reference_audit.EXCEPTIONS_FILENAME,
            command=[
                "plastid",
                "annotate",
                "--audit-only",
                "--reference",
                str(directory),
                *(["--no-confirm"] if args.no_confirm else []),
            ],
        )

    if not args.input:
        print(
            "annotate: --input is required unless --audit-only is asked for — one sample per "
            "-i, repeated for a batch",
            file=sys.stderr,
        )
        return EXIT_USAGE

    # Every sample is read and checked *before* anything is annotated: a batch that mixed two
    # organelles would otherwise have written half its runs before saying so (D5-06 — one run,
    # one organelle, and a mixed input set stops by name).
    samples: list[tuple[str, Genome]] = []
    codes: list[int] = []
    for text in args.input:
        try:
            samples.append((text, _read_input(text)))
        except (OSError, ValueError) as error:
            print(f"{text} — {error}", file=sys.stderr)
            codes.append(EXIT_ERROR)
    if samples:
        try:
            organelle.check_set(
                [genome for _text, genome in samples],
                group=spec.group,
                path="plastid annotate",
            )
        except organelle.OrganelleError as error:
            print(str(error), file=sys.stderr)
            print("annotate: nothing was annotated", file=sys.stderr)
            return EXIT_USAGE

    for text, genome in samples:
        codes.append(
            _annotate_sample(text=text, genome=genome, spec=spec, params=params, args=args)
        )

    failed = [code for code in codes if code != EXIT_OK]
    if failed:
        print(
            f"annotate: {len(failed)} of {len(codes)} sample(s) failed — each one's run directory "
            f"names the reason",
            file=sys.stderr,
        )
        return failed[0]
    print(f"annotate: {len(codes)} sample(s) annotated")
    return EXIT_OK


def _annotate_sample(
    *,
    text: str,
    genome: Genome,
    spec: registry.ModuleSpec,
    params: dict,
    args: argparse.Namespace,
) -> int:
    """One sample of the batch, through the same body every other run uses (seam 3)."""
    return _execute_run(
        genomes=(genome,),
        specs=(spec,),
        params=dict(params),
        out=args.out,
        dataset=args.dataset or genome.accession or genome.name,
        force=args.force,
        docx=args.docx,
        command=[
            "plastid",
            "annotate",
            "-i",
            text,
            *(["--reference", args.reference] if args.reference else []),
            *[part for name in args.exclude for part in ("--exclude", name)],
            *(["--dataset", args.dataset] if args.dataset else []),
            *(["--docx"] if args.docx else []),
        ],
        inputs=[text],
    )


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


def _cmd_reference(args: argparse.Namespace) -> int:
    action = getattr(args, "reference_action", None)
    if action == "audit":
        return _reference_audit(args)
    if action == "verify":
        return _reference_verify(args)
    if action == "loo":
        return _reference_loo(args)
    print("reference: choose an action — verify | loo", file=sys.stderr)
    return EXIT_USAGE


def _reference_verify(args: argparse.Namespace) -> int:
    """Re-hash the pin. Every problem is ``file — reason``; a set that changed is not the set."""
    try:
        problems = plastid_reference.verify_pin(args.panel)
    except (OSError, ValueError, KeyError) as error:
        print(f"{args.panel}: {error} — is this a pinned reference set?", file=sys.stderr)
        return EXIT_ERROR
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        print(f"reference verify: {len(problems)} problem(s) in {args.panel}", file=sys.stderr)
        return EXIT_ERROR
    print(f"reference verify: {args.panel} — every record matches its pin")
    return EXIT_OK


def _reference_audit(args: argparse.Namespace) -> int:
    """Audit a reference set before anything is annotated from it (P3B-03, D5-08).

    A run of its own — the same run directory, manifest and tool record as every other output,
    because the audit *is* something that ran, on inputs that can be named by digest, and because
    the confirmation layer calls an external aligner whose argv and version have to land somewhere.
    Blocking findings stop here with exit 2 and nothing is annotated; reporting findings are named
    and the run proceeds, because the classes they cover are real published data (§8.34).
    """
    panel = Path(args.panel)
    exceptions_path = (
        Path(args.exceptions) if args.exceptions else panel / reference_audit.EXCEPTIONS_FILENAME
    )

    if args.confirm_exception:
        return _reference_confirm_exception(args, panel, exceptions_path)

    return _audit_reference_directory(
        panel=panel,
        out=args.out,
        no_confirm=args.no_confirm,
        exceptions_path=exceptions_path,
        command=[
            "plastid",
            "reference",
            "audit",
            "--panel",
            str(panel),
            *(["--no-confirm"] if args.no_confirm else []),
        ],
    )


def _audit_reference_directory(
    *,
    panel: Path,
    out: str,
    no_confirm: bool,
    exceptions_path: Path | None,
    command: Sequence[str],
) -> int:
    """The audit as a run — its own directory, manifest and tool record (P3B-03, P3B-09).

    Shared by ``plastid reference audit`` and by ``plastid annotate --audit-only``: the flag the
    design gives the user to fix and re-check a set **without annotating anything** (D5-08) is the
    same audit, so it is the same code. Two entry points, one gate.
    """
    from plastid.modules import reference_audit as audit_module

    paths = sorted(panel.glob("*.gb"))
    if not paths:
        print(f"{panel}: no *.gb record to audit", file=sys.stderr)
        return EXIT_ERROR

    try:
        pins = reference_audit.read_exceptions(exceptions_path)
    except (OSError, ValueError, KeyError) as error:
        print(f"{exceptions_path}: {error}", file=sys.stderr)
        return EXIT_ERROR

    digests = "".join(f"{_sha256_file(path)}\x00" for path in paths)
    layout = RunLayout.create(
        out, panel.name, fingerprint=hashlib.sha256(digests.encode()).hexdigest()
    )
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id=audit_module.MODULE_ID, module_order=1)
    manifest = build_manifest(
        run_id=layout.run_id,
        command=list(command),
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params={},
        inputs=[
            {"path": str(path), "sha256": _sha256_file(path), "bytes": path.stat().st_size}
            for path in paths
        ],
        coordinate_system={
            "name": "genbank",
            "rotation_offset": 0,
            "note": "region tables use the normalized frame (LSC → IRa → SSC → IRb)",
        },
        dataset=layout.dataset,
    )

    confirmer, note = _audit_aligner(no_confirm)
    try:
        report = reference_audit.audit_files(
            paths, ctx=ctx, exceptions=pins, confirm=confirmer, confirm_note=note
        )
    except (OSError, ValueError) as error:
        return _fail(manifest, layout, ctx, audit_module.MODULE_VERSION, str(error), EXIT_ERROR)

    for finding in report.findings:
        print(finding.line, file=sys.stderr if finding.blocks else sys.stdout)
    if report.ok:
        print(
            f"reference audit: {len(report.records)} record(s) · 0 blocked · "
            f"{len(report.reported)} reported — {layout.run_dir}"
        )
    else:
        print(
            f"reference audit: {len(report.records)} record(s) · {len(report.blocked)} blocked — "
            f"nothing may be annotated from those files",
            file=sys.stderr,
        )

    files, figures = _write_audit_outputs(report, ctx, layout, panel)
    manifest["files"] = [
        {
            "path": str(path.relative_to(layout.run_dir)),
            "module": audit_module.MODULE_ID,
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
        for path in files
    ]
    manifest["figures"] = figures
    manifest["modules"].append(
        {
            "id": audit_module.MODULE_ID,
            "version": audit_module.MODULE_VERSION,
            "status": "failed" if not report.ok else "ok",
            "error": "\n".join(finding.line for finding in report.blocked) or None,
            "findings": {
                "blocked": len(report.blocked),
                "reported": len(report.reported),
                "confirmed": len(report.confirmed),
            },
            "tools": [asdict(record) for record in ctx.tools],
            "files": [str(path.relative_to(layout.run_dir)) for path in files],
        }
    )
    manifest["status"] = "ok" if report.ok else "failed"
    manifest["finished_at"] = datetime.now().isoformat(timespec="seconds")
    write_manifest(layout.run_dir / "manifest.json", manifest)
    (layout.run_info_dir / "run_info.md").write_text(
        render_run_info(manifest, ctx.tools), encoding="utf-8"
    )
    return EXIT_OK if report.ok else EXIT_USAGE


def _audit_aligner(no_confirm: bool):
    """The second reading, if this machine has an aligner — and the sentence to say if it has not.

    A missing MAFFT is not a failure: the minimum install has no external tools at all, and the
    audit's own rule already reached a verdict. What is forbidden is reporting that verdict as
    *confirmed* when nothing confirmed it.
    """
    if no_confirm:
        return None, "not measured (confirmation not requested)"
    from plastid.external.mafft import MAFFT_VERSION, adjust_direction
    from plastid.run.tools import which

    if which("mafft") is None:
        return None, f"not measured (mafft {MAFFT_VERSION} not found)"
    return adjust_direction, None


def _reference_confirm_exception(
    args: argparse.Namespace, panel: Path, exceptions_path: Path
) -> int:
    """Pin one reported finding to a file's current digest — the person confirms it once."""
    target, separator, check = args.confirm_exception.partition(":")
    if not separator or not check:
        print("--confirm-exception expects FILE:CHECK", file=sys.stderr)
        return EXIT_USAGE
    if not args.confirmed_by or not args.reason:
        print("--confirm-exception needs --confirmed-by and --reason", file=sys.stderr)
        return EXIT_USAGE
    try:
        pin = reference_audit.write_exception(
            exceptions_path,
            file=panel / target,
            check=check,
            reason=args.reason,
            confirmed_by=args.confirmed_by,
            confirmed_on=datetime.now().date().isoformat(),
        )
    except (OSError, ValueError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return EXIT_ERROR
    print(
        f"reference audit: {target} — {check} confirmed by {pin.confirmed_by} "
        f"({pin.confirmed_on}); the pin is void if the file changes"
    )
    return EXIT_OK


def _write_audit_outputs(
    report: reference_audit.AuditReport,
    ctx: RunContext,
    layout: RunLayout,
    panel: Path,
) -> tuple[tuple[Path, ...], dict]:
    """The audit as a result: the tables, the report a person reads, and the dot plot's picture.

    The picture is drawn whenever the rule found an inverted record — **not** whenever an aligner
    was available. It is a picture of the sequences, which the audit read itself, and a machine
    without MAFFT must not lose the one view that lets a person confirm a flip-flop. A machine
    without matplotlib loses only the picture: the matches stay in ``region_matches``, the audit's
    verdict is unchanged, and the manifest says the figure was skipped and why.
    """
    from plastid.modules import reference_audit as audit_module

    dotplot = None
    if audit_module.deviating_files(report):
        names = {path.name: path for path in sorted(panel.glob("*.gb"))}
        dotplot = audit_module.dotplot_for(
            report, reference_audit.region_sequences(names, "SSC"), "SSC"
        )
    result = audit_module.audit_result(report, dotplot=dotplot)
    template = load_template(f"{audit_module.MODULE_ID}.en.md")
    files = _write_result(result, template, ctx)
    figures: dict = {"rendered": [], "skipped": []}
    if dotplot is not None:
        try:
            outcome = render_run_figures(layout.run_dir, figures=["dotplot"])
        except ExportDependencyError as error:
            figures["skipped"] = [{"figure": "dotplot", "reason": str(error)}]
        else:
            files = files + tuple(outcome.files)
            figures["rendered"] = [path.name for path in outcome.files]
    return files, figures


def _write_result(result: Result, template: str, ctx: RunContext) -> tuple[Path, ...]:
    """Write a result's payload for a command that is not the module runner.

    ``plastid reference audit`` builds its own ``Result`` — the audit's input is a set of files,
    not a genome, so it never goes through :func:`execute_step` — and this is the one place that
    turns it into files: the runner's own export and write path, so a table written here and a
    table written by a module are the same file.
    """
    payload = export_files(result, lambda item: render(template, item))
    return write_payload(ctx, payload)


def _reference_loo(args: argparse.Namespace) -> int:
    """Write the leave-one-out manifest — the design, as data a run can read."""
    try:
        reference_set = plastid_reference.build_manifest(
            args.panel, family=args.family, rule=args.rule
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"reference loo: {error}", file=sys.stderr)
        return EXIT_ERROR
    out = Path(args.out) if args.out else Path(args.panel) / plastid_reference.MANIFEST_FILENAME
    if out.is_dir():
        out = out / plastid_reference.MANIFEST_FILENAME
    plastid_reference.write_manifest(reference_set, out)
    n = len(reference_set.records)
    print(f"reference loo: {n} record(s) · family {reference_set.family} → {out}")
    print(f"  every target is annotated by the same {n - 1} record(s), minus itself")
    return EXIT_OK


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

    if args.command == "figures":
        return _cmd_figures(args)

    if args.command == "package":
        return _cmd_package(args)

    if args.command == "qc":
        return _cmd_qc(args)

    if args.command == "annotate":
        return _cmd_annotate(args)

    if args.command == "mito":
        return _cmd_mito(args)

    if args.command == "panel":
        return _cmd_panel(args)

    if args.command == "reference":
        return _cmd_reference(args)

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
