"""`plastid assemble`: reads to a checked plastome, run here, prepared for elsewhere or imported."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Any

from plastid import __version__
from plastid.cli.common import (
    EXIT_ERROR,
    EXIT_OK,
    EXIT_TOOL,
    EXIT_USAGE,
    _cancelled,
    _fail,
    _file_entry,
    _finish,
    _parse_sets,
)
from plastid.cli.run import _read_input, _with_data_label
from plastid.core.assembly import AssemblyError
from plastid.core.digest import sha256_file
from plastid.core.model import Genome
from plastid.modules import assembly as assembly_step
from plastid.modules import registry
from plastid.report.citations import load_citations
from plastid.report.engine import load_template, render
from plastid.run.context import RunContext
from plastid.run.contract import describe_artifact
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.runner import (
    DATA_LABEL_PARAM,
    export_files,
    write_payload,
)
from plastid.run.tools import RunCancelled, ToolError

ASSEMBLY_STEPS = ((assembly_step.MODULE_ID, assembly_step.MODULE_VERSION),)


_IMPORT_LINE = (
    "plastid assemble import --prepared {prepared} --result <dir> --source-name <service> "
    "--retrieved-at <YYYY-MM-DD> -o {out}"
)


def _assembly_params(args: argparse.Namespace) -> dict | None:
    """The step's parameters, or None when a key was refused. A misspelling never takes effect."""
    try:
        params = _parse_sets(args.set or [])
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return None
    unknown = sorted(set(params) - set(assembly_step.ASSEMBLY_PARAMS))
    if unknown:
        accepted = ", ".join(assembly_step.ASSEMBLY_PARAMS)
        for key in unknown:
            print(f"--set {key}: assemble does not take it (accepted: {accepted})", file=sys.stderr)
        return None
    declared = {param.name: [param] for param in assembly_step.ASSEMBLY_PARAM_SPECS}
    problems = registry.value_problems(declared, params)
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return None
    params = registry.anchored_paths(declared, params)
    if args.label:
        params[DATA_LABEL_PARAM] = args.label
    return assembly_step.resolve_params(params)


def _read_pair(
    args: argparse.Namespace,
) -> assembly_step.ReadPair | assembly_step.LongReads | None:
    """The reads the command names, or None after saying what is missing."""
    if args.hifi:
        if args.reads_1 or args.reads_2:
            print("assemble: --hifi replaces -1 and -2: give one or the other", file=sys.stderr)
            return None
        if not args.sample:
            print("assemble: --sample required for this action", file=sys.stderr)
            return None
        if not Path(args.hifi).is_file():
            print(f"assemble: {args.hifi}: no such file", file=sys.stderr)
            return None
        return assembly_step.LongReads(args.sample, Path(args.hifi).resolve())
    missing = [
        flag
        for flag, value in (("-1", args.reads_1), ("-2", args.reads_2), ("--sample", args.sample))
        if not value
    ]
    if missing:
        print(f"assemble: {' '.join(missing)} required for this action", file=sys.stderr)
        return None
    # Absolute, because every tool runs in the run directory and a person types them from elsewhere.
    pair = assembly_step.ReadPair(
        args.sample, Path(args.reads_1).resolve(), Path(args.reads_2).resolve()
    )
    for path in (pair.reads_1, pair.reads_2):
        if not path.is_file():
            print(f"assemble: {path}: no such file", file=sys.stderr)
            return None
    return pair


def _assembly_route(args: argparse.Namespace) -> str:
    """The way forward when this machine cannot run the step: the same commands, spelled out."""
    given = f"--hifi {args.hifi}" if args.hifi else f"-1 {args.reads_1} -2 {args.reads_2}"
    reads = f"{given} --sample {args.sample} -o {args.out}"
    return (
        "Nothing ran, and nothing was sent anywhere. To run the assembly on a service or a "
        "cluster of your choice:\n"
        f"  plastid assemble prepare {reads}\n"
        "  (run the command it prints there, on the reads it wrote)\n"
        f"  {_IMPORT_LINE.format(prepared='<prepared.json>', out=args.out)}"
    )


def _assembly_manifest(
    *,
    layout: RunLayout,
    command: list[str],
    params: Mapping[str, Any],
    inputs: Sequence[Path],
    reference_sets: Sequence[dict],
) -> dict:
    return build_manifest(
        run_id=layout.run_id,
        command=command,
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params={
            key: list(value) if isinstance(value, tuple) else value for key, value in params.items()
        },
        inputs=[_file_entry(path) for path in inputs],
        coordinate_system={
            "name": "assembler",
            "rotation_offset": 0,
            "note": "positions are on the assembled sequence as the assembler wrote it",
        },
        dataset=layout.dataset,
        reference_sets=list(reference_sets),
    )


def _assembly_layout(out: str, dataset: str, digests: Sequence[str]) -> RunLayout:
    fingerprint = hashlib.sha256("".join(digests).encode("utf-8")).hexdigest()
    # Absolute: the assembler is told where to write, and it runs inside the run directory.
    layout = RunLayout.create(Path(out).resolve(), dataset, fingerprint=fingerprint)
    layout.make_dirs()
    return layout


def _cmd_assemble(args: argparse.Namespace) -> int:
    if args.action == "run":
        return _assemble_run(args)
    if args.action == "import":
        return _assemble_import(args)
    return _assemble_prepare(args)


def _assemble_prepare(args: argparse.Namespace) -> int:
    """Quality-control and hash the reads, and print what a run elsewhere has to do."""
    params = _assembly_params(args)
    reads = _read_pair(args)
    if params is None or reads is None:
        return EXIT_USAGE
    layout = _assembly_layout(
        args.out,
        args.dataset or reads.sample,
        [*(sha256_file(file) for file in reads.files), "prepare"],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
        announce_tools=True,
    )
    manifest = _assembly_manifest(
        layout=layout,
        command=args.typed,
        params=params,
        inputs=reads.files,
        reference_sets=(),
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)
    try:
        payload = assembly_step.prepare(reads, params, ctx)
    except RunCancelled:
        return _cancelled(manifest, layout, ctx, ASSEMBLY_STEPS)
    except ToolError as error:
        return _fail(manifest, layout, ctx, assembly_step.MODULE_VERSION, str(error), EXIT_TOOL)
    except (AssemblyError, OSError) as error:
        return _fail(manifest, layout, ctx, assembly_step.MODULE_VERSION, str(error), EXIT_ERROR)
    prepared = ctx.out("prepared", "json")
    prepared.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    assembly_step.clean_work(ctx)
    manifest["status"] = "ok"
    manifest["modules"].append(
        {
            "id": assembly_step.MODULE_ID,
            "version": assembly_step.MODULE_VERSION,
            "status": "prepared",
            "tools": [asdict(record) for record in ctx.tools],
            "files": [str(prepared.relative_to(layout.run_dir))],
        }
    )
    code = _finish(manifest, layout, ctx, EXIT_OK)
    print(f"prepared: {prepared}")
    print("Reads to give the service, hashed in that file:")
    for entry in payload["input_reads"]:
        print(f"  {prepared.parent / entry['path']}")
    print("The assembler's settings, for the service's form or your cluster's command line:")
    print(f"  {' '.join(payload['tool']['command'])}")
    print("Then bring the output directory back:")
    print(f"  {_IMPORT_LINE.format(prepared=prepared, out=args.out)}")
    return code


def _comparison_sequence(args: argparse.Namespace) -> tuple[Genome, str, dict] | None:
    """The record ``--compare-to`` names, its digest and manifest entry; None when unreadable."""
    try:
        reference = _read_input(args.compare_to)
    except (OSError, ValueError) as error:
        print(f"{args.compare_to} — {error}", file=sys.stderr)
        return None
    digest = sha256_file(Path(args.compare_to))
    return (
        reference,
        digest,
        {
            "role": "comparison_reference",
            "path": args.compare_to,
            "accession": reference.accession or reference.name,
            "digest": digest,
        },
    )


def _assemble_import(args: argparse.Namespace) -> int:
    """An assembly made elsewhere: everything after the assembly is the code a local run uses."""
    params = _assembly_params(args)
    if params is None:
        return EXIT_USAGE
    comparison = _comparison_sequence(args) if args.compare_to else None
    if args.compare_to and comparison is None:
        return EXIT_USAGE
    reference, reference_sha256, entry = comparison or (None, None, None)
    try:
        assembled = assembly_step.assemble_imported(
            args.prepared,
            args.result,
            source_name=args.source_name,
            retrieved_at=args.retrieved_at,
            source_url=args.source_url,
        )
    except (AssemblyError, OSError, ValueError) as error:
        print(f"assemble: {error}", file=sys.stderr)
        return EXIT_USAGE
    layout = _assembly_layout(
        args.out,
        args.dataset or assembled.sample,
        [sha256_file(Path(args.prepared)), args.source_name],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
        announce_tools=True,
    )
    manifest = _assembly_manifest(
        layout=layout,
        command=args.typed,
        params=params,
        inputs=(Path(args.prepared), *assembled.mapping_reads),
        reference_sets=[entry] if entry else [],
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)
    return _assembly_outcome(assembled, params, ctx, manifest, layout, reference, reference_sha256)


def _assemble_run(args: argparse.Namespace) -> int:
    """Assemble on this machine, or stop with the reason and the route to a run elsewhere."""
    params = _assembly_params(args)
    if params is None:
        return EXIT_USAGE
    comparison = _comparison_sequence(args) if args.compare_to else None
    if args.compare_to and comparison is None:
        return EXIT_USAGE
    reference, reference_sha256, entry = comparison or (None, None, None)
    reference_sets = [entry] if entry else []

    reads = _read_pair(args)
    if reads is None:
        return EXIT_USAGE
    machine = assembly_step.measure_machine(Path(args.out))
    try:
        tools, database = assembly_step.preflight(
            reads, params, machine=machine, check_resources=not args.skip_resource_check
        )
    except assembly_step.BackendUnavailable as error:
        print(f"assemble: cannot run here: {error}", file=sys.stderr)
        print(_assembly_route(args), file=sys.stderr)
        return EXIT_TOOL
    layout = _assembly_layout(
        args.out,
        args.dataset or reads.sample,
        [*(sha256_file(file) for file in reads.files), str(sorted(params.items()))],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
        announce_tools=True,
    )
    manifest = _assembly_manifest(
        layout=layout,
        command=args.typed,
        params=params,
        inputs=reads.files,
        reference_sets=[database, *reference_sets],
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)
    try:
        assembled = assembly_step.assemble_local(
            reads, params, ctx, tools=tools, database=database, machine=machine
        )
    except RunCancelled:
        return _cancelled(manifest, layout, ctx, ASSEMBLY_STEPS)
    except (assembly_step.BackendUnavailable, ToolError) as error:
        message = f"assemble: {error}\n{_assembly_route(args)}"
        return _fail(manifest, layout, ctx, assembly_step.MODULE_VERSION, message, EXIT_TOOL)
    except (AssemblyError, OSError) as error:
        return _fail(manifest, layout, ctx, assembly_step.MODULE_VERSION, str(error), EXIT_ERROR)
    return _assembly_outcome(assembled, params, ctx, manifest, layout, reference, reference_sha256)


def _assembly_outcome(
    assembled: assembly_step.Assembled,
    params: Mapping[str, Any],
    ctx: RunContext,
    manifest: dict,
    layout: RunLayout,
    reference: Genome | None,
    reference_sha256: str | None,
) -> int:
    """What both backends share: the evaluation, the tables, and the verdict for the next step."""
    version = assembly_step.MODULE_VERSION
    try:
        evaluation = assembly_step.evaluate(
            assembled, params, ctx, reference=reference, reference_sha256=reference_sha256
        )
    except RunCancelled:
        return _cancelled(manifest, layout, ctx, ASSEMBLY_STEPS)
    except ToolError as error:
        return _fail(manifest, layout, ctx, version, str(error), EXIT_TOOL)
    except (AssemblyError, OSError) as error:
        return _fail(manifest, layout, ctx, version, str(error), EXIT_ERROR)
    result = assembly_step.assembly_result(
        assembled, evaluation, data_label=str(params.get(DATA_LABEL_PARAM, ""))
    )
    render_report = _with_data_label(
        partial(
            render,
            load_template(f"{assembly_step.MODULE_ID}.en.md"),
            citations=load_citations(),
        ),
        dict(params),
    )
    files = write_payload(ctx, export_files(result, render_report))
    assembly_step.clean_work(ctx, keep=params["keep_work"])
    manifest["files"] = [
        _file_entry(path, relative_to=layout.run_dir, module=assembly_step.MODULE_ID)
        for path in files
    ]
    usable = evaluation.usable
    if usable:
        manifest["artifacts"] = [
            describe_artifact(
                ctx.out("assembly", "fasta"),
                relative_to=layout.run_dir,
                sample_id=assembled.sample,
                topology="circular",
                status="ok",
                frame="assembler",
                format="fasta",
            )
        ]
    manifest["modules"].append(
        {
            "id": assembly_step.MODULE_ID,
            "version": version,
            "status": "ok" if usable else "incomplete",
            "backend": assembled.backend,
            "where": assembled.where,
            "assembly_status": assembled.status,
            "read_support": evaluation.support,
            "graph": evaluation.graph_call,
            "tools": [asdict(record) for record in ctx.tools],
            "files": [str(path.relative_to(layout.run_dir)) for path in files],
        }
    )
    if not usable:
        why = (
            f"{assembled.sample}: the assembler says {assembled.status_text!r}, read support is "
            f"{evaluation.support}, the graph is {evaluation.graph_call}: no artifact for the "
            "next step."
        )
        if assembled.word_limit_hit:
            why += (
                " It stopped recruiting reads at its word limit. More memory (--set memory_gb), or "
                "a run elsewhere through `plastid assemble prepare`, would recruit more."
            )
        manifest["modules"][-1]["error"] = why
        print(f"assemble: {why} — {layout.run_dir}", file=sys.stderr)
    manifest["status"] = "ok" if usable else "incomplete"
    return _finish(manifest, layout, ctx, EXIT_OK if usable else EXIT_ERROR)
