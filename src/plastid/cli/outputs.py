"""The commands that read a finished run: `figures`, `package`, `qc` and `ir-junction`."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from plastid import __version__
from plastid.cli.catalogue import _describe_param, _param_record
from plastid.cli.common import (
    EXIT_CANCELLED,
    EXIT_ERROR,
    EXIT_OK,
    EXIT_TOOL,
    EXIT_USAGE,
    _batch_provenance,
    _file_entry,
    _finish,
    _manifest_provenance,
    _nothing_finished,
    _parse_sets,
    _print_json,
    _print_table,
    _run_dir_or_error,
    _stored_results,
)
from plastid.core.deps import ExportDependencyError
from plastid.core.genbank import read_genbank
from plastid.core.junctions import JunctionError
from plastid.core.model import Genome
from plastid.modules import registry
from plastid.report.package import build_package
from plastid.report.qc import render_qc_for_run
from plastid.result.figure import FIGURE_FORMATS, FigureRecord
from plastid.result.figures import all_figures, all_r_figures
from plastid.result.result import Result
from plastid.run.context import RunContext
from plastid.run.figures import (
    FigureOutcome,
    read_run_results,
    render_run_figures,
    stale_figure_files,
)
from plastid.run.junction import junction_config, junction_rows, junction_table
from plastid.run.layout import FIGURES_DIRNAME, RunLayout
from plastid.run.manifest import ManifestError, build_manifest, write_manifest
from plastid.run.rfigures import RFigureOutcome, render_run_rfigures
from plastid.run.runner import (
    DATA_LABEL_PARAM,
)
from plastid.run.state import exit_state
from plastid.run.tools import RunCancelled, ToolError


def _cmd_package(args: argparse.Namespace) -> int:
    """The Supplementary package, rebuilt from the run's own stored results.

    A re-package draws figures from the stored tables and copies them, together with the draft
    report, into one folder; it never re-runs an analysis. A missing plotting extra or a missing
    R does not throw the tables away — the absence is written into ``index.md`` and the package
    still exists.
    """
    run_dir = _run_dir_or_error(args.run)
    if run_dir is None or _nothing_finished(run_dir):
        return EXIT_USAGE
    manifest_path = run_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        results = read_run_results(run_dir)
        provenance = _manifest_provenance(manifest)
        outcome = _render_output_figures(
            run_dir,
            sets=manifest.get("params", {}),
            results=results,
            allow_missing_extra=True,
            provenance=provenance,
        )
        reports = _package_reports(run_dir)
        package = build_package(
            run_dir,
            results=results,
            figure_records=outcome.records,
            missing_renderers=outcome.missing,
            reports=reports,
            label=str(manifest.get("params", {}).get(DATA_LABEL_PARAM, "")),
            provenance=provenance,
        )
    except ExportDependencyError as error:
        print(str(error), file=sys.stderr)
        return EXIT_TOOL
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot build the package: {error}", file=sys.stderr)
        return EXIT_ERROR
    for path in package.files:
        print(f"wrote {path.relative_to(run_dir)}")
    for figure_id, reason in outcome.skipped:
        print(f"skipped figure {figure_id}: {reason}", file=sys.stderr)
    for renderer, reason in outcome.missing:
        print(f"{renderer} figures unavailable: {reason}", file=sys.stderr)
    # The package just rendered both catalogues; its status is the run's current figure status,
    # so a later `plastid qc` reports the real reason instead of guessing. A figure the
    # refresh deleted must leave the run's file list, and a rewritten one must carry its new
    # digest — an old path or an old sha256 would pass a stale result off as current.
    manifest["figures"] = outcome.as_manifest()
    refreshed = {
        path.relative_to(run_dir).as_posix(): _file_entry(
            path, relative_to=run_dir, module="figures"
        )
        for path in outcome.files
    }
    manifest["files"] = [
        entry
        for entry in manifest.get("files", [])
        if entry.get("module") != "figures"
        or ((run_dir / entry["path"]).is_file() and entry["path"] not in refreshed)
    ]
    manifest["files"].extend(refreshed.values())
    manifest["supplementary"] = {
        "directory": package.directory.relative_to(run_dir).as_posix(),
        "files": [
            _file_entry(path, relative_to=package.directory)
            for path in package.files
            if path.name != "MANIFEST.json"
        ],
        "figures": [
            {
                "figure": record.figure_id,
                "renderer": record.renderer,
                "status": record.status,
                "reason": record.reason,
            }
            for record in package.figures
        ],
        "skipped_outputs": [
            {"output": output, "reason": reason} for output, reason in package.skipped_outputs
        ],
    }
    try:
        write_manifest(manifest_path, manifest)
    except ManifestError as error:
        # An old run's manifest is not this build's schema. The package is complete and its own
        # MANIFEST.json is authoritative; say the run index was not updated instead of failing.
        print(
            f"warning: the package is complete, but {manifest_path.name} was not updated: {error}",
            file=sys.stderr,
        )
    return EXIT_OK


def _package_reports(run_dir: Path) -> tuple[tuple[str, Path], ...]:
    """The draft report files the run already stored, as (name in the package, source)."""
    candidates = (
        ("run_report.docx", run_dir / "report" / "report.docx"),
        *(
            (f"{directory.name}.report.md", directory / "report.md")
            for directory in sorted(run_dir.glob("[0-9][0-9]_*"))
            if (directory / "report.md").is_file()
        ),
    )
    return tuple((name, path) for name, path in candidates if path.is_file())


def _figure_inventory(
    run_dir: Path, results: Mapping[str, Result], manifest: dict
) -> tuple[tuple[FigureRecord, ...], tuple[tuple[str, str], ...]]:
    """Both catalogues as they stand in the run directory, without drawing anything.

    QC is a text report and must not need a plotting library or R. It still lists every figure
    of both catalogues and says, for each, whether the run holds it or why it does not — the
    reason recorded by the last render when there is one.
    """
    target = run_dir / FIGURES_DIRNAME
    recorded = {
        entry["figure"]: entry
        for entry in (manifest.get("figures") or {}).get("skipped", [])
        if isinstance(entry, dict) and entry.get("figure")
    }
    missing_renderers = tuple(
        (str(entry.get("renderer", "")), str(entry.get("reason", "")))
        for entry in (manifest.get("figures") or {}).get("missing_renderers", [])
        if isinstance(entry, dict)
    )
    records: list[FigureRecord] = []
    for spec in all_figures():
        images = stale_figure_files(target, spec.id, spec.formats)
        caption = target / f"{spec.id}.caption.txt"
        if images:
            files = (*images, caption) if caption.is_file() else images
            provenance = (
                tuple(caption.read_text(encoding="utf-8").splitlines())
                if caption.is_file()
                else tuple(
                    f"{module}.{table} — scope: {results[module].tables[table].scope.label}"
                    for module, table in spec.requires
                    if module in results and table in results[module].tables
                )
            )
            records.append(FigureRecord(spec.id, "matplotlib", files, provenance=provenance))
            continue
        needed = spec.missing(results)
        reason = (
            f"needs {', '.join(needed)} — not in this run"
            if needed
            else recorded.get(spec.id, {}).get("reason")
            or "not rendered — run `plastid figures <run>`"
        )
        records.append(FigureRecord(spec.id, "matplotlib", reason=reason))
    for spec in all_r_figures():
        images = tuple(
            sorted(
                path
                for path in target.glob(f"{spec.id}__*")
                if path.is_file() and path.suffix.lower() in (".png", ".pdf")
            )
            if target.is_dir()
            else ()
        )
        caption = target / f"{spec.id}.caption.txt"
        if images:
            files = (*images, caption) if caption.is_file() else images
            lines = caption.read_text(encoding="utf-8").splitlines() if caption.is_file() else []
            records.append(FigureRecord(spec.id, "R", tuple(files), provenance=tuple(lines)))
            continue
        needed = spec.missing(results)
        reason = (
            f"needs {', '.join(needed)} — not in this run"
            if needed
            else recorded.get(spec.id, {}).get("reason")
            or "not rendered — run `plastid figures <run>`"
        )
        records.append(FigureRecord(spec.id, "R", reason=reason))
    return tuple(records), missing_renderers


def _cmd_qc(args: argparse.Namespace) -> int:
    """The technical QC report — manifest, figure catalogue, and the number audit of reports."""
    run_dir = _run_dir_or_error(args.run)
    if run_dir is None:
        return EXIT_USAGE
    try:
        manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
        results = _stored_results(run_dir)
        records, missing = _figure_inventory(run_dir, results, manifest)
        _, target = render_qc_for_run(run_dir, figure_records=records, missing_renderers=missing)
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot write the QC report: {error}", file=sys.stderr)
        return EXIT_ERROR
    print(f"wrote {target.relative_to(run_dir)}")
    return EXIT_OK


def _cmd_ir_junction(args: argparse.Namespace) -> int:
    """The junction figure: one run directory holding the picture, the table and the provenance.

    The named exception to "a figure reads table cells" (docs/things_we_do_not_do.md): the species
    order comes from a tree and the genes either side of a junction come from the record's own
    annotation, so this figure is drawn from the records themselves.
    """
    records: list[tuple[Path, Genome]] = []
    for item in args.input:
        path = Path(item)
        if not path.is_file():
            print(f"no such record: {path}", file=sys.stderr)
            return EXIT_USAGE
        try:
            records.append((path, read_genbank(path)))
        except (OSError, ValueError) as error:
            print(f"cannot read {path}: {error}", file=sys.stderr)
            return EXIT_USAGE

    try:
        order = _junction_order(records, args.tree, args.order_from)
    except (OSError, ValueError) as error:
        print(f"cannot read the tree: {error}", file=sys.stderr)
        return EXIT_USAGE

    params: dict[str, Any] = {
        "tag": args.tag,
        "width": args.width,
        "tree": args.tree,
        "order_from": args.order_from,
        "drop_pseudo": args.drop_pseudo,
    }
    if args.label:
        params[DATA_LABEL_PARAM] = args.label

    fingerprint = hashlib.sha256(
        (
            "".join(f"{path}:{len(genome)}\x00" for path, genome in records)
            + json.dumps(params, sort_keys=True)
        ).encode("ascii", "replace")
    ).hexdigest()
    layout = RunLayout.create(args.out, args.dataset, fingerprint=fingerprint)
    layout.make_dirs()
    ctx = RunContext(layout=layout, params=params, module_id="ir_junction", module_order=1)
    manifest_path = layout.run_dir / "manifest.json"
    manifest = build_manifest(
        run_id=layout.run_id,
        command=args.typed,
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params=params,
        inputs=[_file_entry(path) for path, _ in records],
        coordinate_system={
            "name": "genbank",
            "rotation_offset": 0,
            "note": "junction positions are the record's own; nothing is rotated into a frame",
        },
        dataset=args.dataset,
    )
    write_manifest(manifest_path, manifest)

    # The drawing is wired here rather than in `run/junction.py`: the `run` layer may not import
    # `external/`, and this is where every other non-module step that calls a tool is wired.
    from plastid.external.cpjsdraw import CPJSDRAW_VERSION, run_cpjsdraw, svg_to_png

    def closed(status: str, code: int, error: str = "") -> int:
        """End the run with its one step: which tools ran, and what stopped them if one did."""
        step: dict[str, Any] = {
            "id": "ir_junction",
            "version": CPJSDRAW_VERSION,
            "status": status,
            "tools": [asdict(record) for record in ctx.tools],
            "files": [entry["path"] for entry in manifest["files"]],
        }
        if error:
            step["error"] = error
        manifest["modules"].append(step)
        manifest["status"] = exit_state(code).status
        return _finish(manifest, layout, ctx, code)

    figures = layout.figures_dir
    figures.mkdir(parents=True, exist_ok=True)
    try:
        rows = junction_rows(records, order)
        config = figures / f"cfg_{args.tag}.cfg"
        config.write_text(junction_config(rows), encoding="utf-8")
        table = figures / f"{args.tag}.csv"
        table.write_text(junction_table(rows), encoding="utf-8")
        svg = run_cpjsdraw(
            ctx,
            config,
            figures / f"{args.tag}.svg",
            script=args.cpjsdraw,
            drop_pseudo=args.drop_pseudo,
        )
        png = svg_to_png(ctx, svg, figures / f"{args.tag}.png", width=args.width)
    except RunCancelled:
        print(f"cancelled — what had finished is in {layout.run_dir}", file=sys.stderr)
        return closed("cancelled", EXIT_CANCELLED)
    except ToolError as error:
        print(str(error), file=sys.stderr)
        return closed("failed", EXIT_TOOL, str(error))
    except JunctionError as error:
        print(f"cannot draw the junction figure: {error}", file=sys.stderr)
        return closed("failed", EXIT_USAGE, str(error))
    except (OSError, ValueError) as error:
        print(f"cannot draw the junction figure: {error}", file=sys.stderr)
        return closed("failed", EXIT_ERROR, str(error))

    for path in (config, table, svg, png):
        relative = path.relative_to(layout.run_dir).as_posix()
        manifest["files"].append(
            _file_entry(path, relative_to=layout.run_dir, module="ir_junction")
        )
        print(f"wrote {relative}")
    for row in rows:
        if row.warning:
            print(f"warning: {row.warning}", file=sys.stderr)
    return closed("ok", EXIT_OK)


def _junction_order(
    records: Sequence[tuple[Path, Genome]], tree: str | None, reference: str | None
) -> list[str] | None:
    """Species names ordered by distance from a reference on the user's own tree.

    The tree is the one the user named, and its tip labels are not this toolkit's species names,
    so a tip is matched to a record by the first two words of the organism name.
    A reference that is not on the tree is said out loud and the input order is kept: silently
    ordering by nothing would look like a tree was used.
    """
    if not tree:
        return None
    from Bio import Phylo

    parsed = Phylo.read(tree, "newick")
    terminals = parsed.get_terminals()
    wanted = _normalised(reference) if reference else None
    root = next(
        (tip for tip in terminals if wanted and _normalised(tip.name).startswith(wanted)),
        None,
    )
    if root is None:
        print(
            f"warning: {reference!r} is not on {tree} — keeping the input order",
            file=sys.stderr,
        )
        return None
    distances = {_normalised(tip.name): parsed.distance(root, tip) for tip in terminals}
    return sorted(
        (genome.name for _, genome in records),
        key=lambda name: _distance_for(name, distances),
    )


def _distance_for(name: str, distances: dict[str, float]) -> float:
    """A record's distance from the reference tip, or infinity when no tip matches it."""
    wanted = _normalised(name)
    for tip, distance in distances.items():
        if tip.startswith(wanted):
            return distance
    return float("inf")


def _normalised(text: str | None) -> str:
    """A name reduced to lowercase letters and digits, for matching a record against a tip."""
    return "".join(character for character in (text or "").lower() if character.isalnum())


def _cmd_figures(args: argparse.Namespace) -> int:
    """Re-render a run's figures from the tables it stored — no module runs again.

    Two catalogues, one pass. The matplotlib figures are this toolkit's own drawing
    (vector only); the R figures are the manuscript's scripts, drawn from the same tables — a
    second source for a figure that already exists, and the only source for the ones the toolkit
    has no drawing for. A failure in one does not stop the other: they are independent, and a run
    with no plotting extra installed can still have its R figures. Every figure the pass did not
    draw is printed with its reason, and recorded in the manifest or ``batch.json``.

    A batch has no single genome to draw the one-genome catalogue for; its sample figures are the
    R ones, read from the aggregate table.
    """
    run_dir = _run_dir_or_error(args.run, allow_batch=True)
    if run_dir is None:
        return EXIT_USAGE
    batch_path = run_dir / "batch.json"
    is_batch = batch_path.is_file()
    manifest_path = batch_path if is_batch else run_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    try:
        sets = _parse_sets(args.sets)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
    declared = {option.name: [option] for option in FIGURE_OPTIONS}
    problems = [
        f"--set {key}: no figure reads it (accepted: {', '.join(sorted(declared))})"
        for key in sorted(sets)
        if key not in declared
    ]
    problems += registry.value_problems(declared, sets)
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return EXIT_USAGE
    sets = registry.anchored_paths(declared, sets)

    records: list[FigureRecord] = []
    batch_results = None
    if is_batch:
        from plastid.run.batch import read_batch_results

        batch_results = read_batch_results(run_dir)
    if args.list:
        return _list_figures(run_dir, batch_results, as_json=args.json)
    if not is_batch and _nothing_finished(run_dir):
        return EXIT_USAGE
    if is_batch:
        chosen = None if args.figure is None else set(args.figure)
        for spec in all_figures():
            if chosen is not None and spec.id not in chosen:
                continue
            reason = (
                "one-genome figure — a batch's sample figures are the R ones, read from "
                "aggregate.csv"
            )
            records.append(FigureRecord(spec.id, "matplotlib", reason=reason))

    try:
        outcome = _render_output_figures(
            run_dir,
            figures=args.figure,
            sets=sets,
            results=batch_results,
            include_matplotlib=not is_batch,
            allow_missing_extra=True,
            provenance=_batch_provenance(manifest) if is_batch else _manifest_provenance(manifest),
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot render figures: {error}", file=sys.stderr)
        return EXIT_ERROR
    records.extend(outcome.records)

    if is_batch:
        payload = manifest
        payload.setdefault("figures", {})["files"] = [
            _file_entry(path, relative_to=run_dir) for path in outcome.files
        ]
        payload["figures"]["skipped"] = [
            {"figure": figure_id, "renderer": _renderer_of(records, figure_id), "reason": reason}
            for figure_id, reason in outcome.skipped
        ]
        payload["figures"]["missing_renderers"] = [
            {"renderer": renderer, "reason": reason} for renderer, reason in outcome.missing
        ]
        payload["figures"]["status"] = outcome.as_manifest()["status"]
        payload["figures"]["failed"] = outcome.as_manifest()["failed"]
        payload["figures"]["policy"] = (
            "a batch has no single genome; the one-genome matplotlib figures are not drawn, "
            "and the sample figures are the R ones read from aggregate.csv"
        )
        payload["figures"]["records"] = [
            {
                "figure": record.figure_id,
                "renderer": record.renderer,
                "status": record.status,
                "files": [path.name for path in record.files],
                "reason": record.reason,
            }
            for record in records
        ]
        if outcome.tools:
            payload["figures"]["tools"] = [asdict(record) for record in outcome.tools]
        manifest_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        for path in outcome.files:
            print(f"wrote {path.relative_to(run_dir).as_posix()}")
        _say_undrawn(outcome)
        return EXIT_OK if outcome.files and not outcome.failed else EXIT_TOOL

    # A refresh can delete a figure whose table disappeared, and it rewrites the ones it does
    # draw. Both cases must reach the run's file list: a deleted path and an old digest are stale
    # results pretending to be current.
    refreshed = {
        path.relative_to(run_dir).as_posix(): _file_entry(
            path, relative_to=run_dir, module="figures"
        )
        for path in outcome.files
    }
    manifest["files"] = [
        entry
        for entry in manifest.get("files", [])
        if entry.get("module") != "figures"
        or ((run_dir / entry["path"]).is_file() and entry["path"] not in refreshed)
    ]
    manifest["files"].extend(refreshed.values())
    for relative in refreshed:
        print(f"wrote {relative}")
    _say_undrawn(outcome)
    # The whole catalogue's status goes into the manifest even when nothing was drawn: a later
    # QC pass and a package read it back instead of guessing why a figure is absent.
    manifest["figures"] = {
        **outcome.as_manifest(),
        "records": [
            {
                "figure": record.figure_id,
                "renderer": record.renderer,
                "status": record.status,
                "files": [path.name for path in record.files],
                "reason": record.reason,
            }
            for record in records
        ],
    }
    try:
        write_manifest(manifest_path, manifest)
    except ManifestError as error:
        print(
            f"warning: the figures are written, but {manifest_path.name} was not updated: {error}",
            file=sys.stderr,
        )
    if outcome.failed:
        return EXIT_TOOL
    if not outcome.files:
        print("no figure could be rendered from this run's tables", file=sys.stderr)
        # Exit 3 is "a tool this run needs is not here" — the plotting extra, or R when the
        # figures asked for are the R ones. A missing *table* stays exit 1: a data problem.
        return EXIT_TOOL if outcome.missing else EXIT_ERROR
    return EXIT_OK


_R_SCRIPT_OPTIONS: dict[str, tuple[str, ...]] = {
    "tree": (
        "tag",
        "root",
        "outgroup",
        "drop-tips",
        "outgroup-label",
        "plot-type",
        "legend",
        "width",
        "height",
        "dpi",
        "title",
    ),
    "rscu": ("species-order", "exclude-species", "width", "height"),
    "amino_acid": ("species-order", "exclude-species", "width", "height"),
    "diversity": ("width", "height", "dpi"),
}


# Every key `plastid figures --set` can reach: a script's option, a file a script reads, or R.
FIGURE_OPTIONS: tuple[registry.ParamSpec, ...] = (
    registry.ParamSpec("rscript", "file", "path to Rscript when it is not on PATH"),
    registry.ParamSpec("tree", "file", "the Newick tree r_tree draws"),
    registry.ParamSpec("tag", "str", "r_tree: base name of its files", default="phylo"),
    registry.ParamSpec(
        "root",
        "str",
        "r_tree: midpoint, outgroup or none",
        default="midpoint",
        choices=("midpoint", "outgroup", "none"),
    ),
    registry.ParamSpec("outgroup", "str", "r_tree: the taxa that root the tree"),
    registry.ParamSpec("drop-tips", "str", "r_tree: the tips to leave out"),
    registry.ParamSpec("outgroup-label", "str", "r_tree: the label of the outgroup clade"),
    registry.ParamSpec("plot-type", "str", "r_tree: which plots to draw", default="both"),
    registry.ParamSpec("legend", "str", "r_tree: where the legend goes", default="bottom"),
    registry.ParamSpec("title", "str", "r_tree: the figure title"),
    registry.ParamSpec("species-order", "str", "RSCU and amino-acid figures: column order"),
    registry.ParamSpec("exclude-species", "str", "RSCU and amino-acid figures: columns to drop"),
    registry.ParamSpec("width", "float", "drawing width, in the script's own unit"),
    registry.ParamSpec("height", "float", "drawing height, in the script's own unit"),
    registry.ParamSpec("dpi", "int", "drawing resolution"),
)


def _list_figures(
    run_dir: Path, batch_results: Mapping[str, Result] | None, *, as_json: bool
) -> int:
    """Both catalogues as this run stands, without drawing anything."""
    manifest_path = run_dir / ("batch.json" if batch_results is not None else "manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = batch_results if batch_results is not None else _stored_results(run_dir)
    records, _missing = _figure_inventory(run_dir, results, manifest)
    if as_json:
        specs = {spec.id: spec for spec in (*all_figures(), *all_r_figures())}
        _print_json(
            {
                "figures": [
                    {
                        "id": record.figure_id,
                        "renderer": record.renderer,
                        "drawn": bool(record.files),
                        "can_draw": not specs[record.figure_id].missing(results),
                        "reason": record.reason,
                    }
                    for record in records
                ],
                "options": [_param_record(option) for option in FIGURE_OPTIONS],
            }
        )
        return EXIT_OK
    _print_table(
        ("FIGURE", "DRAWN WITH", "STATE"),
        [
            (record.figure_id, record.renderer, "drawn" if record.files else record.reason)
            for record in records
        ],
    )
    print("\ndrawing options (--set KEY=VALUE):")
    for option in FIGURE_OPTIONS:
        _describe_param(option)
    return EXIT_OK


@dataclass(frozen=True, slots=True)
class OutputFigures:
    """One figure pass over both catalogues, with every absence accounted for.

    ``run --figures``, ``figures`` and ``package`` all call :func:`_render_output_figures`, so
    the set of figures is the same whichever command the user reaches for. ``missing`` names a
    whole renderer that could not run — the plotting extra, or R — and ``records`` always holds
    one entry per catalogue figure, drawn or not.
    """

    files: tuple[Path, ...]
    skipped: tuple[tuple[str, str], ...]
    records: tuple[FigureRecord, ...]
    missing: tuple[tuple[str, str], ...]
    tools: tuple[Any, ...] = ()
    failed: tuple[tuple[str, str], ...] = ()  # the run could feed these and their renderer did not

    def as_manifest(self) -> dict[str, Any]:
        return {
            "policy": (
                "the toolkit's matplotlib figures and its R figures are collected in one pass; "
                "a missing renderer is a named skip, never a silent absence"
            ),
            "rendered": [path.name for path in self.files],
            "skipped": [
                {"figure": figure_id, "renderer": renderer, "reason": reason}
                for figure_id, reason in self.skipped
                for renderer in [_renderer_of(self.records, figure_id)]
            ],
            "failed": [
                {
                    "figure": figure_id,
                    "renderer": _renderer_of(self.records, figure_id),
                    "reason": why,
                }
                for figure_id, why in self.failed
            ],
            "missing_renderers": [
                {"renderer": renderer, "reason": reason} for renderer, reason in self.missing
            ],
            "status": (
                "failed"
                if self.failed
                else "complete"
                if not self.skipped and not self.missing
                else ("partial" if self.files else "unavailable")
            ),
            "tools": [asdict(record) for record in self.tools],
        }


def _say_undrawn(outcome: OutputFigures) -> None:
    """Every figure a pass did not draw, by name, and a failure in words of its own."""
    failed = dict(outcome.failed)
    for figure_id, reason in outcome.skipped:
        word = "failed" if figure_id in failed else "skipped"
        print(f"{word} {figure_id}: {reason}", file=sys.stderr)
    for renderer, reason in outcome.missing:
        print(f"{renderer} figures unavailable: {reason}", file=sys.stderr)


def _renderer_of(records: Sequence[FigureRecord], figure_id: str) -> str:
    for record in records:
        if record.figure_id == figure_id:
            return record.renderer
    return "unknown"


def _sets_dict(values: Mapping[str, Any] | Sequence[str] | None) -> dict[str, Any]:
    """``--set`` pairs as a mapping, whether they arrive parsed or as ``KEY=VALUE`` strings."""
    if isinstance(values, Mapping):
        return dict(values)
    return _parse_sets(list(values or ()))


def _render_output_figures(
    run_dir: str | Path,
    *,
    figures: Sequence[str] | None = None,
    sets: Mapping[str, Any] | Sequence[str] | None = None,
    results: Mapping[str, Result] | None = None,
    formats: Sequence[str] = FIGURE_FORMATS,
    include_matplotlib: bool = True,
    include_r: bool = True,
    allow_missing_extra: bool = False,
    provenance: Sequence[str] = (),
) -> OutputFigures:
    """The single figure pass every output command uses.

    The matplotlib catalogue is this toolkit's own drawing; the R catalogue is the manuscript's,
    drawn from the same tables, and its captions carry the scope labels. A missing plotting extra
    propagates unless ``allow_missing_extra`` is set, because ``run --figures`` asking for figures
    and getting none is a tool failure (exit 3); ``package`` and ``qc`` still have their tables
    and record the absence instead. A missing R is always a recorded skip: the figures it would
    have drawn are optional, and the run's own figures are not.
    """
    root = Path(run_dir)
    catalogue = read_run_results(root) if results is None else dict(results)
    chosen = None if figures is None else list(figures)
    files: list[Path] = []
    skipped: list[tuple[str, str]] = []
    records: list[FigureRecord] = []
    missing: list[tuple[str, str]] = []
    tools: list[Any] = []
    failed: list[tuple[str, str]] = []

    if include_matplotlib:
        try:
            outcome: FigureOutcome = render_run_figures(
                root,
                figures=chosen,
                formats=formats,
                results=catalogue,
                provenance=provenance,
            )
        except ExportDependencyError as error:
            if not allow_missing_extra:
                raise
            missing.append(("matplotlib", str(error)))
            for spec in all_figures():
                if chosen is not None and spec.id not in chosen:
                    continue
                needed = spec.missing(catalogue)
                reason = f"needs {', '.join(needed)} — not in this run" if needed else str(error)
                skipped.append((spec.id, reason))
                records.append(FigureRecord(spec.id, "matplotlib", reason=reason))
        else:
            files.extend(outcome.files)
            skipped.extend(outcome.skipped)
            records.extend(outcome.records)

    if include_r:
        r_outcome, ctx = _render_r_figures(
            root, sets=sets, figures=chosen, results=catalogue, provenance=provenance
        )
        files.extend(r_outcome.files)
        skipped.extend(r_outcome.skipped)
        failed.extend(r_outcome.failed)
        records.extend(r_outcome.records)
        tools.extend(ctx.tools)
        if r_outcome.blocked is not None:
            missing.append(("R", r_outcome.blocked))

    return OutputFigures(
        files=tuple(files),
        skipped=tuple(skipped),
        records=tuple(records),
        missing=tuple(missing),
        tools=tuple(tools),
        failed=tuple(failed),
    )


def _render_r_figures(
    run_dir: Path,
    *,
    sets: Mapping[str, Any] | Sequence[str] | None = None,
    figures: Sequence[str] | None = None,
    results: Mapping[str, Result] | None = None,
    provenance: Sequence[str] = (),
) -> tuple[RFigureOutcome, RunContext]:
    """The R figure pass, with the tool call wired here — ``run/`` may not import ``external/``."""
    from plastid.external.rscript import missing_packages, rscript_binary, run_figure

    layout = RunLayout(
        out_root=run_dir.parent.parent, dataset=run_dir.parent.name, run_id=run_dir.name
    )
    try:
        params = _sets_dict(sets)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        params = {}
    # `ctx.out` is never used by this pass — it writes into the run's own figures directory — but
    # `ctx.tool` is, and that is what puts R's version and argv in the run's provenance.
    ctx = RunContext(layout=layout, params=params, module_id="r_figures", module_order=0)

    def preflight() -> str | None:
        """Asked once for the whole pass: R itself, then the packages the scripts need."""
        try:
            binary = rscript_binary(params.get("rscript"))
            missing = missing_packages(ctx, binary)
        except ToolError as error:
            return str(error)
        if missing:
            return f"R is missing {', '.join(missing)} — install them (docs/tools/r.md)"
        return None

    def runner(script: str, input_file, input_dir, out_dir: Path) -> Path:
        # The tree script takes its drawing choices as flags, and the run's `--set` values are
        # the only place a user gave them. Forwarding them here is what makes
        # `--set root/outgroup/drop-tips` reach R at all (K-06); the other scripts have no such
        # options, so they are not given them.
        extra: list[str] = []
        # Each script declares the options it reads; forwarding another script's flag would put a
        # value in an argv nothing reads. F3/F4 order and exclude species themselves (`--set
        # species-order`/`exclude-species`), F6 takes the drawing size, the tree takes its own
        # (K-06). Measured: these were reachable only by calling R by hand.
        for key in _R_SCRIPT_OPTIONS.get(script, ()):
            if key in params:
                extra += [f"--{key}", str(params[key])]
        return run_figure(
            ctx,
            script,
            input_file=input_file,
            input_dir=input_dir,
            out_dir=out_dir,
            binary=params.get("rscript"),
            # the preflight above already asked, once for the whole pass
            check_packages=False,
            extra=extra,
        )

    try:
        outcome = render_run_rfigures(
            run_dir,
            runner=runner,
            preflight=preflight,
            figures=figures,
            params=params,
            results=results,
            provenance=provenance,
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"cannot render the R figures: {error}", file=sys.stderr)
        return RFigureOutcome(files=(), skipped=()), ctx
    return outcome, ctx
