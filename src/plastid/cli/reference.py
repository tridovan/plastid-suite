"""`plastid reference` and `plastid panel`: pinning, verifying and auditing reference sets."""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections.abc import Sequence
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

from plastid import __version__
from plastid.cli.common import (
    EXIT_ERROR,
    EXIT_OK,
    EXIT_TOOL,
    EXIT_USAGE,
    _cancelled,
    _fail,
    _file_entry,
)
from plastid.core import plastid_reference, reference_audit
from plastid.core.deps import ExportDependencyError
from plastid.core.digest import sha256_file
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
from plastid.core.samples import load_groups
from plastid.report.engine import load_template, render
from plastid.result.result import Result
from plastid.run.context import RunContext
from plastid.run.figures import (
    render_run_figures,
)
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.runinfo import render_run_info
from plastid.run.runner import (
    export_files,
    write_payload,
)
from plastid.run.state import exit_state
from plastid.run.tools import RunCancelled


def _names_a_record(command: str, given: str | Path) -> bool:
    """Whether a set's flag named one record instead of its folder, having said so if it did.

    Every command that reads ``--reference-set`` answers this mistake in the same words, and as a
    wrong command line (exit 2): the record exists, so "no reference set" would be untrue.
    """
    misnamed = plastid_reference.misnamed_set(given)
    if misnamed is None:
        return False
    print(f"{command}: {misnamed}", file=sys.stderr)
    return True


def _cmd_reference(args: argparse.Namespace) -> int:
    if _names_a_record(f"reference {args.action}", args.reference_set):
        return EXIT_USAGE
    if args.action == "audit":
        return _reference_audit(args)
    if args.action == "verify":
        return _reference_verify(args)
    return _reference_loo(args)


def _reference_verify(args: argparse.Namespace) -> int:
    """Re-hash the pin. Every problem is ``file — reason``; a set that changed is not the set."""
    try:
        problems = plastid_reference.verify_pin(args.reference_set)
    except (OSError, ValueError, KeyError) as error:
        print(f"{args.reference_set}: {error} — is this a pinned reference set?", file=sys.stderr)
        return EXIT_USAGE
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        print(
            f"reference verify: {len(problems)} problem(s) in {args.reference_set}",
            file=sys.stderr,
        )
        return EXIT_ERROR
    print(f"reference verify: {args.reference_set} — every record matches its pin")
    return EXIT_OK


def _reference_audit(args: argparse.Namespace) -> int:
    """Audit a reference set before anything is annotated from it.

    A run of its own — the same run directory, manifest and tool record as every other output,
    because the audit *is* something that ran, on inputs that can be named by digest, and because
    the confirmation layer calls an external aligner whose argv and version have to land somewhere.
    Blocking findings stop here with exit 2 and nothing is annotated; reporting findings are named
    and the run proceeds, because the classes they cover are real published data.
    """
    panel = Path(args.reference_set)
    missing = plastid_reference.missing_set(panel)
    if missing is not None:
        print(f"reference audit: {missing}", file=sys.stderr)
        return EXIT_TOOL
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
        command=args.typed,
    )


def _audit_reference_directory(
    *,
    panel: Path,
    out: str,
    no_confirm: bool,
    exceptions_path: Path | None,
    command: Sequence[str],
) -> int:
    """The audit as a run — its own directory, manifest and tool record."""
    from plastid.modules import reference_audit as audit_module

    paths = sorted(panel.glob("*.gb"))
    if not paths:
        print(f"{panel}: no *.gb record to audit", file=sys.stderr)
        return EXIT_USAGE

    try:
        pins = reference_audit.read_exceptions(exceptions_path)
    except (OSError, ValueError, KeyError) as error:
        print(f"{exceptions_path}: {error}", file=sys.stderr)
        return EXIT_USAGE

    digests = "".join(f"{sha256_file(path)}\x00" for path in paths)
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
        inputs=[_file_entry(path) for path in paths],
        coordinate_system={
            "name": "genbank",
            "rotation_offset": 0,
            "note": "region tables use the normalized frame (LSC → IRa → SSC → IRb)",
        },
        dataset=layout.dataset,
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)

    confirmer, note = _audit_aligner(no_confirm)
    try:
        report = reference_audit.audit_files(
            paths, ctx=ctx, exceptions=pins, confirm=confirmer, confirm_note=note
        )
    except RunCancelled:
        steps = [(audit_module.MODULE_ID, audit_module.MODULE_VERSION)]
        return _cancelled(manifest, layout, ctx, steps)
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
        _file_entry(path, relative_to=layout.run_dir, module=audit_module.MODULE_ID)
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
    code = EXIT_OK if report.ok else EXIT_USAGE
    manifest["status"] = exit_state(code).status
    manifest["exit_code"] = code
    manifest["finished_at"] = datetime.now().isoformat(timespec="seconds")
    write_manifest(layout.run_dir / "manifest.json", manifest)
    (layout.run_info_dir / "run_info.md").write_text(
        render_run_info(manifest, ctx.tools), encoding="utf-8"
    )
    return code


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
            args.reference_set, family=args.family, rule=args.rule
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"reference leave-one-out: {error}", file=sys.stderr)
        return EXIT_USAGE
    out = (
        Path(args.to) if args.to else Path(args.reference_set) / plastid_reference.MANIFEST_FILENAME
    )
    if out.is_dir():
        out = out / plastid_reference.MANIFEST_FILENAME
    plastid_reference.write_manifest(reference_set, out)
    n = len(reference_set.records)
    print(f"reference leave-one-out: {n} record(s) · family {reference_set.family} → {out}")
    print(f"  every target is annotated by the same {n - 1} record(s), minus itself")
    return EXIT_OK


def _cmd_panel(args: argparse.Namespace) -> int:
    if args.action == "refresh":
        return _panel_refresh(args)
    if args.action == "propose":
        return _panel_propose(args)
    return _panel_pin(args)


def _panel_pin_files(args: argparse.Namespace) -> int:
    """Pin records produced outside this machine — the other half of the backend rule."""
    try:
        pinned = pin_files(args.file, args.to, note=args.note)
    except OSError as error:
        print(f"panel pin: {error}", file=sys.stderr)
        return EXIT_USAGE
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
    index = Path(args.index)
    index.parent.mkdir(parents=True, exist_ok=True)
    try:
        ids = client.search_refseq_plastomes(args.taxon, retmax=args.max_records)
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
    save_index(index, snapshot)
    print(f"panel: {len(entries)} RefSeq genome(s) for {args.taxon} → {index}")
    return EXIT_OK


def _panel_propose(args: argparse.Namespace) -> int:
    try:
        snapshot = load_index(args.index)
    except (OSError, KeyError) as error:
        print(
            f"cannot read {args.index}: {error} — run 'plastid panel refresh' first",
            file=sys.stderr,
        )
        return EXIT_USAGE
    try:
        proposal = propose(
            snapshot, args.target, mode=args.mode, level=args.level, per_tier=args.per_tier
        )
    except (ValueError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE

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

    A heavy tool runs where it fits — this machine when it can, a web
    service or a cluster when it cannot — so ``--file`` pins what such a tool produced: the same
    checksums, the same read-back, and the note saying where it came from.

    A RefSeq-only index can never contain a GenBank submission with no RefSeq twin, nor a
    mitogenome — yet those are exactly the genomes a manuscript may need (measured: neither
    *Oldenlandia diffusa* nor *Coleus amboinicus* has a RefSeq plastome).
    """
    if args.file:
        return _panel_pin_files(args)

    if args.groups:
        try:
            wanted = [a for group in load_groups(args.groups) for a in group.accessions]
        except (OSError, ValueError) as error:
            print(f"cannot read {args.groups}: {error}", file=sys.stderr)
            return EXIT_USAGE
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
            return EXIT_USAGE
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
        pinned = pin(entries, client.fetch_genbank, args.to, retrieved_at=retrieved_at)
    except (NcbiError, OSError) as error:
        print(f"NCBI: {error}", file=sys.stderr)
        return EXIT_ERROR
    except ValueError as error:
        print(f"panel pin: {error}", file=sys.stderr)
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
