"""`plastid editing-reference`: rebuild or verify the packaged RNA-editing reference."""

from __future__ import annotations

import argparse
import json
import sys
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
    _finish,
    _sha256_text,
)
from plastid.core.digest import sha256_file
from plastid.core.genbank import read_genbank
from plastid.core.gene_names import (
    load_gene_aliases,
)
from plastid.run.context import RunContext
from plastid.run.layout import RunLayout
from plastid.run.manifest import build_manifest, write_manifest
from plastid.run.tools import RunCancelled, ToolError, tool_version, which


def _cmd_editing_reference(args: argparse.Namespace) -> int:
    if args.action == "build":
        return _editing_reference_build(args)
    return _editing_reference_verify(args)


def _editing_reference_build(args: argparse.Namespace) -> int:
    """Rebuild reference B′ as a run of its own — manifest, tool logs, the JSON written.

    The inputs are pinned twice over: every record by its SHA-256 in ``taxa.yaml``, the edit list
    and its sources by the files packaged beside it. MAFFT is the pinned 7.526 or the build stops,
    because the aligner's version is part of the digest and a reference built with another one
    would carry the same label and different columns.
    """
    from plastid.core import rna_editing
    from plastid.core.resources import data_file
    from plastid.external import mafft

    records_dir = Path(args.records)
    try:
        taxa, sources = rna_editing.load_build_inputs(data_file(rna_editing.REFERENCE_DIRECTORY))
        sites = rna_editing.read_edit_sites(data_file(rna_editing.EDIT_SITES_FILE))
    except (OSError, ValueError, KeyError) as error:
        print(f"editing-reference build: {error}", file=sys.stderr)
        return EXIT_ERROR
    paths = [records_dir / f"{taxon['accession']}.gb" for taxon in taxa]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        print(f"editing-reference build: missing record(s): {', '.join(missing)}", file=sys.stderr)
        return EXIT_USAGE
    pairs = list(zip((str(taxon["accession"]) for taxon in taxa), paths, strict=True))
    records = {accession: read_genbank(path) for accession, path in pairs}
    digests = {accession: sha256_file(path) for accession, path in pairs}

    binary = which("mafft")
    reported = tool_version("mafft", "binary", binary)[0] if binary else None
    if reported is None or mafft.MAFFT_VERSION not in reported:
        print(
            f"editing-reference build: MAFFT {mafft.MAFFT_VERSION} is required (found "
            f"{reported or 'none'}) — the aligner's version is part of B′'s digest",
            file=sys.stderr,
        )
        return EXIT_TOOL

    layout = RunLayout.create(args.out, "editing-reference", fingerprint=_sha256_text(str(digests)))
    layout.make_dirs()
    ctx = RunContext(layout=layout, params={}, module_id="editing_reference", module_order=1)
    manifest = build_manifest(
        run_id=layout.run_id,
        command=args.typed,
        package_version=__version__,
        started_at=datetime.now().isoformat(timespec="seconds"),
        params={},
        inputs=[
            {
                "path": str(path),
                "sha256": digests[str(taxon["accession"])],
                "bytes": path.stat().st_size,
            }
            for taxon, path in zip(taxa, paths, strict=True)
        ],
        coordinate_system={"name": "genbank", "rotation_offset": 0, "note": "records as pinned"},
        dataset="editing-reference",
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)
    try:
        payload = rna_editing.build_reference(
            records=records,
            record_sha256=digests,
            taxa=taxa,
            sites=sites,
            sources=sources,
            aliases=load_gene_aliases(),
            aligner=lambda rows: mafft.align_proteins(rows, ctx),
            aligner_record={
                "name": "mafft",
                "version": mafft.MAFFT_VERSION,
                "reported": reported,
                "argv": list(mafft.PROTEIN_ARGV),
            },
            label=rna_editing.REFERENCE_LABEL,
        )
    except RunCancelled:
        return _cancelled(manifest, layout, ctx, [("editing_reference", "1")])
    except rna_editing.EditingReferenceError as error:
        return _fail(manifest, layout, ctx, "1", str(error), EXIT_ERROR)
    except ToolError as error:
        return _fail(manifest, layout, ctx, "1", str(error), EXIT_TOOL)
    text = json.dumps(payload, ensure_ascii=False, indent=1, sort_keys=True) + "\n"
    out = Path(args.to)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    copy = layout.run_dir / "reference_b_prime.json"
    copy.write_text(text, encoding="utf-8")
    manifest["modules"].append(
        {
            "id": "editing_reference",
            "version": "1",
            "status": "ok",
            "tools": [asdict(record) for record in ctx.tools],
            "files": [str(copy.relative_to(layout.run_dir))],
        }
    )
    manifest["files"] = [_file_entry(copy, relative_to=layout.run_dir, module="editing_reference")]
    manifest["reference_sets"].append(
        {"role": "rna_editing_reference", "digest": payload["digest"], "path": str(out)}
    )
    manifest["status"] = "ok"
    print(f"reference B′ digest {payload['digest']} — {out}")
    return _finish(manifest, layout, ctx, EXIT_OK)


def _editing_reference_verify(args: argparse.Namespace) -> int:
    """B′ must describe itself (its digest) and, when the records are given, its pins."""
    from plastid.core import rna_editing

    try:
        reference = rna_editing.load_reference(args.reference_file)
    except (OSError, ValueError, KeyError) as error:
        print(f"editing-reference verify: {error}", file=sys.stderr)
        return EXIT_USAGE
    problems = []
    if args.records:
        for taxon in reference.payload["taxa"]:
            path = Path(args.records) / f"{taxon['accession']}.gb"
            if not path.is_file():
                problems.append(f"{path} — missing")
            elif sha256_file(path) != taxon["sha256"]:
                problems.append(f"{path} — not the pinned bytes ({taxon['sha256'][:12]}…)")
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return EXIT_ERROR
    print(f"reference B′ {reference.digest} — digest matches its contents")
    return EXIT_OK
