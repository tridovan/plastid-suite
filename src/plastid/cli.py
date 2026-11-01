"""Command-line interface for PlastidToolSuite.

This module is the only place allowed to terminate the process (directly or by returning an
exit code to the console-script wrapper). Library code raises exceptions instead.

Exit codes and what each one means: ``plastid.run.state.EXIT_STATES``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import signal
import sys
import threading
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Any

from plastid import __version__
from plastid.core import organelle, plastid_reference, reference_audit
from plastid.core.assembly import AssemblyError
from plastid.core.codon import UnsupportedTranslationError
from plastid.core.deps import ExportDependencyError
from plastid.core.genbank import read_fasta, read_genbank
from plastid.core.gene_names import (
    NamingConflictError,
    RecordedReference,
    checkout_root,
    load_gene_aliases,
    load_recorded_reference,
)
from plastid.core.ir import (
    IR_SCORING_POLICY,
    IR_SCORING_POLICY_VERSION,
    IRDetectionError,
    IRGateError,
    IRResult,
    detect_ir,
    enforce_ir_gate,
)
from plastid.core.junctions import JunctionError
from plastid.core.model import Genome
from plastid.core.normalize import normalize
from plastid.core.panel import (
    DEFAULT_PROPOSAL_MODE,
    PROPOSAL_LEVELS,
    PROPOSAL_MODES,
    IndexSnapshot,
    PanelEntry,
    entries_from_summaries,
    load_index,
    pin,
    pin_files,
    propose,
    save_index,
)
from plastid.core.rna_editing import EditingReferenceError
from plastid.core.samples import group_of, load_groups
from plastid.modules import assembly as assembly_step
from plastid.modules import registry
from plastid.report.citations import Citation, format_citation, load_citations
from plastid.report.engine import CitationError, load_template, render
from plastid.report.package import build_package
from plastid.report.qc import render_qc_for_run
from plastid.result.docx_export import DocxMeta, export_docx, merge_documents
from plastid.result.figure import FIGURE_FORMATS, FigureRecord
from plastid.result.figures import all_figures, all_r_figures
from plastid.result.result import Result
from plastid.run import batch as batch_run
from plastid.run.cache import Cache
from plastid.run.context import RunContext
from plastid.run.contract import (
    ContractError,
    describe_artifact,
    describe_reference_set,
    sha256_file,
)
from plastid.run.figures import (
    FigureOutcome,
    read_run_results,
    render_run_figures,
    stale_figure_files,
)
from plastid.run.junction import junction_config, junction_rows, junction_table
from plastid.run.layout import FIGURES_DIRNAME, RunLayout
from plastid.run.manifest import ManifestError, build_manifest, write_manifest
from plastid.run.ncbi import NcbiClient, NcbiError
from plastid.run.rfigures import RFigureOutcome, render_run_rfigures
from plastid.run.runinfo import render_run_info
from plastid.run.runner import (
    DATA_LABEL_PARAM,
    DocumentRenderer,
    execute_step,
    export_files,
    write_payload,
)
from plastid.run.state import EXIT_STATES, exit_state, list_runs, read_state
from plastid.run.tools import RunCancelled, ToolError, tool_version, which

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

# The tools no module calls: which command does.
_TOOL_COMMANDS = {"R": "plastid figures", "cpjsdraw": "plastid ir-junction", "perl": "ssr"}

_STOP_SIGNALS = tuple(
    getattr(signal, name) for name in ("SIGINT", "SIGTERM", "SIGHUP") if hasattr(signal, name)
)


def _epilog() -> str:
    states = "\n".join(
        f"  {state.code:<4} {state.status:<11} {state.meaning}" for state in EXIT_STATES
    )
    return (
        "start here:\n"
        "  plastid modules                what can be run\n"
        "  plastid modules ssr            one module and its parameters\n"
        "  plastid tools                  which external tools this machine has\n"
        "  plastid run -m genome_stats -i my_plastome.gb\n"
        "  plastid runs                   what has been run, and how each run ended\n\n"
        f"exit codes:\n{states}\n"
    )


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plastid",
        description=(
            "PlastidToolSuite: plastid genome analysis, from reads or records to tables, figures "
            "and reports."
        ),
        epilog=_epilog(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    modules_parser = subparsers.add_parser(
        "modules", help="list the analysis modules, or one module with its parameters"
    )
    modules_parser.add_argument(
        "module", nargs="?", metavar="MODULE", help="the module to describe"
    )
    modules_parser.add_argument("--json", action="store_true", help=_JSON_HELP)

    tools_parser = subparsers.add_parser(
        "tools", help="list the external tools, and whether this machine has each one"
    )
    tools_parser.add_argument("--json", action="store_true", help=_JSON_HELP)

    runs_parser = subparsers.add_parser("runs", help="list the runs in a directory, newest first")
    runs_parser.add_argument("root", nargs="?", metavar="DIR", help="default: .runs")
    runs_parser.add_argument("--json", action="store_true", help=_JSON_HELP)

    run_parser = subparsers.add_parser(
        "run", help="run one or more analysis modules on one genome, or on several to compare"
    )
    run_parser.add_argument(
        "-m",
        "--module",
        required=True,
        action="append",
        help="module to run, repeatable; `plastid modules` lists them",
    )
    run_parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        help="GenBank or FASTA file holding one record; repeat it for a module that compares "
        "genomes",
    )
    _add_run_directory(run_parser, "the accession")
    run_parser.add_argument(
        "--groups",
        default=None,
        metavar="FILE",
        help="YAML file of sample groups, {group: [accessions]}, labelling every input",
    )
    run_parser.add_argument(
        "--docx",
        action="store_true",
        help="also write a Word report for each module and one for the run (needs the 'docx' "
        "extra)",
    )
    run_parser.add_argument(
        "--figures",
        action="store_true",
        help="also draw every figure this run's tables can feed, into the run's figures directory "
        "(needs the 'figures' extra; an R figure is skipped by name when R is absent)",
    )
    run_parser.add_argument("--no-cache", action="store_true", help=_NO_CACHE_HELP)
    run_parser.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help=_SET_HELP
    )

    batch_parser = subparsers.add_parser(
        "batch",
        help="run single-genome modules on many genomes, one run each, with one table for all",
    )
    batch_parser.add_argument(
        "-m",
        "--module",
        required=True,
        action="append",
        help="single-genome module to run on every sample, repeatable",
    )
    batch_parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        help="GenBank or FASTA file holding one record; repeat it, one run per sample",
    )
    _add_run_directory(batch_parser, "batch")
    batch_parser.add_argument(
        "--docx",
        action="store_true",
        help="also write a Word report in each sample's run (needs the 'docx' extra)",
    )
    batch_parser.add_argument(
        "--figures",
        action="store_true",
        help="after the combined table, draw the multi-sample R figures from it (needs R)",
    )
    batch_parser.add_argument("--no-cache", action="store_true", help=_NO_CACHE_HELP)
    batch_parser.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help=_SET_HELP
    )

    annotate_parser = subparsers.add_parser(
        "annotate", help="annotate plastomes by transfer from a set of reference records"
    )
    annotate_parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        help="FASTA or GenBank file holding one record; repeat it, one run per sample",
    )
    annotate_parser.add_argument(
        "--reference-set",
        default=None,
        metavar="DIR",
        help="directory of GenBank records to transfer from (default: .panel/plastid-reference)",
    )
    annotate_parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="ACCESSION",
        help="leave this record out of the search, repeatable; the run names what it left out",
    )
    _add_run_directory(annotate_parser, "each sample's accession")
    annotate_parser.add_argument(
        "--docx",
        action="store_true",
        help="also write a Word report in each run (needs the 'docx' extra)",
    )
    annotate_parser.add_argument("--no-cache", action="store_true", help=_NO_CACHE_HELP)
    annotate_parser.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help=_SET_HELP
    )

    assemble = subparsers.add_parser(
        "assemble", help="assemble a plastome from short or HiFi reads, here or elsewhere"
    )
    assemble_actions = assemble.add_subparsers(dest="action", metavar="ACTION")
    assemble_run = assemble_actions.add_parser(
        "run", help="assemble on this machine, and check the result against the reads"
    )
    assemble_prepare = assemble_actions.add_parser(
        "prepare",
        help="quality-control and hash the reads for a run elsewhere; nothing leaves this machine",
    )
    assemble_import = assemble_actions.add_parser(
        "import",
        help="read an assembly made elsewhere, and check it against the reads the same way",
    )
    for action in (assemble_run, assemble_prepare):
        action.add_argument("-1", "--reads-1", default=None, help="first reads (FASTQ, gzipped ok)")
        action.add_argument("-2", "--reads-2", default=None, help="second reads")
        action.add_argument(
            "--hifi",
            default=None,
            metavar="FASTQ",
            help="PacBio HiFi reads (gzipped ok), in place of -1 and -2; assembled with Oatk",
        )
        action.add_argument("--sample", default=None, help="sample name for the run and outputs")
        _add_run_directory(action, "the sample")
    _add_run_directory(assemble_import, "the sample")
    for action in (assemble_run, assemble_prepare, assemble_import):
        action.add_argument(
            "--set",
            action="append",
            default=[],
            metavar="KEY=VALUE",
            help="step parameter, repeatable; `plastid modules assembly` lists them",
        )
    assemble_run.add_argument(
        "--skip-resource-check",
        action="store_true",
        help="run even when this machine looks too small for the input. Nothing is ever sent to "
        "a service by this command",
    )
    for action in (assemble_run, assemble_import):
        action.add_argument(
            "--compare-to",
            default=None,
            metavar="FILE",
            help="GenBank or FASTA sequence to compare every assembled path with",
        )
    assemble_import.add_argument(
        "--prepared",
        required=True,
        metavar="FILE",
        help="the prepared.json of `plastid assemble prepare`",
    )
    assemble_import.add_argument(
        "--result",
        required=True,
        metavar="DIR",
        help="the directory holding the assembler's output",
    )
    assemble_import.add_argument(
        "--source-name", required=True, help="where it ran, e.g. usegalaxy.eu"
    )
    assemble_import.add_argument(
        "--retrieved-at", required=True, help="the date it was retrieved, YYYY-MM-DD"
    )
    assemble_import.add_argument(
        "--source-url", default=None, help="a link to the run, such as a history"
    )

    panel = subparsers.add_parser("panel", help="find, propose and fetch reference genomes")
    panel_actions = panel.add_subparsers(dest="action", metavar="ACTION")
    refresh = panel_actions.add_parser(
        "refresh", help="rebuild the list of RefSeq plastomes for a taxon, from NCBI"
    )
    refresh.add_argument("--taxon", required=True, help='taxon name, e.g. "Iridaceae"')
    refresh.add_argument(
        "--index",
        default=".panel/index.json",
        metavar="FILE",
        help="where to write the list (default: .panel/index.json)",
    )
    refresh.add_argument(
        "--max-records",
        type=int,
        default=200,
        help="the most records to ask NCBI for (default: 200)",
    )
    propose_parser = panel_actions.add_parser(
        "propose", help="list candidate genomes for a target from that list, without the network"
    )
    propose_parser.add_argument(
        "--index",
        default=".panel/index.json",
        metavar="FILE",
        help="the list `plastid panel refresh` wrote (default: .panel/index.json)",
    )
    propose_parser.add_argument(
        "--target", required=True, metavar="ACCESSION", help="accession of the target genome"
    )
    propose_parser.add_argument(
        "--mode",
        choices=PROPOSAL_MODES,
        default=DEFAULT_PROPOSAL_MODE,
        help="all: every genome at --level or closer; tiers: up to --per-tier genomes from each "
        "taxonomic tier (default: tiers)",
    )
    propose_parser.add_argument(
        "--level",
        choices=PROPOSAL_LEVELS,
        default="family",
        help="with --mode all: how far from the target to go (default: family)",
    )
    propose_parser.add_argument(
        "--per-tier",
        type=int,
        default=3,
        help="with --mode tiers: how many genomes to take from each tier (default: 3)",
    )
    pin_parser = panel_actions.add_parser(
        "pin", help="download the chosen genomes and record their checksums"
    )
    pin_parser.add_argument(
        "--index",
        default=None,
        metavar="FILE",
        help="list to take entries from; without it every accession is resolved from NCBI",
    )
    wanted = pin_parser.add_mutually_exclusive_group(required=True)
    wanted.add_argument("--accessions", help="comma-separated accessions")
    wanted.add_argument(
        "--groups",
        metavar="FILE",
        help="YAML file of sample groups, {group: [accessions]}; pins every accession in it",
    )
    wanted.add_argument(
        "--file",
        action="append",
        metavar="PATH",
        help="pin a record produced elsewhere, such as a web service or a cluster, repeatable",
    )
    pin_parser.add_argument(
        "--note",
        default=None,
        help="with --file: which service produced it and when, e.g. 'PMGA web, 2026-10-03'",
    )
    pin_parser.add_argument(
        "--to",
        default=".panel/genomes",
        metavar="DIR",
        help="directory to pin into (default: .panel/genomes)",
    )

    reference = subparsers.add_parser(
        "reference", help="check a set of reference records before annotating from it"
    )
    reference_actions = reference.add_subparsers(dest="action", metavar="ACTION")
    ref_verify = reference_actions.add_parser(
        "verify", help="hash every pinned record again and name every file that changed"
    )
    ref_verify.add_argument(
        "--reference-set",
        default=".panel/plastid-reference",
        metavar="DIR",
        help=_REFERENCE_SET_HELP,
    )
    ref_audit = reference_actions.add_parser(
        "audit", help="audit a set: block what cannot be used, report what a person must know"
    )
    ref_audit.add_argument(
        "--reference-set",
        default=".panel/plastid-reference",
        metavar="DIR",
        help=_REFERENCE_SET_HELP,
    )
    ref_audit.add_argument("-o", "--out", default=".runs", help=_OUT_HELP)
    ref_audit.add_argument(
        "--exceptions",
        default=None,
        metavar="FILE",
        help="the pinned exceptions (default: AUDIT_EXCEPTIONS.json in the set); an exception is "
        "bound to the file's digest, so editing the file voids it",
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
        help="pin a confirmed exception for one reported check, then stop; the digest is taken "
        "from the file as it is now",
    )
    ref_audit.add_argument("--confirmed-by", default=None, help="with --confirm-exception: who")
    ref_audit.add_argument(
        "--reason", default=None, help="with --confirm-exception: why the finding is accepted"
    )
    ref_loo = reference_actions.add_parser(
        "leave-one-out",
        help="write which records may annotate which target, each one left out of its own set",
    )
    ref_loo.add_argument(
        "--reference-set",
        default=".panel/plastid-reference",
        metavar="DIR",
        help=_REFERENCE_SET_HELP,
    )
    ref_loo.add_argument(
        "--family",
        required=True,
        help='the family the set was chosen at, e.g. "Iridaceae"',
    )
    ref_loo.add_argument(
        "--rule", default="same family first", help="the sentence the set was chosen by"
    )
    ref_loo.add_argument(
        "--to",
        default=None,
        metavar="PATH",
        help="file or directory to write to (default: REFERENCE_SET.json in the set)",
    )

    editing = subparsers.add_parser(
        "editing-reference", help="the reference of the RNA-editing predictor: rebuild or verify it"
    )
    editing_actions = editing.add_subparsers(dest="action", metavar="ACTION")
    editing_build = editing_actions.add_parser(
        "build", help="rebuild it from its pinned records, the sourced edit list and pinned MAFFT"
    )
    editing_build.add_argument(
        "--records",
        default=".panel/rna_editing",
        metavar="DIR",
        help="directory holding the seven pinned records as <accession>.gb (default: "
        ".panel/rna_editing); `plastid panel pin` fetches them, and any other bytes are refused",
    )
    editing_build.add_argument(
        "--to", required=True, metavar="FILE", help="where to write the reference JSON"
    )
    editing_build.add_argument("-o", "--out", default=".runs", help=_OUT_HELP)
    editing_verify = editing_actions.add_parser(
        "verify", help="check it against its digest, and the records against their pins"
    )
    editing_verify.add_argument(
        "--reference-file",
        default=None,
        metavar="FILE",
        help="a reference JSON (default: the packaged one)",
    )
    editing_verify.add_argument(
        "--records",
        default=None,
        metavar="DIR",
        help="also hash the pinned records in this directory again",
    )

    figures_parser = subparsers.add_parser(
        "figures", help="redraw a run's figures from its stored tables, without running anything"
    )
    figures_parser.add_argument("run", metavar="RUN_DIR", help=_RUN_DIR_HELP)
    figures_parser.add_argument(
        "--figure",
        action="append",
        default=None,
        metavar="ID",
        help="draw only this figure, repeatable (default: every figure the run's tables can feed)",
    )
    figures_parser.add_argument(
        "--list",
        action="store_true",
        help="list every figure, whether this run can draw it, and the drawing options, then stop",
    )
    figures_parser.add_argument("--json", action="store_true", help="with --list: " + _JSON_HELP)
    figures_parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        dest="sets",
        help="drawing option, repeatable; --list names them",
    )

    package_parser = subparsers.add_parser(
        "package", help="rebuild a run's Supplementary package: tables, figures, reports, index"
    )
    package_parser.add_argument("run", metavar="RUN_DIR", help=_RUN_DIR_HELP)

    qc_parser = subparsers.add_parser(
        "qc", help="write a run's technical QC report: what ran, versions, inputs, audits, figures"
    )
    qc_parser.add_argument("run", metavar="RUN_DIR", help=_RUN_DIR_HELP)

    junction_parser = subparsers.add_parser(
        "ir-junction", help="draw the IR junctions (JLB, JSB, JSA, JLA) of records, with CPJSdraw"
    )
    junction_parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        metavar="GENBANK",
        help="a GenBank record with its IR annotated; repeat it, one record per species",
    )
    _add_run_directory(junction_parser, "ir-junction", fixed="ir-junction")
    junction_parser.add_argument(
        "--tag",
        default="ir_junction",
        help="base name of the figure and table (default: ir_junction)",
    )
    junction_parser.add_argument(
        "--width", type=int, default=3200, help="PNG width in pixels (default: 3200)"
    )
    junction_parser.add_argument(
        "--tree",
        default=None,
        metavar="NEWICK",
        help="a tree that sets the species order (default: the order of -i)",
    )
    junction_parser.add_argument(
        "--order-from",
        default=None,
        metavar="NAME",
        help="the species on --tree the others are ordered by distance from",
    )
    junction_parser.add_argument(
        "--cpjsdraw",
        default=None,
        metavar="PATH",
        help="the patched CPJSdraw.pl when it is not on PATH (docs/tools/cpjsdraw.md)",
    )
    junction_parser.add_argument(
        "--drop-pseudo", action="store_true", help="leave pseudogenes out of the figure"
    )

    mito = subparsers.add_parser("mito", help="mitochondrial annotation (prototype)")
    mito_actions = mito.add_subparsers(dest="action", metavar="ACTION")
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
        help="directory of GenBank references to transfer from (default: the packaged mito "
        "reference set)",
    )
    annotate.add_argument("-o", "--out", default=".runs", help=_OUT_HELP)
    annotate.add_argument("--dataset", default=None, help="dataset label (default: accession)")
    annotate.add_argument("--force", action="store_true", help="ignore the cache and re-run")
    annotate.add_argument(
        "--set", action="append", default=[], metavar="KEY=VALUE", help="module parameter"
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
    return None if problems else params


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


def _sha256_file(path: Path) -> str:
    return sha256_file(path)


def _finish(manifest: dict, layout: RunLayout, ctx: RunContext, version: str, code: int) -> int:
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
    ``(module_id, version)`` for each, recorded as ``skipped`` with the reason (M11.16).
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
    return _finish(manifest, layout, ctx, version, code)


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
    return _finish(manifest, layout, ctx, "", EXIT_CANCELLED)


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
            return EXIT_USAGE
    repeated = sorted({spec.id for spec in specs if [s.id for s in specs].count(spec.id) > 1})
    if repeated:
        print(
            f"run: module(s) named twice: {', '.join(repeated)} — one run holds each module once",
            file=sys.stderr,
        )
        return EXIT_USAGE
    if _refused_for_own_command(specs):
        return EXIT_USAGE
    params = _checked_sets(args.set or [], [spec.id for spec in specs])
    if params is None:
        return EXIT_USAGE
    if args.label:
        params[DATA_LABEL_PARAM] = args.label

    try:
        genomes = tuple(_read_input(path) for path in args.input)
    except (OSError, ValueError) as error:
        print(f"cannot read the input: {error}", file=sys.stderr)
        return EXIT_USAGE

    if args.groups:
        # Shared group labels for every comparison output (P2-11). Loaded here because a module
        # never touches the filesystem; the module receives plain data in params.
        try:
            groups = load_groups(args.groups)
        except (OSError, ValueError) as error:
            print(f"cannot read {args.groups}: {error}", file=sys.stderr)
            return EXIT_USAGE
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
            return EXIT_USAGE
        params["groups"] = labels

    try:
        reference, genomes, inputs = _resolve_reference(specs, params, genomes, list(args.input))
    except (OSError, ValueError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
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
        force=args.no_cache,
        docx=args.docx,
        figures=args.figures,
        command=args.typed,
        inputs=inputs,
        references=reference,
    )


def _cmd_batch(args: argparse.Namespace) -> int:
    """Run each arity-1 module once per sample, then copy the stored cells into one table.

    This is the multi-sample orchestrator M11.6 needs, and it is deliberately thin: every sample
    goes through ``_execute_run`` (the same body as ``plastid run``), and the aggregate only reads
    what those runs stored. No value is recomputed, and a failed sample contributes a status, not
    a zero.
    """
    specs: list[registry.ModuleSpec] = []
    try:
        for module_id in args.module:
            spec = registry.get(module_id)
            if spec.id not in {item.id for item in specs}:
                specs.append(spec)
    except registry.UnknownModuleError as error:
        print(str(error), file=sys.stderr)
        return EXIT_USAGE
    if _refused_for_own_command(specs):
        return EXIT_USAGE
    multi = [spec.id for spec in specs if spec.arity != 1]
    if multi:
        print(
            f"batch applies arity-1 modules only; {', '.join(multi)} compares genomes within one "
            f"run — run it once over the whole set with `plastid run`",
            file=sys.stderr,
        )
        return EXIT_USAGE

    params = _checked_sets(args.set or [], [spec.id for spec in specs])
    if params is None:
        return EXIT_USAGE
    if args.label:
        params[DATA_LABEL_PARAM] = args.label

    batch_label = args.dataset or "batch"
    fingerprint = hashlib.sha256(
        json.dumps({"inputs": list(args.input), "params": params}, sort_keys=True).encode()
    ).hexdigest()
    layout = RunLayout.create(args.out, batch_label, fingerprint=fingerprint)
    layout.make_dirs()
    started = batch_run.start_batch(
        run_dir=layout.run_dir,
        dataset=batch_label,
        params=params,
        command=args.typed,
        inputs=list(args.input),
    )
    records: list[batch_run.SampleRecord] = []
    results: dict[str, dict] = {}
    codes: list[int] = []
    try:
        for position, text in enumerate(args.input, start=1):
            print(f"sample {position}/{len(args.input)}: {text}", file=sys.stderr, flush=True)
            record, sample_results = _batch_sample(text, specs, params, records, args)
            records.append(record)
            results[record.sample_id] = sample_results
            codes.append(record.exit_code)
            batch_run.record_sample(layout.run_dir, started, record)
            if record.exit_code == EXIT_CANCELLED:
                break
    except RunCancelled:
        print(f"cancelled — what had finished is in {layout.run_dir}", file=sys.stderr)
        codes.append(EXIT_CANCELLED)
    failed = [code for code in codes if code != EXIT_OK]
    code = EXIT_CANCELLED if EXIT_CANCELLED in failed else (failed[0] if failed else EXIT_OK)
    batch_run.write_batch(
        run_dir=layout.run_dir,
        dataset=batch_label,
        records=records,
        results=results,
        params=params,
        command=args.typed,
        inputs=list(args.input),
        started=started,
        exit_code=code,
    )
    if code == EXIT_CANCELLED:
        return code
    if args.figures:
        # The multi-sample figures read the aggregate, not the per-sample results: the same pass
        # `plastid figures <batch-dir>` runs, kept here so one command can finish the batch.
        from plastid.run.batch import read_batch_results

        payload = json.loads((layout.run_dir / "batch.json").read_text(encoding="utf-8"))
        outcome = _render_output_figures(
            layout.run_dir,
            sets=params,
            results=read_batch_results(layout.run_dir),
            include_matplotlib=False,
            allow_missing_extra=True,
            provenance=_batch_provenance(payload),
        )
        payload.setdefault("figures", {})["files"] = [
            {
                "path": path.relative_to(layout.run_dir).as_posix(),
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in outcome.files
        ]
        payload["figures"]["skipped"] = [
            {
                "figure": figure_id,
                "renderer": _renderer_of(outcome.records, figure_id),
                "reason": reason,
            }
            for figure_id, reason in outcome.skipped
        ]
        payload["figures"]["missing_renderers"] = [
            {"renderer": renderer, "reason": reason} for renderer, reason in outcome.missing
        ]
        payload["figures"]["status"] = outcome.as_manifest()["status"]
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
            for record in outcome.records
        ]
        if outcome.tools:
            payload["figures"]["tools"] = [asdict(record) for record in outcome.tools]
        (layout.run_dir / "batch.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        for path in outcome.files:
            print(f"wrote {path.relative_to(layout.run_dir).as_posix()}")
        for figure_id, reason in outcome.skipped:
            print(f"skipped {figure_id}: {reason}", file=sys.stderr)

    print(f"batch: {len(records) - len(failed)} of {len(records)} sample(s) ok — {layout.run_dir}")
    return code


def _batch_sample(
    text: str,
    specs: Sequence[registry.ModuleSpec],
    params: dict,
    records: Sequence[batch_run.SampleRecord],
    args: argparse.Namespace,
) -> tuple[batch_run.SampleRecord, dict]:
    """One sample through the same body as `plastid run`, read back from what its run wrote."""
    try:
        genome = _read_input(text)
    except (OSError, ValueError) as error:
        print(f"{text} — {error}", file=sys.stderr)
        record = batch_run.SampleRecord(
            sample_id=Path(text).stem,
            species="",
            accession="",
            input=text,
            sha256="",
            status="failed",
            exit_code=EXIT_USAGE,
            error=str(error),
        )
        return record, {}
    accession = genome.accession or genome.name
    label = accession
    taken = {record.sample_id for record in records}
    suffix = 2
    while label in taken:
        label = f"{accession}#{suffix}"
        suffix += 1
    layout = RunLayout.create(args.out, label, fingerprint=_run_fingerprint((genome,), params))
    code = _execute_run(
        genomes=(genome,),
        specs=tuple(specs),
        params=dict(params),
        out=args.out,
        dataset=label,
        force=args.no_cache,
        docx=args.docx,
        command=args.typed,
        inputs=[text],
        layout=layout,
    )
    record = batch_run.read_sample(
        sample_id=label,
        species=genome.name,
        accession=accession,
        input_path=text,
        sha256=_sha256_file(Path(text)),
        exit_code=code,
        run_dir=layout.run_dir,
    )
    return record, dict(batch_run.read_results(record))


def _resolve_reference(
    specs: Sequence[registry.ModuleSpec],
    params: dict,
    genomes: tuple[Genome, ...],
    inputs: list[str],
) -> tuple[dict | None, tuple[Genome, ...], list[str]]:
    """Put a module's declared reference in front of the run's inputs — read, audited, recorded.

    A module that names a ``reference_param`` is working against a *privileged record*, not
    against one of the run's equal inputs (D5-09: the naming reference, chosen by the end user
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
        path = _resolve_recorded_reference(recorded)

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


def _resolve_recorded_reference(recorded: RecordedReference) -> Path:
    """Where the recorded naming standard actually is, or an error saying what to do about it.

    Two shapes have to work. **From a checkout** (the team's runs, CI) the recorded path is
    relative to the checkout root, so a run started in any subdirectory still finds it. **From an
    installed package** there is no checkout: measured 11/10/2026 by the release's Docker check,
    the old code walked up from ``site-packages/plastid/core/gene_names.py``, found no
    ``pyproject.toml``, and ``gene_normalize`` died for anyone who had installed the toolkit rather
    than cloned it. There the path is taken as it comes — absolute, or relative to the working
    directory — and when it is not there the error names the way out: the naming reference is the
    user's own choice (D5-09), not a toolkit default.
    """
    root = checkout_root()
    if root is None:
        return recorded.resolve()
    return recorded.resolve(root=root)


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


def _ir_manifest_fields(ir: IRResult, *, requested_source: str = "annotation") -> dict:
    """The IR evidence a manifest records — success and refusal use the same shape.

    D11-2 requires a refused annotation to stay visible: the source, the measured columns, the
    annotated lengths and the reason the sequence path was not taken all travel in the manifest.

    ``requested_source`` is the operator's choice (``ir_source=annotation|sequence``) and
    ``source`` is the path that produced the result; they are different facts. ``status`` names
    the combination, so "annotation used", "sequence path with an annotation that was not used",
    "sequence path after the annotation was refused" and "no annotation at all" are four states
    rather than one "IR ok". A gate failure is reported by ``gate_passed``/``error``, not by
    collapsing the path state.
    """
    return {
        "requested_source": requested_source,
        "status": _ir_status(ir, requested_source),
        "source": ir.source,
        "annotation_used": ir.source == "annotation",
        "annotation_refused": ir.annotation_note is not None,
        "scoring_policy": IR_SCORING_POLICY,
        "scoring_policy_version": IR_SCORING_POLICY_VERSION,
        "length": ir.length,
        "identity": ir.identity,
        "gate_passed": ir.gate_passed,
        "mirror_constant": ir.mirror_constant,
        # End-to-end alignment evidence (D11-3): the denominator and every column class. The four
        # classes partition `compared_columns`: matched + mismatched + ambiguous + gaps.
        "matched_columns": ir.matched,
        "mismatched_columns": ir.mismatches,
        "ambiguous_columns": ir.ambiguous_columns,
        "gap_columns": ir.internal_gaps + ir.terminal_gaps,
        "compared_columns": ir.compared_length,
        "known_columns": ir.known_columns,
        "internal_gaps": ir.internal_gaps,
        "terminal_gaps": ir.terminal_gaps,
        "alignment_method": ir.method,
        # which end of a length difference the identity was measured at (H-25); "both" means
        # the copies have terminal overhangs at each end and mapping is not supported
        "excess_at": ir.excess_at,
        # The alignment-derived mapping the annotator uses to write a gene between copies of
        # different length (M11.11). None means a consumer that needs one must refuse by name.
        "mapping_supported": ir.mapping is not None,
        "mapping_method": ir.mapping.method if ir.mapping is not None else None,
        # Every copy length the record annotated, even when it was not used, and the reason the
        # annotation path was refused (D11-2). Together with `denovo_length`/`denovo_delta` these
        # are the independent checks a reader can repeat: the de-novo scan checks an annotated
        # result, and the annotation checks a de-novo result. Which one was available, and whether
        # it agreed, is recorded rather than left for the reader to infer from a null.
        "annotated_lengths": list(ir.annotated_lengths),
        "annotation_note": ir.annotation_note,
        "denovo_length": ir.denovo_length,
        "denovo_delta": ir.denovo_delta,
        "independent_check": _independent_check(ir)[0],
        "independent_agreement": _independent_check(ir)[1],
    }


def _independent_check(ir: IRResult) -> tuple[str, str]:
    """Which independent check this result has, and what it said.

    An annotated result is checked by the de-novo scan; a de-novo result is checked by the copy
    lengths the record annotated, when it has any. Neither is "not measured" when no second
    measurement exists — that word is reserved for the absence itself, never written as agreement.
    """
    if ir.denovo_length is not None:
        delta = ir.denovo_delta
        return "de_novo", "agrees" if delta == 0 else "differs"
    if ir.annotated_lengths and ir.source != "annotation":
        return "annotation", ("agrees" if set(ir.annotated_lengths) == {ir.length} else "differs")
    return "none", "not measured"


def _ir_status(ir: IRResult, requested_source: str) -> str:
    """Which IR path produced this result — one state per distinct situation (M11.16).

    ``source`` alone cannot answer this: the de-novo scan runs both when the operator asked for
    the sequence and when the annotation was refused or absent. The requested source separates
    those cases, and whether the refused annotation is named separates "refused" from "not used".
    """
    if ir.source == "manual":
        return "manual"
    if ir.source == "annotation":
        return "annotation"
    # The result came from the de-novo scan, not from the annotation.
    if ir.annotation_note is not None:
        return "sequence_annotation_refused"
    if ir.annotated_lengths:
        return "sequence_annotation_present"
    return "sequence_no_annotation"


def _run_fingerprint(genomes: Sequence[Genome], params: dict) -> str:
    """The run's identity from its ordered sequences and parameters (no clock, no path)."""
    overlay = json.dumps(params, sort_keys=True)
    sequences = "".join(item.sequence + "\x00" for item in genomes)  # ordered: order is meaning
    return hashlib.sha256((sequences + overlay).encode("ascii", "replace")).hexdigest()


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
    reference_sets: Sequence[dict] | None = None,
    layout: RunLayout | None = None,
) -> int:
    """One run's shared body: both ``plastid run`` and ``plastid mito annotate`` end up here.

    The point of sharing it is seam 3 — a second command must not become a second code path.
    Every output, every cache entry and every provenance record is written by the same lines,
    so an annotation produced by the mito command is traceable exactly like any other run.
    """

    label = dataset or genomes[0].accession or genomes[0].name
    if layout is None:
        layout = RunLayout.create(out, label, fingerprint=_run_fingerprint(genomes, params))
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
        reference_sets=list(reference_sets or []),
    )
    write_manifest(layout.run_dir / "manifest.json", manifest)
    contexts = [ctx]
    try:
        return _run_steps(
            genomes=genomes,
            specs=specs,
            params=params,
            manifest=manifest,
            layout=layout,
            contexts=contexts,
            label=label,
            force=force,
            docx=docx,
            figures=figures,
            references=references,
        )
    except RunCancelled:
        return _cancelled(
            manifest, layout, contexts[-1], [(spec.id, spec.version) for spec in specs]
        )


def _run_steps(
    *,
    genomes: Sequence[Genome],
    specs: Sequence[registry.ModuleSpec],
    params: dict,
    manifest: dict,
    layout: RunLayout,
    contexts: list[RunContext],
    label: str,
    force: bool,
    docx: bool,
    figures: bool,
    references: dict | None,
) -> int:
    """The gate, then each module in order. ``contexts`` tells the caller which step a stop reached.

    The manifest on disk moves with the run, so a viewer sees which step it is in and what the
    finished steps wrote.
    """
    ctx = contexts[-1]
    # A module that computes against a packaged reference (M4.5: the RNA-editing reference B′)
    # names it in the manifest by its digest, and the digest enters the step's parameters so the
    # cache key carries it explicitly. A reference that fails its own digest stops the run here.
    for spec in specs:
        if spec.reference_set is None:
            continue
        try:
            entry = spec.reference_set()
        except EditingReferenceError as error:
            return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_ERROR)
        manifest["reference_sets"].append(entry)
        params[f"{spec.id}_reference_digest"] = entry["digest"]
    # The set a GenBank artifact was transferred from; recorded on the artifact itself so a
    # consumer does not have to infer it from the module's parameters (M11.16).
    annotation_reference = next(
        (
            item["sha256"]
            for item in manifest["reference_sets"]
            if item.get("role") == "annotation_db"
        ),
        None,
    )

    if any(spec.requires_ir for spec in specs):
        requested = str(params.get("ir_source", "annotation"))
        gated: list[IRResult] = []
        for candidate in genomes:
            attempted: IRResult | None = None
            try:
                attempted = detect_ir(candidate, source=requested)
                gated.append(enforce_ir_gate(attempted))
            except (IRGateError, IRDetectionError) as error:
                # every genome in a comparison must pass: one diverged IR corrupts them all
                manifest["ir"] = {
                    **(
                        _ir_manifest_fields(attempted, requested_source=requested)
                        if attempted is not None
                        else {"requested_source": requested, "status": "undetected"}
                    ),
                    "gate_passed": False,
                    "error": str(error),
                }
                return _fail(
                    manifest, layout, ctx, specs[0].version, f"IR gate: {error}", EXIT_USAGE
                )
        ctx.ir = gated[0]
        manifest["ir"] = _ir_manifest_fields(gated[0], requested_source=requested)
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
                    "excess_at": result.excess_at,
                }
                for candidate, result in zip(genomes, gated, strict=True)
            ]

    citations = load_citations()
    gated_ir = gated[0] if any(spec.requires_ir for spec in specs) else None
    documents: list[Path] = []
    files: list[dict] = manifest["files"]
    for order, spec in enumerate(specs, start=1):
        ctx = RunContext(layout=layout, params=params, module_id=spec.id, module_order=order)
        ctx.ir = gated_ir
        contexts.append(ctx)
        manifest["current"] = {
            "module": spec.id,
            "order": order,
            "of": len(specs),
            "started_at": datetime.now().isoformat(timespec="seconds"),
        }
        write_manifest(layout.run_dir / "manifest.json", manifest)
        print(f"[{order}/{len(specs)}] {spec.id}", file=sys.stderr, flush=True)
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
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_USAGE,
                downstream=_downstream(specs, order),
            )
        except registry.ParamValueError as error:
            # 2, not 1: the key is declared and the value is one the user can fix on the command
            # line. Raised before the module calls its tool, so nothing is produced from it.
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_USAGE,
                downstream=_downstream(specs, order),
            )
        except UnsupportedTranslationError as error:
            # 2, next to the IR gate: the run was well-formed, and the record's own CDS declare a
            # frame or table this codon reading cannot translate honestly (M11.5).
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_USAGE,
                downstream=_downstream(specs, order),
            )
        except reference_audit.ReferenceAuditError as error:
            # 2, the gate's own number, next to the IR gate: the run was well-formed and it is the
            # reference set that cannot be used. The message is the audit's own `file — reason`
            # lines; the way out (drop the file, or re-check the set) is said here, once.
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                f"{error}\n`plastid reference audit` re-checks the set without annotating "
                f"anything; `--exclude <accession>` drops a record from the search, and the run "
                f"says which records it did not search",
                EXIT_USAGE,
                downstream=_downstream(specs, order),
            )
        except NamingConflictError as error:
            # 4, not 2: the run was well-formed and the data disagrees with the reference — a
            # number of its own keeps it apart from the IR gate and the audit (D5-09).
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_CONFLICT,
                downstream=_downstream(specs, order),
            )
        except ExportDependencyError as error:
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_TOOL,
                downstream=_downstream(specs, order),
            )
        except ToolError as error:
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                str(error),
                EXIT_TOOL,
                downstream=_downstream(specs, order),
            )
        except Exception as error:  # noqa: BLE001 — report any module failure, then exit 1
            message = f"{type(error).__name__}: {error}"
            return _fail(
                manifest,
                layout,
                ctx,
                spec.version,
                message,
                EXIT_ERROR,
                downstream=_downstream(specs, order),
            )

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
            # Current-run tool calls. Empty on `cached`, which is the point: the computation ran
            # in an earlier run, and `provenance` below says which one (M11.12).
            "tools": [asdict(record) for record in ctx.tools],
            "provenance": outcome.provenance,
            "files": [entry["path"] for entry in module_files],
        }
        if spec.reference_param and references is not None:
            # The identity of the privileged record this step worked against — path *and*
            # digest, because a path alone cannot say whether the reference changed (D5-09).
            module_entry["reference"] = references
        manifest["modules"].append(module_entry)
        for entry in module_files:
            if Path(entry["path"]).suffix != ".gb":
                continue
            # A GenBank output is the file the next step reads: say what it is and what it was
            # transferred from, so a consumer can refuse a missing, changed or failed artifact
            # instead of trusting the run directory's shape (M11.16).
            try:
                manifest["artifacts"].append(
                    describe_artifact(
                        layout.run_dir / entry["path"],
                        relative_to=layout.run_dir,
                        sample_id=genomes[-1].accession or genomes[-1].name,
                        topology=genomes[-1].topology,
                        status=outcome.status,
                        reference_sha256=annotation_reference,
                    )
                )
            except ContractError as error:
                return _fail(manifest, layout, ctx, spec.version, str(error), EXIT_ERROR)

    write_manifest(layout.run_dir / "manifest.json", manifest)
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
        # the analysis ran or was restored (P5-01). Both catalogues are collected here, so
        # `run --figures` and `plastid figures` produce the same set (M11.18); a missing plotting
        # extra is still the explicit request's failure, while a missing R is a recorded skip.
        try:
            outcome = _render_output_figures(
                layout.run_dir,
                sets=params,
                allow_missing_extra=False,
                provenance=_manifest_provenance(manifest),
            )
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
        manifest["figures"] = outcome.as_manifest()
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
    """The run-level lines every figure and package must carry (M11.18).

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


def _cmd_package(args: argparse.Namespace) -> int:
    """The Supplementary package, rebuilt from the run's own stored results.

    A re-package draws figures from the stored tables and copies them, together with the draft
    report, into one folder; it never re-runs an analysis. A missing plotting extra or a missing
    R does not throw the tables away — the absence is written into ``index.md`` and the package
    still exists (M11.18).
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
    # so a later `plastid qc` reports the real reason instead of guessing (M11.18). A figure the
    # refresh deleted must leave the run's file list, and a rewritten one must carry its new
    # digest — an old path or an old sha256 would pass a stale result off as current.
    manifest["figures"] = outcome.as_manifest()
    refreshed = {
        path.relative_to(run_dir).as_posix(): {
            "path": path.relative_to(run_dir).as_posix(),
            "module": "figures",
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
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
            {
                "path": path.relative_to(package.directory).as_posix(),
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
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
    reason recorded by the last render when there is one (M11.18).
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
        inputs=[
            {"path": str(path), "sha256": _sha256_file(path), "bytes": path.stat().st_size}
            for path, _ in records
        ],
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
    from plastid.external.cpjsdraw import run_cpjsdraw, svg_to_png

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
        _finish_manifest(manifest, manifest_path, EXIT_CANCELLED)
        return EXIT_CANCELLED
    except ToolError as error:
        print(str(error), file=sys.stderr)
        _finish_manifest(manifest, manifest_path, EXIT_TOOL)
        return EXIT_TOOL
    except JunctionError as error:
        print(f"cannot draw the junction figure: {error}", file=sys.stderr)
        _finish_manifest(manifest, manifest_path, EXIT_USAGE)
        return EXIT_USAGE
    except (OSError, ValueError) as error:
        print(f"cannot draw the junction figure: {error}", file=sys.stderr)
        _finish_manifest(manifest, manifest_path, EXIT_ERROR)
        return EXIT_ERROR

    for path in (config, table, svg, png):
        relative = path.relative_to(layout.run_dir).as_posix()
        manifest["files"].append(
            {
                "path": relative,
                "module": "ir_junction",
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
        print(f"wrote {relative}")
    for row in rows:
        if row.warning:
            print(f"warning: {row.warning}", file=sys.stderr)
    _finish_manifest(manifest, manifest_path, EXIT_OK)
    return EXIT_OK


def _finish_manifest(manifest: dict[str, Any], path: Path, code: int) -> None:
    """Close a manifest — a run that failed still says what it was and what it wrote."""
    manifest["status"] = exit_state(code).status
    manifest["exit_code"] = code
    manifest["finished_at"] = datetime.now().isoformat(timespec="seconds")
    write_manifest(path, manifest)


def _junction_order(
    records: Sequence[tuple[Path, Genome]], tree: str | None, reference: str | None
) -> list[str] | None:
    """Species names ordered by distance from a reference on the user's own tree.

    The tree is the user's (PTS builds none — M10.5), and its tip labels are not this toolkit's
    species names, so a tip is matched to a record by the first two words of the organism name.
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

    Two catalogues, one pass (M11.18). The matplotlib figures are this toolkit's own drawing
    (vector only); the R figures are the manuscript's scripts, drawn from the same tables — a
    second source for a figure that already exists, and the only source for the ones the toolkit
    has no drawing for. A failure in one does not stop the other: they are independent, and a run
    with no plotting extra installed can still have its R figures. Every figure the pass did not
    draw is printed with its reason, and recorded in the manifest or ``batch.json``.

    A batch has no single genome to draw the one-genome catalogue for; its sample figures are the
    R ones, read from the aggregate table (M11.7).
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
            sets=args.sets,
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
            {
                "path": path.relative_to(run_dir).as_posix(),
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in outcome.files
        ]
        payload["figures"]["skipped"] = [
            {"figure": figure_id, "renderer": _renderer_of(records, figure_id), "reason": reason}
            for figure_id, reason in outcome.skipped
        ]
        payload["figures"]["missing_renderers"] = [
            {"renderer": renderer, "reason": reason} for renderer, reason in outcome.missing
        ]
        payload["figures"]["status"] = outcome.as_manifest()["status"]
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
        for figure_id, reason in outcome.skipped:
            print(f"skipped {figure_id}: {reason}", file=sys.stderr)
        for renderer, reason in outcome.missing:
            print(f"{renderer} figures unavailable: {reason}", file=sys.stderr)
        return EXIT_OK if outcome.files else EXIT_TOOL

    # A refresh can delete a figure whose table disappeared, and it rewrites the ones it does
    # draw. Both cases must reach the run's file list: a deleted path and an old digest are stale
    # results pretending to be current (M11.18).
    refreshed = {
        path.relative_to(run_dir).as_posix(): {
            "path": path.relative_to(run_dir).as_posix(),
            "module": "figures",
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
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
    for figure_id, reason in outcome.skipped:
        print(f"skipped {figure_id}: {reason}", file=sys.stderr)
    for renderer, reason in outcome.missing:
        print(f"{renderer} figures unavailable: {reason}", file=sys.stderr)
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
        _print_json(
            {
                "figures": [
                    {
                        "id": record.figure_id,
                        "renderer": record.renderer,
                        "drawn": bool(record.files),
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
    """One figure pass over both catalogues, with every absence accounted for (M11.18).

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
            "missing_renderers": [
                {"renderer": renderer, "reason": reason} for renderer, reason in self.missing
            ],
            "status": (
                "complete"
                if not self.skipped and not self.missing
                else ("partial" if self.files else "unavailable")
            ),
            "tools": [asdict(record) for record in self.tools],
        }


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
        # (K-06). Measured 27/10/2026: these were reachable only by calling R by hand.
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


MITO_MODULE = "mito_annotate"
PLASTID_MODULE = "plastid_annotate"


def _cmd_annotate(args: argparse.Namespace) -> int:
    """Annotate 1..N plastomes from the reference set — one run per sample (P3B-09, M3B.4).

    **The batch is a batch of runs, not one run with N samples.** Each sample gets the run
    directory, manifest, tool record and report that every other output in this toolkit has
    (output is split by run/dataset and never overwritten in place), so a sample's
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

    params = _checked_sets(args.set or [], [spec.id])
    if params is None:
        return EXIT_USAGE
    if args.reference_set:
        params["reference"] = args.reference_set
    if args.label:
        params[DATA_LABEL_PARAM] = args.label
    if args.exclude:
        params["exclude"] = list(args.exclude)

    directory = plastid_reference.resolve_directory(params.get("reference"))
    missing = plastid_reference.missing_set(directory)

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
            codes.append(EXIT_USAGE)
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

    # An external GenBank already carries an annotation; re-annotating it would discard that
    # annotation without saying so. The pipeline has two input paths — a bare assembly goes
    # through this command, an annotated record goes straight to analysis or verify — and the
    # choice is the user's, so it is explicit (M11.16).
    if not _truthy(params.get("reannotate")):
        already = [
            (text, genome)
            for text, genome in samples
            if any(gene.kind in ("CDS", "tRNA", "rRNA") for gene in genome.genes)
        ]
        if already:
            for text, genome in already:
                count = sum(1 for gene in genome.genes if gene.kind in ("CDS", "tRNA", "rRNA"))
                print(
                    f"{text} already carries {count} CDS/tRNA/rRNA feature(s) — "
                    f"re-annotating it would discard them. Analyse it with "
                    f"`plastid run -m <module> -i {text}`, or pass `--set reannotate=yes` to "
                    f"replace its annotation from the reference set.",
                    file=sys.stderr,
                )
            print("annotate: nothing was annotated", file=sys.stderr)
            return EXIT_USAGE

    if samples and missing is not None:
        print(f"annotate: {missing}", file=sys.stderr)
        return EXIT_TOOL

    for text, genome in samples:
        code = _annotate_sample(text=text, genome=genome, spec=spec, params=params, args=args)
        if code == EXIT_CANCELLED:
            return code
        codes.append(code)

    failed = [code for code in codes if code != EXIT_OK]
    if failed:
        print(
            f"annotate: {len(failed)} of {len(codes)} sample(s) failed — the reason is above, "
            f"and a sample that reached its run directory records it there",
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
    directory = plastid_reference.resolve_directory(params.get("reference"))
    reference_sets: list[dict] = []
    if str(directory) not in ("", "."):
        planted = params.get("exclude")
        excluded = (planted,) if isinstance(planted, str) else tuple(planted or ())
        try:
            reference_sets.append(
                describe_reference_set(
                    directory,
                    exclude=excluded,
                    leave_out=genome.accession or genome.name,
                )
            )
        except ContractError as error:
            print(str(error), file=sys.stderr)
            return EXIT_USAGE
    return _execute_run(
        genomes=(genome,),
        specs=(spec,),
        params=dict(params),
        out=args.out,
        dataset=args.dataset or genome.accession or genome.name,
        force=args.no_cache,
        docx=args.docx,
        command=args.typed,
        inputs=[text],
        reference_sets=reference_sets,
    )


def _cmd_mito(args: argparse.Namespace) -> int:
    return _cmd_mito_annotate(args)


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
        return EXIT_USAGE

    params = _checked_sets(args.set or [], [spec.id])
    if params is None:
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
        command=args.typed,
        inputs=[args.input],
    )


def _cmd_reference(args: argparse.Namespace) -> int:
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
    """Audit a reference set before anything is annotated from it (P3B-03, D5-08).

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
    """The audit as a run — its own directory, manifest and tool record (P3B-03, P3B-09)."""
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
    """Pin records produced outside this machine — the other half of the backend rule (P1-12)."""
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

    A heavy tool runs where it fits (P1-12) — this machine when it can, a web
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


def _cmd_editing_reference(args: argparse.Namespace) -> int:
    if args.action == "build":
        return _editing_reference_build(args)
    return _editing_reference_verify(args)


def _editing_reference_build(args: argparse.Namespace) -> int:
    """Rebuild reference B′ (M4.5) as a run of its own — manifest, tool logs, the JSON written.

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
    digests = {accession: _sha256_file(path) for accession, path in pairs}

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
    manifest["files"] = [
        {
            "path": str(copy.relative_to(layout.run_dir)),
            "module": "editing_reference",
            "sha256": _sha256_file(copy),
            "bytes": copy.stat().st_size,
        }
    ]
    manifest["reference_sets"].append(
        {"role": "rna_editing_reference", "digest": payload["digest"], "path": str(out)}
    )
    manifest["status"] = "ok"
    print(f"reference B′ digest {payload['digest']} — {out}")
    return _finish(manifest, layout, ctx, "1", EXIT_OK)


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
            elif _sha256_file(path) != taxon["sha256"]:
                problems.append(f"{path} — not the pinned bytes ({taxon['sha256'][:12]}…)")
    for problem in problems:
        print(problem, file=sys.stderr)
    if problems:
        return EXIT_ERROR
    print(f"reference B′ {reference.digest} — digest matches its contents")
    return EXIT_OK


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


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
        inputs=[
            {"path": str(path), "sha256": _sha256_file(path), "bytes": path.stat().st_size}
            for path in inputs
        ],
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
        [*(_sha256_file(file) for file in reads.files), "prepare"],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
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
    code = _finish(manifest, layout, ctx, assembly_step.MODULE_VERSION, EXIT_OK)
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
    digest = _sha256_file(Path(args.compare_to))
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
        [_sha256_file(Path(args.prepared)), args.source_name],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
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
        [*(_sha256_file(file) for file in reads.files), str(sorted(params.items()))],
    )
    ctx = RunContext(
        layout=layout,
        params=params,
        module_id=assembly_step.MODULE_ID,
        module_order=1,
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
        {
            "path": str(path.relative_to(layout.run_dir)),
            "module": assembly_step.MODULE_ID,
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
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
    manifest["status"] = "ok" if usable else "incomplete"
    code = _finish(manifest, layout, ctx, version, EXIT_OK if usable else EXIT_ERROR)
    if not usable:
        print(
            f"assemble: {assembled.sample}: the assembler says {assembled.status_text!r}, read "
            f"support is {evaluation.support}, the graph is {evaluation.graph_call}: no artifact "
            f"for the next step — {layout.run_dir}",
            file=sys.stderr,
        )
        if assembled.word_limit_hit:
            print(
                "assemble: it stopped recruiting reads at its word limit. More memory (--set "
                "memory_gb), or a run elsewhere through `plastid assemble prepare`, would recruit "
                "more.",
                file=sys.stderr,
            )
    return code


def _print_json(payload: Any) -> None:
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def _print_table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    widths = [max(len(row[i]) for row in (header, *rows)) for i in range(len(header) - 1)]
    for row in (header, *rows):
        cells = [cell.ljust(width) for cell, width in zip(row, widths, strict=False)]
        print("  ".join((*cells, row[-1])).rstrip())


def _param_record(param: registry.ParamSpec) -> dict[str, Any]:
    return {
        "name": param.name,
        "type": param.type,
        "help": param.help,
        "default": param.default,
        "choices": list(param.choices),
        "flag": param.flag or None,
    }


def _module_record(spec: registry.ModuleSpec) -> dict[str, Any]:
    """One module as data: what a person reads, and what a form is built from."""
    return {
        "id": spec.id,
        "version": spec.version,
        "title": spec.title,
        "group": spec.group,
        "genomes": {"min": spec.arity, "max": 1 if spec.arity == 1 else None},
        "ir_gate": spec.requires_ir,
        "tools": list(spec.tools),
        "optional_tools": list(spec.optional_tools),
        "command": spec.command or "plastid run",
        "in_batch": spec.arity == 1 and spec.command is None,
        "reference_param": spec.reference_param,
        "params": [_param_record(param) for param in spec.params],
    }


def _assembly_record() -> dict[str, Any]:
    """The assembly step in the shape of a module: `plastid assemble` is not a registered one."""
    return {
        "id": assembly_step.MODULE_ID,
        "version": assembly_step.MODULE_VERSION,
        "title": "Plastome assembly from short or HiFi reads, checked against the reads",
        "group": "plastid",
        "tools": ["getorganelle", "fastp", "minimap2"],
        "optional_tools": ["oatk"],
        "command": "plastid assemble run",
        "params": [_param_record(param) for param in assembly_step.ASSEMBLY_PARAM_SPECS],
    }


def _describe_param(param: registry.ParamSpec) -> None:
    if param.flag:
        setting = f"set with {param.flag}"
    else:
        default = "not set" if param.default is None else param.default
        choices = f" · one of {', '.join(str(c) for c in param.choices)}" if param.choices else ""
        setting = f"{param.type} · default {default}{choices}"
    print(f"  {param.name}: {setting}")
    print(f"      {param.help}")


def _describe_module(spec: registry.ModuleSpec) -> None:
    usage = (spec.command or f"plastid run -m {spec.id}") + " -i <file>"
    genomes = "one" if spec.arity == 1 else f"at least {spec.arity}"
    print(f"{spec.id} {spec.version}: {spec.title}")
    print(f"  run with: {usage}")
    print(f"  genomes:  {genomes}")
    print(f"  IR gate:  {'yes' if spec.requires_ir else 'no'}")
    print(f"  tools:    {', '.join(spec.tools) or 'none'}")
    if spec.optional_tools:
        print(f"  if asked: {', '.join(spec.optional_tools)}")
    if spec.params:
        print("\nparameters (--set KEY=VALUE):")
        for param in spec.params:
            _describe_param(param)
    print("\nalso accepted on any run:")
    for param in registry.SHARED_PARAMS:
        _describe_param(param)


def _cmd_modules(args: argparse.Namespace) -> int:
    if args.module == assembly_step.MODULE_ID:
        record = _assembly_record()
        if args.json:
            _print_json(record)
            return EXIT_OK
        print(f"{record['id']} {record['version']}: {record['title']}")
        print(f"  run with: {record['command']}")
        print(f"  tools:    {', '.join(record['tools'])}; for HiFi reads, oatk")
        print("\nparameters (--set KEY=VALUE):")
        for param in assembly_step.ASSEMBLY_PARAM_SPECS:
            _describe_param(param)
        return EXIT_OK
    if args.module:
        try:
            spec = registry.get(args.module)
        except registry.UnknownModuleError as error:
            print(str(error), file=sys.stderr)
            return EXIT_USAGE
        if args.json:
            _print_json(_module_record(spec))
        else:
            _describe_module(spec)
        return EXIT_OK
    specs = [spec for _name, spec in sorted(registry.all_modules().items())]
    if args.json:
        _print_json(
            {
                "modules": [_module_record(spec) for spec in specs],
                "assembly": _assembly_record(),
                "shared_params": [_param_record(param) for param in registry.SHARED_PARAMS],
                "exit_states": [asdict(state) for state in EXIT_STATES],
            }
        )
        return EXIT_OK
    _print_table(
        ("MODULE", "GENOMES", "IR GATE", "TOOLS", "WHAT IT DOES"),
        [
            (
                spec.id,
                "1" if spec.arity == 1 else f"{spec.arity}+",
                "yes" if spec.requires_ir else "no",
                ", ".join(spec.tools) or "-",
                spec.title,
            )
            for spec in specs
        ],
    )
    print(
        "\n`plastid modules <module>` shows one module and its parameters; "
        "`plastid modules assembly` shows the assembly step."
    )
    return EXIT_OK


def _cmd_tools(args: argparse.Namespace) -> int:
    """What a run on this machine would find: a report, so a missing tool is not a failure here."""
    from plastid.external.inventory import check_tools

    needed_by: dict[str, list[str]] = {}
    for name, spec in sorted(registry.all_modules().items()):
        for tool in (*spec.tools, *spec.optional_tools):
            needed_by.setdefault(tool, []).append(name)
    statuses = check_tools()
    records = [
        {
            **status.as_dict(),
            "used_by": needed_by.get(status.name)
            or (
                ["plastid assemble"] if status.role == "assembly" else [_TOOL_COMMANDS[status.name]]
            ),
        }
        for status in statuses
    ]
    if args.json:
        _print_json({"tools": records})
        return EXIT_OK
    _print_table(
        ("TOOL", "STATUS", "VERSION", "PINNED", "NEEDED BY"),
        [
            (
                record["name"],
                record["status"],
                record["version"] or "-",
                record["pinned"],
                ", ".join(record["used_by"]),
            )
            for record in records
        ],
    )
    for record in records:
        if record["status"] != "ok":
            print(f"\n{record['message']}")
    return EXIT_OK


def _cmd_runs(args: argparse.Namespace) -> int:
    root = Path(args.root or ".runs")
    if args.root and not root.is_dir():
        print(f"{root}: no such directory", file=sys.stderr)
        return EXIT_USAGE
    runs = list_runs(root) if root.is_dir() else []
    if args.json:
        _print_json({"runs": [state.as_dict() for state in runs]})
        return EXIT_OK
    if not runs:
        print(f"no runs in {root}")
        return EXIT_OK
    _print_table(
        ("STATUS", "DONE", "STARTED", "DATASET", "RUN"),
        [
            (
                state.status,
                f"{state.done}/{state.total}" if state.total is not None else "-",
                (state.started_at or "-").replace("T", " "),
                state.dataset,
                state.run_id,
            )
            for state in runs
        ],
    )
    return EXIT_OK


def _raise_cancelled(signum: int, _frame: object) -> None:
    # Once only: the run is now writing down that it stopped, and must be left to finish doing so.
    for sent in _STOP_SIGNALS:
        signal.signal(sent, signal.SIG_IGN)
    raise RunCancelled(signal.Signals(signum).name)


@contextmanager
def _signals_cancel() -> Iterator[None]:
    """For one command, a signal to stop means the run stops in good order.

    A signal already ignored stays ignored: that is how `nohup` keeps a long run alive after its
    terminal closes, and a run started that way has asked not to be stopped by it. SIGHUP is here
    because every tool runs in its own session, so a closing terminal reaches this process alone.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = {sent: signal.getsignal(sent) for sent in _STOP_SIGNALS}
    for sent, handler in previous.items():
        if handler != signal.SIG_IGN:
            signal.signal(sent, _raise_cancelled)
    try:
        yield
    finally:
        for sent, handler in previous.items():
            if handler is not None:
                signal.signal(sent, handler)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # What was typed, so a run's record can repeat it. Each command takes it from here rather than
    # rebuilding it, because a list rebuilt by hand is a list that loses a flag.
    args.typed = ["plastid", *(sys.argv[1:] if argv is None else argv)]
    with _signals_cancel():
        try:
            return _dispatch(parser, args)
        except RunCancelled:
            print("cancelled", file=sys.stderr)
            return EXIT_CANCELLED


def _subcommands(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    found = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    return dict(found[0].choices) if found else {}


def _dispatch(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    if args.version:
        print(f"plastid {__version__}")
        return EXIT_OK

    if args.command is None:
        parser.print_help()
        return EXIT_OK

    actions = _subcommands(_subcommands(parser)[args.command])
    if actions and args.action is None:
        print(f"{args.command}: choose an action — {' | '.join(actions)}", file=sys.stderr)
        return EXIT_USAGE

    if args.command == "modules":
        return _cmd_modules(args)

    if args.command == "tools":
        return _cmd_tools(args)

    if args.command == "runs":
        return _cmd_runs(args)

    if args.command == "run":
        return _cmd_run(args)

    if args.command == "figures":
        return _cmd_figures(args)

    if args.command == "ir-junction":
        return _cmd_ir_junction(args)

    if args.command == "package":
        return _cmd_package(args)

    if args.command == "qc":
        return _cmd_qc(args)

    if args.command == "batch":
        return _cmd_batch(args)

    if args.command == "annotate":
        return _cmd_annotate(args)

    if args.command == "mito":
        return _cmd_mito(args)

    if args.command == "panel":
        return _cmd_panel(args)

    if args.command == "reference":
        return _cmd_reference(args)

    if args.command == "editing-reference":
        return _cmd_editing_reference(args)

    if args.command == "assemble":
        return _cmd_assemble(args)

    # Unreachable: argparse rejects unknown subcommands before this point.
    return EXIT_USAGE
