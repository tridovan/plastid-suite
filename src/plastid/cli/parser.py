"""The command line itself: every subcommand, flag and help text."""

from __future__ import annotations

import argparse

from plastid.cli.common import (
    _JSON_HELP,
    _NO_CACHE_HELP,
    _OUT_HELP,
    _REFERENCE_SET_HELP,
    _RUN_DIR_HELP,
    _SET_HELP,
    _add_run_directory,
)
from plastid.core.panel import (
    DEFAULT_PROPOSAL_MODE,
    PROPOSAL_LEVELS,
    PROPOSAL_MODES,
)
from plastid.run.state import EXIT_STATES
from plastid.ui.catalogue import DEFAULT_PORT, LOOPBACK


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

    ui_parser = subparsers.add_parser(
        "ui", help="open the local interface in a browser; nothing leaves this computer"
    )
    ui_parser.add_argument(
        "--workspace",
        default=".",
        metavar="DIR",
        help="the folder the interface works in, made if it is not there: records, reads, "
        "references, runs and your own tools live in it",
    )
    ui_parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help="port on this computer; a busy port stops the command, it never moves to another",
    )
    ui_parser.add_argument(
        "--host",
        default=LOOPBACK,
        help="address to listen on; anything but 127.0.0.1 is refused outside the Docker image",
    )
    ui_parser.add_argument(
        "--no-browser", action="store_true", help="print the link and do not open a browser"
    )

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
