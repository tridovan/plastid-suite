"""`plastid annotate` and `plastid mito annotate`: a record annotated from a reference set."""

from __future__ import annotations

import argparse
import sys

from plastid.cli.common import (
    EXIT_CANCELLED,
    EXIT_ERROR,
    EXIT_OK,
    EXIT_TOOL,
    EXIT_USAGE,
    _checked_sets,
    _truthy,
)
from plastid.cli.reference import _names_a_record
from plastid.cli.run import _execute_run, _read_input
from plastid.core import organelle, plastid_reference
from plastid.core.genbank import read_fasta
from plastid.core.model import Genome
from plastid.modules import registry
from plastid.run.contract import ContractError, describe_reference_set
from plastid.run.runner import (
    DATA_LABEL_PARAM,
)

MITO_MODULE = "mito_annotate"
PLASTID_MODULE = "plastid_annotate"


def _cmd_annotate(args: argparse.Namespace) -> int:
    """Annotate 1..N plastomes from the reference set — one run per sample.

    **The batch is a batch of runs, not one run with N samples.** Each sample gets the run
    directory, manifest, tool record and report that every other output in this toolkit has
    (output is split by run/dataset and never overwritten in place), so a sample's
    annotation is traceable on its own — and the command's own numbers are the reference set it
    searched, recorded in each manifest, plus its exit code.

    **A sample that fails does not kill the batch.** It is named on stderr, its own run directory
    records the failure and the reason, and the command exits non-zero so a script cannot mistake
    a partial batch for a whole one ("a failed sample is marked and named without killing
    the run, but the exit code reflects that something failed"). The first failure's code is what
    the command returns: the number is about *why* the run was refused, and the samples that
    failed for the most specific reason are the ones that ran furthest.

    **One organelle, declared at input.** A plastid sample that declares ``/organelle`` otherwise
    is refused by name, as is a reference set of mito records, and so is a batch that mixes the
    two — the path declares for a bare sequence, and where a record declares, the record
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
    if _names_a_record("annotate", directory):
        return EXIT_USAGE
    missing = plastid_reference.missing_set(directory)

    # Every sample is read and checked *before* anything is annotated: a batch that mixed two
    # organelles would otherwise have written half its runs before saying so (one run,
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
    # choice is the user's, so it is explicit.
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
