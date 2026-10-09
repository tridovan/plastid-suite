"""Step execution: cache lookup → run the module → export → store.

The runner is the only place that turns a :class:`Result` into files; modules stay pure.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from plastid.core.ir import IRResult
from plastid.core.model import Genome
from plastid.result.export import Payload, export_csv, export_json
from plastid.result.result import Result
from plastid.run.cache import Cache
from plastid.run.context import RunContext

__all__ = [
    "DATA_LABEL_PARAM",
    "StepOutcome",
    "execute_step",
    "export_files",
    "label_result",
    "write_payload",
]

# The parameter that puts a run's data provenance on every table it exports. A stand-in dataset
# never travels unlabelled (CLAUDE.md: "label every output made from it"), and the label is
# applied here — the one place that turns a Result into files — rather than in twenty modules
# that would each have to remember it.
DATA_LABEL_PARAM = "data_label"

ReportRenderer = Callable[[Result], str]
# A document takes the result *and* the already-rendered report text: the prose is produced once
# and shared, so the Word document and report.md can never disagree about a number.
DocumentRenderer = Callable[[Result, str | None], bytes]

# `run/` may not import `modules/` (layering), so the module callable is typed loosely here and
# the arity check below is what keeps the two contracts honest.
RunFn = Callable[..., Result]


@dataclass(frozen=True, slots=True)
class StepOutcome:
    status: str  # "ok" | "cached"
    key: str
    files: tuple[str, ...]  # absolute paths, deterministic order


def _params_hash(params: dict) -> str:
    """Parameters, with any parameter that names a file replaced by that file's content digest.

    A path is not an identity (D5-09 §7.2): a run told ``--set naming_reference=ref.gb`` must not
    be served a result computed against a *different* ``ref.gb``. Any string parameter that
    resolves to a file on disk therefore contributes its sha256, so editing the file changes the
    key; a parameter that is not a path is hashed as itself.
    """
    expanded = {key: _content_addressed(value) for key, value in params.items()}
    payload = json.dumps(expanded, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _content_addressed(value: object) -> object:
    if isinstance(value, str):
        path = Path(value)
        if path.is_file():
            return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        if path.is_dir():
            # A *directory* of inputs is a parameter the annotate step takes, and the same rule
            # holds: the path is not the identity. Measured 09/10/2026 — `--set
            # reference=.panel/plastid-reference` was hashed as that string, so re-pinning a
            # record, or editing one, served the previous run's tables from the cache with no
            # sign that the reference had moved. The digest is over the file names and the files'
            # own digests, so adding, removing and editing each change the key.
            listing = {
                entry.name: hashlib.sha256(entry.read_bytes()).hexdigest()
                for entry in sorted(path.iterdir())
                if entry.is_file()
            }
            return {"path": str(path), "contents": listing}
    return value


def _input_hash(genome: Genome) -> str:
    """Everything a module can read: the sequence **and** the annotation and the source qualifiers.

    Hashing the sequence alone made two *annotations* of one sequence share a cache key — measured
    05/10/2026 on the corpus's own pair (§8.31): the published and PGA records of *Iris domestica*
    carry the same 153,724 bp and different gene sets (86 vs 85 CDS), so the second run was served
    the first's tables (H-15). A module that reads genes has to be keyed on the genes; one that
    reads ``source`` (submission, for one) has to be keyed on that too.
    """
    payload = json.dumps(
        {
            "accession": genome.accession,
            "name": genome.name,
            "topology": genome.topology,
            "source": genome.source,
            "sequence": genome.sequence,
            "genes": [
                (
                    gene.name,
                    gene.kind,
                    gene.product,
                    gene.pseudo,
                    # The annotator searches the record's own /translation (D5-08): two records
                    # with one sequence and one gene list but different proteins are two
                    # different inputs, and the cache may not conflate them.
                    gene.translation,
                    tuple(
                        (part.start, part.end, part.strand, part.partial_start, part.partial_end)
                        for part in gene.parts
                    ),
                )
                for gene in genome.genes
            ],
            "repeat_regions": [
                (region.start, region.end, region.strand) for region in genome.repeat_regions
            ],
        },
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _module_source_hash(run_fn: RunFn) -> str:
    """The module's own source file, hashed — H-08, and the reason it exists.

    The key used to carry ``module_version`` and nothing else about the code, so a behaviour
    change without a version bump served a **stale payload** (measured 2026-10-04: the CDS-scope
    fix came back ``cached``, four modules had to be bumped by hand). Version bumps stay the
    convention; this is the guard.
    """
    source = inspect.getsourcefile(run_fn) or inspect.getfile(run_fn)
    return hashlib.sha256(Path(source).read_bytes()).hexdigest()


def _ir_fingerprint(ir: IRResult | None) -> str:
    if ir is None:
        return ""
    return f"{ir.length}:{ir.identity:.6f}:{ir.source}"


def label_result(result: Result, label: str | None) -> Result:
    """Prefix every table's scope label with the run's data label (no label, no change).

    Scope labels are what the CSV headers, the JSON export and the report all carry, so one
    prefix here reaches every exported number — a number from a stand-in cannot leave the run
    without saying so. The tables' ids, columns and values are untouched.
    """
    if not label:
        return result
    tables = {
        table_id: replace(table, scope=replace(table.scope, label=f"{label} · {table.scope.label}"))
        for table_id, table in result.tables.items()
    }
    return replace(result, tables=tables)


def export_files(
    result: Result,
    render_report: ReportRenderer | None = None,
    render_document: DocumentRenderer | None = None,
) -> Payload:
    """Logical name → text or bytes. Deterministic: table order, then result.json, then report."""
    files: Payload = {}
    for name, text in result.artifacts.items():
        files[name] = text
    for table in result.tables.values():
        files[f"tables/{table.id}.csv"] = export_csv(table)
    files["result.json"] = export_json(result)
    report = render_report(result) if render_report is not None else None
    if report is not None:
        files["report.md"] = report
    if render_document is not None:
        files["report.docx"] = render_document(result, report)
    return files


def write_payload(ctx: RunContext, files: Payload) -> tuple[Path, ...]:
    """Write a payload into the step's own directory — the one place a Result becomes files.

    Public because a command whose input is not a genome still produces a ``Result`` and still
    writes it the same way (``plastid reference audit``): seam 2 is about *who writes*, not about
    which command is running, and a second writer would be a second set of file semantics.
    """
    written: list[Path] = []
    for logical_with_ext, payload in files.items():
        logical, _, ext = logical_with_ext.rpartition(".")
        path = ctx.out(logical, ext)
        if isinstance(payload, bytes):
            path.write_bytes(payload)
        else:
            path.write_text(payload, encoding="utf-8")
        written.append(path)
    return tuple(written)


def _write(ctx: RunContext, files: Payload) -> tuple[str, ...]:
    return tuple(str(path) for path in write_payload(ctx, files))


def execute_step(
    *,
    module_id: str,
    module_version: str,
    run_fn: RunFn,
    genomes: Sequence[Genome],
    arity: int,
    params: dict,
    ctx: RunContext,
    cache: Cache,
    force: bool = False,
    render_report: ReportRenderer | None = None,
    render_document: DocumentRenderer | None = None,
) -> StepOutcome:
    """Run one module step, or restore it from the cache when nothing relevant changed."""
    if arity == 1:
        if len(genomes) != 1:
            raise ValueError(
                f"module {module_id!r} takes exactly one genome, got {len(genomes)} — "
                f"the v0.1 contract is run(genome, …)"
            )
    elif len(genomes) < arity:
        raise ValueError(
            f"module {module_id!r} needs at least {arity} genomes, got {len(genomes)} — "
            f"the v0.2 contract is run(genomes, …)"
        )
    key = Cache.key(
        module_id=module_id,
        module_version=module_version,
        module_source=_module_source_hash(run_fn),
        params_hash=_params_hash(params),
        input_sha256s=tuple(_input_hash(genome) for genome in genomes),
        ir_fingerprint=_ir_fingerprint(ctx.ir),
        # What the step is asked to *produce* is part of what it is: a run asking for a Word
        # document must not be served a payload that has none (H-14, measured 05/10/2026 — a
        # `--docx` run came back `cached` with 0 .docx files, and the manifest's `command`
        # did not carry the flag either).
        payload_kinds=(
            ("report" if render_report is not None else "")
            + ("+document" if render_document is not None else "")
        ),
    )
    if not force:
        cached = cache.get(key)
        if cached is not None:
            return StepOutcome(status="cached", key=key, files=_write(ctx, cached))

    result = run_fn(genomes[0], params, ctx) if arity == 1 else run_fn(genomes, params, ctx)
    if result.module != module_id:
        raise ValueError(f"module {module_id!r} returned a Result for {result.module!r}")
    result = label_result(result, params.get(DATA_LABEL_PARAM))
    payload = export_files(result, render_report, render_document)
    files = _write(ctx, payload)
    cache.put(key, payload)
    return StepOutcome(status="ok", key=key, files=files)
