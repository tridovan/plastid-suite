# API v0.2 — the multi-genome contract

Extends `docs/API_v0.1.md`. Nothing in v0.1 changes: an arity-one module is untouched.

## What was decided (D2-01, 2026-10-03)

π, SNP and later Ka/Ks take more than one genome. The v0.1 contract `run(genome, params, ctx)`
could not express that, and widening it to `run(genomes, …)` for *every* module would hand each
single-genome module a sequence it must index at `[0]` — a contract wider than the module needs,
which code seam 5 (small explicit contracts) exists to prevent.

So the two contracts stay separate and a module declares which one it follows:

```python
@module("diversity", version="0.1.0", title="Nucleotide diversity (π)", requires_ir=True, arity=2)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result: ...
```

- `arity=1` (the default) — `run(genome, params, ctx) -> Result`, exactly as in v0.1.
- `arity≥2` — `run(genomes, params, ctx) -> Result`, with **at least** `arity` genomes; more is
  allowed because a diversity estimate over N accessions wants all of them.

The runner enforces the count before anything runs: a mismatch raises, and the CLI turns it into
a usage error (exit 2). A module never receives genomes it does not use.

## The command line

`--input` is repeatable, and the order is meaningful (row 0 is the reference):

```
plastid run --module substitutions \
  --input ref.gb --input other.gb --out .runs --dataset iris-pair
```

`arity=1` modules still take exactly one `--input`; passing two is a usage error.

## The IR gate under arity ≥ 2

Every genome must pass the gate, not only the first: one diverged IR corrupts every position-based
comparison that follows. The manifest records the primary genome's IR at `ir` and, when more than
one genome was given, every genome's measured length, identity and source at `ir.genomes`. A
failing gate anywhere stops the run with exit 2 and no module output.

## The alignment

The comparison modules align through MAFFT (`src/plastid/external/mafft.py`, pinned 7.526 per
PROJECT.md §7.8) and never fall back to a hand-rolled aligner — a measured fallback once produced
Ts/Tv 2.667 where MAFFT gave 2.686 (PROJECT.md §9.B). A missing aligner is a named `ToolError`
and exit 3.

All alignment maths lives in `src/plastid/core/align.py` as pure functions, so a number in a
manuscript can be recomputed from the stored alignment without rerunning the aligner.

## Additive changes (D2-04 / P2-09, 2026-10-03)

None of these change an existing contract; they widen what a run can emit.

### A payload may contain bytes

`export_files()` returns `Payload = dict[str, str | bytes]` (declared in `result/export.py`).
Text exports stay text; a `.docx` is a zip and travels as `bytes`. `runner._write()` branches on
the type, and `Cache` stores bytes as base64 under a separate `binary` key, so a cache hit
restores the exact file — including byte-for-byte a Word document.

### `Table.emphasis`

```python
Table(..., key="group", emphasis=("Total",))
```

Key values to render emphasised. Presentation only: it changes no value and no export other than
the Word document's font weight. It requires a key column, and a value that is not one of the
key's values raises at construction — a typo here is invisible in every other export.

### `--docx`

```
plastid run -m gene_table -i genome.gb --docx
```

Writes `tables/*.csv`, `result.json` and `report.md` as before, plus `report.docx` in the same
module directory — explicitly requested, never emitted implicitly. The document carries a title,
a provenance grid (dataset, run id, every input with its sha256, coordinate frame, module, tool
version), the rendered report prose, one captioned table per result table, the module's notes
and the cited references.

Two properties are contractual:

- **Reproducible.** The same result and metadata produce the same bytes. `python-docx` stamps the
  wall clock into the zip, so the package is rebuilt with a fixed timestamp, sorted entries and a
  pinned `create_system`.
- **No number is typed.** Every value on the page is a table cell or a labelled provenance string.
  A number that is not in a table cannot reach the document.

### A run may hold several modules

`-m` is repeatable, and the run — not the module — is the deliverable:

```
plastid run -m genome_stats -m structure -m gene_table -i genome.gb --docx --dataset <label>
```

Each module keeps its own numbered directory (`01_genome_stats`, `02_structure`, …), the manifest
lists them in order, and naming one module twice is a usage error (it would write one directory
twice). The IR gate runs once for the run and applies to every module that needs it.

With `--docx` a multi-module run also writes **`report/report.docx`** — the modules' documents
joined into one. It is assembled by joining rather than re-rendering from the results, because a
cache hit restores *files* and not results: re-rendering would silently produce an empty document
on the second run of an unchanged input. `merge_documents()` keeps every document's page setup and
section breaks (each `sectPr` is re-attached as an inline section break), unions the styles, and
takes the package's plumbing from the first document. A part list that differs raises — that means
the document was not built by this exporter.

The merged document's title and footer label are set from the run, so page 30 of a merged report
does not still claim to be reading the first module.

### The dependency is optional

`python-docx` lives in the `docx` extra (`pip install 'plastid-suite[docx]'`) and is imported
inside the function that needs it. The minimum install stays pure Python: `biopython` + `PyYAML`.
A missing extra is `ExportDependencyError` and exit 3 — a named failure, never an `ImportError`
traceback. `lxml` arrives with `python-docx` and is used for the package pruning.
