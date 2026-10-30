# CLAUDE.md

Rules for changing the code in this repository.

**PlastidToolSuite (PTS)** analyses plastid (chloroplast) genomes from GenBank or FASTA records. It
is academic, MIT-licensed software for researchers on ordinary Linux machines (Windows runs it
through Docker). The mito modules (`mito_annotate`, `mito_genes`) are a frozen prototype: keep them
passing their tests, but do not extend them.

## Commands

- `make test`: fast tests. Needs no external tool, no network and no `.panel/` (~5 min).
- `make golden`: the golden-pair numbers and the pinned-tool checks, the set CI runs (~6 min). A
  test whose tool is missing skips with a printed reason; `PTS_REQUIRE_TOOLS=1` turns that skip
  into a failure. The tools are in `environment.yml` and `docs/tools/`.
- `make acceptance`: every test marked `panel`, which reads the private `.panel/` data that CI never
  has (~40 min). Run it when closing a milestone; a missing tool or panel fails it.
- `make lint` and `make fmt`: ruff. `make check` runs lint, test and golden (~11 min).

## Context and scope

- Start with `HANDOFF.md`, then `docs/pipeline.md` for the canonical product pipeline.
- For maintainer development in this workspace, if `.dev/` is present, read `.dev/TASKS.md` for
  the active task, then its named milestones and Definition of Done in `.dev/ROADMAP.md`. Read
  only design records that task cites. Do not infer the next task from archive documents.
- `.dev/` is optional private planning, excluded from the product repository. In a clone without
  it, use the user's task and public pipeline/known issues; do not invent a replacement TASKS file
  or require private documents to use the toolkit.
- `HANDOFF.md` is the sole current development handoff, not a user guide. Historical session notes
  are archival evidence, not current requirements. Do not load them without a task-specific need.
- `docs/known_issues.md` lists reproduced runtime defects. A documentation change does not fix one.
- Keep shared methods, provenance and requirements in the product repository. Keep personal
  session logs and maintainer planning out of product instructions.
- Do not generalize a manuscript's species list, outgroup or plotting choices into toolkit defaults.

## Architecture

`tests/unit/test_architecture.py` enforces the import layering:

- `core/`: the data model and the algorithms.
- `result/`: tables, figures, exporters.
- `run/`: run layout, manifest, cache, `run_tool()`.
- `external/`: one adapter per tool.
- `modules/`: the analyses.
- `cli.py`: the only place that exits.

A module is a pure function `run(genome, params, ctx) -> Result` (`run(genomes, …)` when its arity is
above 1). Each module is registered in `modules/registry.py` with its group: shared, plastid or mito.
A module never names a file path; it writes through `ctx.out(...)`.

Coordinates are 0-based half-open inside the code and 1-based inclusive in GenBank and in exported
tables. Convert only with `Interval.to_1based` / `Interval.from_1based`.

## Invariants

- **IR gate.** No position-based analysis runs until `core/ir.py` has passed the two
  inverted-repeat copies as ≥ 99.9 % identical. Never bypass, weaken or make it optional.
- **External tools.** Every external tool runs through `run_tool()`, which records name · version
  · argv · exit code. A tool a module uses is declared in the repo: in `environment.yml` with a
  pinned version, or in `docs/tools/<tool>.md`. A tool that exists on only one machine is not part
  of the toolkit.
- **Scope labels.** Every exported number carries a scope label. A value that was not measured is
  written "not measured". Every output made from third-party data carries its provenance label
  (`--label`).
- **Run directories.** Output goes to `<out>/<dataset>/<run_id>/`, and nothing is overwritten in
  place.
- **Minimum install.** biopython + PyYAML alone run every pure-Python analysis. Optional extras
  (docx, figures, stats, excel) are imported inside functions.
- **Gene count.** The count is CDS + tRNA + rRNA. `gene` features are a diagnostic only.

## Code rules

- **Fail loud.** Use specific exceptions with actionable messages. Never fall back silently, and
  never let a gate degrade into a warning.
- **A disagreement has two possible causes: the tool or the data.** Do not assume a published
  record is right. Before blaming the code, test the record's own claim; before accepting the
  record, show why the code is right (an arithmetic construction, an independent measurement, or
  corroborating annotation such as the gene pairs a real IR carries). "Data is wrong" and "tool is
  wrong" both need evidence, and an unexplained difference stays unexplained rather than being
  assigned to either side.
- **Banned:** `sys.exit` outside `cli.py`; `shell=True`; hand-written FASTA/GenBank parsers (use
  BioPython); hardcoded absolute paths; monkey-patching; hidden globals or cwd.
- **Extend by adding.** A new analysis is a new file in `modules/` plus one registry line. Apply DRY
  after the second real copy, never before.
- **Comments.** State units, coordinate conventions and the reason a rule exists. Do not narrate
  history, and do not point to documents outside this repository.
- Report templates reference table cells (`{{table.column}}`), never literal numbers.

## Tests

- A bug fix starts with a failing regression test and an independently justified expected result.
  Never weaken, skip or delete a test to make it pass. Do not close a defect while its regression
  remains xfail. If a golden value changes, explain it with independent evidence and obtain user
  review before accepting a new scientific baseline; do not force an old incorrect value to persist.
- The golden pair (*Iris domestica* NC_050833.1 vs *Iris gatesii* NC_024936.1, `tests/golden/`) is
  one pinned input. It is not evidence about any other genome, so a change must also be tested on the
  other record shapes in `tests/`.
- No test touches the network. Figures are tested through the data they are fed, not their pixels.
- A test that reads `.panel/` carries the `panel` mark. What each kind of golden check proves is in
  `tests/golden/README.md`.
- Report which checks actually ran and every missing tool/panel. A skip is not acceptance evidence.

Commits use the conventional style (`feat:`, `fix:`, `test:`, `docs:`, `chore:`, `refactor:`), in
English.
