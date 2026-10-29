# Development handoff

**State:** 29/10/2026, `main` — phase M11 is merged and acceptance-tested on this tree (host run
set and the full golden). Product documentation, requirements and shared
method live in this repository. Maintainer planning — the active task, the milestones and their
Definition of Done, and the dated evidence records — lives in `.dev/`, which is not part of the
product repository.

## What the product does now

`plastid` reads plastid GenBank or FASTA records and runs the analyses in
[the pipeline](docs/pipeline.md): `structure`, `genome_stats`, `gene_table`, `gene_names`,
`gene_normalize`, `orf`, `codon_usage`, `amino_acid`, `introns`, `ssr`, `repeats`, `substitutions`,
`diversity`, `barcodes`, `kaks`, `popgen`, `stats`, `synteny`, `phylogeny`, `submission`, plus the
reference path `panel`/`reference`/`annotate` and the export path `figures`/`package`/`qc`. The
mito modules are a frozen prototype: keep them passing, do not extend them.

Reproduced defects are in [known issues](docs/known_issues.md). A documentation change does not
fix one.

## Where the work is

- Open work and its order: `.dev/TASKS.md`.
- Milestones and Definition of Done: `.dev/ROADMAP.md`.
- Per-session evidence with the numbers and how they were measured: `.dev/records/`.

A clone without `.dev/` still runs everything; use the user's task, this file and the public
pipeline. Do not invent a replacement task list.

## Verifying a change

- Fast suite: `make test` (no external tool, no network).
- Tool-carrying suite: `make golden`. Missing tools are declared in `environment.yml` and
  `docs/tools/`; a missing one skips with a printed reason, and `PTS_REQUIRE_TOOLS=1` turns that
  skip into a failure. Report which checks actually ran — a skip is not acceptance evidence.
- Lint and format: `make lint`, `make fmt`.

## Environment notes that cost time to find

- `PTS_IQTREE_DIR` points at the directory holding `iqtree2`, not at the shared bin directory.
- MACSE is installed from bioconda (`macse=2.07`) and resolved through `PTS_MACSE_JAR` or `PATH`;
  the bioconda wrapper needs a JDK, so the jar path is the reliable one.
- A run longer than a few minutes goes detached (`nohup … &`) so an interrupted shell does not
  take it down.
