# Development handoff

**State:** 01/11/2026 — `main` holds phases M11, M4.5, the assembly of Phase 7 (M7.1, M7.2, M7.5,
M7.6: short reads and HiFi reads, the graph read, the memory limited) and the review of the command
line (M9.8), merged and acceptance-tested (`make acceptance` 291 passed, 0 skipped). The open work is
the local interface on branch `feat/m9.7-ui` (M9.7: Docker and acceptance remain), then the release
(M9.6). RNA editing is the module `rna_editing` (`docs/tools/prep_cp.md`); its reference B′'s seven
records live in `.panel/rna_editing/`, pinned in `src/plastid/data/rna_editing/taxa.yaml`. Product
documentation,
requirements and shared method live in this repository. Maintainer planning — the active task, the
milestones and their Definition of Done, and the dated evidence records — lives in `.dev/`, which
is not part of the product repository.

## What the product does now

`plastid` reads plastid GenBank or FASTA records and runs the analyses in
[the pipeline](docs/pipeline.md): `structure`, `genome_stats`, `gene_table`, `gene_names`,
`gene_normalize`, `orf`, `codon_usage`, `amino_acid`, `introns`, `ssr`, `repeats`, `substitutions`,
`diversity`, `barcodes`, `kaks`, `popgen`, `stats`, `synteny`, `phylogeny`, `submission`,
`rna_editing` (predicted C-to-U sites, PREP-Cp-compatible), plus the reference path
`panel`/`reference`/`annotate`, `editing-reference` (rebuild or verify the RNA-editing reference) and
the export path `figures`/`package`/`qc`. `plastid assemble run|prepare|import` takes paired short
reads or HiFi reads to a checked plastome through a local or an imported run
([short reads](docs/tools/getorganelle.md), [HiFi](docs/tools/oatk.md)). `plastid modules`, `tools`
and `runs` list what is there, in text or JSON.
The mito modules are a frozen prototype: keep them passing, do not extend them.

Reproduced defects are in [known issues](docs/known_issues.md). A documentation change does not
fix one.

## Where the work is

- Open work and its order: `.dev/TASKS.md`.
- Milestones and Definition of Done: `.dev/ROADMAP.md`.
- Per-session evidence with the numbers and how they were measured: `.dev/records/`.

A clone without `.dev/` still runs everything; use the user's task, this file and the public
pipeline. Do not invent a replacement task list.

## Verifying a change

- Fast suite: `make test` (no external tool, no network, no `.panel/`).
- Golden: `make golden`, the golden-pair numbers and the pinned-tool checks, the same set CI runs.
  Missing tools are declared in `environment.yml` and `docs/tools/`; a missing one skips with a
  printed reason, and `PTS_REQUIRE_TOOLS=1` turns that skip into a failure.
- Panel acceptance: `make acceptance`, the tests marked `panel` that read the gitignored `.panel/`.
  Run it when closing a milestone; it fails instead of skipping. Report which checks actually ran
  — a skip is not acceptance evidence. What each check proves: `tests/golden/README.md`.
- Lint and format: `make lint`, `make fmt`.

## Environment notes that cost time to find

- This machine hides what a new user meets: it has `.panel/`, the extras and the tools. Run a check of
  the install from a clean copy (`docs/release.md`, "A new user's checkout"), not from this tree.
- The Docker image was built on 01/11/2026 through `docker/plastid-ui.bat` (6.5 GB, about 15 minutes; a
  rebuild after a source change takes about 5, because the tool layers stay cached). Check the disk
  first, and the release revision must be built again.
- Docker Desktop mounts a Windows folder writable for the container only when Windows made it. A
  folder made from WSL as root carries root's mode, and the container's user cannot write there.
- The `PTS_*` variables are read by the tests only. The product finds a tool on `PATH` or through
  `--set`, and `plastid tools` says what it found.
- MACSE and IQ-TREE sit in the same conda environment as the other analysis tools, as
  `environment.yml` declares, so `PATH` alone finds them. MACSE's wrapper uses the JDK of its own
  environment; a jar kept elsewhere takes `PTS_MACSE_JAR` in the tests and `--set macse=` in a run.
- A run longer than a few minutes goes detached (`nohup … &`) so an interrupted shell does not
  take it down.
- The assembly tools are in their own environment (`environment-assembly.yml`). The panel test finds
  them through `PTS_ASSEMBLY_DIR`, the GetOrganelle database through `PTS_GETORGANELLE_DB` and Oatk's
  through `PTS_OATK_DB`.
- GetOrganelle exits 0 when it fails. Only the `ERROR` lines of its log say so.
- Its peak memory is read recruitment, and that follows the words it holds, not the input's size. The
  step limits the words from `memory_gb`; on an 8 GB machine two samples that took 18–20 GB in the
  published runs fit in under 6 GB, and one of them stays open for it.
- `.panel/assembly/jin2020/` is what the authors of the published comparison deposited, version 1.6.2.
  `tools/fetch_assembly_panel.py --graphs-only` rebuilds it from a pinned commit.
- ENA answers parallel range requests with 403 or cuts them short. `tools/fetch_assembly_panel.py`
  fetches one file at a time and verifies each.
