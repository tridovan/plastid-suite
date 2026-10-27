# Development handoff

This file is written for one coding session and cleared when that session closes.

## Working tree, uncommitted, on `feat/m11-ir-rework`

Four batches are in the tree. Full fast suite on the final state (tool env,
`PTS_REQUIRE_TOOLS=1`): **1,069 passed, 1 skipped, 66 deselected** (557 s). `ruff check .` and
`ruff format --check .` clean. No tier-B golden has been run on this batch; by the user's rule it
runs once, at the acceptance group (Đợt 4).

1. **Reference-model guard** — `core/caller.py` chooses a rare model only when it is not a 5'-over-
   long annotation that a majority model extends; plastid CDS only. Measured on 39 Crataegus
   records: CDS exact 77.0 % → 84.5 %, ≤1 bp 81.2 % → 88.6 %. Golden annotate 10 passed, pinned
   differences unchanged. Evidence: `.dev/records/reference_model_guard_2026-10-27.md`.
2. **M11.16** — versioned artifact contract (`run/contract.py`, manifest schema 3 with
   `artifacts`/`reference_sets`), the annotate reference-set digest, downstream modules recorded
   `skipped` when an earlier step fails, IR `scoring_policy` + version, a named refusal to
   re-annotate an already-annotated GenBank (`--set reannotate=yes` makes it explicit), and a
   file-routing proof (annotate → `run -m genome_stats` on the artifact, same sha256).
3. **M11.17** — fixed the substitutions raw-vs-normalized region frame, the barcodes raw-vs-
   normalized alignment frame, barcodes pseudogene loci, diversity's raw names/pseudogenes,
   stats' unexported assumption and empty-CDS GC 0.0, the `n_per_group` crash for a one-group
   locus, zero-row crashes in `ssr`/`barcodes`/`diversity`/`genome_stats`, popgen's Gst estimator
   column, and `gene_normalize`'s verify decision.
4. Evidence for 2–3: `.dev/records/M11.16_M11.17_baseline_2026-10-27.md`; plan and progress in
   `.dev/TASKS.md`.

## Đợt 2 (started 27/10/2026): M11.6 batch

- `plastid batch` (`src/plastid/cli.py`, `src/plastid/run/batch.py`): one run per sample via the
  existing `_execute_run`, then `samples.csv` + `aggregate.csv` + `batch.json`; arity-1 only;
  `measured`/`not_measured`; duplicate labels keep `#2`; a failed sample is a status, not a zero.
  `tests/unit/test_cli_batch.py` 5 passed; the three `_execute_run` callers still pass (84 tests).
  Evidence: `.dev/records/M11.6_batch_2026-10-27.md`.
- Slice 2: `read_batch_results` + `batch_frames` (F3/F4/F5 wide from the aggregate),
  `plastid figures <batch>` draws them, and the CLI forwards species-order/exclude-species/width/
  height/dpi. Tests: `test_batch_figures.py` 5, `test_rfigure.py` 16, targeted batch 105 passed.
- Slice 3: panel batch ran — NCBI 75/80 (5 IR-gate, no cells), PGA 30/30, synthetic 27/30
  (3 IR-gate; 216 `not measured`). `batch --figures` and the normalized span split are in.
  Evidence: `.dev/records/M11.6_panel_baseline_2026-10-27.md`.
- Remaining: per-family cell cross-check of batch vs single runs, and the CGAS/paper definition
  comparison before any cell-for-cell claim.

## Still open

- Đợt 1 remainder: the local-first tool-adapter contract (every tool through `run_tool()`;
  assembly is the only local/external split; PREP-Cp is a local adapter for M4.5) and the
  80-panel/30-PGA module cases. Partial-CDS policy and the SSR threshold source are closed.
  (`plastid run -m plastid_annotate` still re-annotates an annotated GenBank — an explicit module
  selection, so the CLI guard does not fire there.)
- Đợt 2: M11.6 (multi-genome aggregation) + the remaining M11.7 (multi-sample figures).
- Then M11.18, then M11.19 → M11.9 → M11.8 + M9.6, where the full golden runs once.
- `rps16` on `MW201730.1` remains 267 bp vs the published 270: an internal junction difference,
  outside the length guard. M4.5 (RNA editing) is independent and still open.
