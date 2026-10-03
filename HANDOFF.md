# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phase 1 — COMPLETE** (2026-10-03). Every `TASKS.md` item is done or explicitly non-blocking:
P1-01…P1-13 delivered, D-01/D-02 done, D-03 (item (c)) documentation only. `make check` green:
**94 fast + 13 golden**, tier B included — tier A asserts every §6 number and tier B reproduces
**SSR 37 = mono 31 · di 4 · tri 2** through the vendored MISA.

**Next action: Phase 2** — genome reading & statistics (ROADMAP M2.x). It slots into the frozen API
(`docs/API_v0.1.md`): a new file in `modules/` + one registry line + a report template. Start with
**M2.1** (gene content & summary statistics).

**Waiting on people (nothing else):** real Illumina reads + the nuclear genome → Phase 6/7 only
(anh Điềm / cô Trang) · item (c) text → documentation only (cô Trang).

**Branches (all local, none pushed):** `main` = init commit · `feat/phase-1` = P1-01 ·
`feat/p1-02-core-model` = P1-02 + P1-03 · `feat/p1-04-result-schema` = P1-04.
Branches stack — consider merging the finished ones into `feat/phase-1` before it grows.

## Log

### 2026-10-03 · Claude · dev WSL

- **Did:** finished Phase 1. P1-10 complete: `core/panel.py` (4 tiers, modes A/B, quality filter,
  `pin()` with an injected fetcher that re-reads the file with our own reader) · `run/ncbi.py`
  (esearch/esummary/taxonomy efetch, 3 req/s, injectable transport) · `plastid panel refresh |
  propose | pin`. Closed P1-12 as a **dev decision** (local for phases 1–5/8, remote for 6–7) and
  added the "Who decides what" rule to CLAUDE.md so engineering calls stop being queued for the PI.
- **Validation:** `make check` (94 fast + 13 golden). Live NCBI run: 76 RefSeq genomes for
  Iridaceae; a re-pinned NC_050833.1 whose **sha256 matches the vendored golden file byte-for-byte**
  (`a946f23e…`), NC_024936.1 → `fd441727…`. Measured and recorded in PROJECT.md §9.C-20: NCBI's
  `chloroplast[filter]` misses "plastid"-worded records (71 vs 76) — including our own comparison
  genome.
- **Residual risks:** the panel's tier model relies on NCBI lineage strings (genus = last component);
  a record with an unusual lineage would land in the wrong tier — the user confirms the final list,
  and `propose` prints the tier for every candidate. `pin()` records unparseable downloads as
  `parse_error` rather than failing the batch — a caller must check it (the CLI prints UNVERIFIED).
- **Blocked:** nothing. Waiting on people only: real reads + nuclear genome (Phase 6/7), item (c).
- **Next:** Phase 2, starting with M2.1 (gene content & summary statistics).

### 2026-10-02 · Claude · dev WSL

- **Did:** created the session set (CLAUDE.md · ROADMAP.md · TASKS.md · HANDOFF.md); squashed git
  history into one `init` commit and force-pushed (amended to remove Claude attribution — repo rule
  + memory); built **P1-01** scaffold, **P1-02** core model + BioPython GenBank reader + vendored
  golden data (sha256-pinned), **P1-03** IR detection + mandatory gate + coordinate normalization.
  §10.6 decided and recorded (GenBank coords per genome · normalized for cross-genome figures).
- **Validation:** `make check` (40 fast + 7 golden, no skips); the 20-genome sweep reproduces
  PROJECT.md §8.3 — 7/9 vs GenBank, gate catches 2/2 diverged genomes with 0 false alarms; golden
  tier-A numbers asserted in tests. Recorded in PROJECT.md §8.5.
- **Did (cont.):** §10.7 decided (reports resolve table cells; an unresolvable reference raises) →
  PROJECT.md §3/§10. **P1-04 green** — `result/result.py` (`Scope` · `Column` · `Table` · `Result`,
  refs `table.column` / `table[row].column`) and `result/export.py` (csv with provenance comments,
  json with raw values, markdown); exports byte-identical across runs.
- **Residual risks:** de-novo detection is still silently short on 1/9 known cases (NC_063920 →
  22,998 bp vs the true 26,161); the annotation cross-check *records* the 3,163 bp disagreement but
  the run does not fail on it — decide the policy later. The de-novo detector is unvalidated on
  mitochondria. `csv` exports carry provenance as `#` comments — readers must strip them (stdlib
  `csv` cannot).
- **Blocked:** §10.8/§10.9 (P1-10) · D-01. Nothing blocks P1-05.
- **Next:** P1-05 test-first — run layout, manifest, `run_tool()` (the single subprocess door), cache.
