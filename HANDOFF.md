# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phase 1 — in progress.** P1-01…P1-09 DONE (all except the blocked ones). `make check` green
(75 fast + 13 golden, tier-B included). **Both golden tiers now run on this machine** — tier A
asserts every §6 number and tier B reproduces **SSR 37 = mono 31 · di 4 · tri 2** via the vendored
MISA (perl 5.38.2).

**Next action:** `P1-13` (hardening: README · environment.yml · Dockerfile · API freeze) and
`P1-11` (sample-group loader — its full Done-when needs Phase 2 comparison tables).

**Blocked on:** §10.8 · §10.9 (Reference Panel → blocks P1-10) · P1-12 (your call: local vs remote
for heavy steps) · D-01/D-02 (doc fixes).

**Branches (all local, none pushed):** `main` = init commit · `feat/phase-1` = P1-01 ·
`feat/p1-02-core-model` = P1-02 + P1-03 · `feat/p1-04-result-schema` = P1-04.
Branches stack — consider merging the finished ones into `feat/phase-1` before it grows.

## Log

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
