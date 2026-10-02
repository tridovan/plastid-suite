# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phase 1 — in progress.** P1-01 · P1-02 · P1-03 DONE. `make check` green (40 fast + 7 golden).

**Next action: `P1-04`** — Result schema + exporters (checklist + Done-when in `TASKS.md`).

**Blocked on:** the §10.7 decision — are reports generated from `Result` objects with cell
references? (PROJECT.md §10.7). Other open decisions: §10.8 · §10.9 (Reference Panel → blocks
P1-10) · D-01 (doc fixes).

**Branches (all local, none pushed):** `main` = init commit · `feat/phase-1` = P1-01 ·
`feat/p1-02-core-model` = P1-02 + P1-03.

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
- **Residual risks:** de-novo detection is still silently short on 1/9 known cases (NC_063920 →
  22,998 bp vs the true 26,161); the annotation cross-check *records* the 3,163 bp disagreement but
  the run does not fail on it — decide the policy later. The de-novo detector is unvalidated on
  mitochondria.
- **Blocked:** §10.7 (blocks P1-04 · P1-07) · §10.8/§10.9 (P1-10) · D-01.
- **Next:** decide §10.7, then P1-04 test-first (write the failing checks before the schema).
