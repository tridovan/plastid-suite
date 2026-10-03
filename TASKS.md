# TASKS.md — active tasks

> Board for the **current phase**. Milestone progress: `ROADMAP.md` · session log: `HANDOFF.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = needs a human decision or data.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> Finished tasks move to the **Done** table at the bottom, with the date.

## Phase 1 — framework & shared data



## Done (Phase 1)

| ID | Task | Milestone | Finished |
|---|---|---|---|
| P1-01 | Repo skeleton + packaging | M1.1 | 2026-10-02 |
| P1-02 | Core model + GenBank reader + vendored golden data | M1.2 | 2026-10-02 |
| P1-03 | IR detection + gate + coordinate normalization | M1.3 | 2026-10-02 |
| P1-04 | Result schema + exporters (refs resolve or raise; byte-identical exports) | M1.4 | 2026-10-02 |
| P1-05 | Run layer: layout · manifest · `run_tool()` · cache · runner | M1.5 | 2026-10-02 |
| P1-06 | Golden harness: tier A complete + tier B probe (`PTS_REQUIRE_TOOLS=1`) | M1.6 | 2026-10-02 |
| P1-07 | Report engine + citations + defaults (templates reference cells only) | M1.7 | 2026-10-02 |
| P1-08 | Module registry + `genome_stats` end-to-end + CLI `plastid run` | M1.8 | 2026-10-02 |
| P1-09 | MISA vendored (md5-pinned) + `ssr` module — tier B: **37 = 31 / 4 / 2** | M1.9 | 2026-10-02 |
| P1-11 | Sample-group metadata (loader; comparison-table use moved to M2.11) | M1.11 | 2026-10-02 |
| P1-10 | Reference Panel: tiers · proposal modes · pinning · NCBI client · CLI (verified live) | M1.10 | 2026-10-03 |
| P1-12 | Execution backend: local for phases 1–5, 8; remote for 6–7 (dev decision) | M1.12 | 2026-10-03 |
| P1-13 | Hardening: README · environment.yml · Dockerfile skeleton · `docs/API_v0.1.md` | M1.13 | 2026-10-02 |
| D-01 | PROJECT.md §12 fixed (mafft + perl added · cpstools → bioconda · MISA download + md5) | — | 2026-10-02 |
| D-02 | PROJECT.md §7.9 — CGAS source review recorded (10 failure modes + what we do instead) | — | 2026-10-02 |

## Housekeeping

| ID | Task | Status | Owner | Done when |
|---|---|---|---|---|
| D-01 | Fix PROJECT.md §12 install lines (cpstools channel · mafft · MISA download) | [ ] | Claude | a fresh machine can follow it end to end |
| D-02 | Record the CGAS source review as PROJECT.md §7.9 | [ ] | Claude | awaiting go-ahead |
| D-03 | Fill Reference Panel item (c) | BLOCKED | Trang | source text was truncated |
