# TASKS.md — active tasks

> Board for the **current phase**. Milestone progress: `ROADMAP.md` · session log: `HANDOFF.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = needs a human decision or data.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> Finished tasks move to the **Done** table at the bottom, with the date.

## Phase 1 — framework & shared data

### P1-10 · Reference Panel (4 tiers) + NCBI fetch/cache · M1.10 · Claude · BLOCKED

**Blocked by:** §10.8 (pinned snapshot vs live NCBI) + §10.9 (item (c) missing).

- [ ] pinned snapshot + checksums; the engine only *proposes* — the user picks the final list

### P1-11 · Sample-group metadata · M1.11 · Claude · [ ]

- [ ] shared group labels for every comparison figure/table

**Done when:** used by ≥ 1 comparison table. *(The comparison tables arrive in Phase 2 — until then
this task can only ship the loader + lookup, and the Done-when cannot be met.)*

### P1-12 · Execution backend decision · M1.12 · Dev + Claude · BLOCKED

**Blocked by:** the team's call — run heavy steps locally or push to a remote server?
Evidence to decide with: PROJECT.md §9.14 (8 GB machines: 1.5–2.5 GB free per job), §8.4.

### P1-13 · Hardening + API freeze v0.1 · M1.13 · Claude · [ ]

- [ ] README quickstart · `environment.yml` · Dockerfile skeleton · `docs/API_v0.1.md`

**Done when:** `make check` passes on a clean clone; Docker `make golden` gives identical numbers.

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

## Housekeeping

| ID | Task | Status | Owner | Done when |
|---|---|---|---|---|
| D-01 | Fix PROJECT.md §12 install lines (cpstools channel · mafft · MISA download) | [ ] | Claude | a fresh machine can follow it end to end |
| D-02 | Record the CGAS source review as PROJECT.md §7.9 | [ ] | Claude | awaiting go-ahead |
| D-03 | Fill Reference Panel item (c) | BLOCKED | Trang | source text was truncated |
