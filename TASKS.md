# TASKS.md — current work and session handoff

> This file holds the **active** work only. Phase milestones and progress live in `ROADMAP.md`;
> when a phase starts, derive its tasks from that phase's milestones and put them here.
> Rules: `CLAUDE.md`. Evidence: `PROJECT.md` (only when needed).
> Statuses: `TODO` · `DOING` · `BLOCKED` · `DONE`. Owners: Dev · Claude · Diem · Trang · Hermes.
> At the end of every session — or as soon as context fills up — update **Current focus** and append
> a handoff entry. Assume the next reader is a fresh agent with zero history.

---

## Current focus

**Status: pre-implementation.** Repo holds docs only (`PROJECT.md`, `CLAUDE.md`, `ROADMAP.md`,
`TASKS.md`). No code yet.

**Open decisions** (need the team's call; details in ROADMAP.md → "Open decisions").
They do **not** block `P1-01` — only later milestones:

- Coordinate system in reports (PROJECT.md §10.6) → blocks P1-02 · P1-03
- Reports generated from `Result` objects (PROJECT.md §10.7) → blocks P1-04 · P1-07
- Reference Panel: snapshot vs live NCBI (§10.8) + item (c) (§10.9) → blocks P1-10
- Doc fixes: §12 install lines (cpstools channel · add mafft · add MISA download) → D-01

**Next action: `P1-01` — scaffold. Unblocked, ready to start.**

---

## Phase 1 tasks (derived from ROADMAP M1.x)

| ID | Task | Milestone | Status | Owner | Done when |
|---|---|---|---|---|---|
| P1-01 | Repo skeleton + packaging: pyproject (hatchling), uv env, Makefile, CLI stub | M1.1 | TODO | Claude | `uv pip install -e '.[dev]'`; `plastid --version`; `make lint` |
| P1-02 | Core model (Genome · Gene · Region) + GenBank reader (BioPython) + vendored golden data | M1.2 | TODO | Claude | unit test reads 153,724 bp · 37.85 % · 86+38+8 · `gene` 133 · 112 names |
| P1-03 | IR detection (annotation-first + de-novo port) + gate ≥ 99.9 % + normalization | M1.3 | TODO | Claude | IR 26,214 bp @ 100.00 %; the 2 known failure modes raise `IRGateError` |
| P1-04 | Result schema (`Table` · `Figure` · `Scope`) + exporters (csv · json · md) | M1.4 | TODO | Claude | cell refs resolve; two exports byte-identical |
| P1-05 | Run layer: layout + manifest + `run_tool()` + content-hash cache | M1.5 | TODO | Claude | manifest validates; re-run `cached`, checksums unchanged |
| P1-06 | Golden tests tier A + `make golden` | M1.6 | TODO | Claude | every §6 tier-A number asserted; tier B prints skip reason |
| P1-07 | Report engine (cell references) + citation registry + `defaults.yaml` | M1.7 | TODO | Claude | every `{{…}}` resolves; no bare digits outside refs |
| P1-08 | Module registry + first module `genome_stats` end-to-end | M1.8 | TODO | Claude | run dir complete: `00_run_info` + tables + report |
| P1-09 | Vendor MISA (md5 pinned) + `ssr` module (tier B) | M1.9 | TODO | Claude | SSR 37 = 31 / 4 / 2 with perl + MISA; missing tool ⇒ exit 3 naming it |
| P1-10 | Reference Panel (4 tiers) + NCBI fetch/cache | M1.10 | TODO | Claude | snapshot pinned + checksummed (§10.8 decided first) |
| P1-11 | Sample-group metadata | M1.11 | TODO | — | used by ≥ 1 comparison table |
| P1-12 | Execution backend decision (local vs remote) | M1.12 | TODO | Dev + Claude | recorded in PROJECT.md §3 |
| P1-13 | Hardening + API freeze v0.1 (README · environment.yml · Dockerfile · `docs/API_v0.1.md`) | M1.13 | TODO | Claude | `make check` on clean clone; Docker `make golden` = same numbers |

## Housekeeping

| ID | Task | Status | Owner | Done when |
|---|---|---|---|---|
| D-01 | Fix PROJECT.md §12 install lines (cpstools channel · add mafft · add MISA download) | TODO | Claude | a fresh machine can follow it end-to-end |
| D-02 | Record CGAS source review as PROJECT.md §7.9 | TODO | Claude | awaiting go-ahead |
| D-03 | Fill Reference Panel item (c) | BLOCKED | Trang | source text truncated |

---

## Handoff log  *(newest first)*

### 2026-10-02 · Claude · dev WSL

- **Did:** created the session set — `CLAUDE.md` (self-contained briefing + rules), `ROADMAP.md`
  (phases 1–9 with milestone checklists), `TASKS.md` (this file). Recorded the 02/10 decisions into
  PROJECT.md (§3 · §10 · §13.1) earlier in the session. CGAS source review done (18,353 lines) —
  lessons prepared for §7.9. Phase 1 thin-slice design drafted; **not approved, no code started**.
- **Blocked:** the 4 decisions in Current focus + doc fix D-01.
- **Next:** get those decisions, then start `P1-01`.
