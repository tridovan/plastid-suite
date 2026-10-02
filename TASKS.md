# TASKS.md — current work and session handoff

> This file holds the **active** work only. Phase milestones and progress live in `ROADMAP.md`;
> when a phase starts, derive its tasks from that phase's milestones and put them here.
> Rules: `CLAUDE.md`. Evidence: `PROJECT.md` (only when needed).
> Statuses: `TODO` · `DOING` · `BLOCKED` · `DONE`. Owners: Dev · Claude · Diem · Trang · Hermes.
> At the end of every session — or as soon as context fills up — update **Current focus** and append
> a handoff entry. Assume the next reader is a fresh agent with zero history.

---

## Current focus

**Phase 1 in progress.** P1-01 (scaffold) · P1-02 (model + reader + data) · P1-03 (IR + gate +
normalization) — **DONE**; `make check` green (**40 passed** fast, **7 passed** golden, no skips).

**Open decisions** (need the team's call; details in ROADMAP.md → "Open decisions"):

- Reports generated from `Result` objects (§10.7) → blocks P1-04 · P1-07 — *next in line*
- Reference Panel: snapshot vs live NCBI (§10.8) + item (c) (§10.9) → blocks P1-10
- Doc fixes: §12 install lines (cpstools channel · add mafft · add MISA download) → D-01

**Next action: `P1-04`** — Result schema + exporters (needs the §10.7 decision: reports generated
from `Result` objects with cell references, or written by hand?).

---

## Phase 1 tasks (derived from ROADMAP M1.x)

| ID | Task | Milestone | Status | Owner | Done when |
|---|---|---|---|---|---|
| P1-01 | Repo skeleton + packaging: pyproject (hatchling), uv env, Makefile, CLI stub | M1.1 | DONE | Claude | `uv sync`; `plastid --version`; `make check` green |
| P1-02 | Core model (Genome · Gene) + GenBank reader (BioPython) + vendored golden data | M1.2 | DONE | Claude | tier A tests read 153,724 bp · 37.8536 % · 86+38+8 · `gene` 133 · 112 names |
| P1-03 | IR detection (annotation-first + de-novo port) + gate ≥ 99.9 % + normalization | M1.3 | DONE | Claude | IR 26,214 bp @ 100.00 %; gate catches 2/2, 0 false alarms; §8.3 reproduced |
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

- **Did:** session working set created (`CLAUDE.md` · `ROADMAP.md` · `TASKS.md`); git history squashed
  into one `init` commit and force-pushed to `origin/main` (amended to remove Claude attribution —
  repo rule + memory); branches `feat/phase-1` (P1-01 scaffold) and `feat/p1-02-core-model`;
  **P1-01 green** (pyproject/`uv.lock`, CLI stub, package tree, import-layering guard, Makefile);
  **P1-02 green** — core model (0-based half-open internal / 1-based inclusive out), BioPython
  GenBank reader, vendored golden data (sha256-pinned); tier A asserts 153,724 bp · GC 37.8536 % ·
  86/38/8 · gene=133 · 112 unique · 21 repeated · truncated ycf1 · trans-spliced rps12.
  `make check`: 22 passed, 1 skipped (IR test → P1-03).
- **Did (cont.):** §10.6 decided (GenBank coords per genome · normalized for cross-genome figures ·
  label every table · `rotation_offset` in run info) → recorded in PROJECT.md §3/§10 and ROADMAP.
  **P1-03 green** — `core/ir.py` (annotation-first + validated de-novo port + gate) and
  `core/normalize.py`; 20-genome sweep reproduces §8.3 (7/9, 2/2 caught, 0 false alarms). Recorded
  in PROJECT.md §8.5: the validated detector is `ir_probe.py strict_best`, **not** `ir_probe_v2.detect_v2`
  (which over-extends: 27,601 vs 26,217 bp on NC_086581); NC_063921 re-measured 98.87 %, not 53.65 %.
- **Branches:** all local, not pushed. `main` = init commit; `feat/phase-1` (P1-01);
  `feat/p1-02-core-model` (P1-02 + P1-03).
- **Blocked:** §10.7 (blocks P1-04 · P1-07), §10.8/§10.9 (P1-10), D-01.
- **Next:** decide §10.7, then P1-04 (Result schema + exporters) test-first.
