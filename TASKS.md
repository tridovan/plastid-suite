# TASKS.md — active tasks

> Board for the **current phase**. Milestone progress: `ROADMAP.md` · session log: `HANDOFF.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = needs a human decision or data.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**

## Phase 1 — framework & shared data

### P1-05 · Run layer: layout · manifest · run_tool · cache · M1.5 · Claude · [ ]

- [ ] `run/layout.py` — `<dataset>/<run_id>/`, `00_run_info/`, numbered module dirs;
      `dataset_slug` strips diacritics (PROJECT.md §9.16)
- [ ] `run/manifest.py` — `schema_version` + `validate()` (no external schema library)
- [ ] `run/tools.py` — the single door: `run_tool()` records name · version · path · argv · time ·
      exit code · stdout/stderr logs; version per kind (binary flag · perl `$^V` · pinned vendored)
- [ ] `run/cache.py` — content-hash key; a cache hit re-exports and marks `cached`

**Done when:** the manifest validates; an identical second run reports `cached` with unchanged file
checksums; deleting one output file and re-running restores it byte-identically.

### P1-06 · Golden harness: tier A complete + tier B skeleton · M1.6 · Claude · [ ]

- [ ] `make golden` runs the whole golden set in one command
- [ ] tier B tool probe (perl · `misa.pl` · mafft); a missing tool prints a reason and skips
- [ ] `PTS_REQUIRE_TOOLS=1` turns tier-B skips into failures (for tool-equipped machines)

**Done when:** every §6 tier-A number is asserted; tier B prints one skip line per missing tool.

### P1-07 · Report engine + citations + defaults · M1.7 · Claude · [ ]

- [ ] `report/engine.py` — `{{table.column}}` / `{{table[row].column}}` resolution; unknown ref raises
- [ ] template scan test: every `{{…}}` resolves; no bare digits outside references
- [ ] `report/citations.py` + `data/citations.yaml`; `data/defaults.yaml`
- [ ] every generated report ends with the draft banner ("Draft — needs human review")

**Done when:** the golden `genome_stats` report renders with every number resolving to a table cell.

### P1-08 · Module registry + `genome_stats` end-to-end · M1.8 · Claude · [ ]

- [ ] `modules/registry.py` — decorator registry with explicit imports (entry points deferred)
- [ ] `modules/genome_stats.py` — pure `run(genome, params, ctx) -> Result`
- [ ] CLI `plastid run --module genome_stats …` (argparse, `main(argv) -> int`, exit codes 0/1/2/3)

**Done when:** one command produces a complete run dir (`00_run_info` + tables + report) and the run
report lists tool name · version · parameters.

### P1-09 · MISA vendored + `ssr` module (tier B) · M1.9 · Claude · [ ]

- [ ] download the IPK zip, verify md5 `cf9d3579…`, vendor `misa.pl` + `misa.ini` + `SOURCE.json`
- [ ] `external/misa.py` — run in an isolated workdir (MISA writes next to its input); parse `.misa`
- [ ] `modules/ssr.py` — by region and motif, scope-labelled

**Done when:** with perl + MISA → SSR 37 = mono 31 · di 4 · tri 2; a missing tool → exit 3 naming it.

### P1-10 · Reference Panel (4 tiers) + NCBI fetch/cache · M1.10 · Claude · BLOCKED

**Blocked by:** §10.8 (pinned snapshot vs live NCBI) + §10.9 (item (c) missing).

- [ ] pinned snapshot + checksums; the engine only *proposes* — the user picks the final list

### P1-11 · Sample-group metadata · M1.11 · Claude · [ ]

- [ ] shared group labels for every comparison figure/table

**Done when:** used by ≥ 1 comparison table.

### P1-12 · Execution backend decision · M1.12 · Dev + Claude · [ ]

- [ ] decide local vs remote for the heavy steps; record the decision in PROJECT.md §3

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

## Housekeeping

| ID | Task | Status | Owner | Done when |
|---|---|---|---|---|
| D-01 | Fix PROJECT.md §12 install lines (cpstools channel · mafft · MISA download) | [ ] | Claude | a fresh machine can follow it end to end |
| D-02 | Record the CGAS source review as PROJECT.md §7.9 | [ ] | Claude | awaiting go-ahead |
| D-03 | Fill Reference Panel item (c) | BLOCKED | Trang | source text was truncated |
