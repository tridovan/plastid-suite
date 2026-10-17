# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

## The R figure engine (M10.2)

**Goal.** The seven R scripts that drew the manuscript's Figures 3–9 become a declared external tool
the toolkit can call, with the scripts it runs carried in the package and the R version pinned.

**Measured before any code (17/10/2026, this machine).** All seven scripts run under **R 4.3.3** and
write all eight figures — the author verified them on 4.1.2, so the engine is not version-fragile.
Two things about the bundle are already known and must not be re-discovered:

- `khoi_tao()` ends with `setwd(OUT_DIR)`, so `--inputdir`/`--input` must be passed **absolute**; a
  relative path makes every script fail to find its own data.
- `00_ham_chung.R` reads `--inputdir`, `--outdir` and per-script flags; the two-panel layout is drawn
  with `grid` because `cowplot` is deliberately not required.

**Steps**
1. `src/plastid/external/rscript.py`: find `Rscript`, read its version, check the nine packages by
   name before anything is attempted, run a script in the run's own directory through `ctx.tool`.
2. `docs/tools/r.md`: the declaration — what R is used for (plotting only: none of these scripts
   computes a number the toolkit does not already hold in a table cell), the install lines, the
   version measured here, and the citation.
3. The seven scripts into `src/plastid/data/r_scripts/`, read through `core.resources.data_file` —
   the same place `misa.ini` lives. They are the team's own code, so they may be carried; the fix
   for the relative-path bug goes in with them.
4. The R pins in `environment.yml` (`r-base` plus the nine packages), and the licence inventory row.
5. `tests/unit/test_rscript.py`: the loud paths without R (missing interpreter and missing package
   are both **named**), and a `tier_b` test that runs all seven scripts over the bundle's sample data
   and counts eight figures.

**Done when:** `make check` is green; with `Rscript` on `PATH` the `tier_b` test runs the seven
scripts and finds eight figures; without it the test skips naming R, and a missing R package is
reported by name before any figure is attempted.

**Traps**
- R is a plotting engine here, never a source of numbers: if a script and its table disagree, the
  table is right.
- The bundle's sample data are the paper's, not PTS's — label anything drawn from them.
- A tool that exists only on this machine is not part of the toolkit: the declaration and the
  `environment.yml` pins ship in the same commit as the first use.

## After it

**M10.3–M10.5** (`ROADMAP.md` → Phase 10): the three figures PTS lacks (amino acid, substitutions,
Ka/Ks), the long→wide pivots for RSCU/SSR/π, and the tree — which stays the user's own script over
the user's own tree.

**H-25** is its own task with its own regression test: the IR gate refuses a record whose two
annotated IR copies differ in length, although the repeat itself is exact.
