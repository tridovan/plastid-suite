# TASKS.md — active tasks

> **Phase 8 is COMPLETE** on `feat/phase-8` — cross-checking, performance, automation (M8.1 · M8.3 ·
> M8.4; **M8.2 is parked with Phase 6** and appears on no board). Every task here is **plastid
> work or infrastructure** — CPStools · CGAS · CPGview · the 34-plastome corpus · CI — and it is
> the last phase whose work needs **no data of its own**. The mito round (Phase 6 entirely, M7.3,
> M8.2) is **parked by the scope decision of 10/10/2026**, not waiting on anybody; what Phase 7
> keeps is assembly, which waits on **real reads**. Phase 9 is packaging, and it wants this
> phase's numbers first.
> Milestones: `ROADMAP.md` → Phase 8 · scope text: PROJECT.md §13.8 · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 8 — Cross-checking, performance, automation (M8.1 · M8.3 · M8.4; M8.2 is parked with Phase 6)

**Goal:** only what cannot be done earlier — every analysis the toolkit promises now exists
(Phase 2B closed the last gap: `repeats`, §8.58), so a number of ours can finally be put beside a
number from another tool.
**In-phase checks:** every disagreement is **explained**, not counted — with the class definition,
the parameter or the coordinate convention that causes it · both tools' **version and argv**
recorded for every cross-check run (a cross-check with an unnamed version is a story) · every pair
of numbers compared on the **same pinned input file** · M8.3's remaining gap either measured or
written as *not measured* · CI demonstrates a failing golden number turning the build red.
*Order: the CPGview decision first (P8-01) — M8.1 names a tool §7.8 has never recorded, and the
coverage rule calls that a gap, not a footnote.*

- [x] **P8-01 · CPGview: declare it or drop it** *(M8.1's first blocker)* — **DONE 10/10/2026:
      surveyed, NOT adopted**; the §7.8 row carries the measurement, §13.8 resolves, §8.60 is the
      entry. The repo GitHub is a landing page (0 code, 0 release), there is no pip/conda/CLI, the
      web app's port does not answer (controlled test), and the only installable form is a 2.87 GB
      container for a figure tool.
  - M8.1 and §13.8 name **CPGview**; §7.8 has no source, no version, no survey line. It is a **web
    service** (a browser tool for chloroplast genome figures), so there is no binary to pin and no
    argv to record — and a record that is not published would **leave the machine** (the §10.12
    question, still undecided).
  - Decide and record: either a §7.8 line with its URL, what it does and the data-handling
    sentence, or **surveyed — not adopted** with the measurement that says why. Either way §13.8
    stops naming a tool the repo cannot account for.
  **Done when:** §7.8 carries the CPGview line (declared, or *not adopted* with its reason), and
  every tool named in §13.8 resolves to a §7.8 row or a `docs/tools/` file.

- [x] **P8-02 · CPStools cross-check** *(M8.1 — the phase's one big cross-check)* — **DONE
      10/10/2026 (§8.61)**: installed pip-only in its own env (3.0), `docs/tools/cpstools.md`
      written in the same session. Five surfaces compared on one pinned input: IR identical,
      SSR 37/37 (2 rows = the documented compound convention), RSCU counts identical on the same
      gene set, π 100/104 genes to 5 dp, Ka/Ks differences named (small-Ks denominators). The
      cross-check found **a real bug in our own `ssr` module** — the `Context` column read a
      trans-spliced gene's outer span, wrong on 11 of 37 golden rows (§9.M-1, fixed with a
      regression test written first).
  - `cpstools` 3.0 (PyPI, §7.8) covers IR · π · RSCU · SSR · KaKs — the overlapping surface. It is
    **not on this machine** (measured 10/10/2026: no import, no directory anywhere): install it
    pip-only in its own env, pin the version in `environment.yml`, and write
    `docs/tools/cpstools.md` **in the same session as the install** (§7.8's rule).
  - Compare on **one pinned input**: `NC_050833.1` for the single-genome numbers (IR · SSR ·
    RSCU), the golden pair for π. Every differing number gets its cause named.
  - The ground already measured is the model: phase 1 found **37 SSRs both ways** and one 1-bp
    disagreement that turned out to be a **compound cluster** `(TA)6 + (A)10` (§8.1/§8.2). That
    is what "differences explained" means here — and cpstools is the tool this project's SSR
    count was first checked against, so this task re-runs that check on the **current** module.
  **Done when:** a new §8 entry holds the side-by-side table (ours · cpstools · the cause of each
  difference), the tool's version and argv are recorded, and each number states its scope.

- [x] **P8-03 · CGAS cross-check** *(M8.1)* — **DONE 10/10/2026 (§8.62)**: installed in its own
      env (v1.0.1), `docs/tools/cgas.md` written. Structure, RSCU and introns match ours exactly
      on the golden; SSR differs by one **threshold** parameter (defaults 54 vs MISA 37 — with
      MISA's thresholds CGAS gives 37, identical). Two findings: CGAS reports `LSC/SSC/IR = 0`
      for `NC_024936.1` with no warning while still counting it "analysed" (our paths find
      IR 26,184 @ 100.00 %), and **our `introns` module publishes no tRNA introns** — the golden
      has 8, now named **H-23** in ROADMAP M2.7.
  - CGAS's source is in the repo's data dir (`.data/cgas-review/`, read 02/10/2026 — §7.9); its
    own environment (R + external tools, its `INSTALL.md`) is what has to be faced first. Decide
    what **can** run here and name what cannot, rather than reporting a check that never ran.
  - Where it runs: the module(s) with a deterministic output — codon usage · IR boundaries · SSR —
    against ours, on the golden. **A difference is neither automatically our bug nor theirs:**
    §7.9 lists CGAS's measured failure modes (substring IR classification, "no IR ⇒ the whole
    genome is LSC", silent success), so its labels are data to explain, never a benchmark.
  **Done when:** a new §8 entry records what ran, on which input, with both sides' parameters —
  and names, for each module that could not run, the dependency that stopped it.

- [x] **P8-04 · M8.3's named gap: `diversity` and `barcodes` at N = 34** *(M8.3)* — **DONE
      10/10/2026 (§8.63)**: `diversity` **13 min 14 s · 138 MB** (π 0.040086 over
      120,964 of 191,707 columns), `barcodes` **12 min 27 s · 138 MB** (859 windows, 13
      hypervariable regions, 10 candidates, 8 found), both `--force`, one at a time. The
      old prediction was wrong in the useful direction: they were left "not measured"
      because `substitutions scope=whole` was *stopped at a cap*, and a cap is not a
      measurement — both finish in ~13 minutes. `docs/M8.3_scale.md` now names no module
      as not measured.
  - `docs/M8.3_scale.md` measured nine single-genome modules plus `substitutions` · `synteny` ·
    `popgen` · `kaks` over the 34-record corpus, and wrote **"not measured"** for `diversity` and
    `barcodes` because both align whole genomes — a reading of the code, not a measurement. That
    honesty is the task: measure them on the same records, `--force`, `/usr/bin/time -v`.
  - The `substitutions` curve predicts they are expensive; a prediction is not an entry. If a run
    does not finish, **that is the finding**, written as what was attempted and for how long.
  **Done when:** `docs/M8.3_scale.md` carries measured wall clock and peak RSS for both modules at
  N = 34 (or the measured timeout with what was tried), and no run module is left "not measured".

- [x] **P8-05 · CI on every change** *(M8.4)* — **DONE 10/10/2026 (§8.64)**: two jobs
      (`fast`: ruff + unit suite + tier-A golden; `golden-tier-b`: the pinned tools via
      micromamba + `PTS_REQUIRE_TOOLS=1` + table2asn checked by md5). The first run was
      red on three real defects — scipy 1.18's silent `nan` exported into the FDR table, a
      test hardcoding `/root/toolbox/bin/table2asn`, and the micromamba env missing from
      `PATH` — all fixed. Proof it catches a broken golden number: run 38024933018, red on
      `37.863555723244254 == 37.8536 ± 0.001`; the revert run 38025149827 is green.
  - The repo is on GitHub (`git@github.com:tridovan/plastid-suite.git`) and has **no CI at all**
    (measured 10/10/2026: no `.github/`, no `.gitlab-ci.yml`). Add a workflow that runs on every
    push: `ruff` lint + format check, the fast unit suite, and the **pure-Python golden tier**
    (`tests/golden/test_tier_a.py` — no tools, no network).
  - The tool-dependent golden (tier B: MISA · MAFFT · repeat-match) needs the **pinned** versions:
    a job that installs the pinned conda env from `environment.yml` — never a bare
    `conda install mafft`, because §9.B-14 measured two mafft builds giving two Ts/Tv numbers from
    one genome.
  - Prove it can fail: one commit that breaks a golden number must turn the build red, recorded in
    the §8 entry. A CI that has never been red is a green light nobody has tested.
  **Done when:** a push shows the workflow running (run URL in the §8 entry), the golden job
  installs the pinned tools, and the deliberately-broken run is recorded as red before the fix
  turns it green.

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Recorded in |
|---|---|---|---|
| D-h ✅ | What CPGview is to this project | **DECIDED 10/10/2026 (P8-01) — surveyed, not adopted.** The recommendation was "checked once, by hand, on an already-public record; not adopted as a dependency"; the check itself turned out to be impossible — the web app does not answer (measured, with a control) — and the survey found no pip/conda/CLI and a 2.87 GB container as the only installable form, for a tool that produces figures rather than numbers. The record is §7.8 · §13.8 · §8.60 | §7.8 · §13.8 · §8.60 |
| D-i ✅ | What CI must run | **DECIDED 10/10/2026 (P8-05) — as recommended.** `.github/workflows/ci.yml`: `fast` (ruff · unit suite · tier-A golden) and `golden-tier-b` (micromamba with every tool pinned to §7.8's versions, `PTS_REQUIRE_TOOLS=1`, table2asn fetched and md5-checked). Watched both ways: the deliberate break turned it red on the golden number (`37.863555723244254 == 37.8536 ± 0.001`), the revert turned it green. README carries the badge | P8-05 · §8.64 · `.github/workflows/ci.yml` |
| D-j ✅ | Where `cpstools` lives | **DECIDED 10/10/2026 (P8-02) — pip into its own env, pinned, declared in `docs/tools/cpstools.md`**, used for the cross-check only. Followed as recommended: `uv pip install cpstools==3.0` into `/root/toolbox/env/cpstools`; `KaKs_Calculator` 2.0.1 (bioconda) alongside it, because cpstools' KaKs route shells out to one; `environment.yml` carries both, annotated as cross-check tools | §7.8 · `docs/tools/cpstools.md` · §8.61 |

## Pitfalls — this phase's own traps

| Trap | What to do instead |
|---|---|
| A cross-check run on a **different input** than ours ("it found 36 SSRs" — on which record?) | one pinned record, named in the §8 entry, the same file for both tools |
| The 26/9 survey's numbers ("9 matches ≥30 bp", "21 blocks", §9.H) | measure here; a survey number is a claim about other data until re-measured — Phase 2B's repeats count was the last one to fall |
| A disagreement read as our bug — or as theirs — with no cause named | name it: class definition (compound · tandem), parameter, coordinate convention (1-based vs 0-based, IR in or out). Phase 1's 1-bp SSR difference is the worked example |
| CGAS's labels taken as truth | §7.9 measured its failure modes (substring IR classification, "no IR ⇒ all LSC", silent success); its output is data to explain, never a benchmark |
| "`diversity` and `barcodes` scale like `substitutions`" | write **not measured** until measured — §6/§9's rule, which the scale doc already follows in those words |
| Timing a cached run instead of a real one (H-08: the cache keys on `module_version`, not on source) | `--force`, and say so in the entry |
| An **unpinned** tool inside CI | install from `environment.yml`; §9.B-14 is the measured reason |

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| A machine with CGAS's full environment (R + its tool chain), if it cannot be built here | anh Điềm | P8-03's second half — the modules that could not run |
| Real samples for the two manuscripts | cô Trang / anh Điềm | M2.10's final version — the stand-in draft re-runs, never gets patched |
| Real Illumina/ONT reads — the project has none; phase 7's smoke tests ran on simulated pairs | anh Điềm / cô Trang | M7.1 · M7.5 (assembly — the only live part of Phase 7) |

*The **mito round** is not in this table on purpose — it is **parked, not waiting**: Phase 6
entirely (M6.1 NUPT/MTPT · M6.2 mito structure · M6.3 plastid–mito comparison), M7.3 and M8.2 are
out of scope for this round by the decision of 10/10/2026 (CLAUDE.md → "Priority today"), and they
are marked parked in `ROADMAP.md`. The structure that makes them addable later is kept
(D5-06/D5-07); the mito code that exists stays and keeps being tested. **M6.4/M6.5** left the plan
on 10/10/2026 — divergence time and ancestral state are the PI's analyses, not modules.*
