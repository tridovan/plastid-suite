# TASKS.md — active tasks

> **Phase 8 is OPEN** on `feat/phase-8` — cross-checking, performance, automation (M8.1 · M8.3 ·
> M8.4; **M8.2 is blocked** and says so below). It is the last phase whose work needs **no data of
> its own**: Phase 6 and 7 come later for three different reasons — a nuclear genome (M6.1, and
> M8.2 through it) · several mitogenomes of one species (M6.2) · real reads (M7.1 · M7.5) — while
> **M6.3 is deferred by decision D5-07, not by data** (the waiting table below names each). Phase 9
> is packaging, and it wants this phase's numbers first.
> Milestones: `ROADMAP.md` → Phase 8 · scope text: PROJECT.md §13.8 · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 8 — Cross-checking, performance, automation (M8.1 · M8.2 · M8.3 · M8.4)

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

- [ ] **P8-01 · CPGview: declare it or drop it** *(M8.1's first blocker)*
  - M8.1 and §13.8 name **CPGview**; §7.8 has no source, no version, no survey line. It is a **web
    service** (a browser tool for chloroplast genome figures), so there is no binary to pin and no
    argv to record — and a record that is not published would **leave the machine** (the §10.12
    question, still undecided).
  - Decide and record: either a §7.8 line with its URL, what it does and the data-handling
    sentence, or **surveyed — not adopted** with the measurement that says why. Either way §13.8
    stops naming a tool the repo cannot account for.
  **Done when:** §7.8 carries the CPGview line (declared, or *not adopted* with its reason), and
  every tool named in §13.8 resolves to a §7.8 row or a `docs/tools/` file.

- [ ] **P8-02 · CPStools cross-check** *(M8.1 — the phase's one big cross-check)*
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

- [ ] **P8-03 · CGAS cross-check** *(M8.1)*
  - CGAS's source is in the repo's data dir (`.data/cgas-review/`, read 02/10/2026 — §7.9); its
    own environment (R + external tools, its `INSTALL.md`) is what has to be faced first. Decide
    what **can** run here and name what cannot, rather than reporting a check that never ran.
  - Where it runs: the module(s) with a deterministic output — codon usage · IR boundaries · SSR —
    against ours, on the golden. **A difference is neither automatically our bug nor theirs:**
    §7.9 lists CGAS's measured failure modes (substring IR classification, "no IR ⇒ the whole
    genome is LSC", silent success), so its labels are data to explain, never a benchmark.
  **Done when:** a new §8 entry records what ran, on which input, with both sides' parameters —
  and names, for each module that could not run, the dependency that stopped it.

- [ ] **P8-04 · M8.3's named gap: `diversity` and `barcodes` at N = 34** *(M8.3)*
  - `docs/M8.3_scale.md` measured nine single-genome modules plus `substitutions` · `synteny` ·
    `popgen` · `kaks` over the 34-record corpus, and wrote **"not measured"** for `diversity` and
    `barcodes` because both align whole genomes — a reading of the code, not a measurement. That
    honesty is the task: measure them on the same records, `--force`, `/usr/bin/time -v`.
  - The `substitutions` curve predicts they are expensive; a prediction is not an entry. If a run
    does not finish, **that is the finding**, written as what was attempted and for how long.
  **Done when:** `docs/M8.3_scale.md` carries measured wall clock and peak RSS for both modules at
  N = 34 (or the measured timeout with what was tried), and no run module is left "not measured".

- [ ] **P8-05 · CI on every change** *(M8.4)*
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
| D-h | What CPGview is to this project | **checked once, by hand, on an already-public record; not adopted as a dependency** — it is a web service (no version to pin, no argv to record), and §10.12 is undecided, so nothing unpublished goes up. M8.1 asks for the comparison; running it on `NC_050833.1` answers that without putting a web service in the install | §7.8 · §13.8 · §8 entry |
| D-i | What CI must run | **lint + fast tests + tier-A golden on every push; tier-B golden in a job that installs the pinned conda env** — tier A is pure Python over vendored data (§6), so it catches a broken golden number with no tool at all; a tier-B job installing *unpinned* tools would report a different number as success (§9.B-14) | P8-05 · `README` badge |
| D-j | Where `cpstools` lives | **pip into its own env, pinned, declared in `docs/tools/cpstools.md`**, used for the cross-check only — the shape D-f gave Vmatch: a comparison tool is not a module dependency, and the minimum install stays pip-only | §7.8 · §8 entry |

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
| A **nuclear genome** — none exists in the repo; every pinned record is an organelle | anh Điềm / cô Trang | **M6.1** (NUPT/MTPT searches nucleus ↔ plastid ↔ mito), and through it **M8.2** — that milestone stays in ROADMAP, marked waiting |
| Several mitogenomes of **one species** (polymorphism is a within-species property) | anh Điềm / cô Trang | **M6.2** — the repo's ten pinned plant mitogenomes are one per species (`mito-reference` ×8, plus `AP045469.1`, `OX459128.1`); phase 4's mito sketch measured what that much data can say |
| A machine with CGAS's full environment (R + its tool chain), if it cannot be built here | anh Điềm | P8-03's second half — the modules that could not run |
| Real samples for the two manuscripts | cô Trang / anh Điềm | M2.10's final version — the stand-in draft re-runs, never gets patched |
| Real Illumina/ONT reads — the project has none; phase 7's smoke tests ran on simulated pairs | anh Điềm / cô Trang | M7.1 · M7.5 — and M7.3 feeds M6.2 |

*Two rows deliberately absent, because they are **not** waiting on data:* **M6.3** (parallel
plastid–mito comparison) has a stand-in pair already — `AP045470.1` (plastid) and `AP045469.1`
(mito, with this project's own prototype annotation), the pair phase 4 sketched; it is deferred by
**D5-07** (plastid first), a decision, not a dataset. And **M6.4/M6.5** left the plan on
10/10/2026 — divergence time and ancestral state are the PI's analyses, not modules.*
