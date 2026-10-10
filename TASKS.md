# TASKS.md — active tasks

> **Phase 9 is OPEN** on `feat/phase-9` — **packaging · licence review · release** (M9.1 … M9.6;
> **M9.7 the Web UI stays deferred** — ROADMAP says "until the core is stable", and that is not a
> task). Phases 1–5, 4a, 3B, 2B, 8 and 8B are closed and squashed on `main`; `phase-8b` is tagged.
> Milestones: `ROADMAP.md` → Phase 9 · full scope: PROJECT.md §13.9 · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

**The phase's own check** (ROADMAP · §13.9): the golden numbers re-run **inside Docker** — one input,
one set of numbers, three environments (WSL · CI · Docker). **Measured 11/10/2026**: `docker` was
absent at the start of the session (*"could not be found… activate the WSL integration"*), then
present — `/usr/bin/docker`, client 29.8.0 · server 29.8.0 (linux) — once the dev turned Docker
Desktop's WSL integration on. **P9-03 is unblocked; nothing on this board waits.**

**Where the phase stands** (11/10/2026): **P9-01 … P9-06 are done**, and the candidate is
**handed over for the team to run** (P9-08) — the licence is MIT and MISA
left the package (§8.74) · the wheel/sdist build and install into a clean venv (§8.75) · the
release image builds and the golden numbers reproduce inside it (§8.76) · README, the module ×
platform × RAM table and the "things we don't do" statement are written and executed (§8.77).
**P9-07 is half done**: `CITATION.cff` exists and validates, `docs/release.md` is the procedure —
the DOI, the PyPI upload, the authors and the video need the team (see "Waiting on data/people").
**Nothing public happens before P9-08 comes back**: the release follows a run on somebody else's
machine, not our own green.

## Phase 9 — the board

- [x] **P9-01 · The licence inventory, and the MISA question** *(M9.3)* — **DONE 11/10/2026
      (§8.74)**: `LICENSE` = **MIT** (D-p) + `license = "MIT"` / `license-files` in pyproject;
      wheel 146 → **145 files, `vendor/misa` gone**; MISA is **obtained** (D-q) with the download
      line and both hashes in `docs/tools/misa.md`, the module taking `--set misa=`, reading the
      release from the script's own header and refusing any but 2.1; `misa.ini` is ours;
      `docs/licensing.md` carries the 14-tool inventory with licences **measured** from conda
      metadata; the guard test was **watched failing on the old tree** before it passed. Golden
      unchanged: 37 = 31 · 4 · 2, `make check` green (856 fast + 85 golden).
- [x] **P9-02 · Package metadata + a clean-venv install** *(M9.2, pip half)* — **DONE 11/10/2026
      (§8.75)**: version **1.0.0** (D-r) + classifiers/keywords/urls (no Windows classifier —
      it does not run natively there); `uv build` → wheel 145 files, sdist 279; installed into a
      **fresh venv** (`uv venv --seed`: this machine's `python3 -m venv` has no `ensurepip`) and
      run from the installed console script — **153,724 bp · 37.85 % · 132 = 86 + 38 + 8**,
      exit 0, on a minimum install of biopython + numpy + PyYAML.

- [x] **P9-03 · The Docker image, and the three-environment check** *(M9.2, Docker half + the
      phase's in-phase check)* — **DONE 11/10/2026 (§8.76)**: `docker build` **exit 0**, image
      1,5 GB, golden suite chạy **ngay lúc build**: **54 passed · 30 skipped · 1 deselected**
      (mọi skip là test cần `.panel/`, **có nêu lý do**; không skip nào vì thiếu công cụ). Chạy
      trong chính image: `genome_stats` → **153.724 bp · 37,85 % · 132 = 86 + 38 + 8 · IR 26.214 @
      100,00 %** và `ssr` → **37 = 31 + 4 + 2** — **y hệt WSL**.
  - The image is real now (`docker/Dockerfile`, replacing the 689 B skeleton): micromamba base,
    **the pinned set** (python 3.12 · perl 5.32.1 · mafft 7.526 · mummer4 4.0.1 · paml 4.10.10 ·
    blast 2.17.0 · aragorn 1.2.41 · vmatch 2.3.1) + **table2asn fetched with its md5** + **MISA
    fetched with its md5** (never baked in from our side, D-q). Versions printed inside; the golden
    suite runs **at build time** under `PTS_REQUIRE_TOOLS=1`, so an image that cannot reproduce the
    numbers **does not build**. `.dockerignore` keeps `.panel/`, `.data/` and `.git` out.
  - **The first build found a real defect** — and only the container could: `gene_normalize` died
    from an installed package (`repo_root()` walked up from `site-packages` looking for
    `pyproject.toml`). Fixed, with two regression tests, §9.N-1. The first build also failed the two
    tests needing Vmatch and table2asn, which were absent by an earlier decision — both are in the
    image now, so nothing skips.
  - **Done when:** `docker build` then `docker run --rm <image>` prints the golden numbers, the
    tier-B run inside is green, versions match the pins.

- [x] **P9-04 · README + quickstart a non-expert can follow** *(M9.1)* — **DONE 11/10/2026
      (§8.77)**: install (Linux · Windows/Docker · from a checkout) → first run → where results
      land → **all five exit codes**, measured: `2` is *both* a wrong command line and the IR gate
      (refused on an amplicon, verified), `3` a missing tool, `4` a naming conflict. Every command
      was run **verbatim in the clean venv** (`plastid --version` → 1.0.0 · one genome · two
      genomes with `scope=cds`), and the run directory matches the README's own picture. The stale
      status line is now a version a test holds to `pyproject.toml`.

- [x] **P9-05 · Module × platform × RAM table + the three support levels** *(M9.4)* — **DONE
      11/10/2026 (§8.77)**: `docs/platforms.md` — all **22 modules**, each with its tool
      requirement, arity, IR gate and measured cost **with the scope named** (34-genome runs vs the
      single-genome `repeats` row measured today, 0.8 s · 87 MB), the machine stated once
      (i7-1165G7 · ~7.8 GB), the three support levels, and the seven modules with **no scale row**
      named rather than implied. Guarded by a test that compares the table's module set to the
      registry's. **Found on the way**: the registry's `citations` field is never read and has
      drifted (the annotators declare no tool there) — so the table's tool column is compiled from
      the adapters each module calls, and the drift is named in the doc.

- [x] **P9-06 · "Things we don't do"** *(M9.5)* — **DONE 11/10/2026 (§8.77)**:
      `docs/things_we_do_not_do.md`, twelve items, each with its **measured** reason and what to do
      instead — no native Windows · no macOS promise · no assembly (Phase 7 waits on reads) · no
      RNA-editing prediction without RNA-seq · no whole-panel alignment (8 genomes 2 min 44 s; 16
      did not finish in 7 min) · one organelle per run · the mito round parked · divergence
      time/ancestral state are the PI's · **no external tool is shipped** · the IR gate refuses
      fragment inputs by design · no web UI · the panel cost, with its machine. `README.md` links
      it from the top.

- [~] **P9-08 · The candidate goes to the team before it goes anywhere public** — **handed over
      11/10/2026**: the artifacts are on the team drive at
      `D:\business\pts\PTS_1.0.0_test_2026-10-11\` (wheel · sdist · the golden record · SHA256SUMS
      · the Vietnamese note `docs/handoff_test_1.0.0.md`), built from `003cfa1` and **re-verified
      after the build** (installed into a clean venv, golden numbers reproduced).
  - **Why it is a task and not a nicety**: every number so far was measured on **one machine**. The
    phase's own check is three environments; a release needs a *fourth* — somebody else's.
  - Three levels, in the note: **mức 1** (5 min, no tools, anyone) · **mức 2** (with MISA/MAFFT,
    anh Điềm) · **mức 3** (`pytest` from the sdist — the sdist carries the tests and the golden
    data, so no repository access is needed).
  - **Done when:** anh Điềm has run mức 2 or 3 on his own machine and reported the outcome — a
    pass, or the verbatim error with its exit code and the run directory — and cô Trang has read a
    `report.md` and said whether the tables and wording are usable. **Anything they find becomes a
    regression test before it is fixed** (CLAUDE.md → testing discipline).

- [~] **P9-07 · Release artifacts: `CITATION.cff` · DOI · the release itself** *(M9.6)* — **half
      done 11/10/2026**: `CITATION.cff` is written for **1.0.0** and **validated** —
      `uv tool run cffconvert --validate` → *"Citation metadata are valid according to schema
      version 1.2.0"*; `docs/release.md` is the ordered procedure with the check between each step.
  - **What is left, and why it is not done here:** the **DOI** needs a Zenodo account and a release
    (and this repository is **private**, which Zenodo's GitHub integration does not archive), the
    **PyPI upload** claims the name `plastid-suite` permanently, the **authors** in `CITATION.cff`
    are the team's to name, and the **how-to video** is recorded by a human. All four are
    outward-facing and effectively permanent — **the dev's explicit go, never a side effect of a
    green build.**
  - **Done when:** `CITATION.cff` carries the DOI and the real authors, the DOI resolves, the tag
    exists on `main`, PyPI serves 1.0.0, and the README shows the DOI.

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Who | Recorded in |
|---|---|---|---|---|
| D-p | **The toolkit's own licence** | **MIT.** The two runtime dependencies are permissive (Biopython's licence, PyYAML MIT), the ecosystem this toolkit cross-checks against is MIT (cpstools, CGAS), and "academic, non-commercial" describes the *project*, not a restriction on the code. GPL-3.0 is the alternative, and it would be the honest choice only if we wanted the toolkit itself to stay copyleft like the tools whose *ideas* it borrows | dev | `LICENSE` + pyproject `license` · PROJECT.md §3 |
| D-q | **Does the released toolkit ship `misa.pl`?** | **No — fetch on first use, md5-verified, and ask IPK in parallel.** Measured: MISA's zip carries no licence and `SOURCE.json` says so; redistribution is the one act a missing licence does not permit. `SOURCE.json` already records the URL, the zip md5 and both file hashes, so the fetch is one named, verifiable step and §6's golden (37 = 31 · 4 · 2) is unchanged. If IPK grants permission in writing, revisit — the third option, replacing MISA with our own finder, is a module job that must reproduce 37 exactly, and is not this phase | dev (PI informed if it changes what runs offline) | PROJECT.md §3 · `docs/licensing.md` |
| D-r | **The first release's version number** | **1.0.0.** A DOI is meant to be cited, the CLI (`plastid run --module … --input … --out … --dataset …`) has been stable across 22 modules and six phases, and the gaps are *named* in M9.5 rather than hidden behind a 0.x. 0.9.0 is the honest alternative if the team reads "assembly is still missing" as "not 1.0 yet" | dev | pyproject + `CITATION.cff` · PROJECT.md §3 |
| D-s | **Where it is published** | **Build and verify here; publish on the dev's explicit go.** PyPI (name `plastid-suite` to be checked as available) + a GitHub release that Zenodo archives. **Note measured 11/10/2026: the repository is private**, and Zenodo's GitHub integration archives public repositories — so either the repo is made public at release time, or the tarball is uploaded to Zenodo by hand. A **conda-forge recipe is a separate, slower path** (a PR to staged-recipes, reviewed by strangers) — `environment.yml` already gives conda users the tool environment, so the recipe is a follow-up unless the team wants it inside this phase | dev | PROJECT.md §3 · `docs/release.md` |
| D-t | **What `gene_normalize` does for a user who installed the package** | **Ask them for the reference, and say so well.** Measured 11/10/2026 (§9.N-1): the recorded naming reference is a **checkout** path (`tests/data/…`), so an installed toolkit has none — it now fails with a message naming `--set naming_reference=<file>`, which is what D5-09 says anyway (the naming record is the end user's choice, not a toolkit default). The alternative is to **ship the record inside the package** (+170 kB) so the pipeline runs out of the box — available, not chosen, because a default nobody chose is the kind of statement this project avoids | dev (PI is the one who picks the ★ record, D5-09) | PROJECT.md §3 · §9.N-1 |

## Pitfalls — this phase's own traps

| Trap | What to do instead |
|---|---|
| Reviewing the **source tree** for licences | Review the **wheel**: measured 11/10/2026 it carries 146 files, `misa.pl` among them — what ships is what is licensed |
| "Not stated" read as "no problem" | No licence = all rights reserved. Record it as a blocker with its remedy, never as a blank cell |
| A Docker image that inherits the host | A bind-mounted tool or an inherited PATH proves the *machine*, not the image. The container's versions are printed from inside it |
| A RAM number without its machine | Quote M8.3's (i7-1165G7 · ~7.8 GB) with every figure, or the figure means nothing |
| macOS ticked because Linux worked | No macOS machine has ever run this (ROADMAP → Waiting). The column reads **unverified** until one does |
| A quickstart written but never executed | Run every command verbatim in a clean venv — the README's stale status line (measured today) is this same defect in miniature |
| A publish run as a side effect | PyPI and Zenodo are permanent and outward-facing: the dev says go, explicitly, each time |

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Written permission from **IPK Gatersleben** (MISA's authors) to redistribute `misa.pl` | cô Trang / anh Điềm — one mail | D-q's preferred path; **the release does not wait on it** if D-q is answered "no, fetch on first use" |
| The **PyPI name** `plastid-suite` availability + the account to publish | dev | P9-07's publish step |
| A **Zenodo** account, and the GitHub release it archives | dev | P9-07's DOI |
| **The team runs the 1.0.0 candidate** — mức 1/2/3 in `docs/handoff_test_1.0.0.md` | anh Điềm (mức 2–3) · cô Trang (đọc `report.md`) | **P9-08**, and therefore the release: one machine's green is not a release's green |
| The **short how-to video** — a human records it | the team | P9-07's last artifact |
| **Who is named as an author** in `CITATION.cff` — it says "PlastidToolSuite contributors" today, and that is a decision, not a placeholder (CLAUDE.md → who decides) | the team | P9-07's DOI: a Zenodo deposit records its authors permanently |
| The **repository's visibility** — it is **private** (measured 11/10/2026), and Zenodo's GitHub integration archives public repositories: make it public at release time, or upload the tarball by hand | dev | P9-07's DOI step |
| A **macOS machine** | the team | P9-05's macOS column — stays "unverified" until one exists |

*Not on this board, so nobody plans around them: **Phase 7's assembly** (waits on the raw Illumina/ONT
reads — the project's oldest blocker: ROADMAP → Waiting) and the **whole mito round** (Phase 6 ·
M7.3 · M8.2 — parked by scope, CLAUDE.md → "Priority today").*
