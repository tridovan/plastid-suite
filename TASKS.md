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

**What exists today, and what this phase owes** — measured 11/10/2026, so nobody re-derives it:
`README.md` (52 lines; its status line still says *"Phase 8 … open on `feat/phase-8`"* — stale twice
over) · `docker/Dockerfile` (689 B **skeleton**) · `environment.yml` (versions pinned) ·
`pyproject.toml` (`version = "0.1.0"`) · 22 modules registered. **P9-01 closed the licence
question**: `LICENSE` is MIT, and the wheel (145 files) carries no third-party program — MISA,
formerly bundled with no stated licence, is obtained instead (§8.74). Still absent: `CITATION.cff`
and any published artifact.

## Phase 9 — the board

- [x] **P9-01 · The licence inventory, and the MISA question** *(M9.3)* — **DONE 11/10/2026
      (§8.74)**: `LICENSE` = **MIT** (D-p) + `license = "MIT"` / `license-files` in pyproject;
      wheel 146 → **145 files, `vendor/misa` gone**; MISA is **obtained** (D-q) with the download
      line and both hashes in `docs/tools/misa.md`, the module taking `--set misa=`, reading the
      release from the script's own header and refusing any but 2.1; `misa.ini` is ours;
      `docs/licensing.md` carries the 14-tool inventory with licences **measured** from conda
      metadata; the guard test was **watched failing on the old tree** before it passed. Golden
      unchanged: 37 = 31 · 4 · 2, `make check` green (856 fast + 85 golden).
  - *(the checklist this closed, kept as the record)*
  - One table over **everything the toolkit runs** (`src/plastid/external/` — aragorn · blast ·
    mafft · misa · mummer · paml · repeatmatch · vmatch — plus table2asn) and **everything it
    ships** (the wheel's own file list, not the source tree), each row: licence · source · is it
    redistributed? · what happens when it is absent.
  - **The sharp one, measured 11/10/2026**: the wheel ships `misa.pl`, and MISA's own
    `SOURCE.json` records `"license": "NOT STATED - the zip ships no licence file and the IPK page
    only requests the citation."` Decide and execute D-q before any artifact leaves this repo.
  - The toolkit's own `LICENSE` file + `pyproject.toml`'s `license` field, on D-p.
  - Fix the inventory defect found today: `environment.yml`'s "NOT here, on purpose" block has an
    entry **with no tool name** (it is PGA's — §7.8 pins its commit `0683a0ec…` and sha256).
  - **Pitfall:** "the zip is public, so it is free" — no licence means all rights reserved, and
    putting the file in a wheel or in git *is* redistribution. A tool we only *call* is not
    redistributed; a tool we *bundle* is, and that is the whole distinction this task draws.
  - **Done when:** `docs/licensing.md` holds the table, and a new
    `tests/unit/test_licensing.py::test_shipped_third_party_code_states_a_licence` is **red on
    today's tree** (MISA) and green after the decision — `make check` green.

- [x] **P9-02 · Package metadata + a clean-venv install** *(M9.2, pip half)* — **DONE 11/10/2026
      (§8.75)**: version **1.0.0** (D-r) + classifiers/keywords/urls (no Windows classifier —
      it does not run natively there); `uv build` → wheel 145 files, sdist 279; installed into a
      **fresh venv** (`uv venv --seed`: this machine's `python3 -m venv` has no `ensurepip`) and
      run from the installed console script — **153,724 bp · 37.85 % · 132 = 86 + 38 + 8**,
      exit 0, on a minimum install of biopython + numpy + PyYAML.

- [ ] **P9-03 · The Docker image, and the three-environment check** *(M9.2, Docker half + the
      phase's in-phase check)*
  - Replace the skeleton with the real image: pinned tool versions **equal to `environment.yml` /
    §7.8**, no host PATH, no bind-mounted tools, `PTS_REQUIRE_TOOLS=1` so a missing tool fails
    instead of skipping.
  - Print the versions **inside** the container — a version nobody can read is not a pin.
  - **Done when:** `docker build` then `docker run --rm <image>` prints the tier-A golden numbers
    (153,724 bp · 37.85 % · 132 = 86 + 38 + 8 · IR 26,214 bp @ 100.00 % · 37 SSR), the tier-B run
    with the tools inside is green, and the versions printed match the pins.

- [ ] **P9-04 · README + quickstart a non-expert can follow** *(M9.1)*
  - The reader is a **researcher in biotechnology/bioinformatics, not a command-line expert**
    (CLAUDE.md → "What this is"): install → first run → where the results land → what exit 2 means.
  - Every command **run verbatim, in order, in a clean venv** — a quickstart that fails at line 3 is
    the first thing a biologist meets.
  - Correct the stale status line, and make it a thing a test checks rather than a thing we remember.
  - **Done when:** the README's own commands were executed verbatim in a clean venv and produced the
    output the README claims, and `tests/unit/test_project_docs.py` fails if the status line drifts
    from the closed phases again.

- [ ] **P9-05 · Module × platform × RAM table + the three support levels** *(M9.4)*
  - All **22 registered modules** × tool requirement (from the registry's own `citations`, not from
    memory) × peak RSS, quoted from `docs/M8.3_scale.md` **with its machine named** (i7-1165G7,
    ~7.8 GB — a RAM number without a machine is not a number, M8.3's own rule).
  - Three support levels: **Linux native · official** · **Windows via Docker · official** ·
    **macOS · experimental, unverified**.
  - **Done when:** the doc lists every module id and a test asserts that set **equals** the
    registry's — so a module added later cannot silently miss the table.

- [ ] **P9-06 · "Things we don't do"** *(M9.5 — block the wrong expectation before it is formed)*
  - No native Windows (Docker only) · assembly needs real reads, which the project does not have —
    Phase 7 · no RNA-editing prediction without RNA-seq data · the whole-genome-alignment patience
    limit with its measured numbers (`substitutions scope=whole`: 8 genomes 2 min 44 s, 16 did not
    finish in 7 min) · one organelle per run, never mixed · the mito analyses parked by scope.
  - **Done when:** every item carries its **measured** reason and a pointer to where the number
    lives, and `README.md` links the statement from the top.

- [ ] **P9-07 · Release artifacts: `CITATION.cff` · DOI · the release itself** *(M9.6)*
  - `CITATION.cff` with the release version, the DOI field, and the team as authors; README shows
    the DOI and the version.
  - `git tag v<version>` on `main` + the GitHub release Zenodo archives.
  - The how-to video is listed as an artifact **the team records** — no code produces it.
  - **Pitfall:** everything here is outward-facing and effectively permanent (a DOI is minted per
    version; a PyPI name is taken the moment it is uploaded). **The dev's explicit go, never a side
    effect of a green build.**
  - **Needs (D-s):** the dev's go, the Zenodo account, and the PyPI name check.
  - **Done when:** `CITATION.cff` parses against the CFF schema with the version and DOI, the DOI
    resolves, the tag exists on `main`, and the README a reader lands on names the released version.

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Who | Recorded in |
|---|---|---|---|---|
| D-p | **The toolkit's own licence** | **MIT.** The two runtime dependencies are permissive (Biopython's licence, PyYAML MIT), the ecosystem this toolkit cross-checks against is MIT (cpstools, CGAS), and "academic, non-commercial" describes the *project*, not a restriction on the code. GPL-3.0 is the alternative, and it would be the honest choice only if we wanted the toolkit itself to stay copyleft like the tools whose *ideas* it borrows | dev | `LICENSE` + pyproject `license` · PROJECT.md §3 |
| D-q | **Does the released toolkit ship `misa.pl`?** | **No — fetch on first use, md5-verified, and ask IPK in parallel.** Measured: MISA's zip carries no licence and `SOURCE.json` says so; redistribution is the one act a missing licence does not permit. `SOURCE.json` already records the URL, the zip md5 and both file hashes, so the fetch is one named, verifiable step and §6's golden (37 = 31 · 4 · 2) is unchanged. If IPK grants permission in writing, revisit — the third option, replacing MISA with our own finder, is a module job that must reproduce 37 exactly, and is not this phase | dev (PI informed if it changes what runs offline) | PROJECT.md §3 · `docs/licensing.md` |
| D-r | **The first release's version number** | **1.0.0.** A DOI is meant to be cited, the CLI (`plastid run --module … --input … --out … --dataset …`) has been stable across 22 modules and six phases, and the gaps are *named* in M9.5 rather than hidden behind a 0.x. 0.9.0 is the honest alternative if the team reads "assembly is still missing" as "not 1.0 yet" | dev | pyproject + `CITATION.cff` · PROJECT.md §3 |
| D-s | **Where it is published** | **Build and verify here; publish on the dev's explicit go.** PyPI (name `plastid-suite` to be checked as available) + a GitHub release that Zenodo archives. A **conda-forge recipe is a separate, slower path** (a PR to staged-recipes, reviewed by strangers) — `environment.yml` already gives conda users the tool environment, so the recipe is a follow-up unless the team wants it inside this phase | dev | PROJECT.md §3 |

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
| The **short how-to video** — a human records it | the team | P9-07's last artifact |
| A **macOS machine** | the team | P9-05's macOS column — stays "unverified" until one exists |

*Not on this board, so nobody plans around them: **Phase 7's assembly** (waits on the raw Illumina/ONT
reads — the project's oldest blocker: ROADMAP → Waiting) and the **whole mito round** (Phase 6 ·
M7.3 · M8.2 — parked by scope, CLAUDE.md → "Priority today").*
