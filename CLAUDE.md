# CLAUDE.md — project briefing & working rules

Injected into every session. Self-contained: reading this file alone should be enough to understand
what the project is and the rules. For where work stands, read `HANDOFF.md` (resume point) plus
`ROADMAP.md` + `TASKS.md` (both small). `PROJECT.md` is the evidence archive — open it only when a
task needs it (see Reading guide).

## What this is

**PlastidToolSuite** — a Python toolkit to analyse **plastid (chloroplast) and mitochondrial genomes**
for a plant research group (PI: cô Trang; bioinformatician: anh Điềm, who runs the real tool chain on
his own machine; dev; Hermes). Academic, non-commercial. End users are biologists who are **not**
command-line experts, working on ordinary computers. ~20 analysis modules: genome structure, gene
content, codon usage/RSCU, SSR, repeats, annotation, assembly, multi-genome comparison, RNA editing,
divergence time, NUPT/MTPT…

Deliverables driving the work:
- **2 species manuscripts** (Oldenlandia, Coleus) — tables, figures, draft Results sections
- **1 tool paper** (Eleutherine bulbosa) — NUPT/MTPT + polymorphic mito structure
- **A released toolkit** — Linux native + Windows via Docker (macOS experimental), DOI

Non-negotiable product requirements: every run emits a **tool + parameters report**, output is
**split by run/dataset** so results are traceable (never overwritten in place), and the minimum
install (pip only) runs the pure-Python analyses without any external tools.

**Two organelles, one per run.** Every run serves exactly one organelle — plastid or mitochondrion —
and the record says which through its own `/organelle` qualifier; a bare sequence declares it by the
path it enters (`plastid mito annotate` is the mito path). A run's input set is single-organelle
too: a run — and a panel — holds plastid records or mito records, **never both**. The one named
exception is the **transfer analysis across compartments** — nucleus ↔ plastid ↔ mito (NUPT/MTPT,
Phase 6): it reads all three deliberately, and it is not a mixed set fed to the single-organelle
modules. Every module belongs to one of three groups — **shared** · **plastid only** · **mito
only** — and a module run against the wrong organelle must stop by name rather than describe the
wrong genome. The group table, the measured reason per module, and what is still to build (the
declaration in the registry, the input-set check, the seven shared analyses not yet opened to mito):
PROJECT.md §13 → "Hai bào quan".

**Priority today: plastid** (team decision 05/10/2026, D5-07). Plant mito data is still scarce, so
new module work serves the plastome; the existing mito path stays as it is, and the structure above
is kept precisely so mito analyses can be added later without rework — the deferred items carry
names in `TASKS.md`.

**The annotate step is ours (D5-08).** The pipeline is assembly → annotate → gene-name unification
→ analyses. Annotation will be built in-house on CGAS's model with our own search-and-decide core;
**Chloe is dropped**, PGA stays only as the measured counter-example, and a reference-set audit
gates every annotation run. Design: `docs/D5-08_plastid_annotator.md`; decision: PROJECT.md §3.

**Current state: Phase 5 in progress.** Progress: `ROADMAP.md` (milestones) · `TASKS.md` (active
checklists with Done-when). Exact resume point: **`HANDOFF.md` → Resume point**.

## Supreme objective — the deliverables decide, not the tools

The project exists to produce the three deliverables above (PROJECT.md §1–§2), for a research group
of biologists who are not command-line experts and work on ordinary computers. Everything else in
this file — the seams, the golden numbers, the layering, the tool declarations — serves those
deliverables. **Assembling a particular set of existing tools is not the goal, and never was.**

So:

- A tool named in `ROADMAP.md` / `TASKS.md` / PROJECT.md is a **hypothesis** until it is measured on
  this project's own data. Measured behaviour beats reputation: PGA was picked on a published
  benchmark and had to be re-examined the moment it was actually run on the golden genome
  (PROJECT.md §8.15: 147 genes where the golden has 132, plus three silent failure modes).
- When a tool is unsuitable, buggy, unmaintained, unlicensable, or simply not there, **replace it —
  or write that part ourselves**. Record what was measured that forced the change (§8/§9) and the
  decision (§3). "The plan says tool X" is never a reason to keep a tool that fails.
- The replacement must meet the same non-negotiables: offline, provenance recorded (name · version ·
  argv), declared in the repo, and no step whose version nobody can name.
- **Never bend a number to fit a tool.** If the tool cannot reproduce the golden numbers, that *is*
  the finding: say it, and change the tool. A pipeline that produces plausible numbers nobody can
  reproduce is the failure this project exists to avoid.

## Golden numbers — the contract

Golden sample: **Iris domestica NC_050833.1** (comparison genome: Iris gatesii NC_024936.1).

| Quantity | Value |
|---|---|
| Length · GC | **153,724 bp** · **37.85 %** |
| Genes | **132 = 86 CDS + 38 tRNA + 8 rRNA** (`gene` features = 133 → diagnostic only; unique names 112) |
| IR | **26,214 bp**, two copies **100.00 % identical** · LSC **83,127** · SSC **18,169** |
| SSR | **37** = mono 31 · di 4 · tri 2 (MISA 2.1 and cpstools agree) |
| Ts/Tv | **2.6917** = 323 Ts / 120 Tv — scope: **CDS only · 79 shared genes · MAFFT · IR excluded** |

The Ts/Tv row was **2.686 = 317 Ts / 118 Tv over 76 genes** until 04/10/2026. The scope definition
never changed; the *set* it ranges over did — three loci (atpI · psaI · psbB) that NC_024936.1
spells atpL · psaL · bspB are the same loci, measured by sequence, position and `/product`, and
were invisible while the shared set was counted by spelling. Re-measured with them in: 79 genes,
68,448 columns, 2.6917 (PROJECT.md §8.27). The old numbers stay in the log as what the narrower
set measured; the whole-alignment scope is still a different number (1.8416).

Changing any of these requires fresh measured evidence. Full provenance: PROJECT.md §6.

## Data — real samples vs development stand-ins

Real datasets arrive from anh Điềm / cô Trang. Until one exists, **development and testing may use
public stand-in data** (GenBank records of the same genus, or the nearest available relative):
choosing a *suitable* stand-in is a dev call, not a PI call.

A stand-in never stops being a stand-in:

- fetch it through `plastid panel pin` — accession · date · sha256 recorded, never by hand;
- label every output made from it; never let a stand-in number reach a manuscript table as if it
  were the team's sample;
- when the real dataset arrives, re-run — do not patch numbers.

## Reading guide — consult, don't consume

| File / section | Read when |
|---|---|
| `HANDOFF.md` | First — the resume point, plus the last sessions' Did / Validation / Risks / Next |
| `TASKS.md` | Always — active tasks as checklists; each ends with a measurable "Done when" |
| `ROADMAP.md` | Always — phase milestones + progress checklist; source for specs/tasks |
| PROJECT.md §3 | Before a design decision not already covered by ROADMAP/CLAUDE |
| PROJECT.md §6 · §9 | Before touching numbers, tests, or claiming "we measured X" |
| PROJECT.md §8 → `docs/log/` | The run log. **One file per phase**; §8 in PROJECT.md is its index, one line per entry. Cite an entry by number (`§8.27`) — never by file and line |
| PROJECT.md §13.x | When a phase needs its full original scope text (ROADMAP.md is the summary) |
| PROJECT.md §7.8 | When installing or pinning external tools (versions, sources, md5) |

## Language rules

| Where | Language |
|---|---|
| Chat replies to the team | Vietnamese |
| `PROJECT.md` (team document) | Vietnamese |
| `CLAUDE.md`, `TASKS.md`, code, comments, commit messages | English |
| Generated tool reports / Results drafts | English |

## Invariants — never violate

- **Gene counting** is always written as **132 = 86 CDS + 38 tRNA + 8 rRNA**. The `gene` feature count
  (133) is a diagnostic only — never write "133 genes".
- **IR gate**: the two IR copies must be ≥ 99.9 % identical before *any* position-based analysis.
  It lives in the framework core; never bypass, disable, or "make optional". (Measured reason:
  naive IR detection failed silently in 2/9 published plastomes — §9A.)
- **Every exported number carries a scope label** (e.g. "CDS only, IR excluded").
- **Evidence rule**: anything a subagent read or measured is a *claim* until verified firsthand.
  Only verified facts go into PROJECT.md.
- **Measurement rule**: if a number is needed and not measured, write "not measured" — never estimate,
  round to something plausible, or silently drop a failing case.
- **A tool that exists only on one machine is not part of the toolkit.** The moment a module uses an
  external tool, that tool must already be declared **in the repo**, in one of three places:
  1. `environment.yml` — conda/bioconda tools, **with the version pinned** to what PROJECT.md §7.8
     records. An unpinned `mafft` is not a shortcut, it is two users getting two different golden
     numbers from the same genome (§9.B-13 measured exactly that: Ts/Tv 2.667 vs 2.686).
  2. `src/plastid/vendor/<tool>/` + `SOURCE.json` — small script tools, vendored byte-for-byte with
     md5/sha256 and the retrieval date. MISA is the worked example.
  3. `docs/` — a tool that fits neither, with its source, its exact install line, and why it is not
     vendored.

  Installing a tool locally to *verify* work is expected and encouraged. What is not allowed is
  stopping there: **the repo change ships in the same session as the install**, never "later".
  Otherwise every number measured here is reproducible by nobody — not by anh Điềm, not by the
  released toolkit on a biologist's computer, not by us in six months. This is the same failure as a
  number that looks plausible and is wrong; it just fails later.

## Code rules (six seams — full detail in PROJECT.md §13.1)

1. A module is a pure function: `run(genome, params, ctx) -> Result`. No cwd, no globals, no `sys.exit`.
2. Results are data (`Result` / `Table` / `Figure`); the export layer writes files. Modules never write
   or name paths — use `ctx.out(...)`.
3. All external tool calls go through `run_tool()` — the only place that calls `subprocess`, and the
   only place that records provenance (name · version · path · argv · time · exit code · logs).
4. The IR gate lives in `core/ir.py`; every position-based module calls `ctx.require_ir()`.
5. One source of truth for packaging: `pyproject.toml`. No `setup.py`. Version declared once.
6. Report templates reference table cells (`{{table.column}}`, `{{table[row].column}}`) — never
   literal numbers in templates.

Banned outright: bare `except:` · `sys.exit` outside `cli.py` · hardcoded absolute paths ·
hand-rolled FASTA/GenBank parsers (use BioPython) · `shell=True` for tool calls.

## Coding practice — principles before mechanics

**This is production code, written once — there is no MVP.** Every line ships as if it will be read,
run and trusted for years: no throwaway stubs, no "fix it later" shortcuts, no placeholder that looks
real but is not. Getting it right the first time is the default, not an aspiration.

Python is not Java: most of SOLID collapses into fewer, simpler rules here. Defaults: KISS · YAGNI ·
explicit over implicit · boring over clever. Prefer the standard library; add a dependency or an
abstraction only when the second real use appears.

**Design principles, translated to this repo**

- **One responsibility per unit** (SRP) — a module does one analysis, a function does one thing.
  A function that parses *and* writes *and* prints is three functions.
- **Extend by adding, not editing** (OCP) — a new analysis = a new file in `modules/` plus one
  registry line. Never grow an existing module to cover a new job.
- **Composition over inheritance** (LSP) — plain functions and frozen dataclasses; a class only
  when state must travel together. A three-level hierarchy means the design is wrong.
- **Small explicit contracts** (ISP) — `run(genome, params, ctx) -> Result`, nothing wider; a
  module never receives parameters it ignores.
- **Depend on the boundary, not the tool** (DIP) — core knows nothing about external programs;
  `run_tool()` and `external/` are the only doors, and the layering test enforces this.
- **DRY after the second real copy** — never before. Counterexample: CGAS, 14 copy-pasted glue
  blocks and ~2,400 duplicated lines. Same logic in two modules → extract to the shared layer.
- **KISS / YAGNI** — build what the current task needs. Snakemake, plugin entry points and report
  i18n are all deliberately deferred for exactly this reason.
- **Fail loud** ("errors should never pass silently" — Zen of Python) — a gate raises, a missing
  tool is named, a disagreement is recorded. Nothing degrades into a warning nobody reads.
- **Explicit over implicit** — no monkey-patching, no `sys.argv` mutation, no hidden cwd or globals
  (each one is a measured CGAS failure mode).

**Mechanics**

- Python 3.10+; `from __future__ import annotations` everywhere. Types on public functions;
  dataclasses for data; `pathlib` over string paths; context managers for resources; no mutable
  default arguments.
- Names: `snake_case` functions/variables, `PascalCase` classes, short concrete modules. Domain
  abbreviations already in use (LSC, IR, SSC, SSR, CDS) are fine; invent no new ones.
- Docstrings: one line on what it does, plus *why* where non-obvious. Always document units and
  coordinate conventions — that is where silent bugs live.
- Small pure functions: data in, data out — no hidden I/O, no globals, no cwd.
- Errors: raise specific exceptions with actionable messages — `ValueError` for bad input, a domain
  exception (`IRGateError`) for a domain failure. Never a bare `except:`.
- Comments explain why, not what; delete commented-out code instead of keeping it.
- Formatting is not a debate: `ruff format` decides (line length 100, imports absolute and sorted).
- Tests are named after the behaviour ("gate flags diverged copies"), one behaviour per test, and
  always cover the negative case — this project's bugs are silent ones.

## Testing discipline — test-driven development

The golden numbers are known **before** the code exists, so this project is a natural fit for
test-first work. Rules:

1. **Core logic is written test-first.** Data model · IR detection & gate · coordinate normalization ·
   result schema · run layout/manifest · report resolution: write the failing check (a §6 golden
   number, a schema assertion, a gate case) **before** the implementation. Red → green.
2. **Every measured failure becomes a regression test before the fix.** The project's own history has
   several silent-wrong modes — naive IR detection (2/9 cases) · pairing constant off by one (K=c−1 →
   1 bp candidate) · hand-written SSR count 36 vs 37 · difflib vs MAFFT (2.667/2.686). Each is a
   permanent test; none may be deleted or weakened.
3. **Negative tests are mandatory for gates and tools.** Bad IR ⇒ `IRGateError` + exit 2 · missing
   tool ⇒ skip with a printed reason (or exit 3 when required) · out-of-scope numbers ⇒ labelled.
   A gate that cannot fail loudly is not done.
4. **Not test-first:** figures and exploratory analysis — test the data fed in, not the pixels.
5. **Never weaken or skip a test to make it pass.** If a golden number disagrees with the code, the
   code (or the measurement) is wrong until proven otherwise — stop, and record the discrepancy in
   PROJECT.md §9.

## Commands

```
make env      # create venv (uv) + install
make test     # fast tests — tier A, no external tools
make golden   # full golden set; prints a skip reason for every missing tool
make check    # lint + test + golden
```

## Workflow for a task

0. Planning or refreshing a phase's tasks? Follow the **`make-tasks`** skill
   (`.claude/skills/make-tasks/SKILL.md`) — it defines what a task must carry.
1. Pick a task in `TASKS.md`; mark it `[~]`. It is done only when its **Done when** check passes.
2. Non-trivial work (multi-file, new module, schema change): plan first and get approval before coding.
3. Verify as you go — write the check first where the Testing discipline applies (below); every step
   has a measurable check (test / golden number). Never defer verification to the end.
4. On completion: run `make check`; record new measured facts in PROJECT.md (§8 log, §9 errors);
   set the task `DONE`.
5. Tick the checkboxes in `TASKS.md`; append a `HANDOFF.md` entry.
6. **Commit and push each task as it completes** — not at the end of the session. The Done-when
   passes, `make check` is green, that task's files are committed and pushed, then the next task
   starts. Session documents that describe several tasks (`TASKS.md` · `PROJECT.md` · `HANDOFF.md`)
   go in one `docs:` commit at the end. Splitting a pile of finished work into per-task commits
   afterwards means reconstructing intermediate states of shared files by hand — and it is the
   step that gets skipped when someone is in a hurry, which is how a session becomes one
   unreviewable commit.

Approved designs that must outlive the session go into PROJECT.md (§3 decisions, §13 roadmap) or
`docs/` — never rely on session-local plan files being readable later.

## Who decides what — so the PI is never the bottleneck

Default: **the dev decides**, with evidence, and records it in the repo. Bring a decision to the team
only when it changes the science or is expensive to reverse.

| Decision | Who decides | Recorded in |
|---|---|---|
| Architecture · formats · libraries · tool choice · local vs remote · test strategy | Dev (+ Claude) | PROJECT.md §3 or `docs/` |
| What to analyse · scope of a scientific claim · what goes in a manuscript · taxonomy choices for a paper | PI (cô Trang) | PROJECT.md §3 |
| Data and tool environments on the bioinformatician's machine | anh Điềm | TASKS.md → "Waiting on data" |
| Big or irreversible work | dev approves the plan first | plan → TASKS.md |

A blocked task must name the **specific** decision or dataset it waits for — "needs the team's call"
is not an acceptable blocker. Cheap-to-reverse decisions are made now and corrected later, not queued.

## Handoff (write this before context runs low)

At the end of every session — or the moment you notice context filling up — update `HANDOFF.md`:
refresh the **Resume point** and append a dated entry with **Did · Validation · Residual risks ·
Blocked · Next**. Assume the next reader is a fresh agent with zero history.

## Commits

Conventional-commit style, English, no diacritics: `feat:` `fix:` `docs:` `test:` `chore:` `refactor:`.
Do **not** add any assistant attribution — no `Co-Authored-By: Claude …`, no "generated with" trailer.
Commits are authored by the team only.

**One commit per task, pushed when the task is done** (workflow step 6). A commit should stand on
its own: check it out, run `make check`, and it passes — which is also what makes `git bisect`
usable. When one file genuinely carries changes for two tasks, say so in the message rather than
leaving the split ambiguous.

**A phase is worked on its own branch** — `feat/phase-N`, opened when the phase's board goes up and
merged to `main` when the phase closes (fast-forward; the pattern every phase so far followed:
`feat/phase-2` · `-3` · `-3a` · `-4` · `-5`). Two reasons, both from this project's own history: a
phase is the unit `git bisect` and a reviewer can hold in their head, and `main` stays a state
somebody can check out and trust. A sub-phase opened mid-flight (4a) belongs on a branch too — the
"it is only a few commits on top" reasoning is exactly how a phase's work ends up interleaved with
the next one's in `main`'s history.

## Don't

- Don't run tests that hit the network — golden data is vendored in the repo. (Fetching a dataset
  once with `plastid panel pin` is a data step, not a test.)
- Don't commit raw data dumps, credentials, or API keys.
- Don't restructure PROJECT.md; append to §9 and update the §3 decisions table.

## The run log — `docs/log/`

§8 (the log of every real run) lives in `docs/log/`, one file per phase, because 33 entries in one
118 KB file made reading §8.27 mean opening a 216 KB document. The rules, and the reasons:

- **Entry numbers are permanent.** `§8.27` is `§8.27` wherever the text sits — **271 references**
  point at them by number (27 in code docstrings, 36 in tests, 23 in `docs/`, 115 across the four
  session files, 70 inside PROJECT.md). Renaming or renumbering an entry breaks citations nobody
  can re-derive. New entry: the next number, in the file of the phase you are in.
- **`PROJECT.md` §8 is the index** — one line per entry, and the only place a citation has to
  resolve from. `tests/unit/test_project_docs.py` fails the build if a cited `§8.N` has no heading,
  if an entry is missing from the index, or if one number appears in two files.
- **Entries are appended, never rewritten.** A correction is a new entry that says what it corrects
  (the Ts/Tv log is the worked example: the old numbers stay, with what they measured).
- **Move text, don't retype it.** The 04/10/2026 split was verified byte-for-byte — the five files
  concatenated are identical to the §8 that was in git. Reorganising this archive is a diff nobody
  reads; making it provably lossless is the only way it can be trusted.

## Subagents

- Use read-only Explore agents for broad reading ("where is X", many files).
- Use a Plan agent to design before large changes; verify its claims before adopting.
- Subagents are for reading and verification — their reports never become evidence without a
  firsthand check.
