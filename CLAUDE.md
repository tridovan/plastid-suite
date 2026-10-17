# CLAUDE.md

Rules for working in this repo. Where the work stands lives in three files — `ROADMAP.md` (the whole
project), `TASKS.md` (the task in flight), `HANDOFF.md` (what the next session needs). Everything else is
archive: open it only when a task names it.

## What this is

**PlastidToolSuite (PTS)** — a Python toolkit that analyses **plastid (chloroplast) genomes** for a plant
research group (PI: cô Trang · bioinformatician: anh Điềm, who runs the real tool chain on a machine of
the lab's own · dev · Hermes). Academic, non-commercial. End users are **researchers in biotechnology /
bioinformatics, not command-line experts**: a command line on Ubuntu (Windows through Docker), on ordinary
computers.

Deliverables: **two species manuscripts** (Oldenlandia, Coleus — tables, figures, draft Results) and a
**released toolkit** (Linux native, Windows via Docker, DOI). The tool paper (Eleutherine bulbosa —
NUPT/MTPT, mito structure) waits with the parked mito round.

Non-negotiable product requirements: every run emits a **tool + parameters report**; output is **split by
run/dataset** (`<dataset>/<run_id>/`), never overwritten in place; the minimum install (pip only) runs the
pure-Python analyses with no external tool.

**This round is plastid only** (D5-07 — a scope rule, not an ordering). A run serves **one organelle** —
declared by the record's `/organelle` qualifier, or by the path a bare sequence enters — and its input set
never mixes plastid and mito. Every module is shared, plastid-only or mito-only (table: PROJECT.md §13 →
"Hai bào quan"); a module run against the wrong organelle must stop by name rather than describe the wrong
genome. The mito milestones are **parked by name** — M6.1 · M6.2 · M6.3 · M7.3 · M8.2 — and never picked up.
The mito code that exists (`mito_annotate`, `mito_genes`) is neither deleted nor extended, and stays under
`make check`. Divergence time and ancestral state are the PI's analyses, not modules.

Two decisions shape the pipeline:

- **Annotation is ours** (D5-08): `plastid annotate`, in-house. Chloe is dropped, PGA survives only as the
  measured counter-example, and a reference-set audit gates every annotation run.
  Design: `docs/D5-08_plastid_annotator.md`.
- **Gene names are unified on the reference DB, before annotation** (D5-09): reference supply → audit → the
  end user picks **one naming record ★** → `gene_normalize` renames the reference DB → assembly → annotate →
  `gene_normalize` in verify mode → analyses. Samples inherit the vocabulary; a change at verify is a
  *finding*, not routine. The ★ record is fixed for the project, never re-picked per run by taxonomy.
  **No table mixes the two vocabularies** (A everywhere; B only in the pairing table). `gene_normalize` is
  its own module — `gene_names` is not grown to cover it. Design: `docs/D5-09_reference_naming.md`.

## The deliverables decide, not the tools

A tool named anywhere in the plan is a **hypothesis** until it is measured on this project's data — PGA was
picked on a published benchmark and gave 147 genes where the golden has 132 (§8.15). When a tool is
unsuitable, buggy, unmaintained, unlicensable or absent, **replace it or write that part ourselves**, and
keep the same non-negotiables: offline, provenance recorded (name · version · argv), declared in the repo.
**Never bend a number to fit a tool** — if a tool cannot reproduce the golden numbers, that *is* the finding.

## Golden numbers — the contract

Golden sample **Iris domestica NC_050833.1** (comparison genome: Iris gatesii NC_024936.1).

| Quantity | Value |
|---|---|
| Length · GC | **153,724 bp** · **37.85 %** |
| Genes | **132 = 86 CDS + 38 tRNA + 8 rRNA** — the `gene` feature count (133) is a diagnostic only; never write "133 genes" |
| IR | **26,214 bp**, two copies **100.00 % identical** · LSC **83,127** · SSC **18,169** |
| SSR | **37** = mono 31 · di 4 · tri 2 (MISA 2.1 and cpstools agree) |
| Ts/Tv | **2.6917** = 323 Ts / 120 Tv — scope: **CDS only · 79 shared genes · MAFFT · IR excluded** |

Changing any of these needs fresh measured evidence (provenance: PROJECT.md §6). Another scope is another
number (whole-alignment Ts/Tv is 1.8416) and must be labelled.

## Data

No real project dataset exists yet (status 08/10/2026): everything in the tree — the pinned panels,
`.data/thanhdiem_data/` — is stand-in or tool-test material. Fetch a stand-in with `plastid panel pin`
(accession · date · sha256 recorded), label every output made from it, and never let its numbers reach a
manuscript table as the team's sample. When the real data arrives, re-run; do not patch numbers.

## Language

Chat with the team: Vietnamese. `PROJECT.md`: Vietnamese. Everything else — this file, `TASKS.md`, code,
comments, commit messages, generated reports: English.

## Invariants — never violate

- **IR gate**: the two IR copies must be ≥ 99.9 % identical before *any* position-based analysis. It lives
  in `core/ir.py`; every position-based module calls `ctx.require_ir()`. Never bypass, disable or make it
  optional — the *gate*, not the detector, is the guard (§8.3, §8.5).
- **Every exported number carries a scope label** (e.g. "CDS only, IR excluded").
- **Measured or not at all.** A number that was not measured is written "not measured" — never estimated,
  rounded to something plausible, or dropped because a case fails. What a subagent read or measured is a
  claim until you verify it firsthand.
- **A tool that exists only on one machine is not part of the toolkit.** Before a module uses an external
  tool, the tool is declared in the repo — in `environment.yml` (conda; **version pinned**: an unpinned
  mafft gave two golden Ts/Tv numbers from one genome, §9.B-14) or in `docs/tools/<tool>.md` (source, exact
  install line and hash, and why it is not in `environment.yml`; MISA is the worked example — obtained,
  never shipped, §8.74). Installing a tool locally to verify work is fine; stopping there is not — the repo
  change ships in the same session, or the number is reproducible by nobody.

## Code rules

The import layering is enforced by `tests/unit/test_architecture.py` and style by `make lint`; the rest is
in the surrounding code — match it.

- A module is a pure function `run(genome, params, ctx) -> Result`: no cwd, no globals. Results are data;
  the export layer writes files, so a module never names a path (`ctx.out(...)`).
- Every external tool goes through `run_tool()` — the only `subprocess` call, and the only place that
  records provenance (name · version · path · argv · time · exit code · logs).
- Report templates reference table cells (`{{table.column}}`), never literal numbers.
- **Extend by adding, not editing**: a new analysis is a new file in `modules/` plus one registry line.
  **DRY after the second real copy** — never before.
- **Fail loud.** A gate raises, a missing tool is named, a disagreement is recorded — nothing degrades into
  a warning nobody reads. Specific exceptions, actionable messages. Docstrings state units and coordinate
  conventions: that is where the silent bugs live.
- **Banned outright:** `sys.exit` outside `cli.py` · hardcoded absolute paths · hand-rolled FASTA/GenBank
  parsers (use BioPython) · `shell=True` · monkey-patching, `sys.argv` mutation, hidden cwd or globals.

This is production code written once: no stubs, no placeholder that looks real. Deferred by decision, so
not to be proposed: Snakemake (PROJECT.md §3), plugin entry points, report i18n, a web UI.

## Testing discipline

1. **Core logic is test-first** — data model · IR gate · coordinates · result schema · run layout · report
   resolution: the failing check (a golden number, a schema assertion, a gate case) comes before the code.
   Figures and exploratory analysis are not; test the data fed in, not the pixels.
2. **Every measured failure becomes a regression test before the fix**, and is never deleted or weakened —
   naive IR detection, an off-by-one pairing constant, SSR 36 vs 37, difflib vs MAFFT are why.
3. **Gates and tools get negative tests**: bad IR ⇒ `IRGateError` + exit 2 · missing tool ⇒ skip with a
   printed reason, or exit 3 when required · out-of-scope numbers ⇒ labelled. A gate that cannot fail
   loudly is not done.
4. **Never weaken or skip a test to make it pass.** If a golden number disagrees with the code, the code or
   the measurement is wrong until proven otherwise — stop and say so. No test touches the network.

## How we work

Three files carry the work, and **git is the history**, so none of them keeps one.

| File | Holds | Rewritten |
|---|---|---|
| `ROADMAP.md` | the whole project: phases, one line per milestone, what each open item waits on | when a milestone closes, or a phase opens or closes |
| `TASKS.md` | the one task in flight: goal · steps · **Done when** (a runnable check) · traps | at the start of each task — the finished one is deleted, not ticked |
| `HANDOFF.md` | what the next session needs and git cannot tell it: state · next · blockers · live risks | at the end of each task — overwritten, never appended |

A session reads `HANDOFF.md`, then `TASKS.md`; `ROADMAP.md` only to pick or close work. A task is done only
when its Done-when passes — no check, no done. A task that spans several files or changes a schema is
**planned and approved before coding**. Then: failing check → code → `make check` green → **one commit**
holding the code, its tests and the three files → pushed to the phase branch. Earlier sessions are
`git log -p -- HANDOFF.md`. A new status, summary or log file is never the answer: a fact belongs in the
code, a test, the commit message, or one of these three files.

**Archive** — read-only; open it only when a task cites it. `PROJECT.md` (the team document, Vietnamese:
decisions §3, golden §6, tool pins §7.8, errors §9, open questions §10) · `docs/log/` (the run log, entries
`§8.1`–`§8.77`) · `docs/` (designs, tool notes). Citations such as `§8.27` and `§9.C-34` are permanent and
test-guarded — never renumber. The log is frozen: a new measurement is pinned by a test and described in
the commit message. Add a log entry only when the team must read it (next number, a row in the `PROJECT.md`
§8 index, the count bumped in `tests/unit/test_project_docs.py`). A decision that changes the science or is
expensive to reverse goes in `PROJECT.md` §3; a design that must outlive the session goes in `docs/` — never
in a session-local plan file. A statement that has become false is corrected in place.

## Who decides

Default: **the dev decides**, with evidence, and records it. Bring a decision to the team only when it
changes the science or is expensive to reverse.

| Decision | Who |
|---|---|
| Architecture · formats · libraries · tool choice · local vs remote · test strategy | dev |
| What to analyse · scope of a scientific claim · what goes in a manuscript · taxonomy choices for a paper | PI (cô Trang) |
| The project's naming reference ★ (D5-09) — picked from the audit report; the PI signs off the names a manuscript prints | the end user |
| The tool environment on the bioinformatician's own machine (the toolkit's own pins stay the dev's) | anh Điềm |
| Big or irreversible work · anything published (PyPI, Zenodo, a public repo) | the dev's explicit go, each time |

A blocked task names the **specific** decision or dataset it waits for — "needs the team's call" is not a
blocker. Cheap-to-reverse decisions are made now and corrected later, not queued.

## Commits

Conventional style, English, no diacritics: `feat:` `fix:` `docs:` `test:` `chore:` `refactor:`. **No
assistant attribution** — no `Co-Authored-By: Claude …`, no "generated with" trailer; commits are authored
by the team only. Never commit raw data dumps, credentials or API keys.

One commit per task, on the phase branch `feat/phase-N`; each stands alone (check it out, `make check`
passes). **The branch holds the detail; `main` holds the phases.** Closing a phase: tag it
(`git tag phase-N feat/phase-N`), then `git merge --squash` the branch into `main` as one commit that lists
its tasks. Branches are never deleted. Between phases `main` takes at most a `docs:` commit.
