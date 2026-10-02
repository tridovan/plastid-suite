# CLAUDE.md — project briefing & working rules

Injected into every session. Self-contained: reading this file alone should be enough to understand
what the project is and the rules. For where work stands, read `ROADMAP.md` + `TASKS.md` (both small).
`PROJECT.md` is the evidence archive — open it only when a task needs it (see Reading guide).

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

**Current state: Phase 1 — framework, pre-implementation. No code yet** (repo holds docs only).
Progress tracker: **`ROADMAP.md`** (milestone checklists). Exact resume point:
**`TASKS.md` → Current focus**.

## Golden numbers — the contract

Golden sample: **Iris domestica NC_050833.1** (comparison genome: Iris gatesii NC_024936.1).

| Quantity | Value |
|---|---|
| Length · GC | **153,724 bp** · **37.85 %** |
| Genes | **132 = 86 CDS + 38 tRNA + 8 rRNA** (`gene` features = 133 → diagnostic only; unique names 112) |
| IR | **26,214 bp**, two copies **100.00 % identical** · LSC **83,127** · SSC **18,169** |
| SSR | **37** = mono 31 · di 4 · tri 2 (MISA 2.1 and cpstools agree) |
| Ts/Tv | **2.686** = 317 Ts / 118 Tv — scope: **CDS only · 76 shared genes · MAFFT · IR excluded** |

Changing any of these requires fresh measured evidence. Full provenance: PROJECT.md §6.

## Reading guide — consult, don't consume

| File / section | Read when |
|---|---|
| `TASKS.md` (small) | Always — it opens with **Current focus** (resume point) + handoff log |
| `ROADMAP.md` | Always — phase milestones + progress checklist; source for specs/tasks |
| PROJECT.md §3 | Before a design decision not already covered by ROADMAP/CLAUDE |
| PROJECT.md §6 · §8 · §9 | Before touching numbers, tests, or claiming "we measured X" |
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
4. **Not test-first:** pure rendering (figures) and exploratory analysis. Test the data fed into a
   figure, not the pixels; snapshot tests optional.
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

1. Pick a task in `TASKS.md`; set it `DOING`; keep **Current focus** current.
2. Non-trivial work (multi-file, new module, schema change): plan first and get approval before coding.
3. Verify as you go — write the check first where the Testing discipline applies (below); every step
   has a measurable check (test / golden number). Never defer verification to the end.
4. On completion: run `make check`; record new measured facts in PROJECT.md (§8 log, §9 errors);
   set the task `DONE`.
5. Update `TASKS.md`: **Current focus** + append a handoff entry (Did / Blocked / Next).
6. Commit only when the user asks.

Approved designs that must outlive the session go into PROJECT.md (§3 decisions, §13 roadmap) or
`docs/` — never rely on session-local plan files being readable later.

## Handoff (write this before context runs low)

At the end of every session — or the moment you notice context filling up — update `TASKS.md`:
**Current focus** (one short block: exact next action) and append a dated handoff entry:
`Did: … / Blocked: … / Next: …`. Assume the next reader is a fresh agent with zero history.

## Commits

Conventional-commit style, English, no diacritics: `feat:` `fix:` `docs:` `test:` `chore:` `refactor:`.
End every commit message with:

```
Co-Authored-By: Claude Code <noreply@anthropic.com>
```

## Don't

- Don't run tests that hit the network — golden data is vendored in the repo.
- Don't commit raw data dumps, credentials, or API keys.
- Don't restructure PROJECT.md; append to §8/§9 and update the §3 decisions table.
- Don't invent benchmark numbers or "expected" values; every threshold comes from §6 or a fresh run.

## Subagents

- Use read-only Explore agents for broad reading ("where is X", many files).
- Use a Plan agent to design before large changes; verify its claims before adopting.
- Subagents are for reading and verification — their reports never become evidence without a
  firsthand check.
