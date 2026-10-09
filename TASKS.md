# TASKS.md — active tasks

> **Phase 2B is OPEN** on `feat/phase-2b` — the plastome analysis set, completed: π per intron/IGS
> (**H-07**, M2.9's remainder), the **repeats** module (**M2.12** · H-16) and the two species'
> manuscript draft re-run (**M2.10**). Why this before Phase 8: Phase 6 and 7 wait on data the
> project does not have, and Phase 8's cross-check against CPStools is not honest while one of
> those analyses does not exist. Milestones: `ROADMAP.md` → Phase 2B (the three milestones
> themselves live in Phase 2's own list).
> Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 2B — the plastome analysis set (M2.9 · M2.10 · M2.12)

**Goal:** the last three unblocked analysis outputs, so that every analysis the toolkit promises
exists before it is cross-checked — and so the two species manuscripts have their remaining tables.
**In-phase checks:** the repeats count on the golden is a **new §8 entry** — measured, never the
26/9 survey's "9 matches" (§9.H's trap) — with the IR pair among its matches · π per intron and per
IGS over the spans the shared layer derives, the **17 introns** `introns` reports as the
cross-check, the IR excluded and the scope labelled · the draft carries its π and Ts/Tv sections
and its figures, every table labelled *stand-in* and scoped.
*Order: the spans come first (P2B-01) — the repeats module and the draft both want the region
machinery that task puts in the shared layer, and DRY's rule fires at its second real copy.*

- [ ] **P2B-01 · The intron and IGS spans, in the shared layer** *(M2.9's remainder · **H-07**)*
  - Extract the intron derivation from `modules/introns.py` into `core/` — it would be DRY's second
    real copy, which is exactly the rule's trigger — and have `diversity`'s region machinery take
    the intron **and** intergenic-spacer spans from it.
  - π per intron and per IGS on the pinned pair; the span count cross-checked against `introns`'
    **17** on the golden; the IR excluded, the scope labelled.
  **Pitfall:** an IGS is a *definition*, not a fact — a gap between two genes on opposite strands,
  or one that runs through the IR, can be read two ways. Say which gaps count in the scope label,
  or two runs will disagree without disagreeing (see D-g).
  **Done when:** `diversity` reports π per intron and per IGS, the intron count equals `introns`'
  17 on the golden, and the numbers are a new §8 entry.

- [ ] **P2B-02 · The repeats module** *(M2.12 · H-16 — the last analysis box with no module)*
  - `repeat-match` (**MUMmer3**, a plain binary, no licence constraint) for the exact matches —
    dispersed and palindromic — with **Vmatch 2.3.1** as the independent second engine, the way
    MISA and cpstools both give 37. Vmatch is declared in `environment.yml` (§7.8) and is installed
    **in the session that first uses it**, never "later".
  - Every match ≥ 30 bp with its coordinates **and its region** (LSC · IRa · SSC · IRb), the IR
    pair included; the count measured on the golden fixture, not inherited.
  **Pitfall:** the 26/9 survey's "**9 matches ≥30 bp**, with the 26,214 bp IR pair" was measured on
  *other* data (§7.4 · §9.H's trap). Inheriting it puts a number in a paper that no run in this
  repo produced — the milestone says "measured, not inherited" for exactly this reason.
  **Pitfall:** Vmatch is **academic-only** (§7.8): a module that needs it cannot be in the minimum
  install, and the release has to say so (M9.3).
  **Done when:** the module runs on the golden, reports every match ≥30 bp with its region, the
  count is a new §8 entry, and a test pins the IR pair among the matches.

- [ ] **P2B-03 · The manuscript draft, re-run** *(M2.10)*
  - Regenerate both species' tables and the Results draft from the pinned stand-in set — now with
    the **π** and **Ts/Tv** sections (MAFFT has been present since 03/10; the draft simply predates
    them) and the figures, plus whatever P2B-01/P2B-02 add.
  - Every table carries its **scope** label and the **stand-in** label; no number is typed into a
    template (the report templates reference table cells, never literals).
  **Done when:** the draft carries π and Ts/Tv, every exported table is labelled *stand-in* and
  scoped, and the numbers agree with §6's golden values (Ts/Tv **2.6917** · SSR **37** ·
  **132 = 86 + 38 + 8**).

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Recorded in |
|---|---|---|---|
| D-f | The repeats engine | **`repeat-match` ships, Vmatch cross-checks** — MUMmer3 is a plain binary with no licence constraint (§7.8) while Vmatch is academic-only, and this project's pattern for a count it publishes is **two independent tools agreeing** (MISA 2.1 + cpstools 3.0 → 37) | §7.8 · PROJECT.md §3 |
| D-g | What an IGS is | **the gap between two called genes inside one region, the IR excluded**, stated in the scope label — the same span set `introns` and `diversity` already use, so the three modules cannot disagree about where a boundary is | the module's scope label · §8 |

## Pitfalls — this phase's own traps

| Trap | What to do instead |
|---|---|
| The 26/9 survey's numbers ("9 matches ≥30 bp", "21 blocks · mean 91.5") were measured on **other data** | measure on the golden fixture and write a **new §8 entry**; never inherit a survey number (§9.H is the worked example) |
| "Fixing" the truncated `ycf1` — the IR copy carries no CDS, so it is not counted | that is right: `count_genes()` counts what the record has (§8.34), and the second `ycf1` is not a second gene |
| π over a region without the gate | the IR gate runs first, the IR is excluded from π, and the scope says so — the golden's 2.6917 is *CDS only · 79 genes · IR excluded* |
| Two engines, one number: repeats counted by tools that define the classes differently | state the class (exact · dispersed · palindromic · tandem) and the ≥30 bp floor **with** every count, and record the second engine's total beside the first |
| A module that needs Vmatch | the minimum install is pip-only and must keep working: the module skips with a named reason when the tool is absent, and `make check` never requires it |

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Real samples for the two manuscripts | cô Trang / anh Điềm | M2.10's **final** version — the stand-in draft can be re-run now, and a real sample means re-running, never patching numbers |
| The ★ naming reference (D5-09) | the end user, from the audit report | nothing in this phase; a run on a real sample would inherit it |
| A macOS machine | the team | M9.4's platform claims (Phase 9) |
| Plant mito data · a nuclear genome · real reads | anh Điềm / cô Trang | M6.1–M6.3 · M7.1 · M7.5 — **Phases 6 and 7 stay closed for now** |
