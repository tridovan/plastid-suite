---
name: make-tasks
description: Build or refresh TASKS.md for a phase in this repo — derive tasks from ROADMAP milestones, give each a measurable Done-when, and surface the hard decisions and known pitfalls. Use when a phase starts or ends, when TASKS.md is stale or full of finished work, or when the user asks to clear/plan tasks.
---

# Making TASKS.md

`TASKS.md` is the **current phase's** board: short enough to read every session, precise enough that
a fresh session can start without asking anything. Finished work does not live here — it lives in
`ROADMAP.md` (ticked milestones) and in git.

## Inputs — read in this order

1. `ROADMAP.md` — the phase's milestones and its **in-phase checks**; milestones are the raw material
   for tasks.
2. `HANDOFF.md` → Resume point — where the last session stopped.
3. `PROJECT.md` §13.x — the phase's full scope text (ROADMAP is the summary; §13 has the detail).
4. `CLAUDE.md` → "Who decides what" — so every open decision lands on the right desk.

## Procedure

1. **Clear the finished phase.** Delete completed blocks and any Done table. The record is the ticked
   milestones in ROADMAP plus the commits; do not keep a cross-phase archive here.
2. **One task per milestone**: id `P<phase>-<nn>`, and state the milestone it serves (`M2.4`).
   Order tasks by the sequence you would actually execute them in, not by milestone number.
3. **Every task carries:**
   - a checklist of concrete sub-items — what will exist when it is done;
   - **Done when:** a runnable, measurable check (a golden number, a test name, a command);
   - **Needs:** the *specific* decision or dataset it waits on, if any — never "the team's call".
4. **A "Decisions to make" section**, one row per open decision: the question, the recommended
   default *with its evidence*, who decides (per CLAUDE.md: dev unless it changes the science), and
   where it gets recorded. Decisions the dev can make with evidence are made — not queued.
5. **A "Pitfalls" section** for the phase: pull from this project's measured failure modes
   (PROJECT.md §9, the CGAS review §7.9) and the phase's own traps. One line each: the trap → what to
   do instead. Add a task-specific **Pitfall:** line inside a task when it is not phase-wide.
6. **A "Waiting on data/people" table**, naming the task each item blocks.
7. **Self-check before finishing:** every task has a Done-when; every blocker names something
   specific; nothing finished is left; the file is still under ~200 lines.

## Quality bar

- A task without a measurable Done-when is not a task — rewrite it or drop it.
- If a Done-when cannot be met in this phase (it needs a later phase's output), move it to that
  phase's milestone list in ROADMAP instead of leaving a permanent `[~]`.
- Prefer a check that already exists (a golden number in §6, a test in `tests/golden/`) over inventing
  a new measurement.
- Pitfalls must be *specific to this repo* — "write clean code" is noise; "the truncated ycf1 has no
  CDS, so `count_genes()` is right and a 'fix' is wrong" is a pitfall.
