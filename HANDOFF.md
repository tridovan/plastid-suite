# HANDOFF.md

What the next session needs and git cannot tell it. Overwritten at the end of every task, never appended;
earlier sessions are `git log -p -- HANDOFF.md`.

## State — 17/10/2026

Branch `feat/phase-9`. Phases 1–5, 4a, 3B, 2B, 8 and 8B are closed and squashed on `main`. Phase 9 is open:
the licence is MIT, the version is 1.0.0, the wheel installs into a clean venv, the Docker image runs the
golden suite when it is built, and the README, the module × platform table and "things we don't do" were
written and executed verbatim. anh Điềm tests the code from git (P9-08); the DOI, PyPI and the video wait
on that and on the team. Milestone detail: `ROADMAP.md`.

The session documents were cut to `CLAUDE.md` (rules) plus the three that carry work — `ROADMAP.md`,
`TASKS.md`, this file. The 39-session log this file used to hold is `git show 4004e1a:HANDOFF.md`.

## Next

**Figure tool — fixes from anh Điềm's package** (`feedback/goi_tool_hinh_IR_2026-10-16/`). Not scoped yet;
`TASKS.md` says how it starts. Then **Phase 7**: pick a published paper with public reads and its assembled
plastome, and judge the assembly on that (`ROADMAP.md` → Phase 7). No team data is needed.

## Waiting on the team — nothing here is ours to do

- **P9-08.** anh Điềm pulls the code from git and tests it by the process that suits the tool: input in,
  output out. No package or guide is handed over (the one on the team drive is stale). Point the pull at
  `feat/phase-9` — `main` holds only the closed phases. Whatever is found becomes a regression test
  before it is fixed.
- **P9-07.** The DOI, the PyPI upload and the video need the authors named in `CITATION.cff`, a Zenodo account,
  the repository made public (it is private today) and the dev's explicit go. `docs/release.md` has the order.
- Nothing else waits on people. The tools are judged on a published paper's public input against the paper's
  own results (`CLAUDE.md` → How the tools are judged): do not ask the team for data, samples or decisions.

## Live risks

- The registry's `citations` field is stored and never read, and it has drifted (the annotators declare no
  tool there). Feed it from each module's `Result`, or drop the field.
- `gene_normalize` from an installed package asks for `--set naming_reference=<file>`, because the naming
  record is the end user's choice and not a default (D5-09). The first-run experience is worth watching.
- The Docker image was built and verified on one machine; macOS has never been run.
- CI's tier B installs from conda-forge and bioconda at run time, so a channel-side change can turn it red
  with nothing changed here; the job prints every pinned version, so the failure is readable.
- Four named-but-unfixed defects (H-21 … H-24) are listed in `ROADMAP.md` → Known defects. Each is a task of
  its own with a regression test, never a silent fix inside another task.
- Tier B needs the pinned tools on `PATH` and the `PTS_*_DIR` variables. Without them the tests skip with a
  named reason and the run looks green — set `PTS_REQUIRE_TOOLS=1` to make a skip a failure.
