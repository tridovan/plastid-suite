# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

## Figure tool — fixes from anh Điềm's package

**Not scoped yet; the first step is the scoping itself.**

**Input.** `feedback/goi_tool_hinh_IR_2026-10-16/` — untracked and in `.gitignore`: it is the team's
material, so read it, don't edit it, don't move it into the repo. In it: the R scripts behind the paper's
Figures 3–9 (`01_R_scripts/`), the IR-boundary tool for Figure 2 (`02_IR_boundary/`, a patched CPJSdraw), a
kit for plugging both into PTS (`03_tich_hop_PTS/`, with a figure → module map) and the author's own
verification record (`04_kiem_chung/KET_QUA_KIEM_CHUNG.md`).

**Steps**
1. Read the package's `README.md` and `03_tich_hop_PTS/README_tich_hop.md`, then its verification record.
2. Compare with what PTS draws today (`src/plastid/result/figures/`): which figures differ, and why.
3. Write the list of fixes here, each with its own Done-when, and get the dev's go.

**Done when (this step):** this file lists the fixes, every one with a runnable Done-when, and the dev has
said go. Nothing is built before that.

**Traps**
- The package's sample data is the paper's, not PTS's stand-ins — label any output made from it.
- A figure's numbers come from a table that carries its scope label; the figure code does not recompute them.
- Any external tool the fixes need (R and its packages, CPJSdraw, `rsvg-convert`) is declared in the repo in
  the same session it is first used — `environment.yml`, or `docs/tools/<tool>.md` (CLAUDE.md → Invariants).
