# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

**Phase 10 is complete** — every figure anh Điềm's feedback bundle
(`feedback/goi_tool_hinh_IR_2026-10-16/`) names is now reachable from a run directory through
`plastid figures`, or is a named skip. Read that bundle — the **whole** of it — before working
anything here: `README.md` (what it holds, and the numbers the team must still settle) ·
`03_tich_hop_PTS/README_tich_hop.md` (how it is meant to be embedded) · `figure_spec_mau.py` (why an
R figure is **not** a `FigureSpec`, in the bundle's own words) · `BANG_MAP_PTS_den_hinh.md` (paper
figure → PTS table — a *claim* to be measured, not a fact) · `04_kiem_chung/KET_QUA_KIEM_CHUNG.md`.
A session that reads only the map designs on a premise the tables do not support; that has happened
once already. The task below is that bundle's sharpest finding, and it is not a figure.

## H-25 — the IR gate refuses a pair whose two annotated copies differ in length

**Goal.** `measure_ir` takes the mirror constant from the two annotated boundaries
(`first.start + second.end - 1`). A copy *d* bases longer therefore shifts every base *d* positions
away from its true partner, and the two copies measure as **unrelated** although the repeat itself
is exact. Measured 17/10/2026 on the bundle's eight Iridaceae records: *d* = 1
(`Chasmanthe_floribunda` NC_072562), 2 (`Iris_hermona` NC_081931) and 4
(`Trimezia_steyermarkii` NC_072662) give identities **0.2924 / 0.2636 / 0.2778**. The gate then
stops `plastid` on three records that are perfectly good plastomes.

**Steps**
1. Reproduce it as a **failing test first**, on those three records — they are in the bundle, under
   `02_IR_boundary/data_mau/`. Do not substitute a synthetic pair: the whole finding is about real
   annotations, and a synthetic one would be built to pass.
2. Decide what the mirror constant *should* be when the two annotated copies differ in length, and
   write down why. The copies are the same sequence; the annotation is what differs. Whatever the
   rule, the golden genome's numbers must not move — NC_050833.1's pair is equal-length, so a
   correct fix is a no-op there, and that is the regression to keep.
3. Re-read what the gate is *for* before loosening anything. It exists because naive detection got
   2 of 9 published plastomes wrong and said nothing (§8.3); a fix that lets a genuinely mismatched
   pair through is worse than the defect. The test says what must **still** be refused.
4. If the rule changes, `core/ir.py`'s docstring and the gate's message change with it, and the IR
   row of `docs/things_we_do_not_do.md` is re-read for anything now false.

**Done when:** the three bundle records pass the gate with an identity ≥ 99.9 %; the golden
genome's IR length and identity are **unchanged**; a pair that genuinely differs is still refused
with `IRGateError` and exit 2; and `make check` is green.

**Traps**
- **Two IR copies of unequal length may be a real annotation error**, not a tool bug — the bundle
  records of `Iris_missouriensis` that a record's own annotation can disagree with itself. Decide
  which this is, per record, and write it down; do not assume it is always ours.
- The identity is a **contract** number (26,214 bp · 100.00 %). A fix that moves it on the golden is
  a wrong fix, not a new number.
- `ROADMAP.md` → Known defects holds H-21 … H-24 beside this one. They are separate tasks, never
  swept up into this one.
