# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

## The `plastid ir-junction` command (M10.1)

**Goal.** Make the junction figure reachable: a user runs one command over the records they have and
gets the figure, its table and the provenance of every tool that touched it.

**Already built and verified** (on this branch): `core/junctions.py` cuts a record into LSC/IRb/SSC/IRa
in the record's own frame, wrap-aware; `external/cpjsdraw.py` runs the patched CPJSdraw and rasterises
with `rsvg-convert`. Measured 17/10/2026 on the bundle's eight Iridaceae records: the derivation
reproduces the bundle's own configuration **8/8**, and the adapter's SVG is **byte-identical** to the
bundle's (`md5 9ae599a9…`).

**Steps**
1. An `ir-junction` subcommand in `cli.py`, following `plastid run`'s shape: `-i/--input` (repeatable),
   `-o/--out`, `--dataset`, and the run layout — one run, one directory, never overwritten in place.
2. Write the configuration, the junction table (CSV) and the run info; every number carries its scope
   label, and `--label` marks output made from someone else's records.
3. `--tree`/`--reference` to order the species by distance, `--width` for the PNG, `--cpjsdraw` for the
   script's path.
4. Report the IR gate's verdict for each record beside the table: H-25 means the gate can refuse a
   record whose annotation this figure still draws, and a reader must not have to discover that.

**Done when:** on a machine with perl, the SVG module and `rsvg-convert`,
`plastid ir-junction -i <the bundle's 8 records> -o /tmp/x` writes the table and an SVG byte-identical
to the bundle's own (`fd0ed35c…`, `9ae599a9…`); without `perl` or `rsvg-convert` it exits **3** naming
the missing tool; and `make check` is green.

**Traps**
- The bundle's records are the paper's, not PTS's: every output made from them is labelled.
- CPJSdraw is **not** bundled, and a pristine clone is refused by name — the compound-IRa patch is a
  correctness fix, not a cosmetic one (`docs/tools/cpjsdraw.md`).
- The figure is drawn in the **GenBank** frame, deliberately: `structure`'s normalized frame moves
  every junction coordinate.

## After it

**M10.2–M10.5** in `ROADMAP.md` → Phase 10: the R engine, the three figures PTS lacks, the pivots, and
the tree. Two facts are already measured, so nobody needs to re-measure them: the seven R scripts run
here under **R 4.3.3** and write all eight figures, and `khoi_tao()` calls `setwd(OUT_DIR)`, so
`--inputdir` must be **absolute** or every script fails to find its own data.

**H-25** (`ROADMAP.md` → Known defects) is its own task with its own regression test: the IR gate
refuses three of those eight records although their inverted repeats are exact.
