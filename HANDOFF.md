# HANDOFF.md

What the next session needs and git cannot tell it. Overwritten at the end of every task, never appended;
earlier sessions are `git log -p -- HANDOFF.md`.

## State — 17/10/2026

Branch `feat/phase-10` (kept). **Phase 9 and Phase 10 are both closed and squashed on `main`**,
tagged `phase-9` and `phase-10`; the Phase 10 squash is `0037379` and the branch still carries the
six commits it summarises.

**Phase 10 is anh Điềm's feedback bundle**, `feedback/goi_tool_hinh_IR_2026-10-16/` — and it is a
*feedback* phase, so the bundle is the specification. `TASKS.md` opens by saying so, because one
session read only `BANG_MAP_PTS_den_hinh.md` first and designed on a premise the tables do not
support. Read the whole bundle.

**H-25 — the IR gate measured an unequal pair as unrelated.** `measure_ir`'s mirror constant
(`a1 + b2 - 1`, §8.3) assumes the two annotated copies are the same length; when they differ by *d*
it pairs every base *d* positions off and an exact repeat comes back at ~25 %. Fixed: the difference
is taken off the **longer copy at whichever end it sits**, both readings are measured, and
`IRResult.excess_at` records which one produced the number. On the bundle's eight records that takes
the gate from 5/8 to **7/8**; Iris hermona stays refused, correctly — its PGA annotation is 48 bp
wider than the true repeat.

**Why "both ends" and not "the outer end":** measured on `tests/data/iridaceae_edge_cases.gb`,
NC_063920's excess sits at the **inner** end — 92.63 % there against 36.46 % at the outer. "Outer
only" would have rewritten §8.3's own finding. That number is now pinned in `test_ir.py`.

**Phase 10 is complete.** `plastid ir-junction` reproduces the bundle's own junction output byte for
byte, and **every** figure the bundle names renders from a run directory through **`plastid
figures`** — or is a named skip. M10.5, the tree, is the last of them: `R_F9_phylogeny_tree.R` was
already packaged but had no door, because the six figures before it all read **tables** and this one
reads a tree.

## What M10.3/M10.4 built

* **The R figures are a second catalogue**, not `FigureSpec`s: `result/rfigure.py` (`RFigureSpec`,
  `RFigureSkip`), the pivots in `result/figures/r_*.py`, and the pass in `run/rfigures.py`. Three
  reasons, all measured: `Figure` accepts only `svg/pdf/eps` and a test pins that; `FigureSpec.build`
  gets no `RunContext`, so R would run with **no provenance**; and one script writes several named
  images. The bundle's own `figure_spec_mau.py` reaches the same conclusion in its own words.
* **The pivots are pure and tested cell-for-cell** (`tests/unit/test_figures_r_tables.py`), and where
  PTS already draws the figure the two must agree (`tests/unit/test_figures_agree.py`). π is the
  stated exception: `pi_window` draws windows, the paper draws features — two pictures of two tables.
* **Figure 7 had no source at all.** `core/align.py` counted only transitions and transversions, so
  the six complementary types the manuscript draws did not exist here. `substitution_counts()` now
  returns the twelve directed pairs and `substitutions` publishes `substitutions_by_type`
  (0.2.0 → 0.3.0). Every existing table and golden number is unchanged; the counters are additive.
* **The tree figure takes the user's own file, not a table.** `RFigureSpec.source_param` names a
  `--set` parameter holding a file the script reads **directly**: `--set tree=<file.treefile>`. The
  file is handed to R as its input and is **never** copied into the run or cached as a result — the
  caption records its path and sha256 instead of a scope label, because it read no table. Measured
  17/10/2026 on the bundle's 34-taxon tree: 32 tips after the script drops two *Asparagus*, rooted on
  two Orchidaceae, which is the bundle's own verification record reproduced.
* **`plastid figures` gained `--set`** (for `rscript=` and `tree=`), runs both passes, and records
  R's version and argv in `manifest.json` under `figures.tools`. A matplotlib failure no longer
  stops the R pass. The R pass is **not** wired into `run --figures` — the golden `test_e2e.py:379`
  asserts that path exits 3 leaving no `figures/`, and that is a deliberate later decision, not an
  oversight. Nor into `package` or `qc`, which still call `render_run_figures` alone: the
  supplementary package and the QC list carry the matplotlib figures only until someone decides
  otherwise.

## Next

**`TASKS.md`** — **M7.1**, batch assembly, judged on a published paper's public reads against its
assembled plastome. Its first step is finding that paper: no answer key, no evaluation, and a
simulated pair is not a substitute. Then M7.2 and M7.6. **H-21 … H-24** stay open in `ROADMAP.md` →
Known defects, each its own task.

**Closing Phase 11** will be the same outward step, when its milestones are done: `git tag phase-N
feat/phase-N`, then `git merge --squash` onto `main` as one commit listing its tasks. Ask the dev
before either — the tag is permanent.

## Waiting on the team — nothing here is ours to do

* **P9-08.** anh Điềm tests the 1.0.0 code; the pull is `main`.
* **P9-07.** DOI, PyPI, video — need the authors named, a Zenodo account, the repo made public, and
  the dev's explicit go (`docs/release.md`).
* A macOS machine, for the one column of `docs/platforms.md` never run.
* **The bundle's own numbers, which are the team's to settle** (`README.md` → "Lưu ý số liệu"):
  Figure 7 was drawn from **20 comparison pairs** while the caption says 33 species; the manuscript's
  five-species compound-IRa list names *Iris hermona* where the data say *Isophysis tasmanica*; the
  CPJSdraw citation is the PeerJ 2023 paper, not "Xu, W. (2024)". None is a toolkit change.

## Live risks

* **A 34-species RSCU / amino-acid / SSR figure is not reachable from one run.** `codon_usage`,
  `amino_acid` and `ssr` are arity-1 — one genome per run — and the run layer concatenates nothing.
  The figures render what the run holds and say so in the caption; merging runs is its own task.
  Recorded in `docs/things_we_do_not_do.md`.
* **F3 and F4 draw only the species on their own hard-coded 34-name list.** A run of another species
  is refused **by name** (`RFigureSkip`) rather than drawn empty — measured on the Coleus stand-in.
  The list is parsed from the packaged script, so a script edit cannot silently change it.
* **F5's `tRNA`/`rRNA` columns are blank**, not zero: `ssrs.context` is Gene·Intron·IGS only and the
  toolkit does not split gene kinds. The caption says so.
* **`R_F8_KaKs_heatmap.R` cannot draw a one-genome matrix** — measured 17/10/2026, isolated on a
  hand-built matrix: its `mat[keep_genes, keep_spp]` drops a single column to a vector and it halts
  with `replacement has 0 rows`. A one-query run is skipped by name (`RFigureSkip`). The script is
  the manuscript's and is not edited; both limits are written down in `docs/tools/r.md`.
* **The IR gate now weighs both ends of a length difference, and that is a rule, not a tweak.** If a
  future record fails with a small `d`, look at `excess_at` before touching the threshold: the gate
  is doing its job and the annotation is what is wrong (Iris hermona's PGA record: IRB and IRA each
  ~50 bp wider than the true 26,132 bp repeat, on the same sequence as the published record that
  passes).
* **Do not "fix" a failing annotation by falling back to de-novo.** Measured 17/10/2026: on
  NC_063920 the annotation is 26,162/26,161 @ 92.63 % (refused) and de-novo finds **22,998 bp** at
  100 % — the exact core. A fallback would accept that record on a 22,998 bp IR and hide the
  3,163 bp of divergence §8.3 exists to catch. `detect_ir` keeps the annotation and records
  `denovo_length` beside it; that is the designed behaviour (§8.5), not an omission.
* **CPJSdraw is not bundled and cannot be** (no licence, measured 17/10/2026). The bundle's README
  claims MIT; that claim is wrong and the repo says so.
* The registry's `citations` field is stored and never read, and it has drifted.
* CI's tier B installs from conda-forge and bioconda at run time, so a channel-side change can turn
  it red with nothing changed here.
* Four named-but-unfixed defects (H-21 … H-24) are in `ROADMAP.md` → Known defects.
* Tier B needs the pinned tools on `PATH` and the `PTS_*` variables, **plus `PTS_RSCRIPT`** for the
  R figure end-to-end; without them the tests skip with a named reason and the run looks green.
