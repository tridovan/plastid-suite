# HANDOFF.md

What the next session needs and git cannot tell it. Overwritten at the end of every task, never appended;
earlier sessions are `git log -p -- HANDOFF.md`.

## State — 17/10/2026

Branch `feat/phase-10`. **Phase 9 is closed and squashed on `main`**, tagged `phase-9`: MIT, the wheel
installed into a clean venv, the Docker image with the golden suite run inside it, and the user docs.
P9-07 (the DOI, PyPI, the video) and P9-08 (anh Điềm's report) stay open on the team's side.

**Phase 10 is the paper's figures.** anh Điềm's bundle (`feedback/goi_tool_hinh_IR_2026-10-16/`) is the
seven R scripts and the patched CPJSdraw that drew the Iridaceae manuscript's Figures 2–9.
**M10.1 and M10.2 are done.** `plastid ir-junction` reads a record's own four regions in the GenBank
frame, drives CPJSdraw, and writes the configuration, the junction table, the SVG and the PNG into a
run directory — on the bundle's eight records the configuration is reproduced **8/8** and the SVG
comes out **byte-identical** (`9ae599a9…`). And the R engine runs: `external/rscript.py` with the
seven scripts as package data, all seven through the adapter in 35 s under R 4.3.3. `TASKS.md` holds
the next task — wiring figures to that engine (M10.3/M10.4); `ROADMAP.md` → Phase 10 has the rest.

**All pushed.** `main` (the Phase 9 squash), `feat/phase-9`, the tag `phase-9` and this branch are on
`origin`; every ref is 0/0 with its remote. anh Điềm's pull for P9-08 is `main`.

## Next

**`TASKS.md`** — wire figures to the R engine (M10.3/M10.4): the three the paper has and PTS lacks
(amino acid, substitutions, Ka/Ks), the long→wide pivots, and — where PTS already draws the figure —
a test that the matplotlib and R versions agree. Then M10.5, the tree, which stays the user's own
script. Then **Phase 7**, assembly, judged on a published paper's public reads against its
assembled plastome.

## Waiting on the team — nothing here is ours to do

- **P9-08.** anh Điềm tests the 1.0.0 code. The pull is now **`main`** — Phase 9 is merged there; the
  `feat/phase-9` pointer in the older notes is out of date. Whatever he finds becomes a regression test
  before it is fixed.
- **P9-07.** The DOI, the PyPI upload and the video need the authors named in `CITATION.cff`, a Zenodo
  account, the repository made public (it is private today) and the dev's explicit go — `docs/release.md`
  has the order.
- A macOS machine, for the one column of `docs/platforms.md` that has never been run.

Nothing else waits on people: the tools are judged on a published paper's public input against the
paper's own results (`CLAUDE.md` → How the tools are judged).

## Live risks

- **H-25** (new, `ROADMAP.md`): the IR gate refuses a record whose two annotated IR copies differ in
  length — three of the eight records in the bundle, whose inverted repeats are in fact exact. It makes
  `plastid` unable to process them at all, and it is why the junction figure reads the annotation
  directly rather than through `detect_ir`. Its own task, with its own regression test.
- **The R half is measured but not yet in the repo.** R 4.3.3, `rsvg-convert` and `libsvg-perl` were
  installed on this machine to run the bundle's scripts, and all seven run and write all eight
  figures. Nothing about R is declared yet, so `docs/tools/r.md` and the `environment.yml` pins must
  ship in the session that first uses them — a tool that exists only on one machine is not part of
  the toolkit. (`rsvg-convert` and `libsvg-perl` are declared: `docs/tools/cpjsdraw.md`.)
- **CPJSdraw refuses a record without a `gene` feature**, and `to_genbank` writes `CDS`/`tRNA`/`rRNA`
  only — so a genome `plastid annotate` produced cannot be drawn as it stands. Measured 17/10/2026,
  written down in `docs/tools/cpjsdraw.md`; the refusal is loud, not silent.
- **CPJSdraw is not bundled and cannot be** (no licence — measured 17/10/2026, GitHub answers
  `license: null`). The bundle's integration README claims it is MIT; that claim is wrong, and the repo
  says so. The patched build is identified by md5, so a pristine clone stops the run.
- The registry's `citations` field is stored and never read, and it has drifted. Feed it from each
  module's `Result`, or drop the field.
- `gene_normalize` from an installed package asks for `--set naming_reference=<file>`, because the naming
  record is the end user's choice (D5-09). The first-run experience is worth watching.
- CI's tier B installs from conda-forge and bioconda at run time, so a channel-side change can turn it
  red with nothing changed here; the job prints every pinned version.
- Four named-but-unfixed defects (H-21 … H-24) are in `ROADMAP.md` → Known defects, now joined by H-25.
  Each is a task of its own with a regression test, never a silent fix inside another task.
- Tier B needs the pinned tools on `PATH` and the `PTS_*_DIR` variables; without them the tests skip
  with a named reason and the run looks green — `PTS_REQUIRE_TOOLS=1` makes a skip a failure.
