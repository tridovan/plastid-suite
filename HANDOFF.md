# Development handoff

Đợt 4 of the batch plan, working tree on `feat/m11-ir-rework`, nothing committed. **M11.19** and
the metadata half of **M11.9** were done first (tool roles, environment split, `4 - Beta`); this
session then ran **M11.8** — the whole-pipeline input → output acceptance — and the verdict is
**not accepted as stable**, for reasons stated below. Read `CLAUDE.md`, then `docs/pipeline.md`; in
the maintainer workspace `.dev/TASKS.md` carries the progress notes and
`.dev/records/M11.8_acceptance_2026-10-28.md` the full ledger with every number and its source.

## M11.8 — what the acceptance found

**Two defects, found by running the pipeline rather than by reading it.**

- **K-08 — fixed.** `00_run_info/run_info.md` printed "(no external tools ran)" for a run whose
  earlier modules had called MISA and `repeat-match`: the page was rendered from the last module's
  `RunContext`, and every module gets its own. The manifest was correct all along. The page now
  reads every module's recorded tools from the manifest and de-duplicates by
  (name, version, path). Held by `tests/unit/test_run_info_tools.py` — red before, green after,
  plus an end-to-end CLI case in the failing module order.
- **K-07 — open.** A single inserted or deleted base inside an IR copy splits the repeat across two
  mirror diagonals; the chain builder keeps one half, and the detector's 10 kb length floor accepts
  a ~12.6 kb "IR" at 100 % that **passes the gate**. A run on such a FASTA exits 0 and every
  position-based table is built on the wrong interval. Measured rate on one record: 5 of 22
  deletion positions. The obvious fix was measured and rejected — it creates a new half-length pass
  on `NC_057683.1` — so the real fix needs a "which overlapping pair is the repeat" rule and
  belongs to the IR milestone. Held by `tests/unit/test_ir_indel_collapse.py` (strict xfail plus
  two passing controls) and named in `docs/things_we_do_not_do.md`.

**What reproduced** (independent expectation → observed): a bare FASTA annotated by transfer gives
exactly the published 132 genes = 86 + 38 + 8 at 153,724 bp; ten modules plus figures and a Word
report on that artifact give IR 26,214 bp @ 100 %, LSC 83,127, SSC 18,169 in 26.3 s; the golden
pair gives Ts/Tv **2.6917** (443 substitutions, 68,382 comparable columns) and π 0.007181 with
MAFFT 7.526 and yn00 4.10.10 recorded. The IR panel is 75 pass / 5 fail / 80-of-80 lengths, with
each of the five named; the rotated FASTA set's five differences are one record whose deposited
annotation is one base short of the maximal exact repeat (PTS right, annotation short); the
synthetic set is 28/30 and the two failures are K-07.

## M11.21 — a phylogeny module, built the CGAS way (28/10/2026)

The maintainer asked for a tree-building module in CGAS's shape. It exists and runs:

- `modules/phylogeny.py` (arity 2 = at least 3 genomes): MAFFT **per CDS**, concatenated in
  gene-name order into a supermatrix, then IQ-TREE. One copy per gene name, pseudo excluded,
  partial excluded unless `include_partial=yes`, a gene needs 2 carriers, and every exclusion is a
  named row with its reason (`matrix_genes`).
- `external/iqtree.py` + `iqtree=2.3.6` in `environment.yml` + `docs/tools/iqtree.md` + inventory
  and licence rows. `-m MFP -bb 1000 -alrt 1000 -bnni -nt AUTO -seed 1`; `-bb` below 1000 is
  refused, and the binary's version is checked against the pin.
- Artifacts beside the tables: `supermatrix.fasta`, `tree.treefile`, `tree.iqtree`.
- Measured on 4 panel records: **78 genes, 70,079 sites**, tree with `99.5/98` support,
  MAFFT 7.526 + IQ-TREE 2.3.6 in the manifest. Two runs give the same topology and support;
  branch lengths differ in the fifth decimal (thread summation order) and the test says so.
- `tests/unit/test_phylogeny_module.py` 9 tests, 9 passed 0 skipped with `PTS_REQUIRE_TOOLS=1`.
- `tests/unit/test_tree_scope.py` was replaced, not deleted: it asserted the module's absence and
  now asserts its declaration, artifacts and the honest limits (no outgroup, no rooting, a
  topology is not a verifiable number).
- Also fixed on the way: `naming_reference` was refused by M11.20's new check (a module's
  `reference_param` was not among the accepted keys) — caught by the full suite and corrected.

Evidence: `.dev/records/M11.21_phylogeny_2026-10-28.md`. Open: partition files, MACSE, a
34-genome scale measurement.

## M11.20 — parameters: no silent keys (started 28/10/2026)

Measured first: `--set this_key_does_not_exist=42` exited 0, and `--set unit_repeats=…` on `ssr`
was recorded in the manifest and then ignored (37 SSRs either way). The toolkit had the parameter
infrastructure wired for part of its surface; the rest accepted keys nobody read.

- `registry.ParamSpec` and `ModuleSpec.params`: every module declares the parameters it reads
  (name, type, description), plus `SHARED_PARAMS` for keys the CLI consumes.
- `registry.unknown_params()` + `cli._accept_params()`: `run`, `batch`, `annotate` and
  `mito annotate` now exit **2** on a key no module takes, naming it and listing what would work.
- Fourteen modules declare their parameters; `tests/unit/test_module_params.py` compares each
  declaration with the module's real `params.get` calls via `ast`, in both directions.
- `ssr` thresholds are **refused by name** rather than accepted and ignored. Honouring them means
  generating `misa.ini` from the chosen values, which may move the golden SSR numbers — its own
  decision, not a side effect of this change.

Evidence: `.dev/records/M11.20_parameters_2026-10-28.md`. M11.20 stays `[~]`: typed parameters for
MAFFT's algorithm and BLAST's E-value/max-HSPs, and the `ssr` threshold path, are not done.

## Also in the tree from this Đợt

- **M11.19**: `docs/tools/README.md` is the tool inventory (32 tools, one role each);
  `environment.yml` is runtime only, with `environment-check.yml` and `environment-future.yml` for
  cross-check and deferred tools; `tests/unit/test_tool_inventory.py` (15 tests, 6 negative
  controls) holds the split from both directions; `tests/golden/test_panel_provenance.py` re-hashes
  all 80 panel records. Minimum install measured from a freshly built wheel.
- **M11.9 (metadata half only)**: `Development Status` moved `5 - Production/Stable` → `4 - Beta`.
  The milestone stays `[~]` because a stability statement waits on M11.8, and M11.8 is not passing.

## Checks actually run

- Full fast suite on the final tree (`PTS_REQUIRE_TOOLS=1`, tool env incl. `PTS_IQTREE_DIR`,
  `-m "not tier_b"`): **1136 passed, 1 skipped, 66 deselected, 1 xfailed, 0 failed** (769 s; the
  phylogeny tests add ~5.5 min of real IQ-TREE runs). The skip is the intentional
  machine-without-table2asn case; the **xfail is K-07's strict regression**, expected to fail
  until the defect is fixed. Targeted guards (`test_module_params`, `test_tree_scope`,
  `test_phylogeny_module`, `test_tool_inventory`, `test_licensing`, `test_product_docs`) pass;
  `ruff check .` and `ruff format --check .` clean.
- Not run: the **full golden** (`make golden`, ~12 min), Docker (this WSL distro has no Docker
  integration), a fresh `conda env create`, macOS, and the nine-family batch-vs-single cross-check.

## Đợt 4 verification — all three unrun checks ran (28/10/2026)

| Check | Result |
|---|---|
| Full golden | **126 passed, 0 failed, 0 skipped** (18 min 29 s) |
| Docker | built; in-container run reproduces 153,724 bp / 132 genes / 37 SSR = 31·4·2 |
| Fresh `conda env create` | exit 0; 9 tools incl. iqtree2; wheel ran 8 modules and `phylogeny` with no `--set iqtree=` |

Three defects found by these runs, all fixed: `iqtree` missing from the Dockerfile; no **R** in
the image (every R check skipped); no optional **extras** in the image (every docx/figure check
skipped — 64 passed / 61 skipped). The image now installs R + the nine R packages and
`.[figures,docx,excel,stats]`; the panel half needs `-v "$PWD/.panel:/work/.panel:ro"` and says so.

Phylogeny scale measured: 34 genomes → **10 min · 179 MB**, 78 blocks, 70,091 sites.

K-07 remains open: two more candidate rules measured and rejected (local maximality; inward
growth, 0 bp both directions). On `NC_057683.1` the truncated pair scores 0.999689 against the
true pair's 0.997342, so no ranking by identity or length can pick correctly — it needs a
structural rule.

## Next

1. **K-07**: design the overlapping-candidate rule for de-novo IR (the fix is not a gap-sign
   change — that was measured and rejected).
2. **The full golden** `make check` — the gate this Đợt still owes.
3. The nine-family batch cross-check from Đợt 2, and the CGAS/paper definition comparison.
4. M4.5 (RNA editing) remains open and must appear in any whole-pipeline verdict; it has no module.

Do not treat this as a release acceptance. M11.8 is `[~]`, not `[x]`.
