# Development handoff

Updated: 2026-10-26. Branch: `main`; reviewed HEAD at session start: `d1ec471`.
This is the sole current development handoff, not a user guide or a record of scientific decisions.
There is no second HANDOFF in the maintainer workspace; earlier versions remain in Git history.

## Current work

The coding session completed **M11.2, M11.3, M11.4, M11.12, M11.13 and M11.14**. M11.3 is closed
with PP078674.1 pinned at `.panel/ir-validation/pp078674/` (LSC 83,820 · IRa 26,331 · SSC 17,195 ·
IRb 26,330) and panel length match 80/80. M11.4 is closed: panel before/after (80 records), scoring
policy, and a final full golden of **86 passed, 0 skipped**. **M11.5 is not started.** The working
tree carries the changes; nothing was committed, pushed or published.

- **M11.2/M11.12/M11.13/M11.14:** regression shapes, cache semantics, rotation annotation
  preservation and the product-documentation checks. K-02 and K-05 are fixed and tested where
  their code lives. Baseline: `.dev/records/M11.2_baseline.md` (private planning).
- **M11.3:** one `RepeatRegion` per GenBank feature, read whole through read/write/rotate. A
  `join()` across the origin is one copy, not several candidates; `core/ir.py` and
  `core/junctions.py` share `RepeatRegion.arc()`. All five wrapped PGA records in the pinned
  corpus now match their annotated lengths; PP078674.1 is pinned locally and matches the DoD
  regions.
- **M11.4:** `align_copies()` measures copy A against the reverse complement of copy B end to end
  with 20-mer anchors and global gap alignment. Unknown bases are ambiguous columns, never matches;
  the gate requires identity ≥ 0.999 and at least half the shorter copy known. Detection is
  annotation-first; a failed gate refuses and names `--set ir_source=sequence`, which is the only
  path that substitutes the sequence's repeat and records the annotation it did not use. De-novo
  detection is circular-aware and chains anchors through SNPs/small indels. Evidence:
  `.dev/records/M11.3_M11.4_baseline.md`.
- **Panel-driven correction (26/10/2026):** the 80-record panel outside Iridaceae found that exact
  runs on either side of a small indel can overlap. The old chain skipped one run and measured
  73–93 % on four records where a global alignment gives 99.73–99.96 %. Anchors may now overlap and
  are trimmed at the next anchor's start; a constructed regression holds the shape. Before/after on
  the same 80 records: HEAD `d1ec471` read 79 records (70 pass · 9 fail · 1 error); the fixed tree
  reads 80 (75 pass · 5 fail). Against the corrected independent circular-arc measurement:
  74 pass + 5 fail agree, 1 case (NC_035154.1) differs only because of gap-scoring (99.9462 % PTS
  vs 99.8923 % cheap-gap parasail) and is flagged for a scoring-policy decision. Two earlier
  "independent" failures on joined records were measured in join-string order, not circular order;
  corrected, they match PTS.

**Scientific baseline changes that need user review before being treated as accepted:**

- NC_063920.1 annotated IR: 92.63 % (refused) → **99.94 % (passes)**. The old value was the
  ungapped mirror pairing; the new alignment reports 10 mismatches + 5 internal gaps over 26,164
  columns.
- *Iris hermona* PGA: 99.86 % (refused) → **99.992 % (passes)**; 0 mismatches + 2 internal gaps.
- NC_024936.1 annotation: 26,220/26,186 bp, **99.80 %, refused by default**; the explicit
  `ir_source=sequence` path gives 26,184 bp at 100 %. This is annotation/sequence disagreement,
  not automatic divergence.
- NC_085783.1 (no annotation): de-novo length 21,148 → **26,184 bp**; the two copies are exact
  reverse complements over that extent and +1 bp breaks the equality.
- Golden NC_050833.1 values are unchanged (26,214 · 100 %).

**Active task:** M11.11 is specified below but **not implemented**. M11.5, M11.10, M11.15 and
M11.7 are not started. Plastid RNA editing (M4.5) remains open and unimplemented; the selected
engine is original PREP-Cp with cutoff 0.8, and no predictor may be substituted. Assembly remains
deferred.

## Next coding session

Read in this order:

1. [CLAUDE.md](CLAUDE.md) for code rules, then [pipeline](docs/pipeline.md) for the product contract.
2. If this maintainer workspace has `.dev/`, open `.dev/TASKS.md` for the active task, then
   `.dev/ROADMAP.md` for the milestone DoD. Read `.dev/design/IR_rework_M11.md` and
   `.dev/records/M11.3_M11.4_baseline.md` before changing IR/alignment behaviour. These private
   files are not required in other clones.
3. Read only the design notes named by the active task.

Order:

1. **M11.11:** map a gene between unequal IR copies through the M11.4 alignment, or refuse by name
   everywhere the equal-copy assumption remains. Do not lower the gate.
2. **M11.5:** make codon/amino-acid analyses depend only on valid CDS, while position-based
   modules stay behind the IR gate.
3. **M11.10/M11.15:** reference flow from a fresh checkout; intron caller coverage.
4. **M11.7:** forward tree root/outgroup/drop-tips to R (K-06 strict xfail will then XPASS).
5. A diverse IR validation panel (outside Iridaceae) is being collected separately; place large
   data under `.panel/ir-validation/`, pin it with `plastid panel pin`, record provenance, and
   classify every record as PASS / refused-by-data / refused-by-algorithm. Do not tune 0.999.

### Active task — M11.11 (annotator output usable downstream)

**Goal:** the annotation artifact is usable by downstream steps, and every place that still
assumes two equal IR copies maps through the M11.4 alignment or refuses by name **before** writing
the artifact.

**Current state:** `plastid_annotate` raises `IRGateError` (exit 2) when `ira.length != irb.length`;
`normalize` records `mirror_constant=0`; gene mapping between unequal copies is not implemented.
The validation panel has 32/80 unequal pairs (1–12 bp), 30 PGA variants and 6 joined copies.
H-24 (`trnV-UAC` missing) must be reproduced before it is fixed.

**Steps:**
1. Failing regression first: one unequal-copy panel record, one boundary gene and `rps12`; reproduce
   H-24 with an answer taken from the locus/annotation, not from the annotator under repair.
2. Map intervals/positions between copies through the M11.4 alignment. If that is not reliable,
   keep the named refusal, but it must happen before the artifact and the reason must be recorded.
   Never lower the IR gate.
3. Replace every remaining equal-copy assumption in `_ir_copies`, `_mirrored`,
   `normalize.mirror_constant` and its consumers with mapping or a named refusal.
4. Fix H-24 from locus evidence, with a negative case and at least one non-Iridaceae family from
   the panel.
5. The artifact must re-read as GenBank with the features/qualifiers downstream needs; gene
   features must not inflate the CDS + tRNA + rRNA count; check exon/kind/coordinates.

**Acceptance:** unequal-copy mapping (or refusal) tested on the panel records; `rps12` and boundary
genes preserved; H-24 positive and negative; artifact re-read; manifest records mapping/refusal;
golden stays at 86 passed unless a baseline change has independent evidence and user review.

## Verification and remaining work

Checks run in this session (final counts are in the session report; see also
`.dev/records/M11.2_baseline.md` and `.dev/records/M11.3_M11.4_baseline.md`):

| Command | Result |
|---|---|
| `uv run pytest -q -m "not tier_b"` | **959 passed**, 3 skipped, 42 deselected, **1 xfailed** (K-06) |
| `uv run pytest -q tests/unit/test_ir.py tests/unit/test_known_issues.py` after the low-complexity guard | 47 passed, 1 xfailed |
| `uv run pytest -q tests/unit/test_gene_normalize.py` with pinned BLAST | 31 passed |
| `uv run pytest -q tests/golden -rs` with `PTS_*` and `PTS_REQUIRE_TOOLS=1`, after the panel-derived anchor-overlap fix | **86 passed, 0 skipped** (8:18) |
| `uv run pytest -q tests/unit/test_ir.py tests/unit/test_known_issues.py tests/unit/test_normalize.py tests/unit/test_junctions.py tests/unit/test_ir_copies.py` after that fix | 79 passed, 1 xfailed |
| Panel: 80 RefSeq records, PTS vs independent circular-arc parasail | 74 pass both · 5 fail both · 1 scoring-only disagreement (NC_035154.1) |
| `uv run python tools/validate_ir_panel.py` | 80 source records · 75 pass · 5 fail · length match 80/80 · PP078674.1 regions match |
| `uv run ruff check src tests` / `ruff format --check src tests` | pass |
| `.dev/tools/ir_truth.py` / `ir_fixed_reading.py` on pinned panels | 5 wrapped PGA records read whole; PTS length ≠ annotation: 0 |

The fast suite predates the low-complexity guard, the `gene_normalize` `ir_source` threading and
the anchor-overlap fix; each has targeted runs above. Known skips are tool/panel skips, not passes.
K-01, K-03 and K-04 are fixed in this working tree; K-06 remains open (M11.7). The unequal-copy
mapping limitation is recorded in [known issues](docs/known_issues.md). No full acceptance run for
the whole pipeline was performed.
