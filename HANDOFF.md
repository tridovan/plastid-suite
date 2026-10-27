# Development handoff

Updated: 2026-10-27. Branch: `main`; HEAD at this session's start: `5735220`.
This is the sole current development handoff, not a user guide or a record of scientific decisions.
There is no second HANDOFF in the maintainer workspace; earlier versions remain in Git history.

## Current work

The previous session closed **M11.2, M11.3, M11.4, M11.12, M11.13 and M11.14**. This session
implements **M11.11**, **M11.10**, **M11.15** and the rest of **M11.7** in the working tree;
nothing was committed, pushed or published. On 2026-10-27 the failures left by that batch's full
golden were resolved (see the verification section).

- **Reference flow and H-21 (M11.10).** A partial CDS whose own `/translation` fits its span is
  **reported**, not blocked; only a model that is genuinely cut blocks. An annotation-first IR
  failure is reported while the sequence path is probed as a named alternative
  (`--set ir_source=sequence`), blocking only when both fail; `RecordFacts.ir_source` and
  `audit_records.ir_source` carry the path. The default 34-record set now audits with **0 blocked
  / 25 reported**, and the golden annotates without `--exclude NC_024936.1`, still 86 + 38 + 8 with
  `rps19` at 83,173–83,451 and the same 35 pinned differences. `pin_files` refuses different bytes
  for an accession already pinned. `tests/golden/test_reference_flow.py` runs pin → audit → choose
  the naming reference → normalize → annotate → verify with each step's real artifact, and
  `tests/golden/test_reference_flow_panel.py` classifies the 80-record/10-family panel (8 blocked
  records, 10 findings) and the 30 PGA variants (all blocked on missing `/translation`), then
  annotates a Solanaceae target from 15 family references.
- **Introns and the caller's range (M11.15).** The `introns` module (0.2.0) publishes CDS and
  tRNA introns with a `kind` column and separate counts, reports "not measured" when the record
  annotates no feature of a class, and names an unclassifiable tRNA gap "not assessed". Golden:
  17 CDS + 8 tRNA = 25 introns, 20,347 bp; the 8 tRNA introns are the record's own multi-part
  tRNAs, not ARAGORN's log. ARAGORN's range is declared as `(0, 3,000)` on every result; a wider
  request is refused by name (measured: 3,000 bp found, 3,001 not, `-i0,4000` no help), and the
  annotator's notes say a locus outside the range was not assessed. Panel: 5 gate failures stop
  by name, 576 tRNA introns across 75 records; the 30 PGA variants give 213 tRNA introns and two
  malformed joins named as not assessed.
- **Figures outside one paper (M11.7).** F3/F4 read the species columns from the input, in input
  order, with `--species-order` and `--exclude-species`; the drawing width scales with the panel
  and the raster DPI is capped. The 80-species panel draws through the packaged scripts, and a
  species outside the paper's list is no longer skipped. F9 writes `<tag>_used.treefile` — the
  tree as drawn — and the tests check the dropped tips, the outgroup root, `--root none` and argv.

- **Unequal-copy mapping (M11.11).** `IRResult.mapping` carries the D11-3 alignment's columns;
  `plastid_annotate` places every exon of the second IR copy through it instead of the equal-copy
  mirror arithmetic. A result without a mapping refuses an unequal pair by name before the
  artifact. The manifest records `mapping_supported` and `mapping_method`. The annotator now
  works in the IR's frame (`frame_offset`) and rotates the finished record back, so a copy that
  crosses the origin is written as a join in transcript order.
- **H-24.** An intron tRNA is named at its locus from the reference set's own tRNA annotation
  (`blastn` of its spliced sequences; strict majority by accession). The golden's `trnC-ACA`
  becomes `trnV-UAC` (32/32 supporting records), `trnE-UUC` becomes `trnI-GAU` (30/32) on both
  copies, and the Crataegus `tRNA-Stop` line is recovered as `trnV-UAC` (4/4) rather than
  dropped. Without reference support, a tie, or a set with no tRNA annotation, ARAGORN's name
  stands.
- **Custom reference directories.** The recorded manifest's allow-list now applies only to the
  recorded directory; any other directory is searched leave-one-out by accession. Before this,
  every record of a user-supplied family set was excluded as "leave-one-out".
- **Modules depend only on the data they read (M11.5).** `codon_usage` and `amino_acid` no longer
  require the IR: they count CDS sequences, so they run on a record whose repeat fails the gate, is
  linear, reduced or N-rich. A CDS whose `/codon_start` is not 1 or whose `/transl_table` is
  neither 1 nor 11 is refused by name (exit 2) rather than translated silently. Position-based
  modules keep the gate. Panel check: 5 gate-fail, 31 ambiguous, 3 linear, 2 reduced records.
- **K-06 (part of M11.7).** `plastid figures` forwards the tree script's `--set
  root/outgroup/drop-tips` (and its drawing options) to R, and the packaged script no longer
  defaults to the manuscript's *Asparagus* tips or Orchidaceae outgroup; the neutral default is a
  midpoint root. The strict xfail is removed and the check now passes.
- **Product language.** The packaged R scripts and the product documentation are being translated
  to English (in progress in this working tree); `.dev/` planning stays private and is not part of
  the product.
- **Evidence:** `.dev/records/M11.11_baseline.md` (private planning) records the panel mapping
  (32 unequal pairs; 272 genes; 267 within 10 bp of the record's own second copy), the golden and
  Crataegus annotate runs, and the three golden differences whose class moved.
- **Earlier session, still current:** the panel-derived anchor-overlap fix, the low-complexity
  index guard and the annotation-preserving rotation are unchanged.

- **psbL 5' start (2026-10-27).** The plastid caller now reads **one** codon before an alignment
  that stops short of the reference's 5' end, and keeps it when it is an accepted initiator
  (`core/caller.py: read_upstream_start`, enabled by `core/plastid_caller.py`; the mito prototype
  keeps its old behaviour). Measured on 39 Crataegus records: psbL was written as 114 bp with no
  start (`no_start`, partial 5'); the codon before the alignment is `ACG`, so the call is now
  117 bp with a complete ORF. Golden annotate is unchanged (10 passed, 35 differences and the two
  partial starts); 76 caller/ORF tests pass. A post-fix re-run of the 39 records raises exact
  agreement 62.9 % → 63.2 % (the 12 records whose NCBI annotation carries psbL at 117 bp); the 27
  records whose NCBI annotation omits psbL remain extras but are now complete ORFs.
- **Crataegus accuracy evaluation (2026-10-27).** The 100-record test DB was screened (50 barcode
  amplicons dropped, 11 plastomes blocked by the audit, 39 kept), re-annotated leave-one-out and
  compared with the NCBI annotation: locus recall 99.9 %, exact coordinates 62.9 % (within 1 bp
  81.4 %), 265 extras and 7 missing. Removing three reference records that carry extended small-CDS
  models raises CDS exact agreement 76.7 % → 85.7 % (within 1 bp 89.9 %), so the largest coordinate
  errors are inherited from the reference set, not from gene finding. The evidence and the psbL fix
  are in `.dev/records/crataegus_accuracy_2026-10-27.md`.

**Scientific baseline changes accepted by the maintainer on 2026-10-27** (reviewed in session; the
evidence below is the record of that review):

- The golden's difference set stays at **35**, but three tRNA rows move from `vocabulary/name` to
  `conclusion/extent`: `trnV-UAC` (was `trnC-ACA`) and `trnI-GAU` at both copies (were
  `trnE-UUC`). The names now agree with the published record; the remaining extent difference is
  ARAGORN's intron boundary (H-22/H-23, M11.15).
- From the previous session, covered by the same review: NC_063920.1 annotated IR 92.63 %
  (refused) → 99.94 % (passes); *Iris hermona* PGA 99.86 % → 99.992 %; NC_024936.1 annotation
  99.80 % refused by default with `ir_source=sequence` giving 26,184 bp at 100 %. Golden
  NC_050833.1 values are unchanged (26,214 · 100 %).
- **Intron count, golden:** 17 → **25** (17 CDS + 8 tRNA). The 8 tRNA introns are the published
  record's own two-part tRNA features and its own exon gaps; the CDS set `diversity` reads is
  unchanged at 17.
- **Default reference audit:** 3 blocked findings on NC_024936.1 → **0** (H-21). The record's
  partial flags sit over a complete 92 aa protein, and its annotated IR (99.80 %) has a named
  sequence-path alternative at 100 %; the report keeps both. Golden calls are unchanged.
- **Figure species selection:** F3/F4 no longer skip a species outside the manuscript's 34-name
  list; they draw the input's columns. The old named-skip tests were replaced with tests that a
  new species draws.

**Not started:** M11.6 (batch join for multi-sample figures — the scripts accept a multi-column
table, the `figures` command still supplies one run's column), M11.16 and later. M11.10, M11.15,
M11.5 and the rest of M11.7 are implemented in this working tree. The Crataegus evaluation opened
two further items that are not yet milestones: a majority/structure guard for the reference-model
arbitration (the 284-error source) and the `MW201730` rps16 exon-loss regression. Plastid RNA
editing (M4.5) remains open and unimplemented; the selected engine is original PREP-Cp with cutoff
0.8, and no predictor may be substituted. Assembly remains deferred.

## Next coding session

Read in this order:

1. [CLAUDE.md](CLAUDE.md) for code rules, then [pipeline](docs/pipeline.md) for the product contract.
2. If this maintainer workspace has `.dev/`, open `.dev/TASKS.md` for the active task, then
   `.dev/ROADMAP.md` for the milestone DoD. Read `.dev/design/IR_rework_M11.md` and the baseline
   record named by the task before changing IR/alignment behaviour. These private files are not
   required in other clones.
3. Read only the design notes named by the active task.

Order:

1. **M11.10/M11.15:** reference flow from a fresh checkout; intron caller coverage (H-22/H-23).
2. **M11.7 remainder:** F3/F4 species/order from the input; multi-sample figures on the 80-record
   panel; check the rendered tree and the caption against the run's actual parameters.
3. Close this batch: one full `make check` (fast + golden) after the 2026-10-27 fixes, and the
   user review of the three tRNA-name baseline rows and the other scientific changes below.
4. A diverse IR validation panel (outside Iridaceae) is being collected separately; place large
   data under `.panel/ir-validation/`, pin it with `plastid panel pin`, record provenance, and
   classify every record as PASS / refused-by-data / refused-by-algorithm. Do not tune 0.999.

## Verification and remaining work

Checks run in this session:

| Command | Result |
|---|---|
| `pytest -q` on the M11.11 affected unit files (IR, IR copies, normalize, junctions, ARAGORN, tRNA naming, annotate module, known issues, architecture) | 111 passed, 1 xfailed (K-06, before its fix) |
| `pytest -q tests/golden -rs` with pinned BLAST+/ARAGORN/MAFFT/PAML/MUMmer/VMATCH/table2asn and the pinned panels, `PTS_REQUIRE_TOOLS=1` | **89 passed, 0 skipped** (11:35; 86 existing + 3 new panel-mapping tests) — M11.11 close-out |
| `pytest -q tests/golden/test_plastid_annotate.py` (also run alone) | 9 passed (100 s) |
| `pytest -q tests/golden/test_ir_panel_mapping.py` (also run alone) | 3 passed (68 s) |
| `pytest -q` on the M11.5 files + substitutions + tier_a + the 9-test module/IR panel | 90 passed |
| `pytest -q` on the final changed surface after the language/R/flag edits (IR, annotate, ARAGORN, tRNA naming, K-06, codon/amino, substitutions, rscript, product docs, gene names, architecture, tier_a) | 180 passed, 1 skipped (the private R bundle is absent) |
| `ruff check src tests` / `ruff format --check src tests` | pass |
| `Rscript -e 'invisible(parse(file=…))'` on all 8 packaged R scripts | OK |
| `plastid run -m plastid_annotate` on the golden (pinned set, `NC_024936.1` excluded) | exit 0; 86 CDS + 38 tRNA + 8 rRNA; artifact re-read; `trnV-UAC`/`trnI-GAU` named |
| `plastid run -m plastid_annotate` on `NC_065154.1` (5 Brassicaceae references) | exit 0; 88 + 39 + 8; artifact re-read |
| `plastid run -m plastid_annotate` on `PX413282.1` (4 Crataegus references) | exit 0; 85 + 37 + 8; `tRNA-Stop` recovered as `trnV-UAC` |
| `plastid run -m codon_usage -m amino_acid` on `NC_047256.1` (IR gate fails) | exit 0 |
| `plastid run -m genome_stats` on the same record | exit 2 (positional analysis blocked) |
| **This batch:** `make test` (`-m "not tier_b"`) with the pinned tools | **1015 passed, 3 failed, 1 skipped** (511 s). All three failures were the pre-H-21 corpus/CLI expectations; those tests were updated and pass alone (9 passed, 1 skipped). |
| **This batch:** `make golden` with all pinned tools and `PTS_REQUIRE_TOOLS=1` | **116 passed, 5 failed** (1229 s). Two failures were the corpus tests collected before their H-21 expectation edit (now pass alone); one is a wrong assertion in the new reference-flow test (fixed after the run); two are open — see below. |
| New targeted evidence this batch | intron limits (4), intron panel (4), reference flow (4), reference-flow panel (3), multi-sample R (4), tree output (3), golden annotate (10), audit unit (32), spans/introns/ARAGORN unit (36), tier A (19) — all passed before the full run. |
| `ruff check src tests` / `ruff format --check src tests` | pass |
| **2026-10-27 fix batch:** `pytest -q tests/golden/test_ir_panel_mapping.py tests/golden/test_reference_audit_corpus.py tests/golden/test_reference_flow.py tests/unit/test_reference_audit.py` with all pinned tools and `PTS_REQUIRE_TOOLS=1` | **43 passed, 0 skipped** (235 s) |
| `pytest -q tests/unit/test_reference_audit.py` after adding the `_ssc_rows` regression test | **33 passed** (58 s) |
| `pytest -q -m "not tier_b"` (full fast suite) after the fixes | **1019 passed, 1 skipped, 0 failed** (490 s); the skip is `test_submission.py`'s deliberate no-table2asn-machine test, and table2asn is installed here |
| `ruff check .` / `ruff format --check .` | pass (262 files formatted) |
| **2026-10-27 psbL fix:** caller/ORF unit tests (`test_plastid_caller.py`, `test_mito_caller.py`, `test_orf.py`) | **76 passed** |
| `tests/golden/test_plastid_annotate.py` after the psbL fix, all pinned tools and `PTS_REQUIRE_TOOLS=1` | **10 passed** (102 s); the 35 differences and the two partial starts are unchanged |
| Post-fix Crataegus re-run (39/39 exit 0) | all 39 psbL calls 117 bp, complete ORF; exact agreement 62.9 % → 63.2 %, within 1 bp 81.4 % → 81.6 %; extras unchanged at 265 |

**Fixed after the full golden, 2026-10-27 (this session):**

1. `reference_audit._ssc_rows` now reads each record's SSC region on the path its own audit
   recorded (`RecordFacts.ir_source`), not annotation-first unconditionally. `NC_024936.1` is
   reported on the sequence path (its annotated pair is 99.80 %), so the second reading gets the
   same 100 % repeat its orientation was measured from instead of raising. A record with no usable
   path is refused by name. A unit regression test builds the shape (annotated IRb 1,000 bp too
   long at 24.26 %, sequence repeat 100 %) and the corpus MAFFT test passes.
2. The panel-mapping distribution did **not** change because of a mapping fix. The 236 exact /
   20 one-off figure came from an ad-hoc probe that compared the mapped copy (IR frame) with
   un-rotated deposited genes; on the two panel records whose `frame_offset` is 1
   (`NC_058280.1`, `NC_058538.1`) that put 9 + 9 genes one base off. Re-running the probe's own
   comparison reproduces `{0: 236, 1: 20, …}` exactly. The test, which rotates both sides, has
   measured **254 exact / 2 one-off / 7 at 4–10 bp / 2 named record disagreements** since its
   first run. The test now pins that whole histogram, and the baseline record carries the
   explanation.
3. `test_reference_flow.py` asserts `reference["accession"] == "NC_050833.1"` (the stem, not the
   file name); the digest assertion was already correct.
4. **psbL 5' start.** `core/caller.py` reads one codon before an alignment that stops short of the
   reference's 5' end, keeping it when it is an accepted initiator; `core/plastid_caller.py`
   enables it. All 39 Crataegus psbL calls are now 117 bp with a complete ORF; the golden is
   unchanged. The mito path keeps the old behavior. Unit tests cover plus strand, minus strand,
   a non-initiator and the mito default.

The fast suite was re-run after the first fixes and is clean in one pass (**1019 passed,
1 skipped**; the skip is the deliberate no-table2asn-machine test and table2asn is installed
here). The full golden was not re-run: TASKS prescribes the three affected golden files, which
pass with every tool and panel present, and no other golden test calls `_ssc_rows`; `_ir_copies`
itself was not changed. The psbL change was verified on the golden annotate file (10 passed); the
annotator-dependent golden files still to re-run before closing the batch are `test_e2e.py`,
`test_tier_b.py` and `test_corpus.py`. A full `make check` remains for the next large close, not
this fix batch.

**Still open:** M11.6 (batch join for multi-sample figures), M11.16–M11.19, M11.8, M11.9 and
plastid RNA editing (M4.5) are not started; assembly is deferred. The Crataegus evaluation added
two items not yet in the roadmap: a majority/structure guard for the reference-model arbitration
and the `MW201730` rps16 exon-loss regression.
