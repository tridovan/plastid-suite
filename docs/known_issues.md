# Known correctness and integration issues

These issues were reproduced during the 26 October 2026 review. **K-01, K-02, K-03, K-04, K-05
and K-06 are fixed.** Documenting a defect did not fix it, and a fixed entry says which test now
holds it. **K-07 is open** (found by the M11.8 acceptance run on 28/10/2026) and its regression is
a strict `xfail`, so the defect is recorded and cannot quietly disappear; **K-08 is fixed** in the
same run. The code version is `1.0.0`; that version label does not establish readiness for
scientific use. The "observed" columns below are the review's before-fix measurements; they are
kept as the reproduction record, not as current behaviour.

Coordinates below are 0-based half-open unless a GenBank position is explicitly stated. The pinned
fixture is [NC_050833.1.gb](../tests/data/NC_050833.1.gb); its provenance and checksum are in
[PROVENANCE.md](../tests/data/PROVENANCE.md) and [SHA256SUMS](../tests/data/SHA256SUMS).

## K-01 — IR extent and identity depend on record shape — **fixed (M11.3/M11.4)**

The fixture has two 26,214 bp IR copies at `[83127,109341)` and `[127510,153724)`. Keeping the
underlying genome but changing its representation produced the following before-fix results:

| Input change | Before the fix | Now |
|---|---|---|
| Remove repeat annotations; change one base at position 140617 in the second copy | De-novo IR 13,107 bp, identity 100 % | 26,214 bp, one mismatch |
| Remove repeat annotations; rotate right by 5,000 bp | De-novo IR 21,214 bp, identity 100 % | 26,214 bp, identity 100 % |
| Delete position 140617; retain annotations, adjusting the second copy's end | 64.68 % identity, refused | 26,213 matches + 1 internal gap, 99.996 %, passes |
| Rotate right by 1 bp; retain a joined annotation across the origin | one copy read as 26,213 bp | the whole joined copy, 26,214 bp |

The fixes are one grouped representation for a joined `repeat_region` (one feature = one copy),
annotation-first detection with an explicit `--set ir_source=sequence` opt-in, and an end-to-end
anchor-chained alignment that measures the full copies rather than the exact core. A passing gate
is still not proof that the annotated *boundaries* are biologically right — it measures whether the
two deposited copies agree — and the manifest now records the annotated lengths, the alignment
columns and any refused annotation.

Held by `tests/unit/test_known_issues.py` (the K-01 cases), `tests/unit/test_ir.py`,
`tests/unit/test_junctions.py` and `tests/golden/test_e2e.py`.

Sources: [GenBank reader](../src/plastid/core/genbank.py), [IR detection](../src/plastid/core/ir.py),
[normalization](../src/plastid/core/normalize.py).

## K-02 — Cache omits annotation fields — **fixed (M11.12)**

`_input_hash` omitted `transl_table`, `codon_start`, `exception` and `trans_splicing`, so two
annotations differing only in those fields shared a cached result even when an analysis read the
field. The cache key now carries every model field an analysis can read, a digest of the shared
package code and packaged resources (not just the module file), and each entry stores the
provenance of the computation that produced it. On a hit the recorded tool versions are revalidated;
a moved, vanished or upgraded tool makes the payload stale and it is recomputed, never served under
the current run's name. Entries written before schema 2 are not read.

Held by `tests/unit/test_run.py::test_a_changed_translation_table_is_not_served_the_cached_orf_result`,
`test_a_changed_tool_version_invalidates_a_cached_payload`,
`test_shared_package_code_and_resources_move_the_key`, and the per-field
`test_each_annotation_qualifier_moves_the_input_hash`.

The original probe used `_input_hash` and `execute_step` directly: the second ORF run returned
`cached` with the first run's table-1 result. The regression now checks the exported status counts
through the runner and separately asserts that the same input is still a cache hit with the
original provenance. The diagnostic used internal helpers and is not a public API contract.

Source: [runner](../src/plastid/run/runner.py), [cache](../src/plastid/run/cache.py).

## K-03 — Junction reader cannot serve as an independent oracle — **fixed (M11.3)**

Before the fix, rotating left by 83,127 or 127,510 bp made one whole IR start at GenBank position 1
and `junction_regions` raised `JunctionError`, mistaking that copy for a wrapped fragment. Rotating
by 88,127 or 96,234 bp and writing the wrapped copy as a `join` made both parts ≥ 5 kb, so the
reader counted three copies and raised again.

A `repeat_region` feature is now one grouped copy, read whole, and the junction reader uses the
same `RepeatRegion.arc()`. The four rotated cases produce the same regions: LSC 83,127 · IRb/IRa
26,214 · SSC 18,169. The reader is still not an oracle: expected boundaries in tests come from the
independently constructed locations, not from this reader.

Held by `tests/unit/test_known_issues.py` (K-03), `tests/unit/test_junctions.py`.

Source: [junction reader](../src/plastid/core/junctions.py).

## K-04 — Unknown bases are counted as matching evidence — **fixed (M11.4)**

Before the fix, replacing every base in both annotated IR copies with `N` returned identity `1.0`
and `gate_passed=True` over 26,214 positions. Unknown bases are now **ambiguous columns**: they
stay in the identity denominator and are never observed matches. The evidence policy requires at
least half of the shorter copy to be known A/C/G/T on both sides (floor 20 columns), so a tiny known
fragment cannot carry a PASS even if it is 100 % identical.

The all-`N` fixture now returns `matched=0`, `known_columns=0`, `ambiguous_columns=26,214`,
identity 0.0, and the gate refuses. The same evidence is recorded in `manifest.json`.

Held by `tests/unit/test_known_issues.py` (K-04), `tests/unit/test_ir.py`
(`test_identity_is_measured_on_full_copies_not_the_exact_core`).

Source: [pair alignment](../src/plastid/core/ir.py).

## K-05 — Rotation drops partial-coordinate flags — **fixed (M11.13)**

`_shift_interval` reconstructed intervals without carrying either partial flag, including when a
feature did not cross the origin. It now keeps each flag on the coordinate it describes and never
marks the artificial split introduced by rotation. Fragments are returned in transcript order:
high then low on the plus strand, low then high on the minus strand, so `coding_sequence` is
unchanged. The read → rotate → write → read path preserves the model's qualifiers as well.

Held by `tests/unit/test_normalize.py::test_rotation_keeps_partial_flags_on_the_endpoints_they_describe`,
`test_a_split_introduced_by_rotation_is_not_partial`,
`test_rotation_preserves_every_coding_sequence`,
`test_rotation_preserves_a_multipart_trans_spliced_gene_across_the_origin`, and
`test_read_rotate_write_read_keeps_the_annotation_model`.

Source: [rotation](../src/plastid/core/normalize.py).

## K-06 — Tree options are not forwarded by the CLI — **fixed (M11.7 forwarding)**

The packaged tree script accepts `--root`, `--outgroup` and `--drop-tips`, and the caption tells a
user to pass them for another family. `plastid figures` forwarded the tree file and R binary but
did not pass these settings through the adapter's `extra` argument, so `--set root=none` did
nothing.

The CLI now forwards the tree script's drawing options from the run's `--set` values
(`root`, `outgroup`, `drop-tips`, `outgroup-label`, `tag`, `plot-type`, `legend`, `width`,
`height`, `dpi`, `title`), and only to the tree script. The script no longer drops the manuscript's
two
*Asparagus* tips or roots on its Orchidaceae outgroup by default: no tip is dropped and no
outgroup is chosen unless the run names them, and the neutral default is a midpoint root
(`--set root=none` keeps the tree's own root). Held by
`tests/unit/test_known_issues.py::test_k06_the_tree_figure_forwards_its_options`.

Sources: [CLI](../src/plastid/cli.py), [adapter](../src/plastid/external/rscript.py),
[script](../src/plastid/data/r_scripts/R_F9_phylogeny_tree.R).

## K-07 — One indel base inside an IR copy can halve the de-novo repeat and still pass the gate — **open (found in M11.8, 28/10/2026)**

Reproduced while running the acceptance of the synthetic IR panel, not by reading code. The inputs
are a record's own annotated copies plus a single-base edit, so the expected answer is arithmetic:

| Input | Detector chose | Expected |
|---|---|---|
| `NC_028027.1` unedited, annotation removed | 25,195 bp @ 100 % | 25,195 bp (correct) |
| the same sequence, one base deleted from copy B's middle | **12,598 bp @ 100 %, gate passed** | 25,194 bp |
| `NC_008641.1`, one base inserted in copy B | **12,796 bp @ 100 %, gate passed** | 25,591 bp |
| `NC_016708.2`, one base deleted from copy B | **12,765 bp @ 100 %, gate passed** | 25,527 bp |
**Mechanism, measured.** The repeat's 20-mer anchors fall on two mirror diagonals (84,124 and
84,125 for `NC_028027.1`). `_mirror_chains` builds a run on each and then chains runs only where
the gap in *both* copies is between 0 and `IR_ALIGN_GAP_BP`; across the indel one copy's gap is
`-1`, so the two halves never join. Each half is a ~12.6 kb run at 100 % identity, which is above
`_DENOVO_LENGTH_WINDOW[0]` = 10,000 and far above `_DENOVO_MIN_IDENTITY` = 0.95, so the detector
returns it as a valid IR and `enforce_ir_gate` accepts it.

**What a user sees.** A run on such a sequence succeeds (exit 0): `genome_stats` reports `IR (bp)
12598 · IR identity 100.00 · IR source denovo_diag`, `structure` splits LSC/IRa/SSC/IRb at the
wrong coordinates, and every position-based module downstream builds on that. Nothing in the
output is marked uncertain. This is the K-01 family ("a record-shape change moves the IR"),
reached from the de-novo path.

**How often.** Not rare enough to dismiss and not universal: deleting one base at 22 positions
spread through `NC_028027.1`'s copy B gave the full length at 17 positions and a halved length at
5. It needs the indel to land where the two diagonals both carry long runs; the exact rate on
other records is not measured.

**Why it is not fixed in M11.8.** M11.8 is acceptance, and the obvious repair was measured before
being rejected: allowing a small negative gap in the chain fixes `NC_028027.1` but makes
`NC_057683.1` (previously no de-novo result at all) return a 12,838 bp half-repeat that passes, and
moves `NC_068894.1` from 22,652 to 25,121 bp. A correct fix has to choose between overlapping
candidate pairs — i.e. it needs a validity rule for "which pair is the repeat", not a gap-sign
change. That is detector design and belongs to the IR milestone (M11.3/M11.4's successor), not to
this one.

**Two more candidate rules were measured on 28/10/2026 and both failed, so the next attempt does not
have to repeat them.** Each was tested on the collapsing edits and on intact records:

| Candidate rule | What was measured | Why it fails |
|---|---|---|
| The detected pair is **locally maximal**; a non-maximal half is truncated | Extending the half inward (A right, B left) by 1–12,000 bp never matched exactly, and extending it to ≥0.999 identity grew it by 0 bp on every case | The half really is a maximal exact pair in its own neighbourhood; this rule cannot see the missing part |
| A truncated pair **grows inward while identity stays ≥ 0.99**; a true pair meets single-copy sequence and stops | Both directions, steps of 100–6,400 bp, threshold 0.999: **inner = 0, outer = 0 on every case**, intact and collapsing alike | Direction and continuity differ per record: `NC_016708.2`'s half is the *inner* half (its complement is outward), so no single growth direction separates them |

The second measurement is the important one: the failure is not "the code looks the wrong way". On
`NC_057683.1` the 12,838 bp half reaches **0.999689** while the true pair is only **0.997342** —
the truncated candidate is *more* repeat-like than the real one, so any rule that ranks candidates
by identity or length picks the wrong pair. A correct rule has to use the record's structure (two
copies of a repeat leave exactly two single-copy arcs of a consistent circlome), which is the
detector redesign this entry already says it needs.

**Held by** `tests/unit/test_ir_indel_collapse.py`: one strict `xfail` that fails if the defect
disappears without the entry being updated, plus two passing tests that (a) confirm the unedited
record still detects its full IR, so the xfail cannot be passing for an unrelated reason, and
(b) confirm the halved pair is accepted by `detect_ir`, which is the entry point the CLI uses.
The test needs the pinned panel and skips with a named reason without it.

**Not claimed.** The annotation path is unaffected — it never runs this search, and every panel
record with an annotation is measured from its annotated copies. The defect reaches users through
a bare FASTA (the route `plastid annotate` and `plastid run` accept) and through any record whose
annotation is absent or unusable.

Sources: [de-novo detector](../src/plastid/core/ir.py) (`_mirror_chains`, `_chain_geometry`,
`_scan_linear`), [regression](../tests/unit/test_ir_indel_collapse.py).

## K-08 — The Methods page named only the last module's tools — **fixed (M11.8, 28/10/2026)**

Found while running the M11.8 acceptance, not by reading code: a run of ten modules ending in a
pure-Python one reported **"(no external tools ran)"** in `00_run_info/run_info.md` while
`misa 2.1` and `repeat-match 4.0.1` had run and written their logs. Reversing the module order
(`-m ssr -m genome_stats`) reproduced it; `-m ssr` alone printed the MISA row correctly, which is
why a single-module test never saw it.

**Cause.** Every module gets its own `RunContext`, and a context owns the tool records of its
module. `_finish` rendered the page from the last context's records, so any module before the last
one lost its rows. The manifest was right all along: `modules[*].tools` recorded MISA and
`repeat-match` in the failing run.

**Fix.** `render_run_info` now reads the tool records from the manifest — every module's — and
additionally accepts the live records, de-duplicating on `(name, version, path)`. A tool called by
five comparison modules appears once, because the table answers "what ran", not "how many times".
A run that genuinely called no tool still says so, which is the honest negative this page has to
keep.

**Why it matters beyond cosmetics.** That page is the Materials & Methods record; an unqualified
"no external tools ran" is a positive false statement about which programs produced the numbers.

Held by `tests/unit/test_run_info_tools.py` — a red-before/green-after unit case, the
de-duplication case, the genuine-negative case, and an end-to-end CLI reproduction in the order
that failed (it skips, with a named reason, where `misa.pl` is absent).

Sources: [page renderer](../src/plastid/run/runinfo.py), [regression](../tests/unit/test_run_info_tools.py).

## Other current limitations

The RSCU and amino-acid R figures draw every species column in their input (M11.7); the tool's
`figures` command still supplies one column per run because `codon_usage` and `amino_acid` take one
genome each, so a combined species panel needs the batch path (M11.6). The reference preparation
and downstream analysis steps are not a single integrated workflow, though the pieces now chain
with each other's real artifacts (`tests/golden/test_reference_flow.py`).

Unequal IR copies are legal and are measured end to end. Mapping a gene between them is now
implemented (M11.11): `plastid_annotate` places the second copy through the gate's own alignment
(`IRResult.mapping`), including the one-base indels the equal-copy constant cannot express. A
result that carries no alignment refuses such a record by name (exit 2) before any artifact is
written. `normalize` still records `mirror_constant=0` for an unequal pair because no equal-copy
constant exists; consumers that need a position map use the alignment, not that provenance field.
Held by `tests/unit/test_ir_copies.py` and `tests/golden/test_ir_panel_mapping.py` (32 unequal
panel pairs, 267 of 272 mapped genes within 10 bp of the record's own second copy; the five
exceptions are deposited-annotation disagreements named in the test). See
[current scope and limits](things_we_do_not_do.md) for details.

## H-21 — The default reference set blocked its own usable record — **fixed (M11.10)**

Reproduced 26/10/2026: `plastid reference audit` on the 34-record set blocked `NC_024936.1` three
times — its two `rps19` CDS copies carry `>83020`/`<153124` over a **complete 92 aa translation**,
and its annotation-first IR measures 99.80 %. The partial rule asked "is any coordinate partial?"
instead of "is the model truncated?", and the IR check stopped at the first path even though the
sequence's own repeat measures 26,184 bp at 100 %. The record is the one that produced a
base-perfect `rps19` call at 83,173–83,451.

Now: a partial CDS whose own `/translation` fits its own span is **reported**; one whose
translation does not fit still **blocks**. The IR check reports the annotation-first refusal and
probes the sequence path as a **named alternative** (`--set ir_source=sequence`, recorded in the
manifest and in `audit_records.ir_source`), blocking only when both fail. The default set audits
with 0 blocked, the annotator runs on it without `--exclude`, and the golden's 35 pinned
differences are unchanged. Held by `tests/unit/test_reference_audit.py`,
`tests/golden/test_plastid_annotate.py` and `tests/golden/test_reference_flow.py`.

## H-22 / H-23 — Intron coverage and the caller's range — **fixed (M11.15)**

H-23: the `introns` table published CDS introns only. It now publishes CDS and tRNA introns with a
`kind` column, counted separately, and checks a class it cannot see as **"not measured"** rather
than printing zero. The golden's 8 tRNA introns are checked against the record's own multi-part
tRNAs (the published count), not against ARAGORN's log line: 17 CDS + 8 tRNA = 25, total 20,347 bp.

H-22: ARAGORN's intron search was passed as a bare `-i` with no stated range. Measured
26/10/2026 on 1.2.41 with a lengthened golden `trnK-UUU`: a 3,000 bp intron is reported, 3,001 bp is
not, and `-i4000`, `-i0,4000` and `-i0,6000` do not change that. The adapter declares the range on
every result, refuses a wider request by name, and the annotator's notes say a tRNA outside the
range was not assessed. A malformed tRNA join (two PGA variants overlap by 5 bp) is named as not
assessed instead of aborting the record's other spans. Held by `tests/unit/test_aragorn.py`,
`tests/unit/test_spans.py`, `tests/unit/test_introns_module.py`,
`tests/golden/test_intron_limits.py` and `tests/golden/test_intron_assessment_panel.py`.

## G-03 — Output commands did not collect one catalogue — **fixed (M11.18)**

Before the fix, `run --figures`, `package` and `qc` called only the matplotlib renderer while
`figures` also ran the manuscript's R scripts. A package could therefore omit an R figure that
the same run directory could draw, and `qc` reported the files it happened to find rather than
the figures that were missing and why.

Now all four commands use one figure pass over both catalogues. `run --figures`, `figures` and
`package` draw the matplotlib figures and the R figures and record one status per catalogue entry
(rendered, or skipped with a reason, plus a whole-renderer status when the plotting extra or R is
absent). `qc` inventories both catalogues without drawing. `package` replaces its folder in one
step, so a table or figure the run no longer holds cannot survive a rebuild, and its index links
every file with its sha256, provenance and the run's data label. The Excel workbook gained a
`plastid.scope` sheet so its numbers carry the same scope as the CSV, JSON and Word exports.

The tool policy per command is explicit: `run --figures` exits 3 when the plotting extra is
missing, because the user asked for figures; `figures` and `package` still draw the R half and
record the matplotlib absence; a missing R is always a recorded skip (the toolkit's own figures
do not depend on it); `qc` needs neither and reports every catalogue entry's status.

`submission` no longer reports a validator that did not run as `0` errors: the counts are blank
and `validator_status` says `not run`; a `.sqn` is described as a conversion, not a validation.
A record with no gene/CDS/tRNA/rRNA features is refused before any file is written — a submission
comes from a checked GenBank annotation, not from a bare sequence or an analysis table.

Held by `tests/unit/test_output_contract.py`, `tests/unit/test_exports.py`,
`tests/unit/test_batch_figures.py`, `tests/unit/test_submission.py` and
`tests/golden/test_e2e.py`.

## Pipeline coverage observations — code inspection, not new runtime reproductions

The [canonical pipeline](pipeline.md) maps every stage to existing components and missing work.
The following require targeted regression/acceptance checks in the next coding sessions:

- `gene_normalize` verify mode returns findings without writing a replacement GenBank. A workflow
  must explicitly decide how differences block or qualify downstream work. Normalized DB artifacts
  now carry their own digest and the run records the raw input's digest and the chosen naming
  reference (`tests/golden/test_reference_flow.py`).
- Alignment coordinates and normalized-reference coordinates in diversity/barcodes still need
  rotation tests on the panel.

Sources for the observations above: [naming](../src/plastid/modules/gene_normalize.py),
[group measurements](../src/plastid/modules/stats.py),
[diversity](../src/plastid/modules/diversity.py), [barcodes](../src/plastid/modules/barcodes.py),
[CLI output commands](../src/plastid/cli.py),
[submission](../src/plastid/modules/submission.py).

## G-04 — a tool in an environment is not a feature — **closed (M11.19)**

`environment.yml` mixed the tools an analysis runs with assemblers, a dating program and engines
that had only been surveyed. The file is now the runtime scope only, and every tool the toolkit
touches has **one role** in [the tool inventory](tools/README.md) — `active-runtime`,
`cross-check`, `future` or `survey-only` — which decides its file:

| File | Contents |
|---|---|
| `environment.yml` | runtime: MAFFT, BLAST+, ARAGORN, MUMmer (`nucmer` and `repeat-match`), PAML `yn00`, R, perl |
| `environment-check.yml` | `cross-check`: Vmatch, CPStools, KaKs_Calculator — no module calls them |
| `environment-future.yml` | `future` (assembly, M7.x) and `survey-only` (Sibelia, BEAST2, mreps, TRF, …) |

`tests/unit/test_tool_inventory.py` holds the split from four directions: an adapter with no row
fails, an `active-runtime` adapter no product file imports fails, a `cross-check` row whose adapter
a product module *does* import fails, and a conda package that no row claims — or that sits in a
file for a different role — fails. Its negative controls fabricate each failure, so the check is
known to bite. Nothing was deleted: a deferred tool keeps its name and its reason, because a later
scope decision must not look like a new discovery.

Measured with a fresh venv, `pip install plastid_suite-1.0.0-py3-none-any.whl`: the minimum
install is **biopython + PyYAML (+ numpy, Biopython's own dependency)**, seven pure-Python modules
run on the golden record in one command, and each missing piece is named with its exit code —
`ssr` without `misa.pl` → exit 3 naming the download line, `--docx` → exit 3 naming the `docx`
extra, `--figures` → exit 3 naming the `figures` extra. `figures` (without `--figures`) still
draws the R half and records the matplotlib absence rather than failing. The panel's own integrity
is now checked too: `tests/golden/test_panel_provenance.py` re-hashes all 80 source records,
matches the 170 manifest rows against the files on disk, and refuses a panel that calls itself
golden.

Sources: [inventory](tools/README.md), [runtime file](../environment.yml),
[cross-check file](../environment-check.yml), [deferred file](../environment-future.yml),
[inventory check](../tests/unit/test_tool_inventory.py),
[panel provenance](../tests/golden/test_panel_provenance.py).
