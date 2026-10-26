# Known correctness and integration issues

These issues were reproduced during the 26 October 2026 review. **K-01, K-02, K-03, K-04 and K-05
are fixed**; K-06 is **open**. Documenting a defect did not fix it, and a fixed entry says which
test now holds it. The code version is `1.0.0`; that version label does not establish readiness for
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

## K-06 — Tree options are not forwarded by the CLI

The packaged tree script accepts `--root`, `--outgroup` and `--drop-tips`. The caption suggests
using these options for another family. However, `plastid figures` forwards the tree file and R
binary but does not pass these settings to the adapter's `extra` argument. Supplying
`--set root=none`, for example, does not override the script's rooting behavior.

The script defaults drop two Asparagus tips and select Dendrobium/Phalaenopsis as outgroups, with
a midpoint fallback. Those are choices from a specific analysis, not general toolkit defaults.
Until the CLI integration is repaired, run the script directly with explicit options when those
choices must be controlled; see [the R adapter documentation](tools/r.md).

Sources: [CLI](../src/plastid/cli.py), [adapter](../src/plastid/external/rscript.py),
[script](../src/plastid/data/r_scripts/R_F9_phylogeny_tree.R).

## Other current limitations

The RSCU and amino-acid R figures accept only the packaged 34-species list. RSCU, amino-acid and
SSR modules each take one genome per run and do not yet produce a combined species panel. The
reference preparation and downstream analysis steps are not a single integrated workflow.

Unequal IR copies are legal and are measured end to end, but **mapping a gene between them** is
not implemented (M11.11): `plastid_annotate` refuses such a record by name (exit 2), and
`normalize` records `mirror_constant=0` because no equal-copy constant exists. See
[current scope and limits](things_we_do_not_do.md) for details.

## Pipeline coverage observations — code inspection, not new runtime reproductions

The [canonical pipeline](pipeline.md) maps every stage to existing components and missing work.
The following require targeted regression/acceptance checks in the next coding sessions:

- `gene_normalize` verify mode returns findings without writing a replacement GenBank. A workflow
  must explicitly decide how differences block or qualify downstream work. Normalized DB artifacts
  also need derived hashes and a link to the raw pin and chosen naming reference.
- `stats._sequence_of` joins forward slices without consulting strand, unlike `coding_sequence`.
  Strand-sensitive GC3/ENC and other CDS-derived results need a controlled comparison. Alignment
  coordinates and normalized-reference coordinates in diversity/barcodes also need rotation tests.
- `run --figures`, `package` and `qc` currently call the matplotlib renderer, while `figures`
  also has an R pass. Rendering an R figure separately does not establish that the package/QC
  includes it or its provenance.
- Submission summary fields use zero for validator counts when no validator result is present.
  Acceptance must distinguish “not run” from “ran and found zero errors”; creating a `.sqn` alone
  is not validation.
- `environment.yml` still mixes runtime tools with survey/deferred tools. Installed packages
  do not establish an assembly, RNA-editing or other unregistered module exists.

Sources: [naming](../src/plastid/modules/gene_normalize.py),
[group measurements](../src/plastid/modules/stats.py),
[diversity](../src/plastid/modules/diversity.py), [barcodes](../src/plastid/modules/barcodes.py),
[CLI output commands](../src/plastid/cli.py),
[submission](../src/plastid/modules/submission.py), [environment](../environment.yml).
