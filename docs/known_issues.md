# Known correctness and integration issues

These issues were reproduced during the 26 October 2026 review. They are **open**: documenting a
defect does not fix it. The code version is `1.0.0`; that version label does not establish readiness
for scientific use. Seventy-four related tests passed during the review while the probes below
still exposed defects. A full tool-equipped acceptance run was not performed.

Coordinates below are 0-based half-open unless a GenBank position is explicitly stated. The pinned
fixture is [NC_050833.1.gb](../tests/data/NC_050833.1.gb); its provenance and checksum are in
[PROVENANCE.md](../tests/data/PROVENANCE.md) and [SHA256SUMS](../tests/data/SHA256SUMS).

## K-01 — IR extent and identity depend on record shape

The fixture has two 26,214 bp IR copies at `[83127,109341)` and `[127510,153724)`. Keeping the
underlying genome but changing its representation produces the following results:

| Input change | Observed | Required behavior for this controlled case |
|---|---|---|
| Remove repeat annotations; change one base at position 140617 in the second copy to a different A/C/G/T base | De-novo IR length 13,107 bp; identity 100%; gate passes | Recover the independently defined full extent and count the mismatch |
| Remove repeat annotations; rotate right by 5,000 bp | De-novo IR length 21,214 bp; identity 100%; gate passes | Changing the origin preserves IR extent |
| Delete position 140617; retain annotations, adjusting the second copy's end to 153723 | Identity about 64.68%; gate refuses | Align the copies with the one-base deletion represented as a gap |
| Rotate right by 1 bp and retain a joined repeat annotation across the origin | One annotated copy is read as 26,213 bp | Read the whole joined copy as 26,214 bp |

The detector can measure an exact-match fragment as if it were the full repeat. The annotation
reader also discards the grouping of compound locations and can silently switch to de-novo when
annotated lengths differ by more than five bases. A passing gate alone is therefore insufficient
evidence that IR boundaries are correct; a refusal alone does not establish that the input is wrong.

Sources: [GenBank reader](../src/plastid/core/genbank.py), [IR detection](../src/plastid/core/ir.py).
The proposed explicit source-selection option is not implemented; do not assume that
`--set ir_source=sequence` currently changes this behavior.

## K-02 — Cache omits annotation fields

`_input_hash` omits `transl_table`, `codon_start`, `exception` and `trans_splicing`. Two annotations
that differ only in those fields can share a cached result even when an analysis reads the field.

A minimal diagnostic, run from the checkout with its Python environment:

```python
from dataclasses import replace
from plastid.core.model import Gene, Genome, Interval
from plastid.core.orf import check_cds
from plastid.run.runner import _input_hash

first = Genome("probe", "GTGAAATAA", genes=(Gene("x", "CDS", (Interval(0, 9),), transl_table=1),))
second = replace(first, genes=(replace(first.genes[0], transl_table=11),))
print(check_cds(first, first.genes[0]).status)  # no_start
print(check_cds(second, second.genes[0]).status)  # alternative_start
print(_input_hash(first) == _input_hash(second))  # True: defective cache identity
```

The review also exercised `execute_step`: the second ORF run returned `cached` with the first
run's table-1 result. Acceptance requires checking exported results through the runner, preserving
valid cache hits, and covering all input fields that affect behavior. The diagnostic uses internal
helpers and is not a public API contract.

Source: [runner](../src/plastid/run/runner.py).

## K-03 — Junction reader cannot serve as an independent oracle

Starting with the same fixture, rotate left by 83,127 or 127,510 bp and move both annotations by
the same offset modulo the sequence length. One whole IR then starts at GenBank position 1.
`junction_regions` raises `JunctionError`, mistaking that whole copy for a wrapped fragment.

Rotate left by 88,127 or 96,234 bp instead and represent the wrapped copy as a `join`. Both joined
parts are at least 5 kb, so the reader counts three copies and raises again. All cases still have
the same two 26,214 bp IR copies, 83,127 bp LSC and 18,169 bp SSC.

Agreement on a small set of original records does not establish correctness under rotation.
Expected boundaries in regression tests must come from independently constructed locations, not
from this reader or the detector being repaired.

Source: [junction reader](../src/plastid/core/junctions.py).

## K-04 — Unknown bases are counted as matching evidence

Replace every base in both annotated IR copies of the fixture with `N`, preserving the sequence
length and locations. `ir_from_annotation` returns identity `1.0` and `gate_passed=True` over
26,214 positions. This is a reproduction of the annotation measurement, not a completed CLI run.

Acceptance must distinguish observed A/C/G/T matches from ambiguous positions, state the identity
denominator and require enough comparable sequence. Merely dropping unknown positions from the
denominator could also make a tiny known fragment look fully supported.

Source: [pair identity](../src/plastid/core/ir.py).

## K-05 — Rotation drops partial-coordinate flags

Set `partial_start=True` on a gene interval and call `rotate_genome(genome, 1)`. The returned
interval loses the flag. `_shift_interval` reconstructs intervals without carrying either partial
flag, including when a feature does not cross the origin.

Acceptance must preserve the original uncertain endpoints on either strand, keep artificial split
points exact, preserve transcript order and coding sequence, and check GenBank write/read behavior.

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
See [current scope and limits](things_we_do_not_do.md) for details.

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
