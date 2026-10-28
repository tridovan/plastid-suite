# Development handoff

Cleared at the start of the 28/10/2026 (Đợt 5) session. Earlier handoffs are in git history
(`git log -p -- HANDOFF.md`); the open work is `.dev/TASKS.md`, the milestones `.dev/ROADMAP.md`.

## This session

Worked `.dev/TASKS.md` in the order agreed with the maintainer (batch the code tasks, then one
gate), then the follow-up items. Each has its evidence record under `.dev/records/`:

- **K-07** — the de-novo IR chain allows the overlap an indel necessarily creates; 132/132 indel
  cases. The 30 source records that disagree with their annotation were re-measured after
  correcting for `IRResult.frame_offset` and adjudicated: 50/80 agree, 20 are boundary-only and
  unexplained, 4 de-novo pairs fail the gate, 2 annotations are short (`NC_068894.1`,
  `NC_053322.1` — the latter overturns the earlier reading), `NC_045310.1` unexplained, and 3 are
  refusals. `docs/known_issues.md` K-07 is updated; user review of the two data-side conclusions
  is pending.
- **M11.20 (SSR)** — `--set unit_repeats` writes MISA's definition line; malformed values exit 2;
  the summary reports the thresholds read back from the ini. Default unchanged (37 = 31·4·2).
- **M11.21** — `partition=none|gene|codon` (`-p`, NEXUS artifact), measured 65 s / 111 s on the
  4-record panel; and **MACSE 2.07** as opt-in `aligner=macse` (bioconda, CeCILL 2.1). Real run:
  78 genes, 70,110 sites, 234 codon partitions, **22:55**, same topology, 100/100 support.
- **M11.7** — `genome_stats.genes` carries a `segments` cell, so an origin-spanning feature is
  drawn as its parts in `ir_junction`/`circular_map`; panel window is 3 kb single-copy + the
  whole IR, caption corrected, figure width 7.4–10 in.
- **CGAS comparison** — modules 10 and 13 measured on the same alignment/scope: Ts/Tv
  **323/120 = 2.692 vs 2.6917**, every type equal; per-gene π identical at 6 dp. **Module 14** was
  then run with MACSE + IQ-TREE on the same 4 records: same topology, **70,101 sites / 78 gene
  partitions** against PTS's 70,110 / 234 codon partitions, support **99.2/97** against **100/100**.
  The 9-site and support differences are unexplained and no PTS number moved; non-comparable cells
  are named in the record.
- **M11.20 (closed)** — the last hard-coded tool parameters are doors: `plastid_annotate.max_hsps`,
  `gene_normalize.evalue/max_hsps`, `kaks.yn00_icode/yn00_weighting` (the ctl is written from the
  values, not a constant). Out-of-range values are refused before the tool; the values used are
  named in the scope or the notes. No default changed.
- **M11.21 (closed for the product)** — `phylogeny --set matrix=cds|introns|igs` builds the CDS
  matrix (default, unchanged) or a matrix from the intron / intergenic-spacer loci
  `core.spans` derives, keyed by canonical gene name; the IR gate runs per genome because a span
  boundary is a coordinate, and `aligner=macse` is refused for a matrix with no reading frame.
  Measured on the 4-report panel: intron matrix + tree in **12.4 s**. CGAS can still concatenate
  the three matrices in one run (`matrix=all` remains open).
- **34-genome MACSE (M11.21, measurement only)** — the first run **failed**: 62 of 78 genes
  finished at ~63 s each (65 min), then gene 63 exceeded the adapter's fixed 600 s and the run
  stopped with `TimeoutExpired`. That is a product limit, not bad data, so `phylogeny` now takes
  `macse_timeout` (seconds per gene, default 600, validated before the jar is resolved; refused
  with `aligner=mafft`). The second run (`macse_timeout=3600`) **completed**: 34 genomes,
  **78 CDS blocks, 70,131 sites**; MACSE 2 h 06, IQ-TREE 4 min 52 s, total **2 h 11**; tree with
  mostly 100/100 support. MACSE at 34 taxa costs ~5.5× the four-genome run for the same number of
  sites, so the cost is taxa, not sites.
- **M11.17** — completed for the requested modules. `tests/golden/test_m11_17_rotated_cases.py`
  runs `structure` on all 30 unannotated rotated FASTA (identical region table across five
  rotations, exact tiling, and de-novo IR/origin reproducing `EXPECTED.tsv` within 2 bp).
  `tests/golden/test_m11_17_rotated_annotation_cases.py` rotates the annotated panel records by the
  same offsets and shows `codon_usage`, `amino_acid`, `stats` (and `barcodes`, `diversity`,
  `substitutions` with the reference rotated) invariant. With the earlier `source/` vs
  `annotation-variants/` cases, that is 11 tests over 9 modules.
  **The panel found two product defects**, both fixed with regressions (see `known_issues.md`):
  **K-10** — `diversity`/`substitutions` aligned the reference as deposited and named its regions
  in the normalized frame (same two genomes, rotated reference: 230,663 vs 154,827 columns, 386 vs
  1,094 substitutions, π `None` on most gene rows); **K-11** — five of six `mafft_strategy` names
  were not options of MAFFT 7.526 (`--fftns`, `--linsi`, …) and the unit test asserted their shape
  instead of running them. A third fix came with K-10: `spans.introns_of` read the two pieces of an
  exon cut by the record's origin as a 153 kb intron.
- **Test cost** — the six real phylogeny tests and the real MACSE test moved from `tests/unit/`
  to `tests/golden/` with `tier_b`, so `make test` no longer pays ~5.5 minutes of IQ-TREE on every
  edit; `make golden` runs them.

**Gate on the final tree, 28/10/2026 (every pinned tool, `PTS_REQUIRE_TOOLS=1`,
`PTS_MACSE_JAR` set):** ruff and `ruff format --check` clean; fast suite **1170 passed, 2 skipped,
91 deselected** (**10:06**; the two skips are the misa.pl and no-table2asn machine checks — the
first passes when `PTS_MISA_DIR` is exported, as in the golden below); **full golden 151 passed,
0 failed, 0 skipped** (**26:18**, log `/tmp/golden_m1117.log`). The 13 new tier-B tests are why
the deselected count rose from 78 to 91; the previous golden on this branch (before M11.17) was
138 passed / 25:32.

**Environment corrections.** `PTS_IQTREE_DIR` is `/root/toolbox/env/iqtree/bin`, not
`/root/toolbox/env/pts/bin`. MACSE is not on the vendor URL any more; the pinned source is
bioconda `macse=2.07` (the jar is at `/root/toolbox/env/macse/share/macse-2.07-0/`). The bioconda
`macse` **wrapper** needs the env's own `java`, which the extracted tarball on this machine lacks:
export `PTS_MACSE_JAR=<jar>` locally, and CI now installs `macse=2.07` + `openjdk=17` in the
tier-B env. The MACSE golden test resolves through `PTS_MACSE_JAR` or `PATH` and no longer probes a
developer's absolute path — the previous version raised `PermissionError` on CI's unreadable
`/root` and failed the build (run 40703129920).

**In flight at handoff:** the MACSE 34-genome scale run started 28/10 21:57
(`/tmp/macse34`, log `/tmp/macse34.log`; ~1.5–2 h expected) and has not finished.

**Not done:** the remaining M11.17 module cases (kaks, popgen, synteny, phylogeny, orf,
gene_table and the 80/10-family PGA pair for the gene-based modules), the MACSE 34-genome result,
separate intron/IGS matrices for `phylogeny` (M11.21), M4.5 RNA editing.
