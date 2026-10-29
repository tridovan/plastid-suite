# Development handoff

Cleared at the start of the 28/10/2026 (Đợt 5) session, with a follow-up on 29/10/2026 (M11.21
`matrix=all`). Earlier handoffs are in git history (`git log -p -- HANDOFF.md`); the open work is
`.dev/TASKS.md`, the milestones `.dev/ROADMAP.md`.

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
  Measured on the 4-report panel: intron matrix + tree in **12.4 s**. **`matrix=all` (29/10/2026)** —
  the complete matrix CGAS module 14 recommends: CDS + intron + intergenic-spacer blocks in **one**
  alignment and **one** tree, block ids prefixed `cds:`/`intron:`/`igs:`, IR gate once per genome;
  `aligner=macse` and `partition=codon` are refused for it. Real run on the 4-record panel: **197
  blocks / 134,962 sites / 60.9 s**. Its internal split differs from the CDS matrix's (different
  data, so a different tree — not recorded as an improvement); no PTS number moved.
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
- **IR gene-mirror rule (29/10/2026)** — the maintainer's criterion: a gene inside one IR copy
  must have its copy inside the other, on the **opposite strand**. `.dev/scripts/ir_gene_mirror.py`
  checks it (read-only), with the golden record as a positive control (**18 loci vs 18, no
  violation**). It confirms the K-07 readings: in `NC_068894.1` `rpl2` (1,491 bp) and `rpl23`
  (282 bp) are annotated twice with equal lengths and opposite strands, both outside the annotated
  repeat and inside the de-novo pair; `NC_053322.1` shows the same for `rrn23`, `rrn4.5`, `rrn5`,
  `trnN-GUU`, `trnR-ACG`; `NC_045310.1` is perfectly mirrored on both pairs, so its +90 bp stays
  `unexplained`. `rps19` is a boundary gene (59–60 bp against 279 bp) and is excluded. A panel
  sweep (77 records, 3 de-novo refusals) found **18 of the 47 length-agreeing records broken in
  their own annotation**. Adjudicating against each record's own sequence (`ir_mirror_adjudicate.py`)
  gives **46 violations in 26 records: 41 are the record's** (same-strand copies whose sequences
  are reverse complements, or a gene annotated once while the copied sequence exists at 100 %),
  **3 more are the record's** once the IR alignment mapping replaces the geometric mirror, and
  **2 stay unexplained** (`NC_057683.1` `rpl2`/`ycf2`, a record whose copies are homologous along
  a shifted alignment). Two mistakes in my own measurement were found and corrected first — the
  annotation pair was compared without `IRResult.frame_offset`, and a geometric mirror was used
  where the alignment mapping is required; they had inflated the count to 100 violations. Third-
  party annotation inconsistencies are not reported by PTS today; adding a diagnostic is a
  maintainer decision, not done. The sweep also found three name-normalization defects, fixed as
  **K-13**. **The rule is now in the reference QA (29/10/2026):** `core/reference_audit.py` reports
  a `mirror` finding per record (REPORT, pinnable by sha256, not BLOCK); the 34-record reference
  set passes it 34/34, the validation panel reports 23/80 with no change to its block counts.
- **Test cost** — the six real phylogeny tests and the real MACSE test moved from `tests/unit/`
  to `tests/golden/` with `tier_b`, so `make test` no longer pays ~5.5 minutes of IQ-TREE on every
  edit; `make golden` runs them.
- **M11.17 rotation, the last six modules (29/10/2026)** —
  `tests/golden/test_m11_17_rotated_remaining_cases.py` **6 passed / 2:31** with the real tools:
  `gene_table` and `orf` are exactly invariant (only `orf`'s three coordinate cells excluded),
  `kaks`/`popgen` tables are identical with each genome rotated by its own offset, `synteny`'s
  reference rotation is exactly invariant (normalisation is origin-independent) while a rotated
  **query** is only claimed for the reference-side measurements (MUMmer splits a block at the
  query's origin), and `phylogeny`'s supermatrix and tree are **byte-identical** on four panel
  genomes. Writing it found **K-12**: `phylogeny` accepted three genomes and then died inside
  IQ-TREE, whose ultrafast bootstrap needs four — the guard is now `< 4`, with the tool's own
  reason, before any alignment; `docs/platforms.md` and `known_issues.md` updated.
- **Roadmap audit (29/10/2026)** — five `[~]` markers were stale and are now closed with evidence:
  **M11.7** (multi-sample figures come from `plastid figures <batch-dir>`; junction work done),
  **M11.18** (full golden 151 passed, nine-family batch check 0 mismatches, CGAS definition
  comparison measured), **M11.19** (Docker build exit 0 and fresh `conda env create` exit 0 in the
  M11.8 acceptance; macOS stays unverified as its DoD allows), **M11.20** (the last hard-coded
  parameters reach argv/ctl), **M11.21** (`matrix=all` done). **M11.8** stays `[~]` with its old
  reasons replaced: its remaining blockers are the open **M11.17** PGA cases and the unexplained
  K-07 residual boundary records, not the full golden or the nine families (both done). **M11.1**
  stays `[~]` — no record shows its cleanup review was completed.

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

**In flight at handoff:** nothing. The MACSE 34-genome scale run (started 28/10 21:57,
`/tmp/macse34`, log `/tmp/macse34.log`) finished at 29/10 01:22 — see the M11.21 bullet above.
Its directory was under `/tmp` and has since been cleared, so only the recorded numbers survive;
a re-run of the artifacts costs ~2 h.

**Not done:** the remaining M11.17 cases (the 80-family/30-PGA pair for `kaks`, `popgen`,
`synteny`, `phylogeny`, `orf`, `gene_table`, and the PGA pair for the gene-based modules),
intron/IGS/all matrices at 34 genomes (measurement only, M11.21), M4.5 RNA editing. M11.16 is
still `[~]` (structured IR validation fields, the 30-PGA artifact cases, the local-tool contract);
M11.8/M11.9 wait on M11.17.
