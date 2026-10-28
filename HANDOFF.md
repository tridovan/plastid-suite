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
  **323/120 = 2.692 vs 2.6917**, every type equal; per-gene π identical at 6 dp. Non-comparable
  cells named (CGAS's 1-dp plotted π, merged 218-region list). Module 14 not yet run.
- **M11.17** — first panel cases: `source/` vs `annotation-variants/` (same sequence) must give
  the same `genome_stats`, `ssr` and `repeats`; 4 passed. The remaining modules and the 30 rotated
  FASTA cases are open.
- **Test cost** — the six real phylogeny tests and the real MACSE test moved from `tests/unit/`
  to `tests/golden/` with `tier_b`, so `make test` no longer pays ~5.5 minutes of IQ-TREE on every
  edit; `make golden` runs them.

**Full gate on the final tree, 28/10/2026 (every pinned tool, `PTS_REQUIRE_TOOLS=1`):** ruff clean;
fast suite **1170 passed, 1 skipped, 78 deselected** (**8:58**; the skip is the intentional
no-table2asn machine test); full golden **138 passed, 0 skipped** (**25:32**). Log:
`/tmp/pts_full_check.log`. The six real phylogeny tests and the real MACSE test now live in
`tests/golden/` as `tier_b`, which cut the fast suite from 12:58 to 8:58 and moved their cost to
the golden (17:07 → 25:32).

**Environment corrections.** `PTS_IQTREE_DIR` is `/root/toolbox/env/iqtree/bin`, not
`/root/toolbox/env/pts/bin`. MACSE is not on the vendor URL any more; the pinned source is
bioconda `macse=2.07` (the jar is at `/root/toolbox/env/macse/share/macse-2.07-0/`).

**Not done:** the rest of M11.17 (per-module cases on the 80/10-family panel and the 30 rotated
FASTA), CGAS module 14 with MACSE/IQ-TREE, MACSE at 34-genome scale, M4.5 RNA editing.
