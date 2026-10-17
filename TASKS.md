# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

## The figures PTS does not have yet (M10.3 + M10.4)

**Goal.** The R engine runs (M10.2); no figure is wired to it. This task connects the toolkit's own
tables to the scripts, for the three the paper has and PTS lacks, and the three whose tables need
reshaping first.

**What the bundle's map says (`feedback/…/03_tich_hop_PTS/BANG_MAP_PTS_den_hinh.md`).**

| Paper figure | PTS table today | What it needs |
|---|---|---|
| 4 amino-acid composition | `amino_acid/amino_acid_composition.csv` | a new figure spec — the table already has the columns |
| 7 substitution types + Ts/Tv | `substitutions/substitutions*.csv` | a new figure spec, and a pivot: the paper's figure is *type × species*, the table is long |
| 8 Ka/Ks heatmap | `kaks/kaks_genes.csv` | a new figure spec, and the same long→wide pivot |
| 3 RSCU · 5 SSR · 6 π | `codon_usage/rscu.csv` · `ssr/*` · `diversity/*` | PTS draws these already (matplotlib); the R versions are a **second source for the same table**, and the two must agree on every number |

**Steps**
1. Decide, per figure, whether it is a `FigureSpec` whose `build` writes the table out and calls R,
   or a step beside the run. The figure contract is "reads cells, never recomputes" — an R-backed
   spec can honour that, but it must not become the place a number is computed.
2. The pivots go in the spec, from table cells, not in the R script.
3. Where PTS already draws the figure, assert the two agree: a disagreement is a defect in one of
   them, and the test says which table each read.

**Done when:** each new figure renders from a run directory with `plastid figures`, its caption
carries the scope label of the table it read, a test asserts its numbers equal that table's cells,
and `make check` is green.

**Traps**
- A figure reads cells and never recomputes — the §8.25 failure in a new costume.
- The R scripts write PNG + PDF into `<outdir>/Figures`; the run's own figures directory is the
  place for them, never the input's.
- The bundle's sample tables are the paper's: label anything drawn from them.

## After it

**M10.5** — the tree (paper figure 9 + S1–S3) stays the user's own script over the user's own tree;
PTS builds no tree. Then **H-25** as its own task with its own regression test: the IR gate refuses a
record whose two annotated IR copies differ in length, although the repeat itself is exact.
