# R — declared here, not bundled

**Why this file exists.** The toolkit never bundles a tool it can call. `Rscript` is the interpreter
the manuscript's seven figure scripts run under, and a missing interpreter must be **named**, never
guessed at.

**What it is.** R (the R Foundation), used for seven scripts that write PNG + PDF. Six consume
tables exported by the toolkit; the tree script instead reads a user-supplied tree and can alter
its root and included tips. Every figure needs traceable inputs and transformations. A table/figure
disagreement needs investigation; neither output is automatically an independent answer key.

**Versions measured:**

| | |
|---|---|
| on this machine, 17/10/2026 | `R version 4.3.3 (2024-02-29)`, `/usr/bin/Rscript` |
| by the scripts' author, 16/10/2026 | `R version 4.1.2 (2021-11-01)` |
| packages required | `ggplot2` `dplyr` `tidyr` `reshape2` `RColorBrewer` `zoo` `scales` `pheatmap` `ape` |

All seven scripts run and write all eight figures under **both** versions — measured, not assumed,
because the bundle was written on 4.1.2 and this toolkit's machine has 4.3.3.

**Packages deliberately NOT required** — `readr` (use `read.delim`), `cowplot` (the two-panel
layout is drawn with `grid`), `gridExtra`. The scripts in the bundle were rewritten so they are not
needed, and a script that needs one of these is a script that fails on the target machine. A missing
package from the nine above is reported **by name** before any figure is attempted.

**Install line** (verified 17/10/2026, Debian/Ubuntu):

```sh
sudo apt-get install -y r-base-core \
  r-cran-ggplot2 r-cran-dplyr r-cran-tidyr r-cran-reshape2 r-cran-rcolorbrewer \
  r-cran-zoo r-cran-scales r-cran-pheatmap r-cran-ape
```

Conda route through `environment.yml`: R's series is pinned; the nine package versions are not:

```yaml
- r-base=4.3
- r-ggplot2
- r-dplyr
- r-tidyr
- r-reshape2
- r-rcolorbrewer
- r-zoo
- r-scales
- r-pheatmap
- r-ape
```

**How it is found.** `Rscript` on `PATH`, or passed explicitly: `--set rscript=/usr/bin/Rscript`.
Nothing else is probed: the interpreter's version is read once, and the packages are checked once.

**The scripts are package data, not this package's code.** They live in `src/plastid/data/r_scripts/`
and are read through `plastid.core.resources.data_file`, the way `data/misa.ini` is. They are the
manuscript's own scripts, kept in the language their authors wrote them in — Vietnamese comments and
all — because they are carried for reproduction, not adopted as this toolkit's source.

**One change was made to them.** `khoi_tao()` ends with `setwd(OUT_DIR)`, so a *relative*
`--inputdir` or `--input` points into the output directory instead of at the caller's data and every
script fails to find its own input. Measured 17/10/2026: `--inputdir data_mau` fails, an absolute
path works. `00_ham_chung.R` now makes those paths absolute before the working directory moves. No
drawing changed.

**How a run wires its tables to these scripts.** `plastid figures <run_dir>` renders them, beside
the toolkit's own figures — `--set rscript=/path/to/Rscript` when R is not on `PATH`. Each R figure
is a pivot in `src/plastid/result/figures/r_*.py` that reads table **cells** and writes the file the
script expects; the script then runs through `external/rscript.py`, so its version and argv land in
the run's provenance. Four things are worth knowing:

* the pivot file the script read is kept in `<run>/figures/<id>__input/`, and the images beside it
  as `<id>__<name>.png` and `.pdf`;
* a **caption** is written as `<run>/figures/<id>.caption.txt`, carrying the scope label of every
  table the figure read — the scripts draw no source note of their own, and an exported number must
  carry its scope;
* a figure is **skipped by name**, never drawn empty: when a source table is absent from the run,
  when R or one of the nine packages is missing, when the run's species is not one the script orders
  by (F3 and F4 intersect their own hard-coded list and would otherwise draw an empty panel), or
  when the run's scope produced no Ts/Tv table (F7);
* the three arity-1 modules — `codon_usage`, `amino_acid`, `ssr` — take **one genome per run**, so
  F3, F4 and F5 carry one species column. `kaks` and `diversity` take the whole panel, so F6 and F8
  cover every genome in the run;
* **F9, the tree, reads no table at all.** This toolkit builds no tree, so the figure takes the
  user's own file — `--set tree=<file.treefile>` — hands it to the script as its input, and never
  copies it into the run. Its caption carries that file's path and sha256 where the others carry a
  scope label. `R_F9` defaults drop two *Asparagus* tips and root on two Orchidaceae, falling back
  to a midpoint root when the outgroup is absent. The CLI currently does not forward overrides
  for these choices; `--set root=none`, `outgroup=...` or `drop-tips=...` does not control them.

**Controlling the tree directly.** Until CLI option forwarding is repaired (K-06 in
[known issues](../known_issues.md)), a source checkout can call the packaged script explicitly:

```sh
Rscript src/plastid/data/r_scripts/R_F9_phylogeny_tree.R \
  --input /absolute/path/to/my.treefile --outdir /absolute/path/to/tree-output \
  --drop-tips "" --root none --og-nhan ""
```

This keeps the existing root and tips. To choose another root, supply the script's explicit
`--root`/`--outgroup` arguments. A direct call is outside the toolkit run recorder; retain the
input checksum, command, R version and output yourself. This example was checked against the
script's argument handling; it was not executed during the documentation-only session.

**Three limits of the scripts themselves**, all measured 17/10/2026 and each turned into a **named
skip** or a pinned test rather than left to a cryptic R error:

* **`R_F8_KaKs_heatmap.R` cannot draw a one-genome matrix.** Its own `mat[keep_genes, keep_spp]`
  drops a single species column to a vector, so `rownames(mat)` is NULL and it halts with
  `replacement has 0 rows`. Isolated on a hand-built matrix: two columns draw, one does not. A run
  comparing the reference against a single genome is skipped with that reason.
* **`R_F3` and `R_F4` draw only the species on their own hard-coded 34-name list** — anything else
  is dropped by their `intersect(SP_ORDER, colnames(data))` and the panel comes out empty, without
  a word. A run of any other species is skipped by name. The list is parsed from the packaged
  script so it cannot drift unnoticed.
* **`R_F6_nucleotide_diversity_pi.R` halts on a repeated `Region` label** — its
  `factor(Region, levels = Region)` refuses a duplicate and R stops with `factor level [6] is
  duplicated`, drawing nothing. Any figure feeding it must therefore label rows by a **unique key**,
  which is why its pivot uses the tables' own `intron_id`/`igs_id` cells rather than a name built
  from parts: on the golden pair a gene-name label collides twice over, and neither collision is
  visible from the name (two IR copies of one intron, and two introns of one gene).

**Current script status.** The absolute-path fix is in `00_ham_chung.R`; the scripts otherwise
retain manuscript-specific behavior. General species selection and configurable tree behavior
are planned work, not a reason to preserve these defaults as permanent toolkit requirements.

**Cite it.** R Core Team. *R: A Language and Environment for Statistical Computing.* R Foundation
for Statistical Computing, Vienna, Austria. https://www.R-project.org/
