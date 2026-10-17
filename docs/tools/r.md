# R — declared here, not bundled

**Why this file exists.** The toolkit never bundles a tool it can call. `Rscript` is the interpreter
the manuscript's seven figure scripts run under, and a missing interpreter must be **named**, never
guessed at.

**What it is.** R (the R Foundation), used here as a **plotting engine and nothing else**: seven
scripts that read a table the toolkit exported and write PNG + PDF. None of them computes a number
the toolkit does not already hold in a table cell — if a script and its table disagree, the table is
right.

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

Conda route, if the toolkit is installed through `environment.yml` — the same nine, pinned there:

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

**Cite it.** R Core Team. *R: A Language and Environment for Statistical Computing.* R Foundation
for Statistical Computing, Vienna, Austria. https://www.R-project.org/
