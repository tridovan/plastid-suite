# Packaged data

Resource files shipped inside the package (not a Python package — no `__init__.py`).
Every path is resolved through `plastid.core.resources.data_file`, so no loader carries its
own "one directory up" guess.

- `defaults.yaml` — central default configuration (arrives in P1-05)
- `citations.yaml` — citation registry: tool → reference (arrives in P1-07)
- `gene_groups.yaml` — gene name → functional group, for the publication gene table (P2-09)
- `misa.ini` — the MISA 2.1 settings the `ssr` module runs with: **ours**, writing down MISA's
  default definitions. The tool's own file is not redistributed (D-q — its zip states no licence,
  `docs/tools/misa.md`), and its script is obtained separately.
- `r_scripts/` — the seven R scripts the manuscript's figures 3–9 were drawn with, plus the helper
  they share. **The team's own code**, carried so a figure can be reproduced, and kept in the
  language its authors wrote it in. R itself is not bundled: it is declared in `docs/tools/r.md`.
  The one change made to them is in `00_common.R` — input paths are made absolute before
  `setwd()`, without which a relative `--inputdir` finds nothing (measured 17/10/2026).

`gene_groups.yaml` is data, not code, on purpose: it holds external knowledge (the grouping
plastome papers use) and the real Oldenlandia/Coleus records will bring gene names it has not
seen yet. A name no group claims is reported as unclassified in the output table — never
dropped, never folded into a neighbour.
