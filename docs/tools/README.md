# Tools — the inventory, and what each one is for (M11.19)

**Why this file exists.** An installed package reads like a capability. It is not: this project
spent a phase installing assemblers it never integrated, and a dependency list is not a feature
list (G-04). Every tool the toolkit touches therefore has **one role**, written here, and the role
decides which environment file carries it:

| Role | Means | Where it is installed |
|---|---|---|
| `active-runtime` | called by a module, a figure or a CLI step during a normal run | `environment.yml` (or obtained separately, with the line in its own page) |
| `cross-check` | **no module calls it**; it exists so a published number can be compared with a second program using the same definition | `environment-check.yml`, or its own page |
| `future` | belongs to a deferred branch (assembly, dating) with no module yet | `environment-future.yml` |
| `survey-only` | tried while choosing an engine and not chosen; kept named so a later scope decision is not mistaken for a discovery | `environment-future.yml` |

`tests/unit/test_tool_inventory.py` reads this table and checks it against the code and the three
environment files, in both directions: an adapter with no row fails; an `active-runtime` adapter no
product file imports fails; a `cross-check` adapter a product file *does* import fails; and a conda
package in any environment file that no row claims — or that sits in a file for a different role —
fails. The table is the claim; the test is what stops it decaying.

## The inventory

Usage is written as `adapter ← caller`. A dash under *Call site* with an `active-runtime` role is a
failure, not a placeholder.

| Tool | Role | Package | Version | Call site | Page |
|---|---|---|---|---|---|
| `perl` | active-runtime | `perl` | 5.32.1 | interpreter for `misa.pl` and `CPJSdraw.pl` — `run/tools.py` (`tool_version`), `external/misa.py`, `external/cpjsdraw.py` | [misa](misa.md), [cpjsdraw](cpjsdraw.md) |
| `mafft` | active-runtime | `mafft` | 7.526 | `external/mafft.py` ← `modules/substitutions.py`, `modules/diversity.py`, `modules/barcodes.py`, `modules/popgen.py`, `modules/kaks.py` | — |
| `aragorn` | active-runtime | `aragorn` | 1.2.41 | `external/aragorn.py` ← `modules/plastid_annotate.py`, `modules/mito_annotate.py` | [aragorn](aragorn.md) |
| `blast` | active-runtime | `blast` | 2.17.0 | `external/blast.py` ← `modules/plastid_annotate.py`, `modules/mito_annotate.py`, `modules/gene_normalize.py` | [blast](blast.md) |
| `mummer4` | active-runtime | `mummer4` | 4.0.1 | `external/mummer.py` ← `modules/synteny.py` | [mummer4](mummer4.md) |
| `repeat-match` | active-runtime | `mummer4` | 4.0.1 | `external/repeatmatch.py` ← `modules/repeats.py` (same package as `mummer4`) | [mummer4](mummer4.md) |
| `paml` | active-runtime | `paml` | 4.10.10 | `external/paml.py` (`yn00`) ← `modules/kaks.py` | [paml](paml.md) |
| `iqtree` | active-runtime | `iqtree` | 2.3.6 | `external/iqtree.py` ← `modules/phylogeny.py` | [iqtree](iqtree.md) |
| `macse` | active-runtime | `macse` | 2.07 | `external/macse.py` ← `modules/phylogeny.py` (`aligner=macse`, opt-in) | [macse](macse.md) |
| `R` | active-runtime | `r-base`, `r-ggplot2`, `r-dplyr`, `r-tidyr`, `r-reshape2`, `r-rcolorbrewer`, `r-zoo`, `r-scales`, `r-pheatmap`, `r-ape` | R 4.3; the nine packages deliberately unpinned (they draw, they do not compute a number) | `external/rscript.py` ← `cli.py` (`figures`, `package`) | [r](r.md) |
| `misa` | active-runtime | — (obtained separately) | 2.1 | `external/misa.py` ← `modules/ssr.py` | [misa](misa.md) |
| `table2asn` | active-runtime | — (NCBI binary) | 1.29.324 | `modules/submission.py` (`_build_sqn`, optional `.sqn` step) | [table2asn](table2asn.md) |
| `cpjsdraw` | active-runtime | — (obtained separately) | 0.0.1 | `external/cpjsdraw.py` ← `cli.py` (`ir-junction`) | [cpjsdraw](cpjsdraw.md) |
| `vmatch` | cross-check | `vmatch` | 2.3.1 | **none** — `external/vmatch.py` is called only by `tests/golden/test_tier_b.py` | [vmatch](vmatch.md) |
| `cpstools` | cross-check | `cpstools` | 3.0 | **none** — compared by hand, recorded in its page | [cpstools](cpstools.md) |
| `kakscalculator` | cross-check | `kakscalculator2` | 2.0.1 | **none** — reached only through cpstools' KaKs sub-command | [cpstools](cpstools.md) |
| `cgas` | cross-check | — (own venv) | 1.0.1 | **none** — compared by hand, a comparison rather than a standard | [cgas](cgas.md) |
| `parasail` | cross-check | — (validation only) | not reported | **none** — `tools/validate_ir_panel.py` compares the IR panel; never a runtime dependency | — |
| `pga` | survey-only | — (retired 10/10/2026) | 2.0 | **none** — its *output* survives as a counter-example fixture in `tests/data/corpus/` | — |
| `pmga` | survey-only | — (never ran here) | — | **none** — the mito annotate step is this project's own caller (P3A) | [pmga](pmga.md) |
| `sibelia` | survey-only | `sibelia` | 3.0.7 | **none** — ran in the 26/9 survey; the synteny route is `nucmer` (D4-02) | — |
| `beast2` | survey-only | `beast2` | 2.7.7 | **none** — divergence time is the PI's analysis, removed from the plan 10/10/2026 | — |
| `mcmctree` | survey-only | — (ships inside the runtime `paml` package) | 4.10.10 | **none** — the program is not a module; the package is runtime because of `yn00` | [paml](paml.md) |
| `mreps` | survey-only | `mreps` | not recorded | **none** — a tandem-repeat finder, a class the toolkit does not publish | — |
| `trf` | survey-only | `trf` | not recorded | **none** — same class as `mreps`; absent from this machine's PATH | — |
| `getorganelle` | future | `getorganelle` | 1.7.7.1 | deferred — Phase 7 assembly, `docs/pipeline.md` box [1] | — |
| `spades` | future | `spades` | 4.3.0 | deferred — Phase 7 assembly | — |
| `bandage` | future | `bandage` | 0.9.0 | deferred — M7.2 graph check; drawing the graph is not checking the assembly | — |
| `oatk` | future | `oatk` | 1.0 | deferred — M7.6 long reads | — |
| `raven` | future | `raven-assembler` | 1.8.3 | deferred — M7.6 long reads | — |
| `novoplasty` | future | `novoplasty` | not recorded | deferred — Phase 7 assembly | — |
| `minimap2` | future | `minimap2` | not recorded | deferred — M7.1 map-back | — |
| `samtools` | future | `samtools` | not recorded | deferred — M7.1 map-back | — |
| `wgsim` | future | `wgsim` | not recorded | deferred — survey of read simulation; no milestone owns it | — |

**`future` is the honest word for assembly, not a promise.** Phase 7 opens on M11.4/M11.8
acceptance; until then these tools are listed so nobody mistakes their presence for a pipeline box
that works.

## The three environments, and what each is for

| File | Contents | Who installs it |
|---|---|---|
| `environment.yml` | the runtime scope: every `active-runtime` tool that conda can supply | anyone running the analyses that need a tool |
| `environment-check.yml` | `cross-check` engines that conda can supply | someone verifying a published number |
| `environment-future.yml` | `future` and `survey-only` tools | nobody, for this release — it is the written state of a deferred branch |

Tools that cannot come from conda — `misa.pl`, `table2asn`, CPJSdraw, CGAS, parasail — have their
exact acquisition line on their own page and are **never** in an environment file; a conda line
that does not work is worse than no line (M11.19's rule for PREP-Cp applies to all of them).

**Minimum install.** `pip install 'plastid-suite'` is biopython + PyYAML. The optional extras
(`docx`, `figures`, `stats`, `excel`) are imported inside the functions that need them, so a
minimum install runs every pure-Python analysis and says which extra to add when one is missing.
Measured with a fresh venv in `docs/platforms.md`.

## Page map

Each page states why the tool is declared, how it is obtained, the version measured, and the
defect a cross-check found. A page is not a promise that the tool is integrated — that is the
role column's job.
