# CGAS 1.0.1 — the second cross-check tool (P8-03, M8.1)

**Why this file exists.** CLAUDE.md's invariant: a tool that is used is declared in the repo. CGAS
is the **other side of M8.1's cross-check** — no module of this toolkit calls it, the minimum
install stays pip-only, and the release notes have to name it. The shape is Vmatch's (D-f) and
cpstools' (D-j).

**What it is.** Chloroplast Genome Analysis Suite (Abdullah, MIT) — 14 modules: assembly ·
annotation · gene comparison · GenBank conversion · comparative analysis · codon usage · amino
acids · SNP · introns · SSR · diversity · phylogeny. Source: `.data/cgas-review/`
(`github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS`, 19 `.py` files, ~21,400 lines —
read and reviewed 02/10/2026, PROJECT.md §7.7 · §7.9).

**Install line** (its own environment; nothing here may move the pinned builds in the shared
toolbox env):

```sh
uv venv /root/toolbox/env/cgas --python 3.12
uv pip install --python /root/toolbox/env/cgas/bin/python -e .data/cgas-review
```

**Measured 10/10/2026 on the dev WSL machine:**

| | |
|---|---|
| Version | `cgas --version` → **CGAS v1.0.1** |
| Python-only modules | 3, 4, 5, 6, 7, 10, 11 — these ran |
| R-dependent modules | 8, 9, 12, 13 **compute their numbers in Python and skip only the plots** when R is absent (`⚠ R is not installed or not in PATH → Skipping visualization`) — measured, not assumed |
| R on this machine | **absent** — so no CGAS figure was produced |
| Not run here | 1 (needs FASTQ reads), 2 (annotation — PGA, retired from this project, §8.56), 14 (needs IQ-TREE and MACSE) |
| License | MIT |

**The argv the cross-check used** (§8.62 is the results):

```sh
cgas -m 7  -i <dir>                              # genome structure: LSC/SSC/IR + GC
cgas -m 8  -i <dir>                              # codon usage / RSCU
cgas -m 11 -i <dir>                              # gene and tRNA introns
cgas -m 12 -i <dir>                              # SSR, its own defaults …
cgas -m 12 -i <dir> -t "10,6,5,5,5,5"            # … and again on MISA's thresholds
```

**What it measured about itself, worth keeping:**

* **It agrees exactly where it detects.** All four golden regions, their GC values, all 61
  sense-codon counts and the 37 SSRs (at equal thresholds) match ours digit for digit.
* **It fails silently where it cannot detect.** On `NC_024936.1` module 7 prints the species name,
  finds no region, prints **no warning**, reports `LSC 0 · SSC 0 · IR 0` in the table and still
  counts the record under *"Successfully analyzed: 2/2 genomes"*. Our own paths find
  26,184 bp @ 100.00 % on that record. This is the §7.9 failure family ("no IR ⇒ silent") seen from
  the other side: the record is not wrong, the zero is.
* **Module 12's defaults are not MISA's.** `{1:10, 2:5, 3:4, 4:3, 5:3, 6:3}` against MISA's
  `{1:10, 2:6, 3:5, 4:5, 5:5, 6:5}` — 54 SSRs versus 37 on the golden, one parameter apart.

**Not vendored, deliberately:** a Python package with its own environment — the same reasoning as
MUMmer, BLAST+ and ARAGORN.
