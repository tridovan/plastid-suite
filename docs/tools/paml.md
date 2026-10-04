# PAML 4.10.10 (`yn00`) — declared, not vendored

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo — `environment.yml` (pinned), `vendor/` + `SOURCE.json` (small
scripts), or this file. PAML is a compiled C package that bioconda and conda-forge ship, so it is
pinned in `environment.yml` and described here. `yn00` is the estimator of Phase 4's Ka/Ks module
(P4-02, decision **D4-01**).

**What it is.** PAML (Yang 2007) is a phylogenetic analysis package; **yn00** is its pairwise
estimator of synonymous and nonsynonymous substitution rates by the Yang & Nielsen (2000)
method. It is the standard pairwise estimator and it was already pinned in `environment.yml` for
Phase 6's `mcmctree` — the same package, a different program in it. Rejected alternatives, with
the reason: a hand-rolled NG86 (biased, and a wrong implementation would be a silent wrong nobody
could check), and KaKs_Calculator2 (a second tool for a job that already has one).

**Versions and provenance, measured 2026-10-04 on the dev WSL machine:**

| | |
|---|---|
| Software | PAML **4.10.10** (`YN00 in paml version 4.10.10, 29 Jan 2026`), build `h87e0c26_2` |
| Package | bioconda / conda-forge |
| Path here | `/root/toolbox/env/pts/bin/yn00` |
| Measured on | the golden pair: **79 shared CDS, 78 scored**; pooled Ka 0.003586 · Ks 0.016550 · Ka/Ks 0.2167 |

**Install line:**

```sh
micromamba install -c conda-forge -c bioconda paml=4.10.10
```

**Not vendored, deliberately:** compiled binaries, resolved per platform by conda — the same
reasoning as BLAST+, MUMmer and ARAGORN.

**How it is invoked, and why.** `yn00` takes **no arguments**: it reads `yn00.ctl` from the
working directory and writes its answer to files named in that ctl. The adapter
(`external/paml.py`) therefore writes both the input and the ctl into a throwaway directory per
pair, and the manifest records the tool with its version; the ctl's settings are:

```
seqfile = input.dat
outfile = yn00.out
verbose = 0
icode = 0     # standard genetic code — NCBI table 11 differs from it only in initiation codons
weighting = 0 # equal weighting of pathways (the YN default)
```

**Three measured facts about the tool** (each of them silent when ignored):

* **the version is not a flag.** Run with no ctl it prints its banner on stdout and exits **1** —
  that is the probe, and a release other than the pin is refused by name.
* **it reads codons and never checks the frame.** A one-base gap in a synthetic pair made it
  report a spurious stop codon and exit 255 — luck, not a check. The frame check is ours
  (`core.align.codon_frame_breaks`), and it runs before yn00 sees anything.
* **stop codons are its poison** — the terminal one must be trimmed by the caller, an internal one
  stops the run. Both are reported, never repaired silently.

`dN` and `dS` are read from the `2YN.dN` / `2YN.dS` matrices and cross-checked against the
Yang & Nielsen row of the run's own output; the two disagreeing is not something this toolkit
averages. `nan` (printed as `-nan` when dS is zero) and `-0.0000` are read as "undefined" and
"zero" — never as a number.

**Licence.** PAML is distributed under a free licence for academic use (the package metadata
states GPL-3.0); flagged for the M9.3 review like every other tool.

**References.** Yang Z. PAML 4: phylogenetic analysis by maximum likelihood. *Molecular Biology
and Evolution* 24:1586–1591 (2007). doi:10.1093/molbev/msm088 · Yang Z, Nielsen R. Estimating
synonymous and nonsynonymous substitution rates under realistic evolutionary models. *Molecular
Biology and Evolution* 17:32–43 (2000). doi:10.1093/oxfordjournals.molbev.a026236
