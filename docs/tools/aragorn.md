# ARAGORN 1.2.41 — declared, not vendored

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo — `environment.yml` (pinned), `vendor/` + `SOURCE.json` (small
scripts), or this file. ARAGORN is a compiled C program (478,176 B, single binary) that bioconda
ships, so it is pinned in `environment.yml` and described here. It is the **tRNA caller** of Phase
3A (P3A-04): the mito annotator gets its CDS from BLAST+ and its tRNAs from this.

**What it is.** ARAGORN (Laslett & Canbäck 2004) detects tRNA and tmRNA genes in nucleotide
sequences by searching both strands for the canonical cloverleaf, scoring each candidate. It is
the tool the heavier annotators run inside — PMGA's own dependency list names it — which is the
reason it is the right choice rather than a hand-written structural caller: an in-house caller
tuned until it returned the expected count would be overfitting with a straight face
(D3A-01).

**Versions and provenance, measured 2026-10-03 on the dev WSL machine:**

| | |
|---|---|
| Software | ARAGORN **1.2.41** (`ARAGORN v1.2.41  Dean Laslett`), build `h7b50bb2_5` |
| Package | bioconda, 147 kB, licence GPLv3 |
| Path here | `/root/toolbox/env/pts/bin/aragorn` |
| Measured on | `OX459128.1`: **17 tRNA genes**, exactly the number the published record carries |

**Install line:**

```sh
micromamba install -c conda-forge -c bioconda aragorn=1.2.41
```

**Not vendored, deliberately:** a compiled binary per platform, and the upstream source
(`ansikte.se`) answers **HTTP 455** from this network — bioconda is the route, exactly as the
decision recorded. 147 kB would be cheap to vendor, but a vendored Linux binary is not what a
Windows/Docker user needs, and conda resolves the right build per platform.

**The invocation the adapter uses, and why:**

* `-w` — batch mode: one line per gene, `N  tRNA-<Aa>  c?[start,end]  <score>  (<anticodon>)`.
  The default output is a secondary-structure drawing, which is for a person to read.
* `-gc1` — the standard genetic code (NCBI `transl_table=1`), which is the code plant mitochondria
  use and the same table the caller stamps on every CDS it calls. Passed explicitly rather than
  left to the default so the report and the ORF checks cannot disagree.
* `-c` for a circular record / `-l` for a linear one — ARAGORN's default assumes circular and
  wraps its search around the ends. Measured on the answer-key record: its first tRNA
  (`tRNA-Phe`, `[258274,73]`) **crosses the origin**, and a wrapped interval is not a GenBank
  coordinate — the adapter splits it into the two parts a `join()` needs.
* no `-i`: introns are off by default, and the measured run reproduces the published count
  without it. A tRNA with an intron would be reported missing rather than silently modelled.

**What the names mean.** ARAGORN prints an amino acid in words and the anticodon as DNA —
`tRNA-Met … (cat)`. The published convention is one letter and RNA (`trnM-CAU`), so the adapter
maps Met→M and `cat`→`CAU` and then hands the result to `core/gene_names.normalize_name`, the same
normalizer a downloaded record goes through. There is no private rename: what the normalizer says
is the name, and anything it cannot resolve is reported as it arrived.

**Licence.** GPLv3 (the package metadata states it; ARAGORN's own documentation asks to cite
Laslett & Canbäck 2004 and 2008).
