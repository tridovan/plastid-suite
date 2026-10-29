# CPStools 3.0 — the cross-check tool for IR · π · RSCU · SSR · KaKs (P8-02, D-j)

**Why this file exists.** CLAUDE.md's invariant: the moment a tool is used, it is declared in the
repo. CPStools is **not a module dependency** — no module of this toolkit calls it, and the minimum
install stays pip-only. It is the other side of M8.1's cross-check, so it is declared here, with the
version and the argv every cross-check number was produced with. The shape is Vmatch's (D-f): a
comparison tool the release notes have to name.

**What it is.** CPStools (Xu Wenbo, China Pharmaceutical University → CACMS) — a package for
analysing chloroplast genomes: IR boundaries, π, RSCU, SSRs, Ka/Ks, GC. MIT licence. Published as
CPStools: a package for analyzing chloroplast genome sequences, doi:10.1002/imo2.25.

**Install line** (its own environment — pip, as D-j decided; it must not move the pinned builds in
the shared toolbox env):

```sh
uv venv /root/toolbox/env/cpstools --python 3.12
uv pip install --python /root/toolbox/env/cpstools/bin/python "cpstools==3.0"
```

**Two more tools came with the cross-check, each in its own environment:**

```sh
# cpstools' KaKs route shells out to KaKs_Calculator (it does not bundle one)
micromamba create -y -p /root/toolbox/env/kakscalculator -c conda-forge -c bioconda kakscalculator2=2.0.1
# MAFFT 7.526 — the project's own pinned build, reused through -m
```

**Measured 10/10/2026 on the dev WSL machine:**

| | |
|---|---|
| Version | `cpstools -v` → **CPStools version 3.0** |
| Package version attribute | `cpstools.__version__` says **1.0.0** — the CLI and the module disagree; the CLI's is the release (PyPI's newest is 3.0) |
| Dependencies it pulls | biopython ≥1.84, numpy, pandas, matplotlib, tqdm |
| KaKs_Calculator | **2.0.1** (bioconda `kakscalculator2`), `KaKs_Calculator` on `PATH` |
| `environment-check.yml` | `cpstools=3.0` (bioconda) — the same release this pip install is. It left `environment.yml` in M11.19: the runtime file is the analyses, and no module calls this |

**The argv the cross-check used:**

```sh
cpstools IR   -i NC_050833.1.gb
cpstools SSRs -i NC_050833.1.gb
cpstools RSCU -d <dir>            # default -l 300;  and again with -l 0
cpstools Pi   -d <absolute dir> -m /root/toolbox/env/pts/bin/mafft
cpstools KaKs -i <reference.gb> -m YN
```

**Defects measured in the tool itself** — recorded because a cross-check tool's behaviour is data to
explain, never a benchmark (a lesson from reading CGAS's source):

1. **`Pi` only runs when `-d` is absolute.** The work dir is resolved with `abspath()` when
   `gene_cp_sort.txt` is *written*, and with the raw argument when it is *read* — so `-d .` writes
   the file one directory up and then crashes reading it (`FileNotFoundError`). Measured twice.
2. **It refuses to overwrite.** A second run stops with `The save directory has been existed,
   please delete the directory!` — outputs land **next to** the input directory
   (`common_gene/`, `IGS/`, `info/`), not inside it.
3. **`KaKs` reports success when it failed.** With an invalid `-m` it prints one
   `File not found error for …KaKS_results/gene_KaKs.txt` per gene and still ends
   `All KaKs calculations are completed and combined results have been saved.` Only `-m YN`
   (a real KaKs_Calculator method) produces values.
4. **`RSCU` counts the three stop codons** in its RSCU table (76 of them on the golden), and its
   default `-l 300` filter drops 28 of the golden's 79 loci. It does print what it dropped.
5. **Its common-gene set is an intersection by name**, so two records that spell one locus
   differently lose it from both sides — measured: `atpI`, `psaI`, `psbB` are in neither the gene
   table nor any spacer name, because NC_024936.1 spells them `atpL`, `psaL`, `bspB`. This is
   **the gene-set-by-name trap the Ka/Ks module hit**, reproduced independently by another tool.

**Not vendored, deliberately:** a pip/conda package with its own dependency tree, resolved per
platform — the same reasoning as MUMmer, BLAST+ and ARAGORN.
