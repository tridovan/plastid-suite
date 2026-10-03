# PMGA (MGAVAS) — declared, and measured as far as it can be without running it

**Why this file exists.** CLAUDE.md: a tool that exists only on one machine is not part of the
toolkit, and a tool a module uses must be declared in the repo — in `environment.yml`, in
`vendor/`, or here. PMGA is neither on conda nor a small script: it ships as a **Singularity
container**, so it is declared here with its source, its real command line, and what it costs.

**What it is.** PMGA (**P**lant **M**itochondrial **G**enome **A**nnotator) — internally called
**MGAVAS** (Mitogenome Annotation, Visualization, Analysis and Submission). Li et al.,
*Plant Communications* 2025, 6:101191, <https://doi.org/10.1016/j.xplc.2024.101191>. It annotates
a plant mitogenome: protein-coding genes, tRNA, rRNA, repeats, and RNA-editing sites. It is the
tool D3-01 chose for the mito path, after MITOS2 was rejected — MITOS2's reference set
(`refseq81m`) is metazoan, so on a plant mitogenome it would return *wrong* genes, not none.

**Artifact, measured 2026-10-03:**

| | |
|---|---|
| Record | Figshare 27201798, "Source code of PMGA", published 2024-10-11, **GPL-3.0+** |
| Download | `https://ndownloader.figshare.com/files/49723380` — **3,176,338,860 B (3.18 GB)**, sha256 `81a23743eb6786d9…` |
| What it unpacks to | `PMGA_tools/` — a **Singularity container filesystem** (`.singularity.d`, `bin`, `lib`, `usr`, `apps`), not source code and not a `.sif` image |
| Unpacked size | **more than 6 GB**, 148,650 entries |
| Inside | a conda environment (`apps/miniconda3`, env `cpgavas2`), the pipeline's own scripts (`apps/mgavas`, 78,816 entries), and its reference data (`modules/mganno`, 64,833 files) |
| Inner git | one commit, `464166cf…`, 2022-01-27 — the pipeline itself has no released version number |

**The real command line**, quoted from the container's own POD (`/usr/local/bin/mgavas_n`):

```
mgavas -pid <id> -in <fasta> -db <1|2|4> [-conf <file>] [-ref <genbank>] [-outdir <dir>]
```

* `-pid` — required, alphanumeric; names the output files;
* `-in` — required, a FASTA file. Multi-contig input is sorted, split and **concatenated into one
  sequence** before annotation (their help page warns that the concatenated file in the results
  folder must then not be submitted to NCBI);
* `-db` — required: which reference dataset to use (`mgavas_n` accepts 1–6 and its error message
  names 1, 2 and 4; `mgavas_m` accepts 1, 2, 4, 5, 6);
* `-conf` / `-ref` — optional parameter file and reference GenBank; both default to the bundled
  `sample_files/`;
* `-outdir` — optional, defaults to `/tmp/dir_<pid>`.

Entry-point chain: `/usr/local/bin/mgavas` → `mgavas_n` (Perl) → `mgavas_s` (bash, sources
`/apps/mgavas/.rc_bash`) → `perl /apps/mgavas/bin/mgavas_m`. **Call `mgavas` or `mgavas.sh` from
inside the container**, not `mgavas_m`: the wrapper is what sets the environment the tools need.

**What it runs underneath** (counted in its own scripts): **MAKER** (the gene-model engine, ~50
call sites), vmatch, ARAGORN, tRNAscan-SE, MISA, TRF, BLASTn/BLASTx, EMBOSS, HMMER, circos, and
Deepred-mt for RNA-editing prediction. Sample inputs are bundled: *Arabidopsis thaliana*
mitochondrion (NC_037304) and plastome (NC_000932).

## Can it go in the app, and can an ordinary computer run it?

**No, and not everywhere — both measured, not guessed.**

* It is a **container**, so it needs Singularity/Apptainer (Linux; on Windows that means WSL2 or
  Docker, i.e. a container inside a container), plus **>6 GB of disk** and a working directory.
* It **bundles ten tools including MAKER**, which is the heaviest step in organelle annotation.
  The project's target machines are 4–8 GB of RAM (anh Điềm's is 3.7 GB).
* It therefore belongs to the **declared-external-tool tier**, exactly like Chloe and table2asn —
  never to the pip-only tier, which must keep working with no external tool at all.

**Two ways to use it, and the second is the product answer** (the same reasoning D3-01 recorded
for GeSeq):

1. **Container, locally** — for a machine that can afford it, declared here and called through
   `run_tool()`; the toolkit records name · version · path · argv like any other tool.
2. **The web service, once, pinned** — `http://www.1kmpg.cn/pmga/` annotates the same inputs. Run
   it once per species, then **pin the output** (`plastid panel pin`: accession · date · sha256)
   and let the toolkit read it with its own reader. The manual step happens once; every later run
   is offline, and the biologist's computer never has to hold the container.

## What is NOT yet measured — and what it needs

**Nobody has run it here.** The container is downloaded and unpacked, its command line is read
from the container itself, but `apptainer exec PMGA_tools/ …` was **refused by the permission
classifier**: it is third-party code that the user has not named in this session, which is the
same boundary that stopped PGA earlier (and was cleared then by naming it). Until it runs, these
are open:

* whether it runs at all on a 7 GB machine (and, in particular, whether its conda environment
  survives WSL2);
* **peak RAM and wall-clock time** for one mitogenome — the number that decides question 1 above;
* whether the gene set it returns agrees with `modules/mito_genes.py` on the stand-in
  `OX459128.1` (51 = 31 CDS + 17 tRNA + 3 rRNA, PROJECT.md §8.14).

`external/pmga.py` is **not written yet, deliberately**: an adapter is written against observed
argv and observed output file names, and neither has been observed. The moment the container runs,
the adapter is a small job — the shape of `external/pga.py`, with `-pid`/`-in`/`-db` and the
output directory read back through our own reader.
