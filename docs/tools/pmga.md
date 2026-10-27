# PMGA (MGAVAS) — declared, and measured as far as it can be without running it

> **Pipeline decision changed 05/10/2026.** The mito annotate step is this project's **own**
> reference-transfer caller (`core/mito_caller.py` + `modules/mito_annotate.py`, Phase 3A) —
> PMGA never ran here and **nothing waits on it** (H-05: its web service still refuses
> connections). "Route 2" below stays as the record of the 03/10 decision and as the fallback if a
> machine that can run the container ever appears; `external/pmga.py` remains unwritten on purpose.

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
| Unpacked size | **11 GB**, 138,045 files |
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
  The project's target machines are 4–8 GB of RAM (Điềm's is 3.7 GB).
* It therefore belongs to the **declared-external-tool tier**, exactly like Chloe and table2asn —
  never to the pip-only tier, which must keep working with no external tool at all.

**This is not a new question: the plan already answered it** (P1-12 — the backend
rule): *a heavy step runs where it fits*. Locally on an ordinary machine when it can, on a web
service or a bigger server when it cannot. The same rule D3-01 applied to GeSeq applies here:

1. **Container, locally** — for a machine that can afford it, declared here and called through
   `run_tool()`; the toolkit records name · version · path · argv like any other tool. It is the
   route for whoever has the RAM; the paper's own authors ship the container for exactly that.
2. **The web service, once, pinned** — `http://www.1kmpg.cn/pmga/` annotates the same inputs. Run
   it once per species, save the GenBank it returns, and **pin that file**:

   ```
   plastid panel pin --file pmga_output.gb --note "PMGA web, 2026-10-03"
   ```

   The pin records what an NCBI download records — bytes · sha256 · md5 · the record re-read with
   our own reader — plus where it came from. Every later run is offline, and the biologist's
   computer never holds the container. This is the route for the team's own machines, and the
   `--file` half was built on 2026-10-03 precisely because the rule requires it.

**Decided 03/10 (dev, confirmed by Trí): route 2 is what this project uses, and
`external/pmga.py` is not written.** The reason is measured, not a preference: the container is
11 GB unpacked, it drives MAKER, the WSL VM here is capped at 8 GB and the C: drive had 42 GB
free — **this machine is not strong enough to run the annotator.** The adapter stays open with a
named trigger: the first machine that *can* run the container. It is a small job then, because the
command line above is already read from the container itself.

## What is NOT yet measured — and what it needs

**Nobody has run it here, and the container is no longer on the machine.** It was downloaded and
unpacked (3.18 GB → 11 GB, 138,045 files), its command line was read out of it, and it was then
**deleted**: measured the same day, the WSL virtual disk reports ~1 TB free while the real C:
drive had **42 GB left (89 % full)**, and `ext4.vhdx` does not shrink when files are deleted
inside the guest — WSL refuses its own sparse mode ("currently disabled due to potential data
corruption"). Fourteen gigabytes was a third of what was actually free, for a tool this machine
cannot comfortably run. Re-downloading it is one `curl` away (the URL and sha256 are above), so
nothing is lost — and the whole point of route 2 below is that this machine never needs it.

`apptainer exec PMGA_tools/ …` was also **refused by the permission classifier** at the time:
third-party code the user had not named in-session, the same boundary that stopped PGA earlier.
Until it runs somewhere, these stay open:

* whether it runs at all on a 7 GB machine (and, in particular, whether its conda environment
  survives WSL2);
* **peak RAM and wall-clock time** for one mitogenome — the number that decides question 1 above;
* whether the gene set it returns agrees with `modules/mito_genes.py` on the stand-in
  `OX459128.1` (51 = 31 CDS + 17 tRNA + 3 rRNA).

`external/pmga.py` is **not written yet, deliberately**: an adapter is written against observed
argv and observed output file names, and neither has been observed. The moment the container runs,
the adapter is a small job — the shape of `external/pga.py`, with `-pid`/`-in`/`-db` and the
output directory read back through our own reader.
