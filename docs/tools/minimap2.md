# minimap2 2.31 — the map-back of an assembly

**Why this file exists.** `plastid assemble` tests an assembly by mapping the reads back onto it,
and a module's tool is declared in the repository. It lives in `environment-assembly.yml`.

**Why this aligner.** GetOrganelle recruits its reads with Bowtie2, so a check built on the same
aligner and the same kind of index would share its blind spots. The check has to be independent of
the thing it checks, and for the same reason it reads the raw reads, not fastp's output.

**Install line:**

```sh
conda env create -f environment-assembly.yml
```

**The invocation the adapter uses:**

```
minimap2 -ax sr --eqx --secondary=no --sam-hit-only -t N -o mapped.sam target.fasta R1 R2
```

| Option | Reason |
|---|---|
| `-ax sr` | the short-read preset, the same for every sample |
| `--eqx` | writes `=`/`X` CIGAR operators, so mismatches are counted from the CIGAR alone. The reader refuses a plain `M` |
| `--secondary=no` | one placement per read. A read inside an inverted repeat has two equal placements; its primary falls on either copy, which spreads the repeat's reads over both and leaves the depth level |
| `--sam-hit-only` | leaves unmapped reads out. Most reads of a genome skim are not plastid, and writing them would make the SAM as large as the input |

**The target is the sequence plus its own first 1,000 bases**, because a circle has no end and a read
across the origin needs somewhere to land. Positions are folded back onto the circle, and each read
counts once because only primary alignments are read. The SAM is read directly, with no `samtools`
step, and is deleted with the run's working files.

**The version is checked** (`MINIMAP2_VERSION` in `src/plastid/external/minimap2.py`).

Citation: Li H. *Bioinformatics* 34:3094–3100 (2018).
