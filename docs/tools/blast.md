# BLAST+ — declared, not vendored

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo — `environment.yml` (pinned), `vendor/` + `SOURCE.json` (small
scripts), or this file. BLAST+ is a compiled NCBI tool set that conda ships, so it is pinned in
`environment.yml` and described here. Until P3A-03 it was only a transitive requirement of PGA;
now `modules/mito_annotate.py` runs it directly, which makes its version part of every annotation
the toolkit produces.

**What it is.** NCBI BLAST+ — the reference implementation of sequence similarity search. The mito
caller uses three programs:

* `makeblastdb` — builds the searchable database from the committed reference FASTA;
* `blastx` — translates the query genome in all six frames and aligns it to the reference
  **proteins**, which is how a CDS's coordinates are transferred;
* `blastn` — aligns nucleotides to nucleotides, used for the three rRNA genes.

**Versions and provenance, measured 2026-10-03 on the dev WSL machine:**

| | |
|---|---|
| Software | BLAST+ **2.17.0** (`blastx: 2.17.0+`, build 2026-08-22) |
| Binaries | `blastx` 260,632 B · `blastn` 264,728 B · `makeblastdb` 207,384 B (the `bin/` tree is 52 MB) |
| Path here | `/root/toolbox/env/pts/bin/` — **not on `PATH`**, so a run is pointed at it with `--set blast=<directory>` |
| Checked by | `external/blast.py` refuses a different release by name: a different BLAST is a different annotation |

**Install line** (the route this machine used; also the route `environment.yml` encodes):

```sh
micromamba create -n pts -c conda-forge -c bioconda blast=2.17.0
# then either put its bin/ on PATH, or pass --set blast=<env>/bin to a run
```

**Not vendored, deliberately:** 52 MB of platform-specific NCBI binaries against a repository whose
whole test dataset is under a megabyte — and conda already resolves it per platform (Linux,
Windows-via-Docker), which a vendored Linux binary would not.

**Three decisions in the adapter were measured, not assumed.** Each of them fails silently when it
is wrong, so each is a test in `tests/unit/test_blast.py`:

1. **`-seg no`.** SEG masks low-complexity stretches of the query — sensible when hunting an
   unknown protein, wrong when transferring an annotation to a plant mitogenome, which is wall to
   wall hydrophobic membrane proteins. Measured on the answer-key record: with the default mask,
   `atp9`'s best alignment drops from 137 to 67 bits and 21 other loci lose score too. The
   identity floor, not the mask, is what keeps junk out.
2. **No `-parse_seqids`.** The obvious flag for a FASTA whose ids are `accession|gene|…` is the one
   that breaks: measured 03/10/2026, `makeblastdb -parse_seqids` fails outright with
   `Error: Duplicate seq_ids are found: LCL|NC_006581.1|SDH3|-|NICOTIANA` — the reference set
   legitimately holds repeat copies of a gene — and for the ids it does parse it **upper-cases**
   them, which would break every lowercase gene symbol. Without it, `sseqid` is the defline's first
   token, which is exactly the subject key the caller wants.
3. **`-task blastn` for rRNA.** Megablast and blastn found the same three rRNAs with the same spans
   on the answer-key record, but a real run's target may be a divergent sample and the database is
   25 sequences: the sensitive task costs almost nothing here. DUST masking on/off was measured to
   give bit-identical output, so the default is kept.

**One behaviour of the tool's own output that everything depends on:** BLAST truncates `sseqid` at
the first whitespace. The reference headers carry the organism name, so
`NC_006581.1|matR|matR|Nicotiana tabacum` arrives as `NC_006581.1|matR|matR|Nicotiana` — a subject
key must be built from the first three fields, and `external/blast.py` refuses an id with fewer.

**Licence.** Public domain (NCBI, U.S. Government work).
