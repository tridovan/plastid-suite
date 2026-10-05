# Chloe — declared, not vendored

**Why this file exists.** CLAUDE.md: a tool that exists only on one machine is not part of the
toolkit. The moment a module uses an external tool, that tool must be declared in the repo — in
`environment.yml`, in `vendor/` with a `SOURCE.json`, or here. Chloe fits neither of the first two:
it is a Julia package plus a 1.7 MB reference artifact that Julia fetches on first use, so it is
declared here, with its source and its install line, and every run records name · version · path ·
argv (see `src/plastid/external/chloe.py`).

> **Pipeline decision changed 05/10/2026 (D5-08):** Chloe is dropped from the pipeline — the
> annotate step becomes this project's own tool (`docs/D5-08_plastid_annotator.md`, phase 3B in
> `ROADMAP.md`). This file stays as the tool's record while `modules/annotation.py` still offers
> `--set annotator=chloe` / `pga`, until the new engine passes its Done-when.

**What it is.** An offline chloroplast annotator from Ian Small's group (UWA / Curtin), written in
Julia. It aligns a target plastome against its own curated reference set and projects gene models,
which is why it does not need the caller to supply references the way PGA does.

**Versions and provenance, measured 2026-10-03 on the dev WSL machine:**

| | |
|---|---|
| Software release | **0.1.18** — the checkout's `Project.toml`; this is the number the toolkit pins |
| GFF's own declaration | `##source-version Chloe 1.0` — a second, different number; also pinned |
| Citable artifact | Zenodo software deposit, **doi:10.5281/zenodo.18194040**, MIT licence |
| Checkout | `git clone https://github.com/ian-small/chloe`, commit `1d7b710` (2025-06-09) |
| Julia | **1.13.0**, run through `micromamba` |
| Reference set | `Chloe_refs` **1.0.1**, 23 chloroplast references + 6 rDNA ones, 1.7 MB |
| Runtime | ~10 s of annotation + ~12 s of Julia start-up and compilation — **22 s** wall, once per genome |

**Install line** (the route this machine used, 2026-10-03):

```sh
micromamba create -y -n chloe -c conda-forge julia
git clone https://github.com/ian-small/chloe
cd chloe && julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` is where the reference set is downloaded (once, from
`https://github.com/ian-small/chloe_references/archive/refs/tags/v1.0.1.tar.gz`); after that the
annotator runs **offline**, which is the non-negotiable the toolkit is built around. Julia's own
documented route is [juliaup](https://julialang.org/downloads/) if you prefer it to conda.

**Not vendored, deliberately:** Julia itself (a 500 MB toolchain) plus up to 750 MB of Julia
artifacts on this machine. The clone is 6 MB but useless without them.

**The reference set is pinned by its tree hash**, `9163915d99fdfbc7d8510f0479b3d08905682ca5` —
also the name of its directory under `~/.julia/artifacts/`. A different hash is a different
reference set, so the adapter refuses the run rather than annotating from references nobody can
name. **Re-measuring after an update means updating two constants** (`CHLOE_VERSION`,
`CHLOE_REFS_TREE` in `external/chloe.py`) against the golden genome, and recording the result in
PROJECT.md §8.

**Scope, for the papers.** Chloe annotates **angiosperm chloroplast genomes only** — its README
says so and offers to contact the authors about other plants. For a gymnosperm, fern or
lycophyte plastome, use the `pga` annotator instead (the same module, `--set annotator=pga`).

**Licence.** The repository carries **no LICENSE file**; the Zenodo deposit for version 0.1.18
declares **MIT**. Both facts are recorded here for the M9.3 licence review — a declaration in a
deposit is not the same artifact as a licence file in the tree we clone, and saying which is which
is the point of writing it down.

**Three of its habits are handled in the adapter**, each measured rather than assumed:

1. **It transforms the sequence by default.** `--no-transform` is documented as "do not flip and
   orient sequence to standard configuration", and with the transform on, the GFF's coordinates
   belong to a sequence the caller never had. Measured on the golden genome: the single-copy
   region (18,169 bp) is reverse-complemented — 13,001 of its positions differ — and the gene
   set is *identical* (same 132 genes, same names). A gene count cannot see this; the ORF checks
   can (12 of 86 CDS stop translating). The adapter passes `--no-transform` and refuses any run
   that writes a transformed FASTA.
2. **It names itself in its output, not on the command line.** A script is never probed with a
   flag here (PGA's lesson: an unknown flag made it run itself with its defaults), so the release
   is read from `Project.toml` and the format from the GFF's `##source-version` line.
3. **Its output carries no `/transl_table`.** The caller declares the code (plastids: 11) and the
   adapter stamps it on the CDS that declare none; the note in every annotation report says so.
