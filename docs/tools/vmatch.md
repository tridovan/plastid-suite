# Vmatch 2.3.1 — declared, not vendored, and deliberately not a module's engine

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo. Vmatch is used by the **repeats** analysis (M2.12, decision **D-f**)
— but as the *independent second engine* the published count is checked against, never as the
engine a module runs on:

* **`repeat-match` (MUMmer 4.0.1) ships.** It is a plain binary with no licence constraint, and
  it is in `environment.yml` with the rest of MUMmer. The minimum install keeps working.
* **Vmatch is academic-only.** Its licence covers academic use; anything else needs a commercial
  licence from the authors. A toolkit module that *required* it could not be released to everyone
  (M9.3's licence review), so the module does not: `external/vmatch.py` exists for the
  cross-check, and the release note has to say so. Since M11.19 it lives in
  `environment-check.yml`, not in the runtime file, and `tests/unit/test_tool_inventory.py`
  fails if any product module imports it.

**What it is.** Vmatch (Stefan Kurtz) finds maximal exact matches — `mkvtree` builds the index,
`vmatch` queries it. `-d` reports **direct** repeats and `-p` **palindromic** (inverted) ones,
which are exactly the two classes `repeat-match` reports, so the two engines' counts are
comparable at all. It is the same software family as REPuter.

**Measured 10/10/2026 on the dev WSL machine:**

| | |
|---|---|
| Software | Vmatch **2.3.1**, build 2020-01-06 (`vmatch -version`, `mkvtree -version`) |
| Package | bioconda `vmatch=2.3.1` — `http://www.vmatch.de` is the author's page |
| Path here | `/root/toolbox/env/vmatch/bin/` — `mkvtree`, `vmatch` (its own env: installing it into the shared toolbox env risked moving the pinned `mafft`/`mummer4` builds) |
| Licence | academic only — see above |
| Measured on | the golden: **9** matches at ≥30 bp, the same nine `repeat-match` gives, coordinates and classes agreeing row for row |

**Install line** (its own environment, so the pinned builds in the toolbox env cannot move):

```sh
micromamba create -y -p /root/toolbox/env/vmatch -c conda-forge -c bioconda vmatch=2.3.1
```

The test that compares the two engines finds it through `PTS_VMATCH_DIR` (like `PTS_MUMMER_DIR`,
`PTS_BLAST_DIR` and the others), and **skips with a named reason** when it is absent: this machine
is the only one in the project that has it, and a check that quietly turned into "one engine
agreed with itself" would be worse than no check.

**The invocation the adapter uses:**

```
mkvtree -db <genome.fa> -dna -pl -allout -indexname idx
vmatch -l <floor> -d -p -nodist idx
```

`-nodist` prints one line per match instead of an alignment dump. **Coordinates are 0-based
starts**, unlike `repeat-match`'s 1-based ones; the conversion lives in
`external/vmatch.py` and nowhere else, which is what makes the two engines' rows comparable.

**Not vendored, deliberately:** compiled binaries, resolved per platform by bioconda — the same
reasoning as MUMmer, BLAST+ and ARAGORN.

**Reference.** Kurtz S, Choudhuri JV, Ohlebusch E, Schleiermacher C, Stoye J, Giegerich R.
REPuter: the manifold applications of repeat analysis on a genomic scale. Nucleic Acids Research
29:4633–4642 (2001). doi:10.1093/nar/29.22.4633
