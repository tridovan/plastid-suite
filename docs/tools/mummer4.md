# MUMmer 4.0.1 — declared, not vendored

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo — `environment.yml` (pinned), `vendor/` + `SOURCE.json` (small
scripts), or this file. MUMmer 4 is a set of compiled C++ binaries that bioconda ships, so it is
pinned in `environment.yml` and described here. It is the aligner of Phase 4's synteny module
(P4-01, decision **D4-02**): `nucmer` → `delta-filter` → `show-coords` is the pipeline PROJECT.md
§7 records for the %identity-by-coordinate (mVISTA-style) result, and another tool would make the
phase check unreproducible.

**What it is.** MUMmer 4 (Marçais et al. 2018) aligns whole genomes by exact-match anchoring
followed by extension. `nucmer` writes a **delta** file — the coordinates of every alignment plus
each alignment's internal gaps — and the sibling utilities filter (`delta-filter`) and tabulate
(`show-coords`) it. The toolkit reads both the tabulated blocks (the tool's own %IDY) and the
delta itself (the per-window series), and refuses a run where the two disagree
(`external/mummer.py`).

**Versions and provenance, measured 2026-10-04 on the dev WSL machine:**

| | |
|---|---|
| Software | MUMmer **4.0.1**, build `pl5321h9948957_0` |
| Package | bioconda, 877 kB (plus perl and ruby, which its scripts use) |
| Path here | `/root/toolbox/env/pts/bin/` — `nucmer`, `delta-filter`, `show-coords` |
| Measured on | the golden pair (*I. domestica* × *I. gatesii*): **7 blocks, %IDY 89.66–99.11 (mean 97.37)**, 152,332 paired bases, coverage 99.09 % |

**Install line:**

```sh
micromamba install -c conda-forge -c bioconda mummer4=4.0.1
```

**Not vendored, deliberately:** compiled binaries per platform, and bioconda resolves the right
build per platform — the same reasoning as ARAGORN and BLAST+. The Windows/Docker route
(ROADMAP M9.2) gets it from the same channel.

**The invocation the adapter uses, and why:**

```
nucmer --maxmatch -l 20 -c 65 -t 1 -p out <target.fa> <query.fa>
delta-filter -1 out.delta > out.filter.delta
show-coords -rcl out.filter.delta
```

* `--maxmatch` — every maximal exact match is an anchor, not only matches unique in both
  sequences. It is the flag PROJECT.md §7 records, and it is what a plastome comparison wants:
  the inverted repeat is by definition a repeat, and uniqueness filtering would drop half of it.
* `-l 20 -c 65` — MUMmer 4's own defaults, **passed explicitly**: they live in
  `data/defaults.yaml` (`synteny:`) because the block count moves with them (measured: 12 raw
  alignments at the defaults, 22 at `-c 40`) and a default is a number nobody wrote down.
* `-t 1` — threads must not be able to move a number a manuscript rests on.
* `delta-filter -1` — 1-to-1: each reference position maps to its best query position. The other
  modes answer different questions and can list one target region in several rows, which would
  make a coverage number lie.
* `show-coords -rcl` — sorted by reference, with coverage and sequence lengths. `-rcl` is the
  exact flag set §7 recorded.

**A version mismatch is refused, not noted** (`nucmer_version`): the three binaries are also
checked to come from one directory, because a mixed installation would be a silent inconsistency
nobody could name afterwards.

**What was measured about the format** (recorded because it is undocumented in the pipeline's own
files, and a wrong reading is silent):

* MUMmer 4 writes **one `>` header per sequence pair**, with every alignment of that pair as a
  7-number record terminated by a `0` line — a reader that expects a header per record returns
  one block for a comparison that has seven.
* In the gap list, a value `d` consumes `|d| − 1` *paired* columns, then `d > 0` advances the
  target alone and `d < 0` the query alone. Verified on all twelve raw alignments of the golden
  pair: every block lands on its own end coordinates and `mismatches + unpaired bases` equals the
  header's error count.

**Licence.** MUMmer 4 is released under the Artistic License 2.0 (the bioconda package metadata
states it; the upstream repository carries the LICENSE file). Flagged for the M9.3 review like
every other tool.

**Reference.** Marçais G, Delcher AL, Phillippy AM, Coston R, Salzberg SL, Zimin A. MUMmer4: a
fast and versatile genome alignment system. PLOS Computational Biology 14:e1005944 (2018).
doi:10.1371/journal.pcbi.1005944
