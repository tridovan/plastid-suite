# IQ-TREE 2.3.6 — declared, not vendored (M11.21)

**Why this file exists.** CLAUDE.md's invariant: the moment a module uses an external tool, that
tool is declared in the repo — `environment.yml` (pinned), `vendor/` + `SOURCE.json` (small
scripts), or this file. IQ-TREE is a compiled binary set that bioconda ships, so it is pinned in
`environment.yml` and described here. It is the **inference step of the `phylogeny` module**
(M11.21), and the version is part of every tree: ModelFinder's model choice and the ultrafast
bootstrap implementation have both changed between 2.x builds.

**What it is.** IQ-TREE 2 (Bui Quang Minh et al., 2020) — maximum-likelihood phylogenetic
inference with ModelFinder model selection and ultrafast bootstrap. The `phylogeny` module follows
the route CGAS's module 14 uses: MAFFT per protein-coding gene → concatenated supermatrix →
IQ-TREE.

**Install line** (its own environment; nothing here may move the pinned builds in the shared
toolbox env):

```sh
micromamba create -y -p /root/toolbox/env/iqtree -c conda-forge -c bioconda iqtree=2.3.6
# or, if using environment.yml:
conda env create -f environment.yml
```

**Measured 28/10/2026 on the dev WSL machine:**

| | |
|---|---|
| Version | `iqtree2 --version` → **IQ-TREE multicore version 2.3.6 for Linux x86 64-bit built Dec 11 2024** |
| Package | bioconda `iqtree=2.3.6`; the executables are `iqtree2` and `iqtree` (both present) |
| Path here | `/root/toolbox/env/iqtree/bin/iqtree2` |
| Licence | GPL-2.0-or-later (the IQ-TREE repository's own licence) |
| Measured on | 4 panel records (`NC_008641.1`, `NC_016708.2`, `NC_028027.1`, `NC_016736.1`): 78 CDS blocks, 70,079 sites, tree with UFBoot/SH-aLRT support values |

**The invocation the adapter uses:**

```
iqtree2 -s alignment.fasta --prefix tree -m MFP -bb 1000 -alrt 1000 -bnni -nt AUTO -seed 1 --redo
```

| Flag | Why |
|---|---|
| `-m MFP` | ModelFinder Plus: the model is chosen from the data rather than fixed. The choice is written into `tree.iqtree` and recorded |
| `-bb 1000` | 1000 ultrafast bootstrap replicates. Below 1000 the adapter **refuses** — a support value from fewer replicates is not one this toolkit publishes |
| `-alrt 1000` | 1000 SH-aLRT replicates, shown beside UFBoot as `UFBoot/SH-aLRT` |
| `-bnni` | optimise the bootstrap trees by NNI, the setting IQ-TREE recommends with `-bb` |
| `-nt AUTO` | threads; recorded in argv, and it does not change the topology |
| `-seed 1` | a fixed seed, so two runs of one matrix give one tree. Changing it is a deliberate act and shows in argv |
| `--redo` | an existing `--prefix` output is overwritten rather than silently reused |

**What the adapter does not do.** It does not choose an outgroup or root the tree — a tree drawn
from an unrooted ML result needs the user's own outgroup, and `docs/things_we_do_not_do.md`
records that the toolkit does not infer one from taxonomy. The `r_tree` figure takes the
`.treefile` this module writes and still asks for `--set root/outgroup` explicitly.

**Not vendored, deliberately:** a compiled binary, resolved per platform by bioconda — the same
reasoning as MUMmer, BLAST+ and ARAGORN.

**Reference.** Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A,
Lanfear R. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic
era. Molecular Biology and Evolution 37:1530–1534 (2020). doi:10.1093/molbev/msaa015
