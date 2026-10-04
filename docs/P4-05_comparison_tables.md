# P4-05 — the comparison tables for the two manuscripts

**What this is.** The phase-4 deliverable: one run per manuscript species over the **pinned
stand-in panel**, plus the mito sketch, each writing the table set the manuscripts' comparison
sections are built from. The runs are the deliverable; this file records how to reproduce them
and what they contain.

**Data, stated once.** Everything here comes from `.panel/` — public GenBank records pinned with
their sha256 (PROJECT.md §8.6), **not the team's samples**. The one exception is the Coleus run's
reference, `PX923121.1`: it is the group's own GenBank record of *C. amboinicus* (HANDOFF session
5), and its run label says so.

Every table carries its run's label — `stand-in (public records, not the team's samples)` for
Oldenlandia, the same plus the reference caveat for Coleus, and a cross-genus-sketch wording for
the mito run — because `--label` writes it into each table's scope line (the line the CSV header,
the JSON export and the report all read) and into every report's header. When the team's own
samples arrive, the numbers change; the label does not until the data does.

## The runs

The commands are exact — copy them, or read them out of each run's `manifest.json`, which records
the inputs' sha256, the parameters and every tool's argv.

```sh
TOOLS=<directory with mafft, nucmer, delta-filter, show-coords, yn00>
OL="stand-in (public records, not the team's samples)"
CO="stand-in panel (public records) — the reference PX923121.1 is the group's own GenBank record of C. amboinicus"
MT="stand-in (public records) — cross-genus mito sketch; AP045469.1 carries this project's own prototype annotation (P3A-07)"
P=.panel/oldenlandia-coleus
M=.panel/oldenlandia-coleus-mito

# mito — reference OX459128.1, query this project's own prototype annotation of AP045469.1
plastid run -m kaks -m popgen \
  -i $M/OX459128.1.gb -i .panel/mito-prototype/AP045469.1.gb \
  -o .runs/phase4 --dataset standin-mito --label "$MT" \
  --groups data/oldenlandia-coleus-mito.yaml --set paml=$TOOLS --force

# Coleus — reference PX923121.1 (C. amboinicus, the group's own record), the other ten as queries
plastid run -m substitutions -m kaks -m synteny -m barcodes -m popgen \
  -i $P/PX923121.1.gb -i $P/PZ428230.1.gb -i $P/NC_069583.1.gb -i $P/AP045470.1.gb -i $P/NC_058328.1.gb -i $P/NC_057679.1.gb \
  -i $P/NC_057985.1.gb -i $P/PP078674.1.gb -i $P/PQ325515.1.gb -i $P/NC_057983.1.gb -i $P/PX977793.1.gb \
  -o .runs/phase4 --dataset standin-coleus --label "$CO" \
  --groups data/oldenlandia-coleus.yaml \
  --set scope=cds --set mummer=$TOOLS --set paml=$TOOLS --force

# Oldenlandia — reference NC_057985.1 (O. diffusa, RefSeq stand-in), the other ten as queries
plastid run -m substitutions -m kaks -m synteny -m barcodes -m popgen \
  -i $P/NC_057985.1.gb -i $P/PP078674.1.gb -i $P/PQ325515.1.gb -i $P/NC_057983.1.gb -i $P/PX977793.1.gb \
  -i $P/PX923121.1.gb -i $P/PZ428230.1.gb -i $P/NC_069583.1.gb -i $P/AP045470.1.gb -i $P/NC_058328.1.gb -i $P/NC_057679.1.gb \
  -o .runs/phase4 --dataset standin-oldenlandia --label "$OL" \
  --groups data/oldenlandia-coleus.yaml \
  --set scope=cds --set mummer=$TOOLS --set paml=$TOOLS --force
```

`--force` is deliberate in the commands above: a cached module restores its files without running
its tools, so a cached run has no tool logs — and the logs are what makes the numbers checkable
from the run directory alone.

`--set scope=cds` matters: the substitutions module's default is the **whole-alignment** scope,
and the manuscripts quote the **CDS scope** (§6/§8.25/§8.27) — the tables carry whichever scope
they were run with, in their own scope line, and the two are never mixed.

## What each run directory holds

| Path | What it is |
|---|---|
| `01_substitutions/tables/substitutions_summary.csv`, `substitutions_by_genome.csv`, `substitutions_by_gene.csv` | Ts/Tv, CDS scope: pooled, per genome compared, and per gene — each with the scope label |
| `02_kaks/tables/kaks_summary.csv`, `kaks_pooled.csv`, `kaks_genes.csv` | Ka/Ks per gene and per genome, every gene keeping a row (scored, or with the reason it was not) |
| `03_synteny/tables/synteny_summary.csv`, `synteny_genomes.csv`, `synteny_blocks.csv`, `synteny_identity.csv` | nucmer blocks and the % identity series |
| `03_synteny/mummer_*.coords`, `mummer_*.delta` | the tool's own outputs, so the tables can be checked without re-running nucmer |
| `04_barcodes/tables/*.csv` | ranked windows, hypervariable regions, candidate loci |
| `05_popgen/tables/*.csv` | Hd · θ(S) · θ(π) · Tajima's D · Fu & Li · Fst per gene and group |
| `00_run_info/` | the manifest (input sha256, tool name · version · argv, IR evidence) and every tool's stdout/stderr log |
| `report/report.md` | the run's prose, every number a table cell |

**Re-derivable from the run directory alone.** The inputs' sha256 are in the manifest; every
tool's argv, version and **full output** are in `00_run_info/` (one numbered log pair per call —
after the fix of 04/10 that made each call keep its own; §9.C-22); the MAFFT logs hold the
per-gene alignments the Ka/Ks and popgen numbers were computed from, and the mummer artifacts hold
the alignment the synteny tables describe. Nothing here needs the network. *A run whose manifest
says a module is `cached` did not re-run its tools, so that run has no logs for it — read the
status column.*

## Measured results

**Mito — the cross-genus sketch (D4-05), run `standin-mito/20261004-143824-29db3666`.** Reference
`OX459128.1` (*Oldenlandia corymbosa*), query the project's own prototype annotation of
`AP045469.1` (*Coleus hadiensis*, P3A-07). 19 genes are carried by both; of those, **12 produced a Ka/Ks** (pooled Ka 0.028443 ·
Ks 0.112924 · ratio 0.2519), 6 are reported as `frameshift` and 1 has an internal stop codon —
and the other 12 names are listed too: 8 the prototype annotation does not carry, 4 the reference
does not. The popgen tables are a **demonstration, not a result**: one sequence per group, so
every statistic is blank with *"nothing is estimable"* and every Fst row says a group has fewer
than two sequences. That is the honest output of running the machinery on two mitogenomes from
different genera, one of them a prototype annotation.

**Oldenlandia — run `standin-oldenlandia/20261004-142600-a99f235f`.** Reference NC_057985.1 (*Scleromitrion diffusum*), 8 annotated genomes compared (the three records with no CDS are named in the notes).

*Ts/Tv, CDS scope:* **1.5157** pooled = 16959 Ts / 11189 Tv over **77 shared genes** (68743 columns). Per pair — the manuscript's own table:

| Genome | Group | Shared genes | Comparable columns | Ts | Tv | Ts/Tv |
|---|---|---|---|---|---|---|
| PP078674.1 | Oldenlandia | 77 | 63810 | 161 | 168 | 0.9583 |
| NC_057983.1 | Oldenlandia | 77 | 62552 | 921 | 587 | 1.5690 |
| PX977793.1 | Oldenlandia | 77 | 62610 | 925 | 597 | 1.5494 |
| PX923121.1 | Coleus | 77 | 62771 | 3683 | 2392 | 1.5397 |
| NC_069583.1 | Coleus | 77 | 63424 | 3769 | 2491 | 1.5130 |
| NC_058328.1 | Coleus | 77 | 62697 | 3684 | 2439 | 1.5105 |
| NC_057679.1 | Coleus | 77 | 63478 | 3816 | 2515 | 1.5173 |

*Ka/Ks* (PAML yn00 4.10.10, per shared CDS): **80 genes considered**, 77 carried by every annotated genome, 7 genomes produced pooled values —

| Genome | Group | Scored | Without a value | Pooled Ka | Pooled Ks | Pooled Ka/Ks |
|---|---|---|---|---|---|---|
| PP078674.1 | Oldenlandia | 72/77 | 5 | 0.002471 | 0.015586 | 0.1585 |
| NC_057983.1 | Oldenlandia | 69/78 | 9 | 0.008712 | 0.086485 | 0.1007 |
| PX977793.1 | Oldenlandia | 67/78 | 11 | 0.009031 | 0.091612 | 0.0986 |
| PX923121.1 | Coleus | 57/78 | 21 | 0.033998 | 0.398480 | 0.0853 |
| NC_069583.1 | Coleus | 57/78 | 21 | 0.034192 | 0.399155 | 0.0857 |
| NC_058328.1 | Coleus | 59/78 | 19 | 0.037611 | 0.392385 | 0.0959 |
| NC_057679.1 | Coleus | 57/78 | 21 | 0.034475 | 0.397337 | 0.0868 |

*Synteny* (nucmer, 1-to-1): 10 genomes; the closest PP078674.1 gives 1 block(s) at 98.39 %IDY, the furthest NC_057679.1 61 block(s) at 89.07 %IDY.

*Barcodes*: 29 hypervariable regions at 600/200 bp; 9/10 candidate loci found, the most variable being trnH-psbA (π 0.251894), psbK-psbI (π 0.212766), ndhF (π 0.152529).

*Population genetics*: 80 genes over Coleus: 6, Oldenlandia: 5; 19 of 160 group rows carry no neutrality test (n < 4 or no segregating site) and 2 of 80 Fst rows are blank — every blank names its reason.

**Coleus — run `standin-coleus/20261004-142600-a614528f`.** Reference PX923121.1 (*Coleus amboinicus*), 8 annotated genomes compared (the three records with no CDS are named in the notes).

*Ts/Tv, CDS scope:* **1.5395** pooled = 15235 Ts / 9896 Tv over **77 shared genes** (68743 columns). Per pair — the manuscript's own table:

| Genome | Group | Shared genes | Comparable columns | Ts | Tv | Ts/Tv |
|---|---|---|---|---|---|---|
| NC_069583.1 | Coleus | 77 | 63252 | 81 | 62 | 1.3065 |
| NC_058328.1 | Coleus | 77 | 63147 | 328 | 218 | 1.5046 |
| NC_057679.1 | Coleus | 77 | 63222 | 198 | 124 | 1.5968 |
| NC_057985.1 | Oldenlandia | 77 | 62771 | 3683 | 2393 | 1.5391 |
| PP078674.1 | Oldenlandia | 77 | 62807 | 3691 | 2409 | 1.5322 |
| NC_057983.1 | Oldenlandia | 77 | 62267 | 3617 | 2332 | 1.5510 |
| PX977793.1 | Oldenlandia | 77 | 62475 | 3637 | 2358 | 1.5424 |

*Ka/Ks* (PAML yn00 4.10.10, per shared CDS): **80 genes considered**, 77 carried by every annotated genome, 7 genomes produced pooled values —

| Genome | Group | Scored | Without a value | Pooled Ka | Pooled Ks | Pooled Ka/Ks |
|---|---|---|---|---|---|---|
| NC_069583.1 | Coleus | 78/80 | 2 | 0.000976 | 0.005868 | 0.1663 |
| NC_058328.1 | Coleus | 76/80 | 4 | 0.005002 | 0.020268 | 0.2468 |
| NC_057679.1 | Coleus | 76/79 | 3 | 0.002562 | 0.012393 | 0.2067 |
| NC_057985.1 | Oldenlandia | 57/78 | 21 | 0.033998 | 0.398480 | 0.0853 |
| PP078674.1 | Oldenlandia | 59/79 | 20 | 0.030475 | 0.366815 | 0.0831 |
| NC_057983.1 | Oldenlandia | 54/78 | 24 | 0.033635 | 0.393052 | 0.0856 |
| PX977793.1 | Oldenlandia | 56/78 | 22 | 0.034156 | 0.370289 | 0.0922 |

*Synteny* (nucmer, 1-to-1): 10 genomes; the closest PZ428230.1 gives 1 block(s) at 99.70 %IDY, the furthest PP078674.1 63 block(s) at 89.10 %IDY.

*Barcodes*: 20 hypervariable regions at 600/200 bp; 9/10 candidate loci found, the most variable being trnL-trnF (π 0.260571), atpF-atpH (π 0.250487), trnH-psbA (π 0.243906).

*Population genetics*: 80 genes over Coleus: 6, Oldenlandia: 5; 19 of 160 group rows carry no neutrality test (n < 4 or no segregating site) and 2 of 80 Fst rows are blank — every blank names its reason.

## Known limits of these tables

- The comparison groups are **genera, not populations**: the Fst columns measure divergence
  between taxa (stated in the popgen tables and report).
- A **CDS comparison skips records with no annotation** — three of the eleven panel records —
  and names them in the notes rather than emptying the table.
- Ka/Ks across genera is saturated at synonymous sites; the within-genus rows are the readable
  ones.
- The mito run is a **cross-genus sketch** (D4-05): two mitogenomes from different genera, one of
  them this project's own prototype annotation, so its Ka/Ks rows are a demonstration of the
  machinery, not a result about a species.
