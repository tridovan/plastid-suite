# Oatk 1.0 — assembling a plastome from HiFi reads

**Why this file exists.** `plastid assemble --hifi` runs Oatk, and a module's tool is declared in the
repository. This page is the route's specification. What goes in, what the statuses mean and how
reads and graph are read are the same as for short reads, and [that page](getorganelle.md) has them.

**Why this assembler.** It is built for plant organelles from PacBio HiFi reads, and it is the
assembler of the 195 species in Zhou et al. 2025 (*Genome Biology* 26:235). It takes HiFi reads only:
its documentation names no other technology, so nanopore reads have no assembly route here.

## Install

```sh
conda env create -f environment-assembly.yml
git clone https://github.com/c-zhou/OatkDB && git -C OatkDB checkout 75e8db0
```

Oatk 1.0 and HMMER 3.4 are in the assembly environment. The profiles are not a package. Point the run
at `OatkDB/v20230921` with `--set oatk_db=<dir>` or `OATK_DB`. Both `.fam` files must be there with
their four pressed files (`.h3f`, `.h3i`, `.h3m`, `.h3p`). OatkDB states no checksums, so the run
records each file's digest in the manifest under `reference_sets`, role `assembly_seed`. Measured:

| File | sha256 |
|---|---|
| `embryophyta_pltd.fam` | `4c8010f9562556ab069a4e445bfeb5b61ea4f5729ac555853b022c957c3c87cd` |
| `embryophyta_mito.fam` | `3b63cc59e110e50b18e4d4f0a6e4fff53aa9b23c9582266a9c9dfdeb41e8b35e` |

## The step

```sh
plastid assemble run --hifi HIFI.fastq.gz --sample NAME --out .runs [--reference REF.gb]
```

| Stage | Tool | Output |
|---|---|---|
| Assembly | Oatk, driving HMMER's `nhmmscan` | paths, graph, log |
| Read support | minimap2 (`map-hifi`) on every read given | `junction_support`, `weak_regions` |
| The graph | `core/assembly_graph.py`, no external tool | `assembly_graph`, `graph_segments` |
| Comparison, with `--reference` | `core/assembly.py`, no external tool | `reference_comparison`, `reference_differences` |

There is no quality control stage: a HiFi read is already a consensus. The invocation, with the
defaults from `src/plastid/data/defaults.yaml`:

```
oatk -k 1001 -c 30 -D 1500000000 -t 4 --nhmmscan BIN/nhmmscan
  -m DB/embryophyta_mito.fam -p DB/embryophyta_pltd.fam -o out/oatk HIFI.fastq.gz
```

`-k` and `-c` are Oatk's own defaults, `--set syncmer` and `--set min_coverage`. Oatk's documentation
advises `-c` at 5 to 10 times the nuclear coverage, to keep nuclear sequence out, which a run cannot
know ahead of a sample. For a plastome it did not matter: see below.

**One isomer in, two out.** Oatk writes one contig, whose header states its walk through the graph.
A plastome with an inverted repeat exists in two forms, so the other is spelled here from the same
graph, and `assembly.paths.fasta` holds the pair, as a short-read run does. With more structures than
that there is no pair to complete, and the graph reads `unresolved`. Oatk's overlaps differ link by
link, and its segment for the repeat stops short of the repeat, so the repeat's ends are read from the
bases: on *Arabidopsis* the segment is 25,987 bp and the repeat 26,264 bp, the sequence's own reading.

**How it fails.** Unlike GetOrganelle it exits non-zero, and its reason is the last lines of stderr,
which the step reports. Measured: `empty syncmer graph` when nothing reached `-c`, at `-c 100000`, and
`no organelle component found` on 20 Mb of reads. No sequence is written then, and raising `max_bases`
or lowering `min_coverage` is the way forward.

## What limits the memory

The reads it takes. `-D` stops it reading after that many bases, and `--set max_bases` changes it:

| Bases taken | Peak memory | Time, 4 threads | Plastome |
|---|---|---|---|
| 150 Mb | 1.12 GB | 31 s | closed |
| 300 Mb | 1.12 GB | 34 s | closed, the same circle |
| 700 Mb | 1.24 GB | 44 s | closed, the same circle |
| 1.45 Gb, the whole file | 2.19 GB | 1 min | closed, the same circle |

The default is the most that was measured, so nothing here is an extrapolation: 2.2 GB at any input
size, on this sample's mix. A plastome is reached long before a sample's reads run out, since 150 Mb
of 1.45 Gb was enough. A sample with far fewer plastid reads than this one may need more, and it says
so by stopping. Oatk writes megabytes, so the step asks nothing of the disk. Mapping reads every base
given, so a very large file costs time there. The word limit of [the short-read route](getorganelle.md#what-decides-the-memory)
has no part in this one.

## Measured on an answer from other data

*Arabidopsis thaliana* Col-0, `ERR13987671`, Sequel II, 120,477 reads, 1.45 Gb, N50 14 kb, Q35 by its
own qualities, about 1,450× over the plastome. The answer is `AP000423.1`, assembled from BACs and
PCR with Sanger reads (Sato et al. 1999), so it owes nothing to these reads. On the development
machine (8.33 GB, 5 threads), 2 min 45 s for the whole step, 2.23 GB at the peak:

| | |
|---|---|
| Status, read support, graph | `complete`, `supported`, `resolved` |
| Length | 154,480 bp; `AP000423.1` 154,478 bp; the short-read assembly of `DRR173054` 154,479 bp |
| Against `AP000423.1` | 0 mismatches, 2 single-base differences |
| Against the short-read assembly | 0 mismatches, 1 single-base difference |
| Read support | depth 1,453 ± 105, fewest crossing any boundary 1,294, no weak region |
| Inverted repeat | 26,264 bp in the graph, in the sequence and in the short-read assembly |

The result did not move with the settings: `-c` from 30 to 400 and 150 Mb to 1.45 Gb of reads all gave
one circle, base for base.

**Where it differs from the record.** Each place was read straight from the reads: exact flanks of 21
bases, and the length of the run between them, by technology.

| Place in `AP000423.1` | Record | Short-read assembly | HiFi assembly | Reads say |
|---|---|---|---|---|
| T run after 28,685 | 13 | 14 | 14 | HiFi 14: 999, 13: 67, 15: 276. Short 14: 234, 13: 37 |
| G run after 67,176 | 13 | 13 | 14 | HiFi 14: 732, 13: 575, 15: 132. Short 12: 24, 11: 10, 13: 8 |

The first is the base short reads left open: a second technology, and another study's plants, read 14
where the record has 13. Whether 1999's 13 was an error or a difference between stocks of Col-0 these
reads still cannot say. The second is not settled. HiFi reads divide on a run that long, Illumina
reads of it are poor, and Oatk does not polish, so it keeps whatever the majority of its path holds.
**A one-base difference in a long homopolymer from this route is not a finding** until something
reads that run better.

## Reads longer than the repeat

Reads shorter than a repeat cannot tell the isomers apart, which is why two paths are reported and
neither is called the molecule. Nanopore reads can. The first 1.5 GB of `ERR13987449`, the same study,
holds 53,207 reads, N50 40 kb, 28,940 of them longer than the repeat. Both isomers were spelled from
the graph and mapped together, and a read counted when one alignment covered a whole repeat copy and
500 bp beyond each end. **305 reads fit one isomer and 307 the other.** That is the flip-flop seen in
single molecules, in equal parts, where everything else here infers it from a graph. The HiFi reads
give 2 such reads, so this took length and not depth.

This is a measurement, in `tests/golden/test_assembly_long_panel.py`, and no part of the step. It is
what could choose among the structures of an unresolved graph, with reads long enough for the repeats
of that graph, and no sample with both has been measured.

## Local first, import second

The route is the one on [the short-read page](getorganelle.md#local-first-import-second), with `--hifi`
in place of `-1` and `-2`. `prepare` hashes the reads and prints Oatk's command. Keep its stderr beside
the output as `oatk.log`: that is what states the version, and the import refuses a directory without
it. The directory needs `*.pltd.ctg.fasta` and `*.pltd.gfa`. usegalaxy.eu offers Oatk 1.0+galaxy0 and
usegalaxy.org.au does not (queried 31 October 2026). No run on a service has been measured for this
route, and the import has only read this machine's own output.

Citation: Zhou C, Brown M, Blaxter M, Darwin Tree of Life Project Consortium, McCarthy SA, Durbin R.
*Genome Biology* 26:235 (2025).
