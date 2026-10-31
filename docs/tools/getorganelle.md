# GetOrganelle 1.7.7.1 — the assembly step, and how its result is tested

**Why this file exists.** `plastid assemble` runs GetOrganelle, and a module's tool is declared in
the repository. This page is also the step's specification: what goes in, what comes out, what the
statuses mean and which published numbers it is judged against.

**Why this assembler.** It finished the most plastomes in an independent benchmark of seven tools
(Freudenthal et al. 2020, *Genome Biology* 21:254), and its authors published the raw comparison
this step is judged against: 50 public read sets with their published plastomes, commands, logs
and results, including the three it could not close (Jin et al. 2020, *Genome Biology* 21:241,
Additional file 2; `github.com/Kinggerm/GetOrganelleComparison`).

## Install

```sh
conda env create -f environment-assembly.yml
conda run -n plastid-suite-assembly get_organelle_config.py -a embplant_pt,embplant_mt --config-dir <dir>
```

The environment is about 3 GB. The database is 77 MB and is fetched once; point the run at it with
`--set database=<dir>` or `GETORG_PATH`. Measured here: `embplant_pt` and `embplant_mt` 0.0.1 in
both the seed and the label database. The run checks each FASTA against the sha256 in its own
`VERSION` file and records them in the manifest under `reference_sets`, role `assembly_seed`.

## The step

```sh
plastid assemble run -1 R1.fastq.gz -2 R2.fastq.gz --sample NAME --out .runs [--reference REF.gb]
```

| Stage | Tool | Output |
|---|---|---|
| Quality control | fastp ([page](fastp.md)); `--set trim=no` skips it | `read_qc` |
| Assembly | GetOrganelle, driving Bowtie2, SPAdes and BLAST+ | paths, graph, log |
| Read support | minimap2 ([page](minimap2.md)) on the raw reads, read by `core/assembly.py` | `junction_support`, `weak_regions` |
| The graph | `core/assembly_graph.py`, no external tool | `assembly_graph`, `graph_segments` |
| Comparison, with `--reference` | `core/assembly.py`, no external tool | `reference_comparison`, `reference_differences` |

The assembler's invocation, with the defaults from `src/plastid/data/defaults.yaml`:

```
get_organelle_from_reads.py -1 R1 -2 R2 -F embplant_pt -o out -R 10 -k 21,55,91,127 -t 4
  --random-seed 12345 --config-dir DB --which-spades BIN --which-bowtie2 BIN --which-blast BIN
  --spades-options "-m GB"
```

`-R 10 -k 21,55,91,127` are the settings of the published comparison. With them 39 of the 50 closed;
the authors closed 8 more with settings chosen for each sample, and 3 stayed open, which gives the
paper's 47. A k-mer longer than the reads is dropped by the assembler.

### What the statuses mean

| Status | The assembler's own words | Written | Artifact for the next step |
|---|---|---|---|
| `complete` | `circular genome` | `assembly.fasta`, and `assembly.paths.fasta` when there are several paths | only with read support `supported` and a `resolved` graph |
| `gapped` | `circular genome with gaps` | `assembly.fasta` | no: an unknown base is not a finished sequence |
| `fragmented` | `N scaffold(s)` | `scaffolds.fasta` | no |
| `no_result` | no status line | the log only | no; this is not a measured zero |

The log and the names of the path files state the result independently, and a directory where they
disagree is refused. Exit `0` means an artifact was declared in the manifest; exit `1` means the
step ran and its sequence may not be used downstream. The sequence stays in the run directory for a
person to read.

**Two paths are the normal result.** A plastome with an inverted repeat has two isomers, because the
small single-copy region recombines between the copies, and reads shorter than the repeat cannot
tell them apart. Both are kept. `assembly.fasta` is the first, and when `--reference` is given
`reference_comparison` says which path is that record. More than two is a different state, and is
what [the graph](#the-graph) is read for.

### Read support

The raw reads, not the trimmed ones, are mapped to the sequence followed by its own first 1,000
bases, so a read can cross the origin, and positions are folded back onto the circle. Quality control
serves the assembler. A check that read its output would lose the evidence it removed: on
*Selaginella*, trimming took one spot from five crossing reads to one with the sequence unchanged. The
aligner clips adapter and poor ends by itself, and the depth stays comparable with published mapping,
which is done on raw reads. A read **crosses** a boundary when it has at
least 30 aligned bases on both sides of it with no indel of 10 bp or more between them: an aligner
bridges a stretch of the target that the molecule lacks with a deletion, and a read does not vouch
for what it jumps over.

A boundary is `supported` when at least **2 reads** cross it, and otherwise `undetermined`: absence
of evidence is never reported as a verified join. The four inverted-repeat junctions (JLB, JSB, JSA,
JLA) and the origin are lookups in one profile, so a plastome with no inverted repeat is read by the
same rule. `weak_regions` applies it to every boundary of the circle, so a wrong join nobody thought
to ask about still shows, and the assembly is `supported` only with none.

Reads that **stop** at a boundary are counted and reported and do not decide. On correct sequence
half the reads can stop at a motif the sequencer reads badly, and at a wrong join inside an inverted
repeat they stop on one side only, as they do at those motifs, so stopping tells the two apart in
neither direction. Only crossing is evidence for a join.

The rule was measured on sequences of known state, each with its own reads:

| Sequence | State | Regions under 2 crossing reads | Fewest crossing |
|---|---|---|---|
| *Magnolia officinalis* `KY085916.1` | correct: three assemblers agree | 0 | 62 |
| *Laurus nobilis*, reassembled | correct: the paper's reassembly | 0 | 87 |
| *Selaginella kraussiana* `MH549643.1`, depth 784 | correct: reassembled identically | 0 | 5 |
| *Eleutherococcus senticosus* `KY085901.1`, depth 22 | correct: reassembled identically | 0 | 2 |
| *Arabidopsis thaliana* `AP000423.1`, depth 1,442 | one base from the assembly (below) | 0 | 182 |
| *Laurus nobilis* `KY085912.1` | wrong joins, named by Jin et al. | **2** | 0 |
| *Musa balbisiana* `KT595228.1`, depth 2,216 | the assembler cannot close it either | **8** | 0 |

Requiring 3 reads already flags a correct plastome at depth 22; an absolute 10 flags it 216 times; a
share of the depth flags *Selaginella*, where 5 reads cross one spot at depth 784 on a sequence that
two independent assemblies give identically. So the requirement is 2, and the margin is one read:
more low-depth sets would test it harder than these did. The values are `flank`, `min_spanning`, `min_clip` and
`clip_window` in `defaults.yaml`, and each can be set with `--set`.

**When an assembly and a reference disagree, the reads are asked about both.** If no path equals the
reference, or the assembly stayed open, the reference is mapped and read by the same rule and
`reference_read_support` puts the two side by side. Neither is taken as the truth.

### The graph

Reads test joins. They cannot test the arrangement: two arrangements of the same pieces cross the same
joins and hold the same k-mers, so a read shorter than a repeat fits both. `core/assembly_graph.py`
reads the selected graph, with no external tool and none of the assembler's words:

| Question | A failure |
|---|---|
| Is each path a walk of the graph: every step a link, and the walk spelling the path as a circle? | refused: the path file and the graph are not from one assembly |
| Is the graph closed: one component, no free end, every segment and every link used? | `unresolved` |
| How many structures close through it? Every circular walk that uses each segment as often as the path does is enumerated here. | `unresolved`, unless one, or two that are isomers |

Two structures are isomers when one inversion turns the first into the second and both cross the same
links. The first is what recombination between inverted copies does. The second is what lets the
reads mapped to the first path speak for both.

The call is `resolved`, `unresolved` or `not measured`, in `assembly_summary` beside what the graph
leaves open. A result with no selected graph, or with a path whose title does not state its walk, is
`not measured` and gives no artifact, because the output contract includes the graph.

Measured on the 58 graphs deposited for the published comparison (GetOrganelle 1.6.2):

| The 47 closed runs | |
|---|---|
| Each path a walk that spells its sequence, the graph closed and wholly used, and the walks enumerated here exactly the paths the assembler wrote | all 47 |
| 2 structures, isomers | 42 |
| 1 structure | 2: *Selaginella kraussiana*, whose repeats are direct, and *Cytinus hypocistis*, 19,400 bp, one segment closing on itself |
| **12 or 15 structures: `unresolved`** | 3: *Picea abies* 15, *P. sitchensis* 15, *Jasminum tortuosum* 12 |

**Why those three are not passed.** The order of the pieces is all that separates their structures:
the number of 41-mers each path of *P. abies* lacks from `HF937082.1` is 1,102 on all 15. Read along
the records by exact anchors, the published order is path 3 of 15 on *P. abies* and path 1 on *P.
sitchensis*, so the order the paths are written in means nothing. The assembler says as much: "Please
check the existence of those isomers by using reads mapping (library information) or longer reads."
Reads longer than a repeat, or pairs across it, could separate them, and this step does not try. All
the paths stay in `assembly.paths.fasta`, and none is offered as the answer.

**Depth is reported and decides nothing.** `graph_segments` gives each segment's depth over one copy's
depth. It is a guide, not a rule: on 3 of the 47 a segment used twice reads 1.32 (*Juniperus cedrus*),
1.47 (*Mankyua chejuensis*) and 2.82 (*P. sitchensis*).

**The inverted repeat, read a second way.** The longest stretch a walk holds twice is the repeat
according to the graph, independent of the sequence: its ends are where the copies stop being one
segment. It equals the sequence's reading on all 41 of the 47 where that finds an inverted repeat. On
*Mankyua chejuensis* the graph holds one of 9,756 bp that the sequence's reading does not report,
because its validated window starts at 10 kb. `ir_bp` then says "no inverted repeat found", and a note
names the graph's reading.

## Measured against the published comparison

GetOrganelle 1.7.7.1 + SPAdes 4.3.0 here, against Jin et al. 2020 (version 1.6.2). The sets are in
`.panel/assembly/`, and `tests/golden/test_assembly_panel.py` holds the expectations.

| Read set | Jin Table S1 | Here | Peak memory · time, 5 threads |
|---|---|---|---|
| `SRR5602589` *Magnolia officinalis* | closed, 160,136 bp, identical to `KY085916.1` | closed, 160,136 bp, identical | 2.34 GB · 14 min |
| `SRR5602602` *Laurus nobilis* | closed, 153,040 bp; `KY085912.1` 152,750 bp differs by "1 sites, 429 bp" | closed, 153,040 bp; the same 429 bp base for base, and 139 bp of the record's own tail | 2.06 GB · 12 min |
| `SRR2037123` *Selaginella kraussiana* | closed, 129,971 bp, identical to `MH549643.1`, two direct repeats | closed, 129,971 bp, identical | 1.18 GB · 5 min |
| `SRR2057084` *Musa balbisiana* | not closed, 14 fragments | not closed, 15 scaffolds | 2.45 GB · 6 min |
| `DRR173054` *Arabidopsis thaliana* Col-0, not in the paper | — | closed, 154,479 bp, one base from the Sanger record `AP000423.1` | 3.54 GB · 11 min |

**The depth is a second engine's check.** minimap2 here and Bowtie2 in Table S2, on the same
sequences: *Magnolia* 150.62 ± 29.51 against 151.30 ± 29.22, *Laurus* reassembled 221.07 ± 31.82
against 221.60 ± 31.43, its record 221.35 ± 38.10 against 221.65 ± 37.88.

**When the reads judge a record.** The reassembled *Laurus* is supported at every boundary, and
`KY085912.1` under the same reads is undetermined at exactly the two places the sequences differ. No
read crosses 113,552–113,602, where the record lacks 429 bp from one repeat copy, or the record's
origin, whose last 139 bases are the reverse complement of positions 93,382–93,520: one copy runs
past its true end into the LSC. Those 429 bp are, base for base, the reverse complement of the insert
Jin et al. print.

**One base from a Sanger record.** *Arabidopsis* is the one set whose answer comes from other data
and another technology (Sato et al. 1999). The assembly differs from it by one T in an intergenic
run at 28,673–28,685: the record has 13, the assembly 14, and of the reads that hold the whole run,
308 say 14 and 50 say 13. Whether that is a 1999 error or a difference between two Col-0 stocks
cannot be told from here, and it is recorded as unexplained.

**Why *Musa* stays open.** Three tandem arrays longer than its 100 bp reads. In the wider graph they
are loops of short edges at 3–5× the depth of one copy, and the contigs stop at them:

| Array | On `KT595228.1` | Reads wholly inside it | Pairs: reads minus record |
|---|---|---|---|
| 24 bp unit | 61,229–61,432 | 1,842 | +75 bp |
| 17–23 bp units | 90,604–90,791, and its mirror | 9,303 | +167 bp |
| 48 bp unit | 118,894–118,994, and its mirror | 5,846 | +268 bp |

A read made wholly of an array's 31-mers can only exist if the array is at least a read long, so no
read crosses it, and a graph's loop cannot be walked a measured number of times. The last column uses
no aligner and no assembler. This library's insert is tight (509 ± 14 bp), each mate is placed by its
first 31 bases, and where the reads hold more bases than the record a pair spans fewer. Of 2,627
windows along the record 2,485 are within 10 bp of 509, so the three departures stand out. On the
record the first array is 84 bp, shorter than a read. These are the record's own reads, so these are
places where it does not agree with them.

At over 2,000× the reads still leave eight stretches of the record that fewer than two reads cross.
They are six places, two of them seen through both repeat copies. Two fall where pairs show the record
departing from the reads (−103 bp at 91,250–91,750, +48 bp at 114,750–115,200). Two sit in sequence of
3–4 % GC, where too few pairs can be placed to measure anything. Two (124,289 and 124,575), which one
read crosses and where pairs show nothing, are not explained.

**fastp does not change the answer on the two sets tried.** On *Magnolia* it removed 9.2 % of the
reads as low quality, and the assembly is still the same molecule as `KY085916.1`; the two isomers
swapped order, so the published orientation is path 2 there, and the comparison reads both paths. On
*Selaginella* it removed 6.6 % and the assembly is still identical to `MH549643.1`.

## Local first, import second

`--backend auto` runs here when the tools, the database and the machine allow it. Otherwise it stops
with the reason and the route, and nothing has run or been sent anywhere:

```sh
plastid assemble prepare -1 R1.fastq.gz -2 R2.fastq.gz --sample NAME --out .runs
# run the printed command on a service or a cluster, on the reads it wrote, then:
plastid assemble run --backend import --prepared prepared.json --result DIR \
    --source-name usegalaxy.eu --retrieved-at YYYY-MM-DD --out .runs
```

`prepare` runs fastp here and hashes what it wrote. The import re-hashes those reads, refuses a log
that used more reads than were prepared, and accepts the result only through
`validate_external_result`; the description is kept as `external_result.json`. The evaluation that
follows is the same code as for a local run, on the raw reads, which the import re-hashes too. So the
import needs the raw reads still in place and minimap2 here, and nothing else from the assembly
environment.

**Services that run this tool** (queried through each server's public API on 31 October 2026):

| Server | GetOrganelle |
|---|---|
| usegalaxy.eu | 1.7.7.1+galaxy1 |
| usegalaxy.org.au | 1.7.7.1+galaxy1 |
| usegalaxy.org | not offered |
| usegalaxy.fr | not offered |

Putting reads on a public service is the user's decision; the toolkit never does it.

**Measured on usegalaxy.eu, 31 October 2026.** The server fetched the public *Magnolia* reads from ENA
and ran `-R 10 -k 21,55,91,127` in 15 minutes. Through the import, the result is `supported`, and both
paths are base for base the ones assembled on the development machine, under the same names, although
the server ran SPAdes 3.15.5, Bowtie2 2.4.5 and BLAST+ 2.15.0 where this machine runs 4.3.0, 2.5.5 and
2.5.0. That is one sample: it shows the route works and that this sample does not depend on those
versions, and says nothing about samples with a harder graph. The server's log calls its database
`customized`, though its path names `getorganelledb/0.0.1`, so the seed cannot be called the same one
from here. `assembly_summary` names the programs and the database each run states, and the notes name
any version that differs from the local pins.

**The result directory.** Galaxy names what it hands back after its own history. Give the files these
names in one directory before the import:

| Galaxy output | Name |
|---|---|
| Logfile | `get_org.log.txt` |
| Path sequences, each element | `<element>.path_sequence.fasta` |
| Selected graph | `<element without its last number>.selected_graph.gfa` |
| Simplified graph, optional | `extended_K<k>.assembly_graph.fastg.extend-embplant_pt-embplant_mt.fastg` |
| Contig labels, optional | the same name with `.csv` |

The element names read like `embplant_pt.K127.complete.graph1.1`, and the path files' names state the
result, so they are not free to change. The selected graph is not optional: without it the graph is
`not measured` and there is no artifact. A path paired with another run's graph is refused, which is
what a hand-assembled directory is exposed to.

## What decides the memory

The size of the input does not. Across the 50 published runs (version 1.6.2, one thread), the
logarithm of peak memory correlates at r = −0.70 with the share of plastid reads and at r = 0.26
with the bases given. A 1.14 Gb sample took 60.9 GB and a 1.04 Gb sample took 2.5 GB. 31 of the 50
stayed under 7 GB; the median was 4.4 GB and the worst 65.1 GB. A file's size therefore cannot say
what a run will use, and `--backend auto` does not pretend to: it checks a floor (`min_ram_gb`) and
disk (`disk_factor`), and gives SPAdes a limit, so a run that needs more stops with a message rather
than meeting the kernel's out-of-memory killer. An estimate from the sample's plastid share belongs
to the resource measurement still open.

Citation: Jin J-J, Yu W-B, Yang J-B, Song Y, dePamphilis CW, Yi T-S, Li D-Z. *Genome Biology* 21:241
(2020).
