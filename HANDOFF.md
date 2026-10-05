# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phases 1, 2, 3, 3A, 4 — COMPLETE and merged to `main`.** **Phase 4a (the Điềm corpus as the test
dataset) — COMPLETE, all four tasks done and pushed** (P4a-01 pinned · P4a-02 vendored + Phase-1
checks · P4a-03 the scale numbers · P4a-04 the annotator at n = 34). Three real bugs fell out of it,
all of the same family — a *record* that killed a whole run (§9.I). **Immediate next action:**
Phase 5 (figures & reports) — **its board is built** (`TASKS.md`, five tasks with Done-whens).
Start with **P5-05** (the queryable-number check over generated reports — the phase's own check, and
the one that can be written first) or **P5-01** (the figure system + the first figure). Of M5.1's
eleven figure types, **seven have a source table and four do not** — the four are named in the board
and in ROADMAP, not scheduled.

**Pending, not yet scheduled — the self-built plastid annotator (D5-08).** The annotate step
becomes this project's own tool (CGAS's model, our search-and-decide core; **Chloe is dropped**,
PGA stays only as the measured counter-example), with a reference-set audit gating every annotation
run. Proposed as a new phase **3B after Phase 5**; the four open decisions (D-a tRNA caller ·
D-b reference set · D-c slot · D-e shared core) and the full design are in
`docs/D5-08_plastid_annotator.md`, and the slot question is in `ROADMAP.md` → Open decisions.
Nothing in Phase 5 waits on it.

**The corpus is now a permanent part of the test set.** 34 Iridaceae plastomes × 2 annotations of
the same sequences (published NCBI · PGA) + 34 bare FASTAs, in **`.data/thanhdiem_data/`** (working
tree, **gitignored**, 136 files · 28 MB, every sha256 matching `D:\`); both annotation sets pinned as
panels and the FASTAs as a manifest (§8.31) · **five files vendored** into `tests/data/corpus/` with
twelve tier-A tests (§8.31/P4a-02) · every file under `tests/data` now guarded by a `SHA256SUMS`
completeness test.

**The corpus in one table (§8.31):** the golden genome's own file is *in* it (sha256 `a946f23e…`,
byte-identical to the vendored fixture) · set 02 is the same sequence annotated by PGA — **85 CDS
where the contract says 86**, and **3/34 records fail the IR gate at 26–29 %** while de-novo finds
100.00 % on the very same sequence · **34/34 records differ in gene spelling**, and the differences
split into vocabulary (`clpP1`/`clpP` 23 · `pbf1`/`psbN` 21 · `ndhL`/`ndhI` 8 · the §8.27 trio 7)
and **calls** (`ycf1` has no CDS in PGA's annotation in **24/34**; `rps12` loses a copy in 12) ·
`ndhL` and `ndhI` are the same locus (same coordinates, identical sequence — the golden pair is
unaffected, both NCBI records say `ndhL`).

**What Phase 4 moved** — read `TASKS.md` (the board) and PROJECT.md §8.26–§8.31 first:

* **§6's golden Ts/Tv is now 2.6917** = 323 Ts / 120 Tv over **79** shared genes, 68,448 columns
  (was 2.686 over 76). The scope never changed; the shared *set* did — atpI/psaI/psbB are spelled
  atpL/psaL/bspB in NC_024936.1 and are the same loci (§8.27). CLAUDE.md's golden table carries
  the new numbers and the reason.
* **Synteny's in-phase check is 7 blocks, %IDY 89.66–99.11 (mean 97.37), coverage 99.09 %** on the
  golden pair — the board's earlier "21 blocks / mean 91.5" does **not** reproduce from the
  probe's own inputs under any parameter set; the sweep is §8.26 and the finding §9.H.
* New tools declared and pinned: **MUMmer 4.0.1** (`docs/tools/mummer4.md`) and **PAML 4.10.10**'s
  `yn00` (`docs/tools/paml.md`); MAFFT 7.526 and BLAST+ 2.17.0 were already pinned.

**The deliverable** — three runs over the pinned stand-in panel, recorded in
`docs/P4-05_comparison_tables.md` with their run directories. Every table carries the stand-in
label (`--label` writes it into each table's scope line and each report's header).

**Merged to `main` 04/10** (fast-forward; `feat/phase-4` was stacked on nothing), per the
established pattern for a closed phase. Phase 5 (figures & reports — the plots the two manuscripts
need) is **still the phase after Phase 4a**.

**Open for the PI (wording only, nothing blocks):** `trnG-GCC` vs `trnG-UCC` (§10.10, H-04) and
which gene name the manuscript/NCBI submission prints (D3-06).

**Waiting on data / people:**

| What | Who | Blocks |
|---|---|---|
| One RNA-seq sample (any tissue) | team | M4.5 RNA editing — nothing else |
| The team's own Oldenlandia/Coleus samples | anh Điềm / cô Trang | the manuscript tables being *the team's* — the stand-in tables exist today |
| A second mitogenome of the same genus as a manuscript species | team / public | a mito comparison worth a manuscript (D4-05 ran as a cross-genus sketch) |
| A curated PGA reference set (P3-03) · the team's `.sbt` (P3-05 `.sqn` half) | anh Điềm / cô Trang | nothing in Phase 4 |

## Log

### 2026-10-04 · Claude · dev WSL · session 19 (Phase 4a: the corpus pinned, the fixtures vendored, and two runs that died on real data)

**Did.** P4a-01 and P4a-02 DONE and pushed; P4a-03 (scale) started, and it has already paid for
itself.

* **P4a-01 · the corpus in the tree and pinned.** `.data/thanhdiem_data/` (gitignored, 136 files,
  28 MB, every sha256 matching `D:\`), both annotation sets pinned as panels, the FASTAs as
  `fasta.json`. Recorded in §8.31 with the measurements that came out of pinning.
* **P4a-02 · fixtures + Phase-1 checks.** Five files · 1.3 MB in `tests/data/corpus/` and twelve
  tier-A tests in `tests/golden/test_corpus.py`: the golden re-annotated is one CDS short, *Iris
  hermona* passes the gate as published and fails it at 26.36 % under PGA **on the same sequence**
  (exit 0 / exit 2 through the CLI), the group's own record 87 → 85, the `ndhI`/`ndhL` alias
  resolving both spellings to one key. Plus a new guard: **every file under `tests/data` must be in
  `SHA256SUMS` and match it** — the corpus files were not, and the test found that before commit.
* **Two corrections to my own board, both worth the ink.** (1) D4a-05 said "add the `ndhI`↔`ndhL`
  alias" — it has been in `data/gene_aliases.yaml` since 03/10 (§8.17); today only re-measured it at
  8/34 pairs. (2) The corpus's tool is identified by the **files themselves**
  (`/note="Annotation Method :: PGA-Plastid Genome Annotator"` in 34/34 of set 02, 0/34 of set 01),
  not by the naming style I had inferred it from.

**P4a-03 is where the corpus earned its place — two records killed two whole runs, and neither
failure was a wrong number.**

1. **`structure`, *Chasmanthe floribunda* NC_072562.1.** LSC 81,003 bp ⇒ the 1,000 bp windowing
   ends in a **3 bp** window (`TTA`). No G, no C, `gc_skew` undefined, module raised ⇒ **run dead**,
   1 of 34. The golden genome survives only because its LSC ends in a 127 bp window that happens to
   hold G and C. Fixed: an undefined window is blank (cumulative blank from there on), counted in
   `n_undefined_windows` and named in a note; a *region* without G and C still raises.
2. **`orf`, *Libertia pulchella* NC_072575.1.** The record annotates `rps16` **twice**
   (`5741-5780,4684-4895` and `5737-5780,4684-4891`); both collapse to span `4684-5780`, so both
   took the id `rps16:4684-5780`, and `Table` refuses a repeated key ⇒ **run dead**, 1 of 34. Fixed:
   ids unique per feature (`#2`), both rows checked, `n_duplicated_cds` counted and named — a
   duplicated annotation is a finding about the record, not rubbish to drop.

Both modules → **0.2.0** (H-08: the cache keys on `module_version`). Regression tests for both are
offline and synthetic; PROJECT.md **§9.I** carries the lesson — *neither fixture could have come
from data this project picked for itself, because it picks records that work.*

**P4a-04 · Chloe over the 34 bare sequences — DONE (§8.33).** Running the annotator on a bare FASTA
needed a **FASTA door in `plastid run`** first (`-i` read GenBank only), and opening it exposed a
third bug of the same family: `orf` could not analyse a record with no CDS, though its own comment
says that is a legitimate input — the per-CDS table was built *before* the branch that guards it.

34/34 annotated, none failed. Totals: published 2,924 CDS · PGA 2,887 · **Chloe 2,922**, with
**0 CDS failing the ORF checks** and 201 alternative starts. **No record has all three agreeing on
the gene set (0/34)**, and the differences are not scattered — they cluster on the two hardest loci
of a plastome: **`ycf1`** (a CDS in 24/34 published records that neither other annotator calls — the
IR/SSC boundary gene §8.11 already flagged) and the **trans-spliced `rps12`** (Chloe calls 2 copies
where 22/34 published records carry 1 — the number of IR copies is an annotator's decision, not the
sequence's). The rest is vocabulary (4/34 differ in tRNA spelling only). And **132 = 86 + 38 + 8 is
not a constant**: published and Chloe both range 84–87 CDS here — the P2-10 trap.

**P4a-03 · scale (M8.3) — DONE, `docs/M8.3_scale.md` + §8.32.** Machine stated with every number
(i7-1165G7, 6 vCPU, 8,135,488 kB — the ordinary machine the toolkit targets).

* Nine single-genome modules: **13–26 s for all 34 genomes each**, peak **87 MB**.
* The whole-genome MAFFT (`substitutions`' default scope): 4 genomes 72 s · 8 → 2 min 44 s ·
  **16 → over 7 min** · **34 → over 11 min**. Peak stays under 100 MB throughout.
* The **CDS scope** — the scope the manuscripts quote — does all 34 in **5 min 24 s** at 179 MB, and
  finishes. On a machine with 8 GB, memory is the cheaper currency.
* Per-pair modules: **`synteny` 19.6 s** (99 tool calls) · **`popgen` 6 min 13 s** (80) · **`kaks`
  8 min 26 s** — 5,202 MAFFT + 4,766 `yn00` = **9,968 tool calls**.
* **The finding: call count is not what costs.** 9,968 two-sequence alignments finish in 8.5 min; ONE
  34-plastome alignment does not finish in 11. A module's price cannot be read off its loop — the
  estimate written into the file first said "about 2,600 of each", and the 33 `nucmer` calls it
  expected to be the slow half turned out to be the fastest thing measured.
* `diversity` and `barcodes` at 34 are **named as not measured**, never estimated.

**Late in the session: the run log moved out of `PROJECT.md`** — §8 was **1,531 of 2,379 lines
(118 of 216 KB)**, so reading one entry meant opening a document that no longer even fits a default
Read. It is now five files under `docs/log/` (`phase-1` · `phase-2-3` · `phase-3a` · `phase-4` ·
`phase-4a`), with a 33-line index left in §8. **Entry numbers did not move**, because **271
references** point at them by number (27 in code docstrings, 36 in tests, 23 in `docs/`, 115 across
the four session files, 70 inside PROJECT.md) — renumbering would have broken citations a docstring
cannot re-derive. The split was proven a **pure move**: the five files concatenated are
byte-identical to the §8 that was in git, and §1–§7 / §9–§13 are untouched. New guard:
`tests/unit/test_project_docs.py` fails the build if a cited `§8.N` has no heading, if the index
misses an entry, if one number appears in two files, or if a cited `§9.x` item is gone — with a
negative test, because a guard nobody has watched fail is not a guard. Rule written where it is
read: `CLAUDE.md` → "The run log".

**Validation.** `make check` green after each commit — **564 fast + 42 golden** at the last run
(tier-B skips name their missing tool; this shell has no tool env unless `PATH` carries
`/root/toolbox/env/pts/bin`).

**Residual risks.** (1) The corpus is gitignored and re-fetchable only from anh Điềm's disk; the
durable record is §8.31 plus `PINNED.json`. (2) The two panels are keyed differently by
construction — anything joining them must join on the **file stem**, never the accession. (3) The
scale numbers are this machine's (WSL, 8 GB cap): they bound a *shape*, not a promise for another
machine.

**Blocked.** Nothing.

**Next.** Finish P4a-03 (§8.32 + board), then **P4a-04** — Chloe over the 34 bare sequences, which
needs a FASTA door in `plastid run` (`-i` reads GenBank only, measured) or the de-annotated-GenBank
route P3A-07 used.

### 2026-10-04 · Claude · dev WSL · session 18 (Phase 4a opened — the Điềm corpus becomes the test dataset)

**Did.** Opened **Phase 4a** (board + ROADMAP milestones M4a.0–M4a.4) and closed **P4a-01**.

* **The corpus, in the tree.** 34 Iridaceae plastomes × 2 annotations of the same sequences
  (published NCBI · PGA) + 34 bare FASTAs = 68 `.gb` + 68 `.fasta`, 136 files, 28 MB. Copied from
  `D:\business\pts\thanhdiem_data` into **`.data/thanhdiem_data/`** (working tree, **`.data/`
  gitignored**), 136/136 sha256 verified. Measured why it belongs in WSL: 300 small files
  7 ms vs 700 ms on `/mnt/d` (**~100×**), while one cold 28 MB bulk read is only 0.38 s vs 0.55 s
  (**1.4×**) — a run is thousands of small writes, so that is the number that matters.
* **Both sets pinned, verified.** `.panel/iridaceae-34/` and `.panel/iridaceae-34-pga/`, 34 records
  each, every sha256 · md5 · byte count and measured value re-read (`PINNED.json` is the table).
  Set 01: **0** gate failures, 2,924 CDS. Set 02: **3** — *Chasmanthe* 29.24 % · *Iris hermona*
  26.36 % · *Trimezia* 27.78 %, each named — 2,887 CDS.
* **The FASTA half needed its own route.** Measured: `panel pin --file` on a `.fasta` pins it under
  the file stem with `parse_error: no GenBank records found`, prints UNVERIFIED and **exits 1** — the
  pin path is GenBank-only. Pinned instead as `fasta.json` (68 entries: sha256 · bytes · length ·
  the paired record's sha256), **and every one of the 68 sequences was proven identical to its
  pinned record**, so the manifest records bytes while the sequence's provenance is the panel.
* **Pinning the two sets into one directory would have been wrong** — `panel pin --file` keys a
  record by its own id, and PGA keys one `Aristea_ecklonii_voucher_KK2022_1_NC_072561` where NCBI
  keys it `NC_072561.1`: one genome, two names, nothing pairing them. A directory per set; the
  pairing is the **shared file stem**.

**What the corpus measured (all new, §8.31).**

* **The golden genome is in it**: `01_…/Iris_domestica_NC_050833.1.gb` is byte-identical to the
  vendored fixture (`a946f23e…`) — the corpus contains the contract it is testing.
* **Set 02 = PGA on the same sequence**: the golden reads **85 CDS where the contract says 86**
  (the second `rps12` copy is gone), and the group's own *E. bulbosa* goes 133 → 131 the same way.
* **3/34 PGA records fail the IR gate** at 26–29 % — and de-novo detection on the **very same
  sequence** finds **100.00 %** every time. On *Chasmanthe*, de-novo returns PGA's own 26,440 bp
  span at 100 % ⇒ the span is right and the pairing is not. 8.8 % of the corpus would be
  silently wrong if anyone trusted the annotation's IR.
* **34/34 records spell at least one locus differently** — and the split matters: **vocabulary**
  (`clpP1`/`clpP` 23 · `pbf1`/`psbN` 21 · `ndhL`/`ndhI` 8 · the §8.27 trio 7 · `lhbA`/`psbZ` 1 ·
  bare tRNA 8 · the rRNA dialect 1) versus **a different call** (`ycf1` has no CDS in PGA's
  annotation in **24/34**; `rps12` loses a copy in 12/34). §8.27's lesson at 34× scale.
* **`ndhL` (published) and `ndhI` (PGA) are the same locus** — same coordinates
  (`118156-118698(-)`), identical 543 bp sequence, on the golden genome itself. D4a-05; the golden
  pair is unaffected (both NCBI records say `ndhL`), so **no golden number moves**.
* **A gap that belongs to P4a-04**: `plastid run -i` reads **GenBank only** (`read_fasta` exists
  solely on the `mito annotate` path) — running Chloe over 34 bare sequences needs a FASTA door, or
  the de-annotated-GenBank route P3A-07 already used (the two carry the same sequence, 68/68
  verified).

**Validation.** `make check` green (**533 fast + 30 golden**; tier-B skips name their missing tool —
this shell has no tool env). Both panels re-read clean; `fasta.json` re-reads clean; the corpus
measurements above are all scripts run in this session, not readings.

**Residual risks.** (1) The corpus is **gitignored and re-fetchable only from anh Điềm's disk** —
the durable record is §8.31 plus `PINNED.json`; a re-fetch must re-produce the same sha256s.
(2) The two panels are keyed differently by construction, so anything joining them must join on the
**file stem**, never on the accession. (3) PGA is present only as pinned output, never as a runnable
tool here — unchanged from Phase 3.

**Blocked.** Nothing in Phase 4a.

**Next.** P4a-02 — vendor the three pairs into `tests/data/corpus/` with `SHA256SUMS` + PROVENANCE,
the tier-A gate checks, and the measured `ndhI`↔`ndhL` alias.

### 2026-10-04 · Claude · dev WSL · session 17 (Phase 4 complete — five tasks, and the numbers they moved)

**Did.** The whole Phase 4 board, task by task, each committed and pushed on `feat/phase-4`:

* **P4-01 · synteny.** `external/mummer.py` (nucmer --maxmatch -l 20 -c 65 -t 1 -> delta-filter -1
  -> show-coords -rcl, MUMmer 4.0.1 pinned), `core/delta.py` (the delta format read as measured:
  one header per pair, `|d|-1` paired columns then one side advances), `modules/synteny.py` (block
  table, per-genome summary, 600/200 % identity series, raw coords/delta as artifacts).
  **The board's "21 blocks, %IDY 83.7–98.2 (mean 91.5)" does not reproduce** from the probe's own
  inputs under any parameter set (§8.26 sweep); the check now holds the measured **7 blocks,
  89.66–99.11, mean 97.37, coverage 99.09 %** and the finding is §9.H.
* **P4-02 · Ka/Ks.** PAML `yn00` through `run_tool()` on a per-gene MAFFT codon alignment, with our
  own frame check (yn00 reads codons and never checks, measured) — and the "list the genes one
  genome does not carry" clause found **three loci spelled differently** (atpL/psaL/bspB). Aliases
  added; the shared set of the golden pair went 76 -> 79, so **§6's golden Ts/Tv was re-measured:
  2.6917 = 323 Ts / 120 Tv over 79 genes, 68,448 columns** (CLAUDE.md updated too).
* **P4-03 · population genetics**, pure Python (D4-03): Hd · θ(S) · π · Tajima's D · Fu & Li D*/F*
  · Nei's Gst, each naming its formula. Tajima's D and π match **scikit-allel 1.3.13 to six
  decimals**; Fu & Li are pinned by hand-computed vectors (EggLib's sdist build returns scrambled
  site views here — recorded); Hudson's Fst was implemented first and **rejected with a measured
  reason** (it returns −1 for identical composition at n=2).
* **P4-04 · barcodes**: π reused from the same core function the diversity module calls (a tier-B
  check compares **773/773 windows**), regions by an exported rule, candidate loci from
  `data/barcode_loci.yaml`. `trnH-psbA` is found in the panel's Oldenlandia records and **not** in
  the golden pair — rps19 sits between the two genes there, and the row says so. Two silent naming
  failures fell out (`trnk-UUU`, `trnl-GUA`) — fixed by a syntax rule and one alias.
* **P4-05 · the deliverable**: three runs over the pinned panel (`standin-oldenlandia`,
  `standin-coleus`, `standin-mito`) with `--label` putting the stand-in provenance into every
  table's scope line and every report header; `docs/P4-05_comparison_tables.md` records how to
  reproduce them and what they hold.

**Three more real defects, all found by running the deliverable rather than by a test:**

1. **The CDS scope emptied the whole table when a record had no annotation** (3 of the panel's 11
   records) — it now leaves them out and **names** them; fewer than one annotated genome is a
   raised error, not an empty CSV.
2. **No shared gene between two mitogenomes** — `OX459128.1` carries no `/gene`, so every CDS
   "name" is a `/product` sentence, and the comparison modules were passing the *plastid* alias
   table explicitly, so the mito vocabulary (`data/mito_gene_names.yaml`, present since P3-04) was
   never asked. `core.gene_names.aliases_for` now picks the table by the record's own
   `/organelle`: 0 shared genes -> **19**.
3. **Every tool call overwrote the previous call's logs.** The log index came from
   `glob("*__{name}__*.stdout.log")`, which does not match the names written
   (`001__mafft__stdout.log`) — 616 MAFFT calls kept **one** log, so "the raw output travels with
   the run" was false for every tool but the last, since Phase 1. Fixed, regression test added
   (three calls -> three pairs, first content intact); the mito run went from 1 to **74 MAFFT ·
   24 yn00 logs**.

Also recorded: the run cache keys on `module_version` but not on the module's source, so a
behaviour change without a version bump served a **stale payload** (the CDS-scope fix came back
`cached`); four modules were bumped to 0.2.0 and the deeper fix is **H-08**.

**Validation.** `make check` green — **488 fast + 30 golden**; the full tier-B set with the tool
environment (MAFFT 7.526 · MUMmer 4.0.1 · PAML 4.10.10 · ARAGORN 1.2.41 · BLAST+ 2.17.0 ·
table2asn 1.29.324 · Chloe 0.1.18): **48 passed, 1 skipped** (PGA, not installed). The three
deliverable runs are in `.runs/phase4/` with their table sets, tool logs and manifests.

**Residual risks.** (1) The phase's checks are new golden numbers measured today (§8.26–§8.30) —
as reproducible as the pinned tools and the pinned panel (gitignored; re-fetch per §8.6).
(2) Ka/Ks on **cross-genus** pairs hits frameshifts in ~20/80 genes (MAFFT is not a codon
aligner); every such gene keeps a row with that status. (3) The mito half is a cross-genus sketch
(D4-05): one sequence per group, so every popgen statistic is blank with its reason.

**Blocked.** M4.5 (RNA editing) — needs one RNA-seq sample, unchanged.

**Next.** Merge `feat/phase-4` to `main`, then Phase 5 (figures & reports — the plots the two
manuscripts need) or the team's own data, whichever comes first. Open for the PI (wording only):
`trnG-GCC` vs `trnG-UCC` (H-04) and which gene name the manuscript prints (D3-06).

### 2026-10-03 · Claude · dev WSL · session 16 (P2-07/08 run for real — and the bug that was hiding there)

**Did.** Merged the three-deep branch stack into `main` (fast-forward, pushed). Then closed the
last Phase-2 leftovers, which had been waiting on MAFFT since Phase 2: **P2-08 (Ts/Tv)** and
**P2-07 (π)** now run on real data.

**The run found a real bug, and it is this project's favourite kind.** MAFFT 7.526 prints its
alignment in **lowercase**; the counters classify a column by `_DNA = {"A","C","G","T"}`, so every
column was silently "skipped": the first real Ts/Tv run reported **zero comparable columns and a
blank ratio** — not a wrong number, no number — while every table still rendered. It survived
because the unit tests build `Alignment` objects by hand (uppercase) and the **MAFFT adapter had no
test at all**. Fixed at both ends: the adapter upper-cases on read (with the measurement in a
comment), and `Alignment` now **refuses** a lowercase row — fail loud instead of counting zeros.

**Validation.** `make check` green: **443 fast + 41 golden** (with the tool env; the 4 remaining
skips are tools this shell lacks). New tier-B checks: the CDS-scope Ts/Tv and the whole-alignment
one, plus π.

**The golden number reproduces exactly, through the toolkit's own path** — §6 is 2.686 = 317 Ts /
118 Tv over 76 shared CDS and 66,066 columns, and that is what the run printed, to the base. The
whole-alignment default is a different number (**1.8416**) over 152,338 columns, and every table
labels its scope. π: **0.007181** whole; **LSC 0.009185 · IRa 0.001682 · IRb 0.001681 · SSC
0.014066**; per gene and a 600/200 sliding window.

**Residual risks.** (1) **π per intron/IGS is still missing** (M2.9's last piece) — named as
**H-07**: intron spans are derived inside `modules/introns.py` and IGS spans nowhere, so the
shared layer needs that derivation extracted first (DRY's second real copy). (2) The
whole-alignment Ts/Tv aligns two circular genomes without normalizing their origins first; the
number is honest for the alignment it describes and its label says so, but a per-region view of a
wider scope would be the more useful shape. (3) M2.10 (manuscript draft) still waits on the team's
real samples — unchanged.

**Blocked.** Nothing.

**Next.** **Phase 4 is open** — branch `feat/phase-4`, board built (`TASKS.md`, five tasks with
measured Done-whens). Start with **P4-01 (identity + synteny)**: install and declare `mummer4`, then
reproduce the in-phase check — **21 blocks, %IDY 83.7–98.2 (mean 91.5)** — whose pipeline §7 records
(`nucmer --maxmatch → delta-filter → show-coords`) and whose inputs are still on disk at
`/mnt/d/business/pts/probe/run_local/{target,compare}.fa`. M4.5 (RNA editing) stays blocked on RNA
data; MAFFT 7.526 is installed and pinned, so the alignment path under every one of these tasks is
already proven.

### 2026-10-03 · Claude · dev WSL · session 15 (Phase 3A closed — the CLI, the golden check, the first real sample)

**Did.** P3A-05 (the `plastid mito annotate` CLI), P3A-06 (the itemized golden check) and P3A-07
(the unannotated sample, pinned). `core/genbank.py` gained `read_fasta` — an assembled sequence is
the input door this phase needed; `cli.py` grew a `mito annotate` command that shares
`_execute_run` with `plastid run`, so the mito path inherits the manifest, the cache and the
provenance instead of re-implementing them (seam 3). The boards and ROADMAP now show Phase 3A
complete.

**Validation.** `make check` green (441 fast + 30 golden, tier-B skips named); with `PTS_BLAST_DIR`
all **8** golden checks pass. Measured: the de-annotated `OX459128.1` (102 features stripped) →
`mito_annotation.gb` → re-read **51 = 31 CDS + 17 tRNA + 3 rRNA** → `orf` + `mito_genes` both run
unchanged. The consolidated count check prints `CDS 31/31 · tRNA 17/17 · rRNA 3/3` beside the
published numbers and lists every difference.

**P3A-07 is the honest one.** `AP045469.1` (*Coleus hadiensis*, 313,585 bp, 0 features) →
**45 genes = 23 CDS + 19 tRNA + 3 rRNA**, pinned `3e7c7114ae5d…`, read back offline from the pin.
**8 vocabulary loci are not called**, and a permissive re-search (E-value 10 against 1e-5) shows
four of them — `mttB` · `rps1` · `rps13` · `rps4` — have **no alignment at all**, so the floors are
not what excluded them. Whether those genes are missing from the genome or from this assembly is
not a question a reference-transfer annotator can answer; that is recorded rather than papered
over by loosening a threshold.

**Residual risks.** (1) The Coleus run is a prototype annotation with 13 of 23 CDS failing the ORF
checks (9 `no_start` — alignments that stop short of the reference's N-terminus); it is **not** a
submission and every table says so. (2) The three tRNA differences on the answer key (one missing,
one extra, one anticodon disagreement) are still awaiting a curator's ruling. (3) The floors are
calibrated on one record; the Coleus run is the first evidence of how they behave on a genome with
no answer key, and it says: they behave fine, the genes really are not findable.

**Blocked.** Nothing.

**Next.** Phase 4 — M4.1 Ka/Ks · M4.2 %identity + synteny · M4.3 barcodes · M4.4 population
genetics · M4.6 Ts/Tv (**scope fixed: CDS only · 76 shared genes · MAFFT · IR excluded**). MAFFT
7.526 is installed and pinned, so the Phase-2 numbers that were waiting on it (P2-07/P2-08) can
also be measured now.

### 2026-10-03 · Claude · dev WSL · session 14 (P3A-04 — the tRNA caller, and 17 is not one set)

**Did.** Installed and declared **ARAGORN 1.2.41** (`environment.yml` + `docs/tools/aragorn.md`,
same session as the install), wrote `external/aragorn.py`, wired it into `modules/mito_annotate.py`
(a `trna_calls` table, tRNA columns in the summary, tRNA genes in the artifact), and added the
golden comparison. 12 new unit tests; the parameter that points at the binaries is now one name
(`--set tools=<dir>`) for both tools.

**Validation.** `make check` green (438 fast + 30 golden, 9 tier-B skips with named reasons); with
`PTS_BLAST_DIR` the five golden tests pass, including the new one. Measured on `OX459128.1`:
**ARAGORN calls 17 tRNA, the record carries 17, and they are not the same 17** — `trnA-UGC
13244-13348` is in the record and not called (105 bp, longer than a tRNA), `trnL-AAG 36355-36427`
is called and not in the record, and the tRNA-His at 11741-11819 is **`trnH-AUG` to ARAGORN and
`trnH-GUG` to the record** — an anticodon disagreement, printed rather than resolved. 9 spans
differ by 1–6 bp; the largest is the one tRNA that crosses the origin (`[258274,73]`, written as
the two parts a `join()` needs).

**Two traps that changed the code.** (1) `aragorn -V` is not a version flag: it exits **0** and
prints "No sequence file specified", which `run_tool` was recording as the version of the tool
that calls tRNA genes — the probe now requires a digit in the answer, and the adapter verifies the
release from `-h` itself. Same family as the BLAST+ `-version` bug from the previous session: the
`--version` convention does not exist. (2) A name is not an identity: the record's 17 tRNAs carry
only 14 distinct names, so any comparison that deduplicates by name reports 14.

**Residual risks.** (1) The three tRNA differences above are **findings, not errors** — nobody has
ruled on whether ARAGORN's `trnH-AUG` or the record's `trnH-GUG` is right for that locus; it is
recorded in PROJECT.md §8.23/§9.E and flagged for review. (2) `-i` (intron search) is off, so an
intron-containing tRNA would be reported missing rather than modelled. (3) ARAGORN's
`trnL-AAG` extra call may be a genuine plastid-derived tRNA or a false positive — not yet
adjudicated.

**Blocked.** Nothing.

**Next.** P3A-05 — `plastid mito annotate --input <fasta> --reference <dir> --out <dir>`, with the
output going through `orf` + `mito_genes` unchanged; then P3A-06's itemized golden check
(31 · 17 · 3, differences listed) and P3A-07 (`AP045469.1`, then pin it).

### 2026-10-03 · Claude · dev WSL · session 13 (P3A-03 — the mito gene caller, measured on the answer key)

**Did.** Built the caller that transfers an annotation from the pinned mito reference set:
`core/mito_caller.py` (pure: chaining, exon merging, floors, refinement), `external/blast.py`
(BLAST+ through `run_tool`, with the outfmt-6 parser beside the argv that produces it),
`modules/mito_annotate.py` (wiring + six tables + the GenBank artifact), `defaults.yaml`'s `mito:`
section, the report template, `docs/tools/blast.md`, and `environment.yml` repinned to
`blast=2.17.0`. 47 unit + 4 tier-B tests, written test-first.

**Validation.** `make check` green (424 fast + 30 golden, 8 tier-B skips with named reasons). Real
run on `OX459128.1` with BLAST+ 2.17.0: **31 CDS + 3 rRNA called, nothing not called**; the golden
test pairs every call with the published CDS its exons overlap (bijection over 31) and prints
**18/31 spans byte-identical, 13 differing by ≤ 75 bp, strand 31/31**, ORF profile
**21 ok · 3 no_start · 7 no_stop**.

**Five measurements changed the design** (all in PROJECT.md §8.22 / §9.D, each with a test):
BLAST truncates `sseqid` at the first whitespace; the committed reference FASTA holds 8 duplicate
headers — 6 identical repeat copies and **2 symbols carrying two different sequences**
(Arabidopsis `atp6` 385/349 aa), so an HSP belongs to the smallest reference length containing it
and two sequences under one key are never chained; `-seg no` (SEG halves atp9's alignment,
137 → 67 bits); `-parse_seqids` **fails** on our reference set and upper-cases ids; and ranking
candidates by bitscore **builds chimeras** (`atp6` pulled 381 bp upstream onto a 75.6 % alignment)
⇒ the winner is the highest-identity candidate that passes the floors.

**A bug outside this module:** `run_tool`'s version probe tried `--version`/`-v`/`-V`, but BLAST+
accepts **only** `-version` — so the tool the module depends on was recorded as "unknown" while
every probe "succeeded" in failing. `-version` added, with a regression test.

**Residual risks.** (1) Two published exons of nad1/nad5 (19 and 7 codons) are not recoverable from
any alignment; the caller reports 4 exons and says so rather than inventing a fifth — P3A-06 owns
listing that difference. (2) The caller cannot know a gene is a pseudogene (the answer key marks
4), so its ORF profile reads 21 ok · 3 no_start · 7 no_stop against the record's 21 · 4 · 3 · 3.
(3) The floors are calibrated on one record; a divergent sample is exactly what P3A-07
(`AP045469.1`) will test. (4) BLAST version mismatch is a hard refusal (the Chloe precedent) — if
that proves too strict for anh Điềm's machine it must be decided in PROJECT.md §3, not loosened
quietly.

**Blocked.** Nothing.

**Next.** P3A-04 — ARAGORN 1.2.41 for tRNA (bioconda; `environment.yml` + `docs/tools/aragorn.md`
in the same session), then P3A-05's CLI (`plastid mito annotate`), P3A-06's itemized golden check
(31 · 17 · 3) and P3A-07 (`AP045469.1`, then pin it).

### 2026-10-03 · Claude · dev WSL · session 12 (Phase 3A opened — a mito annotator we build ourselves)

**Did.** Opened **Phase 3A** (`ROADMAP.md`, `TASKS.md`, branch **`feat/phase-3a`** off
`feat/phase-3`) after measuring what a self-built mito annotator would cost. Then **P3A-01** and
**P3A-02**.

* **Why the phase exists.** P3-04's count check cannot close by either PMGA route: the container is
  11 GB unpacked driving MAKER against an 8 GB WSL cap, and PMGA's web service refuses connections.
  CLAUDE.md's rule for a tool that cannot run is *write that part ourselves* — and the estimate is
  **under 10 MB on disk and seconds of runtime**, against PMGA's 14 GB.
* **P3A-01 · GenBank writer — DONE (§8.20).** `to_genbank()` beside the reader. Round trip on the
  golden: **265/265** features, parts in transcript order, sequence and source intact; 324 unit
  tests green, so the `.tbl` path verified against NCBI's own tables is unchanged. Two findings:
  `/pseudo` and `/trans_splicing` are **bare** qualifiers, and BioPython writes `=""` for an empty
  string while reading its own output back happily (§9.B-25 — a round trip through the same library
  is not evidence); and the qualifier builder moved to **`model.gene_qualifiers()`** so the table
  and GenBank writers cannot disagree.
* **P3A-02 · Reference set — DONE (§8.21).** Eight RefSeq mitogenomes pinned by the panel's own
  taxonomy rule (three Lamiaceae for *Coleus*, three Rubiaceae for *Oldenlandia*, plus Arabidopsis
  · Nicotiana · Oryza), 5.2 MB **not** committed; **287 proteins + 25 rRNAs (149 KB)** are, and a
  test rebuilds them from the pinned panel and compares. Coverage **31/31**.
  Two measurements changed the design: **self-translating is impossible** — all 530 reference CDS
  publish a `/translation`, none declares a `/transl_table`, and plant mito CDS begin with `ACG`
  which no table accepts as an initiator, so CDS-mode translation raises; and **the references are
  polluted** — 236 anonymous ORFs annotated as CDS (121 in *Nicotiana* alone), plastid-derived
  `rbcL`/`petG`/`atpB` inside the mitogenome, and a *Psychotria* "5S rRNA" that is 722 bp.

**Validation.** `make check` green after each task (362 unit + 30 golden; tier-B skips name their
missing tool). Every claim above is a measurement on the pinned records, not a reading of a paper.

**Residual risks.** (1) **The vocabulary is narrower than the gene set** — 23 reference genes have
no entry, ten of them standard (`rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` ·
`rps14` · `rps19` · `sdh3`). They are kept in the reference with `-` so a run can say it found one;
named as **H-06** and deliberately *not* fixed here. (2) The writer is **not** a lossless GenBank
editor — `locus_tag`, `note`, `protein_id`, `translation` and `misc_feature` are outside the model
and stated as such in its docstring. (3) The golden round trip has **no partial coordinates** to
exercise; the partial path is covered by a synthetic test only.

**Blocked.** Nothing in Phase 3A. P3-04's count check stays open (H-05) and is now *also* the thing
the phase is trying to make unnecessary.

**Next.** **P3A-03 — the gene caller** (`blastx` → best hit per gene → codon boundaries → `Gene`;
`blastn` for rRNA), then P3A-04 (ARAGORN), P3A-05 (CLI), P3A-06 (the 31·17·3 golden check),
P3A-07 (`AP045469.1`, then pin it). All measured facts and the pitfalls are in `TASKS.md`; the
design decision that matters most is **D3A-02 — boundaries come from the reference alignment, never
from ORF-walking**, because the answer key's own 31 CDS score 21 ok · 4 pseudogene · 3 no_stop ·
3 no_start.

### 2026-10-03 · Claude · dev WSL · session 11 (Phase 3 closed)

- **Did:** closed P3-04 by decision — the mito annotator is PMGA on its web service, pinned
  (`plastid panel pin --file`, §8.19), not a local adapter, because the machine cannot run the
  annotator (11 GB unpacked + MAKER vs an 8 GB WSL cap and 42 GB free on C:). Deleted the 14 GB
  container after reading its real CLI out of it. Board, ROADMAP M3.1, PROJECT.md §3 and
  `docs/tools/pmga.md` all say so, with the trigger that would reopen the adapter.
- **Validation:** `make check` green — 345 fast + 30 golden. Phase 3's own check (§13.3): every
  run's output re-reads with our reader, the ORF checks pass, the gene recount uses the fixed rule
  — all three hold for both the plastome and the mito path.
- **Risks / notes:** (1) the mito **count check** against a real PMGA run has still never been
  made — it needs one web-service run by the team; until then `mito_genes` is verified against a
  downloaded record, not against an annotator's output; (2) the disk: deleting inside WSL returns
  nothing to Windows (`ext4.vhdx` does not shrink; WSL refuses its own sparse mode) — a manual
  `diskpart compact vdisk` from Windows is what gives the 14 GB back.
- **Blocked:** nothing. Open for the PI: D3-06 (which gene name is printed in the manuscript) and
  the `trnG` naming question (§10.10).
- **Next:** Phase 4 — M4.1 Ka/Ks · M4.2 %identity + synteny · M4.3 barcodes · M4.4 population
  genetics · M4.6 Ts/Tv (**scope fixed: CDS only · 76 shared genes · MAFFT · IR excluded**).
  M4.5 (RNA editing) needs RNA data nobody has yet. Note MAFFT is not on this machine — Phase 4's
  golden checks need it installed and pinned (`docs`/§7.8, measured earlier at 7.526).

### 2026-10-03 · Claude · dev WSL · session 10 (PMGA downloaded, and the backend rule completed)

- **Did:** downloaded and unpacked the PMGA container (Figshare 27201798), read its real command
  line out of its own POD, and closed the half of the plan's backend rule that was missing:
  `plastid panel pin --file/--note` (§8.19), tested, so a tool that ran on a web service or
  another machine can be pinned and consumed offline. `docs/tools/pmga.md`, TASKS.md P3-04,
  PROJECT.md §3 + §8.19.
- **Validation:** `make check` green — **345 fast + 30 golden**. `panel pin --file` measured on the
  stand-in mitogenome: 381,931 bytes · sha256 `e31da1d3…` · 258,274 bp · 31 CDS + 17 tRNA +
  3 rRNA · circular (the §8.14 numbers).
- **Risks / notes:** (1) the container is 11 GB unpacked and made the WSL VM balloon toward its
  8 GB cap — the user saw `vmmemWSL` at 4 GB after the extraction; dropping the page cache inside
  the guest does **not** shrink the VM (that takes `wsl --shutdown` or idle time with
  `autoMemoryReclaim=gradual`). (2) PMGA is still unrun: peak RAM and wall clock are unmeasured,
  and that is the number that decides whether it can ever be a requirement rather than an option.
- **Blocked:** nothing. The container is on this machine and can be deleted whenever the disk or
  the VM footprint is unwanted — the web-service route does not need it.
- **Next:** P3-04 via the pin route (pin a PMGA/GeSeq mito annotation, run `mito_genes`), or
  Phase 4's comparison modules.

### 2026-10-03 · Claude · dev WSL · session 9 (P3-05's `.sqn` half — run against NCBI's own converter)

- **Did:** ran `table2asn` for real (a stand-in `.sbt` is vendored at
  `tests/data/standin_template.sbt`), which found that our table was **not submittable**: fixed the
  missing `source` feature, the dropped `/exception`, and the partial-start rule for a CDS whose
  codon 1 is not an initiator. `Gene` gained `/exception`, `Genome` gained the record's source
  qualifiers, the module reads NCBI's `.val` into its tables, and a tier-B test runs the whole
  `.sqn` path on a machine that has the tool.
- **Validation:** `make check` green — **343 fast + 30 golden**, plus the real-tool tier-B tests
  (`PTS_TABLE2ASN`, `PTS_CHLOE`, `PTS_JULIA`). Measured: NCBI's validator on the golden table
  **13 errors → 7**, the remainder being the class a real accepted plastome submission also has
  (control run: *Dioscorea alata*'s published table → 9 errors, 7 of them the same class).
- **Risks / notes:** (1) the partial-start rule is applied by the writer and named in the report —
  it adds a marker, never moves a coordinate, but it *is* a claim, so a genome whose CDS start is
  wrong for another reason will look "partial" rather than "broken" until someone reads the note;
  (2) the 7 remaining errors are inherent to duplicated IR genes in this format and cannot be
  removed from our side; (3) the stand-in template must never be used for a real submission.
- **Blocked:** nothing. The **team's own `.sbt`** is still needed *to submit* (NCBI's form, cô
  Trang's details) — that is now a submission prerequisite, not a code one.
- **Next:** P3-04's mito annotator (PMGA container — record its real entry point), or Phase 4's
  comparison modules, which can key on the canonical names.

### 2026-10-03 · Claude · dev WSL · session 8 (Chloe adapter — P3-03 closes)

- **Did:** `core/gff.py` (GFF3 reader), `external/chloe.py`, `modules/annotation.py` (annotator
  selection: `chloe` default · `pga` for non-angiosperms), `tests/unit/test_gff.py`,
  `tests/unit/test_chloe.py`, two tier-B tests, a real GFF fixture + provenance, `docs/tools/chloe.md`,
  D3-01 updated, PROJECT.md §8.17 · §9.B.20-21 · §10.10 · §3.
- **Validation:** `make check` green — **333 fast + 30 golden** + the two Chloe tier-B tests run for
  real (`PTS_CHLOE`/`PTS_JULIA`). Measured through our reader on Chloe's own file:
  **132 = 86 + 38 + 8**, ORF **80 `ok` + 6 `alternative_start`, 0 failures**; end-to-end CLI 22 s,
  manifest carries `chloe 0.1.18` · `julia version 1.13.0` · the full argv with `--no-transform`.
- **Risks / notes:** (1) the transform guard is a FASTA comparison — a future Chloe that transforms
  without writing that file would slip past it (the release pin and the reference-set pin are the
  second line); (2) `--no-transform` is measured on one genome; (3) `core/gff.py` is ours because
  BioPython has no GFF3 parser — strict, but it has met one producer's output shape.
- **Blocked:** nothing. Open for the PI: the `trnG` name (§10.10). Licence review (M9.3): the
  checkout has no LICENSE file, the Zenodo deposit declares MIT.
- **Next:** P3-04's annotator (PMGA container — its entry point must be recorded, not invented) ·
  P3-05's `.sqn` (needs the team's `.sbt`) · Phase 4 comparison modules, which can now key on the
  canonical names.

### 2026-10-03 · Claude · dev WSL · session 7 (Phase 3: P3-01 · P3-02 · P3-05 · P3-03)

- **Did:** P3-01 (ORF checks), P3-02 (gene-name normalization + the whole printed form),
  P3-05's Python half (`.fsa` + 5-column `.tbl`, artifacts through the Result/export/cache path),
  P3-03 (PGA adapter + annotation module, gated). P3-04's tool decided (PMGA) but not implemented.
- **Validation:** `make check` green — **276 fast + 30 golden + 1 skipped** (the PGA test, with a
  named reason). Real CLI runs produced, for the golden genome: the ORF report (0 failures, 5
  alternative starts), the 11-record panel name report (155 → 115, 0 unresolved, both bare tRNA
  names resolved), and `submission.fsa` + `submission.tbl` (265/265 features round-trip, cache
  restores them byte-identically).
- **Six bugs found by the work, each now a regression test:** (1) a `gene` feature named `trnK`
  spans the intron, so resolving bare names per *feature* compared an intron against a tRNA and
  reported a false "unresolved" — resolution now runs from the tRNA feature only; (2) the CLI
  refused more than `arity` genomes, contradicting API v0.2 and making an 11-record panel
  impossible; (3) `annotated-not-coding` matched by position alone and `ndhF`'s CDS "covered" the
  truncated ycf1; (4) the TBL writer put `>` after the number (`1050>`) where NCBI writes it as a
  prefix (`>1050`); (5) a table key not present in its own columns (caught by the template fixture);
  (6) `name_mappings` keyed on a column that did not exist.
- **Residual risks:** PGA has never run here, so `external/pga.py`'s argv is from PGA's README and
  unverified; the `.tbl`'s partial/qualifier layout is per NCBI's documentation and their verbatim
  examples, not against a real submitted `.tbl`; the `.sqn` path has never run (no `table2asn`, no
  `.sbt`); template prose cannot contain digits at all, so gene names with digits must come from
  table cells.
- **Blocked:** P3-04 on the PMGA container; P3-03's measured half on a machine with PGA; the `.sqn`
  half on a submitter template (`.sbt`) — `table2asn` itself is installed and declared.
- **Late in the session, three more things landed:** the **mito path** (`modules/mito_genes.py`:
  a plant mitogenome described with no IR gate, 51 = 31 + 17 + 3 on the stand-in, every CDS named
  by a symbol); **`table2asn` 1.29.324 installed and declared** (`docs/tools/table2asn.md`, md5
  recorded) with both of its error paths tested; and **three real NCBI-produced `.tbl` files
  vendored** (`tests/data/real_*.tbl`), which closed the "field for field" clause for everything
  but joins — and which corrected the writer (three fields per interval line, not five) and the
  parser (a table file may hold several records; one VADR file carried eight and they had been
  merged silently).
- **Next:** see the Resume point.

### 2026-10-03 · Claude · dev WSL · session 6 (blocker clearance)

- **Did:** cleared what could be cleared before Phase 3 starts. Installed **MAFFT 7.526** from
  conda-forge and verified it byte-for-byte against §7.8's pin, plus BLAST+ 2.17.0. Found that two
  named blockers were not blockers, and that one tool in the plan is the wrong kingdom. Wrote
  `docs/PHASE3_BLOCKERS.md`.
- **Validation:** the MAFFT tarball is 2,610,460 B with md5 `cf3b2acf649bc31e576ed43a2fa4aaa3` —
  both numbers are §7.8's, so this is the pinned build and not merely a 7.526. `environment.yml`
  dry-runs clean and the solver picks `mafft 7.526 h4bc722e_0` from conda-forge rather than
  bioconda's 7.525.
- **A mistake I made and fixed:** I installed tools on this machine and changed nothing in the repo.
  Every number measured here would have been reproducible by nobody. Worse, the repo's own
  `environment.yml` listed a **bare `mafft`** while §9.B-13 had already measured that the aligner
  changes Ts/Tv (2.667 vs 2.686) — so following our own install file produced a different golden
  number. Both fixed; the rule is now an invariant in CLAUDE.md.
- **Residual risks:** `table2asn` is only partly downloaded (NCBI's FTP timed out twice). PGA is
  cloned nowhere useful — running freshly-cloned third-party code was refused by the permission
  classifier, correctly, since no one had named the source in-session; one manual `perl PGA.pl`
  clears it. PMGA's benchmark numbers come from a web summary, **not the paper**.
- **Blocked:** `docs/PHASE3_BLOCKERS.md` — four asks, each with a named owner.
- **Next:** Phase 3 A→Z. Start with P3-01 and P3-02 (both pure Python, nothing blocks them).

### 2026-10-03 · Claude · dev WSL · session 5 (Phase 2 closed, Phase 3 board)

- **Did:** pushed `feat/phase-2`; opened branch `feat/phase-3` (stacked — Phase 2 is not merged to
  `main` yet); cleared `TASKS.md` of Phase 2 and built the Phase 3 board from the milestones
  (M3.1–M3.4) with the `make-tasks` skill: 5 tasks, each with a measurable Done-when.
- **Validation:** the board's numbers are measured, not guessed — the whole dataset was swept for
  start codons (`/transl_table=11` and `/codon_start=1` on all 1,099 CDS; starts ATG 746 · GTG 11 ·
  ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1) and for gene names the grouping file does not claim
  (exactly four, all in the two manuscript-species records). `TASKS.md` is 108 lines.
- **Corrected a wrong phase check:** ROADMAP's Phase 3 check read "gene recount uses the fixed
  counting rule (**132**)". 132 is *Iris domestica*'s value, not a constant — P2-10 measured 131 and
  133 for two other species. The *rule* is fixed; the *number* is not. Left uncorrected this is the
  "133 genes" trap in a new costume.
- **Did (cont.):** **D3-02 decided** — the four unusual gene names are the same loci as
  `ycf3` · `ycf4` · `psbN` · `clpP`, proven three ways (syntenic neighbours, identical length — `pbf1`
  and `psbN` both exactly 132 bp — and `/product`). PROJECT.md §8.9.
- **Did (cont.):** the user identified `Vu,H.-T.` as **cô Trang**, which made §8.9's conclusion
  incomplete and flipped its main recommendation. Measured afterwards: **two of the eleven records are
  the group's own** (`PX923121.1` *C. amboinicus* — the manuscript species — and `PX977793.1`); the
  naming difference is a **whole dialect** (tRNA hyphen/bracket/**bare**, rRNA `rrn16`/`rrn16S`, protein
  classic/functional), not four genes; and a bare tRNA name resolves by **sequence identity**, not by
  reading an anticodon at a fixed offset (offset 34 matched only 106/293). The comparison key moved to
  the **team's own names**, because normalising to the classic ones would make the team's paper print a
  name that disagrees with the team's own submitted record. PROJECT.md §8.10.
- **Residual risks:** P3-03/P3-04 (PGA, MITOS2) have no tool on this machine, so they can only be
  implemented and `PTS_REQUIRE_TOOLS=1`-gated — the same honest position P2-07/P2-08 are in.
- **Process note worth keeping:** I first routed D3-02 to the PI *because I had no evidence*, then
  kept it there after measuring — which is the wrong half to keep. The mapping is a dev call with
  evidence; only the name **printed in the manuscript / submitted to NCBI** (D3-06) is the PI's.
  A dev-decidable item parked on the PI blocks silently until someone asks.
- **Blocked:** PGA/MITOS2 (anh Điềm's toolbox) block the measured halves of P3-03/P3-04. A
  confirmation of whether `Vu,H.-T.` — submitter of the two records that use the new names — is this
  group would settle whether those names are ours to explain.
- **Next:** P3-01.

### 2026-10-03 · Claude · dev WSL · session 4 (P2-09, D2-04, P2-10)

- **Did:** finished Phase 2. **P2-09** — a categorised publication gene table
  (`data/gene_groups.yaml` is the grouping, `modules/gene_table.py` builds it) plus a Word export
  (`result/docx_export.py`, `plastid run … --docx`) that is **byte-reproducible** and pruned of the
  template's dead weight. **D2-04** — `python-docx` as an optional extra, imported inside the
  function, so the minimum install stays biopython + PyYAML. **P2-10** — `-m` is repeatable, so a
  run holds several modules and writes one combined `report/report.docx`; draft produced for both
  manuscript species from the stand-in set.
- **Validation:** `make check` green — **208 fast + 28 golden**. Every document was rendered
  (LibreOffice → PDF → PNG) and **looked at**, not just asserted on. Measured: two builds of one run
  are byte-identical, and so is a fully cached re-run (including the merged document).
- **Four real bugs found by this work, each now a regression test:**
  1. the gene table's own reconciliation reported **0** on its Total row — `by_kind.get("Total")`
     on a Counter holding only CDS/tRNA/rRNA returns the default. The column caught it on the
     first run, which is the whole reason it exists (PROJECT.md §9.B-19).
  2. the OOXML style pruner read styles as **attributes** (`element.get("w:pStyle")`) when they are
     **child elements** (`<w:pStyle w:val="…"/>`), so it removed nearly every style the document
     used — including `Title` and `Heading 1`. **LibreOffice still rendered something plausible**;
     only a structural test ("every style the document references must still be defined") caught it
     (PROJECT.md §9.C-20).
  3. a table's own NUL-looking allocation bug: scaling every column by one factor starved the short
     ones, so a header broke mid-word ("Gene / s"); columns now get a floor before the rest is
     shared out.
  4. `Table.emphasis` (a bold totals row) accepted any string — a typo would have been invisible in
     every export except the Word document's font weight. It now must name an existing key.
- **Residual risks:** the Word layout is verified by eye against LibreOffice; Word itself is not
  installed here, so a Word-only rendering difference would not have been seen. The stand-in drafts
  are labelled but are **public records of the same two species** — same species is not the same
  specimen.
- **Blocked:** P2-07/P2-08 numbers still need MAFFT 7.526 (absent).
- **Next:** Phase 3, or collect the three open sub-items in the resume point.

### 2026-10-03 · Claude · dev WSL · session 3 (P2-03 … P2-08, P2-11, D2-01/02/03)

- **Did:** completed seven Phase-2 tasks. `core/codon.py` + `modules/codon_usage.py` (RSCU · ENC ·
  CAI · CBI · Fop), `core/protein.py` + `modules/amino_acid.py` (GRAVY · aromaticity),
  `modules/introns.py`, SSR by region + motif, `modules/diversity.py` (π) and
  `modules/substitutions.py` (Ts/Tv), plus group labels on comparison output.
  D2-01 (multi-genome contract, `docs/API_v0.2.md`), D2-02 (π 600/200), D2-03 (NCBI table 1).
- **Validation:** `make check` green — **169 fast + 22 golden**. Every new number is either a
  §6 golden number or an analytic check (ENC = 61 unbiased / 20 at the extreme; CAI = Fop = CBI = 1
  for a self-referenced sequence; π of a known alignment; the two IR copies' skew mirrors exactly).
- **Six real bugs found by the new work, each now a regression test:**
  1. `_parse_gff` read only `monomeric_repeat` — MISA writes every longer unit as `microsatellite`
     and a compound as one `repeat_region` line, so **7 of 37 SSRs were silently dropped**. The old
     test fixture invented a type MISA never writes, which is why it passed. P2-06's reconciliation
     column caught it.
  2. `coding_sequence` — re-sorting a trans-spliced gene's parts or reading them all forward gives
     a wrong ORF; only per-part strand + annotation order yields rps12's `ATG…TAA`.
  3. Two independent copies of a `Table` key were allowed (`rpl2.1` twice) — references resolved
     silently to one row. `Table` now rejects a repeated key.
  4. `load_groups` let `yaml.YAMLError` through (not a `ValueError`) — a malformed file crashed
     with a traceback instead of a message.
  5. `NcbiClient.resolve` could leak a `KeyError` where `NcbiError` belongs.
  6. `entry_from_summary` claimed `annotated=True` for every record; GenBank submissions are not
     annotated by default — 3 of the 11 stand-in plastomes have zero gene features.
- **Residual risks:** the comparison modules are unit-verified but **no real-data number was
  measured** — MAFFT is absent and installing a different version than the pinned 7.526 would
  measure the wrong thing. π per intron/IGS is still missing from P2-08.
- **Blocked:** nothing. P2-10 still needs the team's real samples.
- **Next:** P2-09 (see the binary-export note in the resume point), then P2-10.

### 2026-10-03 · Claude · dev WSL · session 2 (stand-in data, P2-01, P2-02)

- **Did:** **D-04** — a public stand-in dataset for Oldenlandia/Coleus: `data/*.yaml` definitions
  (13 records), `panel pin --groups` (new) + `--index` now optional + live accession resolution via
  the new `NcbiClient.resolve`; `PinnedGenome` measures gene counts · topology · IR outcome at pin
  time, so `PINNED.json` *is* the verification table. Hardened what the first cut got wrong
  (all found by a review pass): malformed YAML raised a traceback (`yaml.YAMLError` is not a
  `ValueError`), `resolve` could leak a `KeyError`, `annotated=True` was assumed for GenBank
  records. **P2-01** done (e2e tests: golden numbers in report *and* CSV · 2nd golden genome ·
  both diverged-IR genomes exit 2). **P2-02** done — `modules/structure.py` + `gc_skew()` +
  `structure.en.md` + `tests/unit/test_structure.py`.
- **Validation:** `make check` green — **124 fast + 18 golden**. All 13 sha256 re-verified
  byte-identical on a second pin. Structure on the golden genome reproduces LSC 83,127 ·
  IRa = IRb = 26,214 · SSC 18,169 · tiling 153,724, and the two IR copies skew exactly opposite
  (reverse complements) — a real correctness signal for the sign convention.
- **Residual risks:** `PZ428230.1` is a bad record (15,756 bp IR, `linear`) — excluded from
  position-based analysis, and nothing in the tool would have caught it automatically (§9.A-8).
  The API doc said `table[rowkey].column` but the implementation parses `table.rowkey.column`;
  no template had used a row key until now — docs corrected to the tested syntax, code unchanged.
  3 of 11 stand-in plastomes have no annotation (0 CDS) — sequence-level analyses only.
- **Blocked:** nothing. P2-10 still needs the team's real samples.
- **Next:** P2-03 (decide D2-03 first), then P2-04 … P2-11.

### 2026-10-03 · Claude · dev WSL

- **Did:** finished Phase 1. P1-10 complete: `core/panel.py` (4 tiers, modes A/B, quality filter,
  `pin()` with an injected fetcher that re-reads the file with our own reader) · `run/ncbi.py`
  (esearch/esummary/taxonomy efetch, 3 req/s, injectable transport) · `plastid panel refresh |
  propose | pin`. Closed P1-12 as a **dev decision** (local for phases 1–5/8, remote for 6–7) and
  added the "Who decides what" rule to CLAUDE.md so engineering calls stop being queued for the PI.
- **Validation:** `make check` (94 fast + 13 golden). Live NCBI run: 76 RefSeq genomes for
  Iridaceae; a re-pinned NC_050833.1 whose **sha256 matches the vendored golden file byte-for-byte**
  (`a946f23e…`), NC_024936.1 → `fd441727…`. Measured and recorded in PROJECT.md §9.C-20: NCBI's
  `chloroplast[filter]` misses "plastid"-worded records (71 vs 76) — including our own comparison
  genome.
- **Residual risks:** the panel's tier model relies on NCBI lineage strings (genus = last component);
  a record with an unusual lineage would land in the wrong tier — the user confirms the final list,
  and `propose` prints the tier for every candidate. `pin()` records unparseable downloads as
  `parse_error` rather than failing the batch — a caller must check it (the CLI prints UNVERIFIED).
- **Blocked:** nothing. Waiting on people only: real reads + nuclear genome (Phase 6/7), item (c).
- **Next:** Phase 2, starting with M2.1 (gene content & summary statistics).

### 2026-10-02 · Claude · dev WSL

- **Did:** created the session set (CLAUDE.md · ROADMAP.md · TASKS.md · HANDOFF.md); squashed git
  history into one `init` commit and force-pushed (amended to remove Claude attribution — repo rule
  + memory); built **P1-01** scaffold, **P1-02** core model + BioPython GenBank reader + vendored
  golden data (sha256-pinned), **P1-03** IR detection + mandatory gate + coordinate normalization.
  §10.6 decided and recorded (GenBank coords per genome · normalized for cross-genome figures).
- **Validation:** `make check` (40 fast + 7 golden, no skips); the 20-genome sweep reproduces
  PROJECT.md §8.3 — 7/9 vs GenBank, gate catches 2/2 diverged genomes with 0 false alarms; golden
  tier-A numbers asserted in tests. Recorded in PROJECT.md §8.5.
- **Did (cont.):** §10.7 decided (reports resolve table cells; an unresolvable reference raises) →
  PROJECT.md §3/§10. **P1-04 green** — `result/result.py` (`Scope` · `Column` · `Table` · `Result`,
  refs `table.column` / `table[row].column`) and `result/export.py` (csv with provenance comments,
  json with raw values, markdown); exports byte-identical across runs.
- **Residual risks:** de-novo detection is still silently short on 1/9 known cases (NC_063920 →
  22,998 bp vs the true 26,161); the annotation cross-check *records* the 3,163 bp disagreement but
  the run does not fail on it — decide the policy later. The de-novo detector is unvalidated on
  mitochondria. `csv` exports carry provenance as `#` comments — readers must strip them (stdlib
  `csv` cannot).
- **Blocked:** §10.8/§10.9 (P1-10) · D-01. Nothing blocks P1-05.
- **Next:** P1-05 test-first — run layout, manifest, `run_tool()` (the single subprocess door), cache.
