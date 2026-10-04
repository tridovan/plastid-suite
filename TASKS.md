# TASKS.md — active tasks

> **Phase 4a is active; Phase 4 is kept above it as instructed (all five tasks DONE).** Progress:
> `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 4 — Multi-genome comparison  *(ROADMAP.md → Phase 4 · PROJECT.md §13.4)*

**Goal:** the comparison content of the two species manuscripts — Oldenlandia and Coleus — for
**plastid and mitochondrion**. **Input:** the pinned stand-in panel (`.panel/oldenlandia-coleus/`,
11 plastomes · `.panel/oldenlandia-coleus-mito/`, 2 mitogenomes). Public records, **not** the team's
samples: every output from them carries the stand-in label.

**In-phase checks:** Ts/Tv = **2.6917** for *I. domestica* × *I. gatesii* — **CDS only · 79 shared
genes · MAFFT · IR excluded** (re-measured 04/10, §8.27: three differently-spelled loci are shared);
the module reproduces it today (§8.25,
`tests/golden/test_tier_b.py`) · the identity plot on the same pair must return **7 blocks**,
%IDY **89.66–99.11** (mean **97.37**), coverage **99.09 %** — `nucmer --maxmatch -l 20 -c 65 -t 1
→ delta-filter -1 → show-coords -rcl`, reference in the normalized frame (§8.26).
**The board's earlier "21 blocks, %IDY 83.7–98.2 (mean 91.5)" could not be reproduced from the
probe's own inputs** — a full parameter sweep is in §8.26 and the finding in §9.H; the check now
carries the measured numbers, never a tuned one. A mismatch ⇒ stop, re-check inputs/parameters.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D4-01 | The Ka/Ks estimator | **PAML `yn00`** — already pinned in `environment.yml` (`paml=4.10.10`, §7.8) and the standard pairwise estimator (Yang & Nielsen 2000). Rejected: a hand-rolled NG86 (biased, and a wrong implementation would be a silent-wrong nobody could check); KaKs_Calculator2 (a second tool for the same job) | Dev | PROJECT.md §3 |
| D4-02 | The synteny / %identity tool | **nucmer (MUMmer4)** — the pipeline §7 records for the identity plot is `nucmer --maxmatch → delta-filter → show-coords`; another tool makes the phase check unreproducible. §7.8 records no version ⇒ pinned to the build measured here (4.0.1) | Dev | `docs/tools/mummer4.md` + `environment.yml` |
| D4-03 | Population genetics | **Written here, pure Python** — DnaSP is a Windows GUI (unscriptable, unusable in Docker), pixy/libsequence want VCF or read data we do not have; the estimators are published closed-form formulas, testable against analytic vectors, and this keeps the minimum install pure-Python (CLAUDE.md) | Dev | PROJECT.md §3 |
| D4-04 | What "21 blocks" counts, and with which parameters | **Answered by measurement 04/10**: the probe's inputs are on disk and are the Iris pair, but **no** parameter set reproduces 21 blocks (full sweep in §8.26) — the number belongs to a different, more diverged pair. The parameters the module uses are written down (`-l 20 -c 65 -t 1`, `delta-filter -1`, `show-coords -rcl`); the finding is §9.H and the check carries the measured numbers | Dev | PROJECT.md §8.26, §9.H |
| D4-05 | The **mito** half of M4.1–M4.4 | Two mitogenomes from **different genera** (*O. corymbosa*, *C. hadiensis*), one of them **our own prototype annotation** (§8.24). Run it, label it, read it as a cross-genus sketch; a mito comparison worth a manuscript waits for a second record of the same genus | Dev | PROJECT.md §3 |

### Pitfalls for this phase

| Trap | Instead |
|---|---|
| "the Ts/Tv" quoted without its scope | **2.6917 is CDS only · 79 shared genes · IR excluded** (§8.27; the earlier 2.686 was the same scope over 76 names — three shared loci were spelled differently). The whole-alignment scope is a different number (1.8416, §8.25), and every table labels which one it is |
| MAFFT's output feeding a counter | MAFFT 7.526 prints **lowercase** while the counters read `ACGT` — adapters upper-case at the boundary, and `Alignment` refuses lowercase (§9.G) |
| aligning two circular genomes as they arrive | they start at different origins, so a raw whole-genome alignment is mostly gaps — normalize to the IR frame (Phase 1) or align per gene |
| a stand-in number reaching a manuscript table | every comparison output from `.panel/` says `stand-in`; the team's samples change the number, not the label |
| Tajima's D on n = 3, or Fst between species | report **n** and what the statistic means here; something needing more data is *not computable*, never a number |
| a truncated gene in Ka/Ks | the golden ycf1 has no CDS (P3-01) — excluded, and its absence reported, never repaired |
| a tool installed but not declared | `environment.yml` (pinned) + `docs/tools/` **in the same session** as the install (CLAUDE.md invariant) |

### P4-01 · Identity + synteny, mVISTA-style · M4.2 · [x]

- [x] install `mummer4` and declare it (pinned `environment.yml` + `docs/tools/mummer4.md`)
- [x] `external/mummer.py` through `run_tool()`: `nucmer --maxmatch` → `delta-filter` →
      `show-coords`, with the exact parameters written down (D4-04 — §8.26)
- [x] `modules/synteny.py`: one row per alignment block (query/target spans, length, %IDY,
      orientation) + the %identity-by-coordinate series with its window and step stated
- [x] the table carries the block count and the **mean** %IDY the in-phase check quotes

**Done when:** `plastid run -m synteny` on the Iris pair returns **7 blocks** with %IDY
**89.66–99.11** (mean **97.37**) and coverage **99.09 %**, and the manifest names nucmer's version
and argv. *The board first carried "21 blocks · 83.7–98.2 · mean 91.5 — the probe's own result";
that number could not be reproduced from the probe's inputs under any parameter set (§8.26) —
the finding is §9.H and the check now holds the measured numbers.* **Needs:** D4-02, D4-04 (both
decided and recorded, §3).

### P4-02 · Ka/Ks selection pressure · M4.1 · [x]

- [x] PAML `yn00` through `run_tool()`, declared with its measured version
- [x] a **codon-aware** alignment per shared CDS pair (MAFFT on the CDS itself, never a protein
      alignment translated back by hand) — plus a frame check yn00 does not do itself
- [x] `modules/kaks.py`: Ka, Ks and Ka/Ks per gene plus the pooled row, with the estimator named
- [x] genes one genome does not carry are **listed**, not silently dropped — and listing them
      found three loci that were only spelled differently (atpL/psaL/bspB), which re-measured the
      golden Ts/Tv (§8.27)

**Done when:** on the golden pair the module reports per-gene Ka/Ks with the aligner and PAML's
version in the manifest, and the analytic checks hold — identical CDS give Ka = Ks = 0 with the
ratio reported **undefined** (never 0 or 1), and a Ks = 0 column is labelled rather than divided by.
**No published Ka/Ks exists for this pair** (measured 03/10), so the run's values are recorded in
PROJECT.md §8.27 as measured, with n and the method named: **79 shared genes, 78 scored**, pooled
Ka 0.003586 · Ks 0.016550 · Ka/Ks 0.2167; `accD` reported as a frameshift, 25 genes with Ks = 0.
**Needs:** D4-01 (decided, §3).

### P4-03 · Population genetics, DnaSP-style · M4.4 · [x]

- [x] `core/popgen.py`, pure Python, one function per statistic — Hd · θ(S) · θ(π) · Tajima's D ·
      Fu & Li's D* and F* · Fst — each docstring naming its formula and its source
- [x] `modules/popgen.py` over the panel's groups: per-group values, Fst between groups
- [x] every statistic travels with its **n**, and one that needs more sequences than exist is
      reported *not computable* rather than as a number

**Done when:** the estimators pass analytic vectors (a hand-built alignment with a known Hd and θ; a
constant alignment ⇒ D undefined and said so; two identical groups ⇒ Fst = 0) **and** a real run
over the pinned panel writes the tables with n = 3 / n = 2 labelled as thin. **Needs:** D4-03.
**Done 04/10** — vectors pass (Tajima's D and π also cross-checked against scikit-allel to six
decimals; Fu & Li by hand computation), and the panel run (11 records, Coleus 6 · Oldenlandia 5)
writes 80 genes with n = 4 for most, the n = 3 / n = 1 rows blank with their reason, and the
three records with no CDS named in the notes (§8.28).

### P4-04 · Barcode candidates & hypervariable regions · M4.3 · [x]

- [x] `modules/barcodes.py`: rank the panel's windows and regions by variability, **reusing**
      `modules/diversity.py`'s π — no second computation of the same quantity
- [x] report the hypervariable regions as coordinates and span, with the window/step they were
      found at
- [x] name the candidate barcode loci the plastome carries (`matK` · `rbcL` · `trnH-psbA` …) and say
      plainly which ones the run did **not** find

**Done when:** for one window chosen at random the module's π equals `diversity`'s π for the same
window (a cross-check, not a re-implementation), and the run prints its candidate loci with the
reference genome its coordinates are in. **Needs:** nothing — π exists.
**Done 04/10** — the tier-B check runs `-m barcodes -m diversity` on the golden pair and compares
**every** window: 773/773 identical (§8.29). 26 hypervariable regions (threshold exported),
8 of 10 candidate loci found; `trnH-psbA` is not one of them and the row says why (rps19 lies
between the two genes in both golden records) — plus two silent naming failures found on the way
(`trnk-UUU`, `trnl-GUA`) and fixed (§8.29, §9.H).

### P4-05 · The comparison tables for the two manuscripts · M4.6 + the phase's deliverable · [x]

- [x] one run per manuscript species over the **pinned** panel, offline: Ts/Tv (CDS scope), Ka/Ks,
      synteny, barcodes (+ popgen, which the panel's two groups make meaningful)
- [x] the Ts/Tv tables carry the **CDS-scope** numbers and the fixed scope label — never the
      whole-alignment ones (§8.25/§8.27); the scope also gained a **per-genome** table, because a
      manuscript table is per pair, not one pooled number
- [x] a mito run over the two mitogenomes, labelled a cross-genus sketch (D4-05) — reference
      `OX459128.1`, query this project's own prototype annotation (`.panel/mito-prototype/`)
- [x] every table produced from `.panel/` carries the **stand-in** label — `--label` writes it into
      each table's scope line *and* each report's header, from one place (the runner)

**Done when:** `plastid run` over the pinned panel writes the table set for both species with **no
network**, every scope label names its scope, and the numbers are re-derivable from the run
directory alone. **Needs:** the stand-in panel (pinned ✓); the team's samples change what the tables
are *called*, not whether they run.
**Done 04/10 — see `docs/P4-05_comparison_tables.md` for the three run directories, the exact
commands and the measured highlights (§8.30).** Found and fixed on the way: the CDS scope emptied
the whole table when a record had no annotation (3 of the panel's 11) — it now names them; the
mito run found **no** shared gene because one record names its CDS by `/product` sentence and the
comparison modules were passing the plastid vocabulary explicitly (`core.gene_names.aliases_for`
picks by `/organelle` now); and the run cache served a **stale payload** after a behaviour change
whose module version had not moved (H-08, §9.C-21).

## Phase 4a — The Điềm corpus as the test dataset  *(opened 04/10/2026)*

**Goal:** Phases 1–4 were built and tested on data this project picked for itself. The corpus anh
Điềm provided — **34 Iridaceae plastomes, each in two annotations, plus the bare sequences** — is
the first dataset nobody here chose: it is big enough to fail, it carries its own answer key in
three forms (published · PGA · the sequence itself), and one of its records is the group's own.
**This phase makes it the test corpus; it produces no manuscript tables** — Iridaceae is not
Oldenlandia or Coleus.

**What the corpus is — every row measured 04/10/2026, none read off a website:**

| | |
|---|---|
| Files | 68 `.gb` + 68 unannotated `.fasta` · 28 MB · 34 genomes × 2 annotation sets |
| Where | `.data/thanhdiem_data/` — **in the working tree, gitignored**; copied from `D:\business\pts\thanhdiem_data`, 136/136 sha256 identical to the source |
| Set 01 | the published NCBI annotation |
| Set 02 | **PGA** re-annotation of the *same* sequences (34/34 identical) — the two sets cannot differ by sequence, only by annotation |
| The golden genome is in it | `01_…/Iris_domestica_NC_050833.1.gb` is byte-identical to `tests/data/NC_050833.1.gb` (sha256 `a946f23e…`) |
| The group's own record is in it | `PX516262.1` *Eleutherine bulbosa* (Vu,H.-T.) — the tool-paper species; set 01 **133 = 87+38+8**, set 02 **131 = 85+38+8** |
| IR gate | set 01 **0/34 fail** · set 02 **3/34 fail** — *Chasmanthe* 29.24 % · *Iris hermona* 26.36 % · *Trimezia* 27.78 %, all `ir_source: annotation`, while **de-novo on the same sequence finds 100.00 %** every time |
| Gene totals | set 01 **2,924** CDS + 1,292 tRNA + 272 rRNA · set 02 **2,887** CDS (37 fewer) |
| Gene names | **34/34** records spell at least one locus differently between the annotations |
| The differences, categorised (04/10) | **vocabulary** — `clpP1`→`clpP` 23/34 · `pbf1`→`psbN` 21/34 · `ndhL`→`ndhI` 8/34 · `atpL`/`psaL`/`bspB`→`atpI`/`psaI`/`psbB` 7/34 (§8.27's trio) · `lhbA`→`psbZ` 1/34 · bare tRNA names 8/34 · the rRNA dialect 1/34. **A call difference** — `ycf1` is a CDS in the published record and has none in PGA's in **24/34** (71 % of the corpus); `rps12` loses a copy in 12/34. Spelling aside, the two annotations disagree about whether a gene *exists* |
| Read speed (why it lives in WSL) | 300 small files: **7 ms** here vs **700 ms** on `/mnt/d` (~100×); one cold 28 MB parse only 0.38 s vs 0.55 s (~1.4×) — a run is thousands of small writes, so the copy pays for itself there, not in the bulk read |

**In-phase checks:** the vendored fixtures decide Phase 1's gate offline · both panels pin with
every sha256 re-readable · M8.3 carries measured numbers, not extrapolations · Chloe annotates
**34/34** bare sequences and every count difference from the published record is itemized.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D4a-01 | Is this corpus a stand-in, or the team's data? | **A test corpus.** 34 public Iridaceae records, none a manuscript species — but one (`PX516262.1`) *is* the group's own, so the label is **per record**, and no output from here is ever a manuscript table | Dev | PROJECT.md §3 |
| D4a-02 | Which records are vendored into `tests/data/` | **Three pairs ≈ 1.8 MB** — the golden genome (set 01 costs nothing, it is already vendored; set 02 is the golden re-annotated by PGA: **85 CDS where the contract says 86**), *Iris hermona* (gate PASS / gate FAIL on one sequence), the group's own *E. bulbosa*. The smallest subset where every file carries a written-down expected answer | Dev | PROJECT.md §3 |
| D4a-03 | How the 34 `.fasta` get pinned | **A sha256 manifest**, not `panel pin` — measured: a FASTA pins under its file stem with `parse_error: no GenBank records found` and the CLI exits 1 (UNVERIFIED). A FASTA door in `pin_files` is the alternative, and waits for a second caller (YAGNI) | Dev | PROJECT.md §3 |
| D4a-04 | Which annotation is the answer key for the annotator test | **Set 01, the published record.** Set 02 is a third opinion to compare against, never truth — §8.15 measured PGA at 147 genes where the golden has 132 | Dev | PROJECT.md §3 |
| D4a-05 | `ndhL` (published) vs `ndhI` (PGA) | **The same locus** — on the golden genome both annotations place it at `118156-118698(-)` with the identical 543 bp `ATGTTCCCT…TAA`. Fourth of §8.27's class (a locus hidden by spelling); add the alias with that evidence | Dev | `data/gene_aliases.yaml` + §8 |

### Pitfalls for this phase

| Trap | Instead |
|---|---|
| two files, one genome | pin each annotation set into its **own** directory: a pinned file is keyed by the record's own id, and PGA keys a record `Aristea_ecklonii_voucher_KK2022_1_NC_072561` where NCBI keys it `NC_072561.1` — one directory holds one genome twice under two names with nothing pairing them. The pairing is the **shared file stem** |
| a name difference read as a gene difference | 34/34 records differ in spelling (measured) — an alias needs **coordinates and sequence**, or it is a guess |
| trusting an annotation's IR | 3/34 PGA records carry IR copies at 26–29 % while the sequence is 100 %; on *Chasmanthe* de-novo returns PGA's own 26,440 bp span at 100.00 %, so the span is right and the pairing is not — exactly what the gate exists for |
| a `.fasta` through `panel pin` | measured: UNVERIFIED + exit 1; the sequence half needs its own route (D4a-03) |
| extrapolating the time budget from 11 records | the 11-record barcodes run took ~15 min (measured 04/10) — measure at 34, never multiply |
| the run cache serving a stale payload | H-08: a behaviour change without a version bump is invisible; `--force` or bump the version |
| the corpus drifting into a manuscript table | it is Iridaceae — every run from it is a test run and labelled as one |

### P4a-01 · The corpus in the tree, and pinned · M4a.1 · [x]

- [x] copy into the working tree, gitignored — `.data/thanhdiem_data/`, 136 files · 28 MB · every sha256 matching the source (04/10; a `tar` stream, under a second)
- [x] pin set 01 → `.panel/iridaceae-34/` and set 02 → `.panel/iridaceae-34-pga/`, one `--file` per record, `--note` naming the source directory and the date
- [x] pin the 34 `.fasta` as a sha256 manifest (D4a-03), not through `panel pin` — `fasta.json`, 68 entries
- [x] `PINNED.json` is the verification table: sha256 · length · gene counts · topology · IR outcome, measured at pin time for every record

**Done when:** both `PINNED.json` files hold 34 records, every pinned file's sha256 re-reads as
recorded, and set 02's file shows the three gate failures as `ir_identity` < 0.999 with each record
named. **Needs:** D4a-01, D4a-03.
**Done 04/10 — §8.31.** Both panels re-read clean (sha256 · md5 · bytes · the measured values), 34
records each; set 01 **0** gate failures · 2,924 CDS, set 02 **3** — *Chasmanthe* 29.24 % · *Iris
hermona* 26.36 % · *Trimezia* 27.78 %, each named — · 2,887 CDS. The `.fasta` half: every one of the
68 sha256s re-reads, and **68/68 sequences are identical to the pinned record's** — so the manifest
records bytes, while the sequence's provenance is the panel itself.

### P4a-02 · Vendored fixtures + the Phase-1 checks · M4a.2 · [ ]

- [ ] vendor the three pairs (D4a-02) under `tests/data/corpus/` + `SHA256SUMS` + PROVENANCE
- [ ] tier-A checks: `NC_050833.1.gb`'s sha256 is still `a946f23e…` · set 02's golden measures **85 CDS** with `rps12`'s second copy absent · *Iris hermona* set 01 normalizes while set 02 raises `IRGateError`
- [ ] the `ndhI` ↔ `ndhL` alias into `data/gene_aliases.yaml` with the measured evidence (D4a-05)
- [ ] every vendored record through the CLI's gate + normalize, exit code asserted per record (0 / 2)

**Done when:** `make test` covers the corpus offline with no network, the golden fixture's sha256 is
unchanged, and the hermona pair gives PASS / `IRGateError` (exit 2). **Needs:** D4a-02, D4a-05.

### P4a-03 · Scale — what a 34-genome run costs · M4a.3 (M8.3) · [ ]

- [ ] wall clock and **peak RSS** per module over the pinned panel, one module at a time (`/usr/bin/time -v`)
- [ ] the two heavy costs named separately: whole-genome MAFFT, and Ka/Ks's per-gene MAFFT + yn00 (≈80 + 80 calls **per pair**, × 33 pairs)
- [ ] the batching rule for an ordinary machine — the largest N actually measured, and what fits the 8 GB cap
- [ ] the numbers into PROJECT.md §8; M8.3 ticks only once they exist

**Done when:** a module × N × wall clock × peak RSS table exists with the largest N **measured** for
each module, and a stated batch size for a 34-genome run. **Needs:** the pinned panel (P4a-01).
**Pitfall:** H-08 — run every timed module with `--force`, or the clock measures nothing.

### P4a-04 · Chloe on the 34 bare sequences — three annotations of one sequence · M4a.4 · [ ]

- [ ] `-m annotation` (Chloe) over the 34 `.fasta` — the annotator's first real n = 34
- [ ] every output re-read with our own reader · ORF checks · the fixed counting rule (never "133 genes")
- [ ] the three annotations side by side per record — published · PGA · Chloe
- [ ] every difference categorised: **vocabulary** (`clpP1`/`clpP` · `pbf1`/`psbN` · `lhbA`/`psbZ` · `ndhL`/`ndhI` · the rRNA dialects · bare vs anticodon tRNA) versus a **different call** (the `ycf1` CDS, the lost `rps12` copy) versus a **different span** — the first two categories already have counts (see the corpus table)
- [ ] §8.15's PGA verdict re-measured at n = 34 instead of n = 1

**Done when:** 34/34 annotate and re-read, every count difference from the published record is
itemized per record, and the vocabulary-versus-call split is a number. **Needs:** P4a-01 (the panel)
and P4a-02's answer key. **Pitfall:** a name is not a gene — this corpus is where that distinction
*is* the result.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| The team's own Oldenlandia/Coleus samples | anh Điềm / cô Trang | the manuscript tables being *the team's* — P4-05 runs on the stand-in panel either way |
| RNA-seq (one sample, any tissue) | team | **M4.5 RNA editing** — nothing else in this phase |
| A second mitogenome of the **same genus** as a manuscript species | team / public | a mito comparison whose numbers are readable (P4-01…04 are plastid; D4-05) |
| `trnG-GCC` vs `trnG-UCC` — one name for two tRNA sequences (H-04) | cô Trang | the gene-table wording in the manuscript |
| Which gene name the manuscript and the NCBI submission print (D3-06) | cô Trang | wording only; the pipeline prints `canonical (synonym)` today |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| H-04 | `trnG-GCC` vs `trnG-UCC` (PROJECT.md §10.10) | OPEN | cô Trang decides the name printed in the manuscript table |
| H-05 | PMGA web service re-test (re-tested 03/10: still refuses connections) | OPEN | the `curl` in ROADMAP answers; then one web run + `panel pin --file` — nothing here waits on it |
| H-06 | `data/mito_gene_names.yaml` holds **31** loci; the references carry **10 standard genes more** (`rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` · `rps14` · `rps19` · `sdh3`) | OPEN, measured 03/10 | each new locus carries the reference `/product` string that justifies it, and a test proves the vocabulary covers the set — it changes naming for **every** mito record, so it is its own task |
| H-08 | The run cache keys on `module_version` but not on the module's source: a behaviour change without a version bump serves a stale payload (measured 04/10, §9.C-21) | OPEN | either the cache key carries a source hash, or CI fails a changed module whose version did not move |
| H-07 | π per **intron/IGS** — the last piece of M2.9 | OPEN | intron spans are derived inside `modules/introns.py` and IGS spans nowhere: extract that derivation to the shared layer first (DRY's second copy) and feed `diversity`'s region machinery; a run whose regions cover LSC/IR/SSC **plus** every intron and IGS, cross-checked against `introns`' own count (17 on the golden genome) |
