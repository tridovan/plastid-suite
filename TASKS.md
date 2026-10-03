# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 4 — Multi-genome comparison  *(ROADMAP.md → Phase 4 · PROJECT.md §13.4)*

**Goal:** the comparison content of the two species manuscripts — Oldenlandia and Coleus — for
**plastid and mitochondrion**. **Input:** the pinned stand-in panel (`.panel/oldenlandia-coleus/`,
11 plastomes · `.panel/oldenlandia-coleus-mito/`, 2 mitogenomes). Public records, **not** the team's
samples: every output from them carries the stand-in label.

**In-phase checks:** Ts/Tv = **2.686** for *I. domestica* × *I. gatesii* — **CDS only · 76 shared
genes · MAFFT · IR excluded**; the module reproduces it today (§8.25,
`tests/golden/test_tier_b.py`) · the identity plot must show **21 blocks**, %IDY **83.7–98.2**
(mean **91.5**). A mismatch ⇒ stop, re-check the inputs and the parameters — never adjust the
number.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D4-01 | The Ka/Ks estimator | **PAML `yn00`** — already pinned in `environment.yml` (`paml=4.10.10`, §7.8) and the standard pairwise estimator (Yang & Nielsen 2000). Rejected: a hand-rolled NG86 (biased, and a wrong implementation would be a silent-wrong nobody could check); KaKs_Calculator2 (a second tool for the same job) | Dev | PROJECT.md §3 |
| D4-02 | The synteny / %identity tool | **nucmer (MUMmer4)** — the in-phase check's 21 blocks came from `nucmer --maxmatch → delta-filter → show-coords` (§7); another tool makes the check unreproducible. §7.8 records no version ⇒ pin the build installed here | Dev | `docs/tools/mummer4.md` + `environment.yml` |
| D4-03 | Population genetics | **Written here, pure Python** — DnaSP is a Windows GUI (unscriptable, unusable in Docker), pixy/libsequence want VCF or read data we do not have; the estimators are published closed-form formulas, testable against analytic vectors, and this keeps the minimum install pure-Python (CLAUDE.md) | Dev | PROJECT.md §3 |
| D4-04 | What "21 blocks" counts, and with which parameters | The probe's **inputs are on disk** (`/mnt/d/business/pts/probe/run_local/{target,compare}.fa` = the Iris pair), but its exact nucmer command is **not** in the probe tree (searched 03/10). Re-establish the parameters that give 21 blocks and write them down; if they cannot be reproduced, that is a finding for §9 — not a number to tune | Dev | PROJECT.md §8 |
| D4-05 | The **mito** half of M4.1–M4.4 | Two mitogenomes from **different genera** (*O. corymbosa*, *C. hadiensis*), one of them **our own prototype annotation** (§8.24). Run it, label it, read it as a cross-genus sketch; a mito comparison worth a manuscript waits for a second record of the same genus | Dev | PROJECT.md §3 |

### Pitfalls for this phase

| Trap | Instead |
|---|---|
| "the Ts/Tv" quoted without its scope | **2.686 is CDS only · 76 shared genes · IR excluded**; the whole-alignment scope is a different number (1.8416, §8.25), and every table labels which one it is |
| MAFFT's output feeding a counter | MAFFT 7.526 prints **lowercase** while the counters read `ACGT` — adapters upper-case at the boundary, and `Alignment` refuses lowercase (§9.G) |
| aligning two circular genomes as they arrive | they start at different origins, so a raw whole-genome alignment is mostly gaps — normalize to the IR frame (Phase 1) or align per gene |
| a stand-in number reaching a manuscript table | every comparison output from `.panel/` says `stand-in`; the team's samples change the number, not the label |
| Tajima's D on n = 3, or Fst between species | report **n** and what the statistic means here; something needing more data is *not computable*, never a number |
| a truncated gene in Ka/Ks | the golden ycf1 has no CDS (P3-01) — excluded, and its absence reported, never repaired |
| a tool installed but not declared | `environment.yml` (pinned) + `docs/tools/` **in the same session** as the install (CLAUDE.md invariant) |

### P4-01 · Identity + synteny, mVISTA-style · M4.2 · [ ]

- [ ] install `mummer4` and declare it (pinned `environment.yml` + `docs/tools/mummer4.md`)
- [ ] `external/mummer.py` through `run_tool()`: `nucmer --maxmatch` → `delta-filter` →
      `show-coords`, with the exact parameters written down (D4-04)
- [ ] `modules/synteny.py`: one row per alignment block (query/target spans, length, %IDY,
      orientation) + the %identity-by-coordinate series with its window and step stated
- [ ] the table carries the block count and the **mean** %IDY the in-phase check quotes

**Done when:** `plastid run -m synteny` on the Iris pair returns **21 blocks** with %IDY
**83.7–98.2** (mean **91.5**) — the probe's own result — and the manifest names nucmer's version and
argv. **Needs:** D4-02, D4-04.

### P4-02 · Ka/Ks selection pressure · M4.1 · [ ]

- [ ] PAML `yn00` through `run_tool()`, declared with its measured version
- [ ] a **codon-aware** alignment per shared CDS pair (MAFFT on the CDS itself, never a protein
      alignment translated back by hand)
- [ ] `modules/kaks.py`: Ka, Ks and Ka/Ks per gene plus the pooled row, with the estimator named
- [ ] genes one genome does not carry are **listed**, not silently dropped

**Done when:** on the golden pair the module reports per-gene Ka/Ks with the aligner and PAML's
version in the manifest, and the analytic checks hold — identical CDS give Ka = Ks = 0 with the
ratio reported **undefined** (never 0 or 1), and a Ks = 0 column is labelled rather than divided by.
**No published Ka/Ks exists for this pair** (measured 03/10), so the run's values are recorded in
PROJECT.md §8 as measured, with n and the method named. **Needs:** D4-01.

### P4-03 · Population genetics, DnaSP-style · M4.4 · [ ]

- [ ] `core/popgen.py`, pure Python, one function per statistic — Hd · θ(S) · θ(π) · Tajima's D ·
      Fu & Li's D* and F* · Fst — each docstring naming its formula and its source
- [ ] `modules/popgen.py` over the panel's groups: per-group values, Fst between groups
- [ ] every statistic travels with its **n**, and one that needs more sequences than exist is
      reported *not computable* rather than as a number

**Done when:** the estimators pass analytic vectors (a hand-built alignment with a known Hd and θ; a
constant alignment ⇒ D undefined and said so; two identical groups ⇒ Fst = 0) **and** a real run
over the pinned panel writes the tables with n = 3 / n = 2 labelled as thin. **Needs:** D4-03.

### P4-04 · Barcode candidates & hypervariable regions · M4.3 · [ ]

- [ ] `modules/barcodes.py`: rank the panel's windows and regions by variability, **reusing**
      `modules/diversity.py`'s π — no second computation of the same quantity
- [ ] report the hypervariable regions as coordinates and span, with the window/step they were
      found at
- [ ] name the candidate barcode loci the plastome carries (`matK` · `rbcL` · `trnH-psbA` …) and say
      plainly which ones the run did **not** find

**Done when:** for one window chosen at random the module's π equals `diversity`'s π for the same
window (a cross-check, not a re-implementation), and the run prints its candidate loci with the
reference genome its coordinates are in. **Needs:** nothing — π exists.

### P4-05 · The comparison tables for the two manuscripts · M4.6 + the phase's deliverable · [ ]

- [ ] one run per manuscript species over the **pinned** panel, offline: Ts/Tv (CDS scope), Ka/Ks,
      synteny, barcodes
- [ ] the Ts/Tv tables carry the **CDS-scope** numbers and the fixed scope label — never the
      whole-alignment ones (§8.25)
- [ ] a mito run over the two mitogenomes, labelled a cross-genus sketch (D4-05)
- [ ] every table produced from `.panel/` carries the **stand-in** label

**Done when:** `plastid run` over the pinned panel writes the table set for both species with **no
network**, every scope label names its scope, and the numbers are re-derivable from the run
directory alone. **Needs:** the stand-in panel (pinned ✓); the team's samples change what the tables
are *called*, not whether they run.

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
| H-07 | π per **intron/IGS** — the last piece of M2.9 | OPEN | intron spans are derived inside `modules/introns.py` and IGS spans nowhere: extract that derivation to the shared layer first (DRY's second copy) and feed `diversity`'s region machinery; a run whose regions cover LSC/IR/SSC **plus** every intron and IGS, cross-checked against `introns`' own count (17 on the golden genome) |
