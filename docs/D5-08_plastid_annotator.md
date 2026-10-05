# D5-08 — a self-built plastid annotator (proposal)

**Status: proposal.** This shapes the tool; the open decisions at the end need a call before any
code. When they are settled, the decision row goes into `PROJECT.md` §3 and this file stands as its
design record. It extends `D3-01` ("which annotator") — D3-01 picked an offline tool from what
existed; this builds ours.

## What it is, and where it sits

The **annotate** step of the pipeline the team described:

```
assembly (Phase 7: GetOrganelle · OATK)  →  ANNOTATE (this tool)  →  gene-name unification
(modules/gene_names.py + the alias tables)  →  the analyses (Phases 2/4 modules)
```

- **In:** one assembled plastome as a single-record FASTA (no features — the assembler's output),
  plus a **reference set** of annotated GenBank records pinned by the project's taxonomy rule
  (same family first — the Reference Panel / `mito-reference` pattern).
- **Out:** a GenBank artifact that **re-reads with our own reader** (P3A-01's rule: a file that
  does not re-read is not a result), tables, and a report. The next steps — name unification and
  every analysis module — consume it unchanged; no second code path.
- **Plastid first** (D5-07). The mitochondrial counterpart already exists
  (`core/mito_caller.py` + `modules/mito_annotate.py`); this is its plastid sibling.

## Why build it at all

- **PGA is measurably not the answer**: 147 genes where the golden has 132, plus three silent
  failure modes (§8.15). The reference-transfer idea is sound; that implementation is not.
- **Chloe is the current pipeline annotator** (§8.17: 132 = 86 + 38 + 8, 0 ORF errors, 34/34 on
  the corpus) — it is the bar, not the target. Building ours adds: the arbitration layer below,
  no Julia in the annotation path, and one engine whose every decision is ours to record.
- **The pieces exist**: `external/blast.py` (run_blastx · run_blastn · HSP parsing, BLAST+ 2.17.0
  pinned), `core/orf.py`, the GenBank reader/writer, `run_tool()`, the IR gate, the run manifest.

## The engine — four layers

### 1. Search — orthogonal queries, each labeled with its program

| What | Query | Program | Why this program |
|---|---|---|---|
| CDS | target nucleotide vs **the reference set's own `/translation`** | `blastx` | never self-translate a CDS — §8.21 measured that plant mito CDS cannot be translated at all, and a curator's protein also carries trans-splicing and edited ends correctly |
| rRNA | target vs reference rRNA sequences | `blastn` | the mito caller's proven pair (P3A-03) |
| IR | the sequence against itself | `blastn` (`-perc_identity 99`) | PGA's approach, kept — but behind our structural detection and the ≥ 99.9 % gate (§9.A: naive IR detection fails silently 2/9) |

`blastp` joins only when both sides already carry proteins (naming/verification), and `tblastx`
is not used here: with a curated protein on one side it adds cost without sensitivity.

### 2. Candidates — every call carries its evidence

The six measured traps of §8.22 are rules here, not folklore: `-seg no` (measured: SEG masking
cost `atp9` 137 → 67 bits); **never stitch HSPs of two different sequences behind one key**
(that is how a chimeric gene gets built); the minimum-length subject that contains the hit's span
when keys collide; strand read from the subject side for rRNA-class calls; a version probe that
accepts `-version`; identity **and** coverage floors with the measured separation on the answer
key. Codon/frame checks come from the sequence via `core/orf.py` — the same function that scores
every downloaded record.

### 3. Arbitration — one decision per locus, and it is reported

- **Structure first, score second.** Frame integrity, start/stop codons, intron phase, and
  syntenic neighbours decide before any score does; a score only tie-breaks **within the same
  program** — `blastn` % is nucleotide, `tblastn`/`blastx`/`blastp` % is amino acid, and the two
  are not one scale.
- **Overlap is not error by default.** Measured on the golden (2026-10-05): four overlapping
  pairs are *in the contract* — `psbC`×`psbD` 53 bp, `ndhC`×`ndhK` 10 bp, `atpB`×`atpE` 4 bp,
  `rps12`×`rps12` 114 bp. A rule that drops the weaker overlap would delete real genes. Default:
  **keep both, report the overlap** (the mito caller's measured stance, §8.22-era).
- **Multiple references disagreeing at one locus get one call, not two.** PGA pools every
  reference's sequences and annotates each independently — measured: 147 genes from 7 references,
  four loci called twice under two spellings (§8.15); CGAS's later modules either merge silently
  or show both as different genes (source review, 2026-10-05). Here: one call per locus plus the
  losing candidates **listed with their reasons**.
- **Ambiguity is reported, never guessed** — the rule the naming layer already follows.

### 4. Export — nothing leaves without its reason

GenBank (re-read before it is returned), tables for called / not-called / unnamed hits /
**rejected hits each with a reason** (the §8.22 standard: 68 unused hits, one reason each), the
`source` declares **`/organelle`** (H-11), the run declares its organelle (H-09), and every
external call goes through `run_tool()`.

## Validation — the proposed Done-when

1. **Golden**: annotate the golden *sequence* (features stripped) → exactly **132 = 86 CDS +
   38 tRNA + 8 rRNA**, and **every difference from the published record listed** — the P3A-06
   pattern, kinds and spans, nothing summarized away.
2. **Corpus, leave-one-out**: all 34 bare FASTAs, each with a reference set that **excludes the
   target**; diff each result against the published record **and** Chloe; classify every
   difference as vocabulary or conclusion, the way §8.33 did.
3. **No silent failure**: a search that finds nothing is a failed run; every below-floor hit is
   reported with its reason.
4. **Provenance**: every tool call recorded (name · version · path · argv · exit code); the
   artifact re-reads; `/organelle` present.

**The bar: equal-or-better than Chloe on (1) and (2).** "Better than PGA" is not a claim worth
making — PGA is measurably broken here.

## Open decisions (a call is needed before coding)

| # | Question | Options (recommendation in bold) |
|---|---|---|
| D-a | tRNA calling | **ARAGORN de novo** (already declared, `aragorn=1.2.41`; matches the mito path; no reference needed) vs `blastn` transfer from reference tRNAs (PGA's way) |
| D-b | The plastid reference set | built and pinned by the taxonomy rule (same family first, like `mito-reference`); for validation, leave-one-out on the corpus. Which records — a dev call with the rule applied |
| D-c | Phase slot | **a separate phase (3B, mirroring 3A) after Phase 5** vs inside Phase 6 |
| D-d | Replacement policy | if it passes the Done-when: **replace Chloe as the pipeline annotator** (D3-01 gains an update line) vs stay an alternative |
| D-e | Shared core | **extract the calling core from `core/mito_caller.py`** (the second real copy is arriving — DRY's own rule) with the mito tests as the regression guard |

## Non-goals (deliberately kept out)

No de-novo gene discovery (this is reference **transfer**), no intron model beyond alignment
evidence, no RNA-editing prediction, no repeat discovery beyond the gate's need, no assembly
(Phase 7).

## Evidence this rests on

§8.15 · §8.17 · §8.20–§8.25 · §8.27 · §8.31 · §8.33; the 2026-10-05 source review of PGA
(`PGA.pl:1337–1393` query sets: blastn for tRNA/rRNA and intron-bearing genes, tblastn for CDS
and per-exon; `:1518–1560` the four OS branches; `grep -i overlap` = 0; README: "One genome
should be chosen as the reference") and of CGAS (`cgas_module2.py` passes the whole reference
directory to PGA; `cgas_module3.py` renames to one reference; `cgas_module5.py` keys everything on
post-normalization names). The one measurement taken this session — the golden's four overlap
pairs, 2026-10-05 — is recorded here with its date; it belongs in the §8 log once this becomes a
phase.
