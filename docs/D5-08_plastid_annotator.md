# D5-08 — a self-built plastid annotator (proposal)

**Status: the decision is recorded — D5-08 in `PROJECT.md` §3 — and this file is its design
record.** The decision recorded in `PROJECT.md` §3 (D5-08, dev): CGAS's flow with our
search-and-decide core, a reference-set audit in front of every run, **Chloe dropped from the
plan**. The open decisions at the end still need a call before coding. It extends `D3-01` ("which annotator") — D3-01 picked an offline tool
from what existed and now carries an update pointing here; this builds ours.

## What it is, and where it sits

The **annotate** step of the pipeline the team described:

```
gene_normalize on the reference set (D5-09)  →  assembly (Phase 7: GetOrganelle · OATK)  →
ANNOTATE (this tool)  →  gene_normalize, verify mode (D5-09)  →  the analyses (Phases 2/4 modules)
```

- **In:** assembled plastomes as single-record FASTA files (no features — the assembler's output),
  **one or many per run** (batch is part of the CGAS model the team wants), plus a **reference
  set** of annotated GenBank records pinned by the project's taxonomy rule (same family first —
  the Reference Panel / `mito-reference` pattern), and **name-normalized before it is used**
  (D5-09 — the reference set is itself normalized to the naming reference first).
- **Out:** a GenBank artifact per sample that **re-reads with our own reader** (P3A-01's rule: a
  file that does not re-read is not a result), tables, and one report. The next steps — the name
  **verify** pass (D5-09) and every analysis module — consume it unchanged; no second code path.
- **Plastid first** (D5-07). The mitochondrial counterpart already exists
  (`core/mito_caller.py` + `modules/mito_annotate.py`); this is its plastid sibling.
- **Chloe is dropped from the plan** (D5-08, 2026-10-05). When this engine passes the Done-when it
  becomes the annotate step, and `modules/annotation.py`'s engine set (`chloe`, `pga`) changes
  with it — until then that flag set *is* the code.

## Why build it at all

- **CGAS's model is the right shape; its search-and-decide core is measurably broken.** The team's
  call (2026-10-05): the annotate step follows CGAS's flow — a reference set, per-sample
  annotation, a summary, then normalization and comparison (this toolkit already has those as
  `gene_names` + the alias tables + the analysis modules) — but not its core. CGAS module 2 is a
  PGA wrapper, and PGA returned 147 genes where the golden has 132 with three silent failure modes
  (§8.15); CGAS adds no arbitration of its own (source review 2026-10-05: `grep -i overlap` = 0 in
  `PGA.pl`; references pooled into one query set, `:1106-1141`).
- **The current annotators are not wanted.** Chloe is third-party and Julia-dependent — dropped by
  D5-08; PGA is measurably not the answer. The annotate step becomes ours.
- **The pieces exist**: `external/blast.py` (run_blastx · run_blastn · HSP parsing, BLAST+ 2.17.0
  pinned), `core/orf.py`, the GenBank reader/writer, `run_tool()`, the IR gate, the run manifest.

## The model — CGAS's flow, with its measured defects fixed

| CGAS (source review, 2026-10-05) | Here |
|---|---|
| module 2 passes the reference directory straight to PGA; each reference's sequences pool into one query set | our own search layer below; the reference set is pinned, versioned and recorded like every other input |
| no arbitration: a locus called by two references stands twice (measured: 147 genes from 7 references, §8.15) | one call per locus; losing candidates listed with their reasons |
| no overlap logic anywhere in **PGA.pl** (`grep -i overlap` = 0 in its 16,320 lines — commit `0683a0ec…`); CGAS adds no arbitration of its own, and its own `overlap` matches are coordinate-range merging (`cgas_module5.py:545`), not gene-call arbitration | overlap kept and reported by default — the golden itself carries four overlapping pairs |
| "total genes" counts `gene` **features**, and that same number is used as the reference denominator (annotation rate ≡ 100 %) | the counting convention is 132 = 86 + 38 + 8, reconciled by the project's own table — never a self-congratulating rate |
| IR features taken from PGA's `repeat_region`, regions labeled by order and by size; fewer than two IR features → debug-level silent skip | structural IR detection + the ≥ 99.9 % gate before anything position-based |
| module 3 renames genes to one reference, name-based only, original spelling overwritten in the file | arbitration by structure first; name unification stays a **separate** step with its evidence (`gene_names` + the alias tables) — **[07/10/2026, D5-09: that step now runs before annotation, on the reference DB, and writes the normalized GenBank; the reference set must itself be name-normalized before the annotator uses it — `docs/D5-09_reference_naming.md`]** |
| no tool provenance; module outputs to fixed folder names | `run_tool()` records name · version · argv · exit code; every run is run-scoped |

## The reference-set audit — settled 2026-10-05

A gate in front of every annotation run, and available stand-alone (`--audit-only`) so the user can
fix and re-check without annotating anything: every audit failure stops the run — **exit 2, nothing
annotated, the manifest records the error** — and the loop is *fix or drop the bad file → run
again → the audit runs from scratch*.

| Verdict | Checks | Why it is the right verdict |
|---|---|---|
| **BLOCK** — this file cannot be used | unreadable / no GenBank record · no gene features · a CDS missing `/translation` · a gene feature with partial coordinates (`<`, `>`) · IR copies not disjoint/ordered, or the ≥ 99.9 % gate fails · a CDS whose own sequence fails the ORF check | blastx needs the published translation (§8.21 — self-translation is impossible for plant mito, and not trusted for plastid); a partial reference stamps partial flags on the target (§8.15 trap 1); the IR structure is what positions are normalized against |
| **REPORT, do not block** | SSC inverted relative to the majority (flip-flop) · gene-set deviations from the majority — compared as **name sets**, and only with N ≥ 3 (with N = 2 the report says "the two references disagree") | measured 2026-10-05 (§8.34): 32/34 published records share one SSC orientation and **2/34 are inverted — with all 12 SSC CDS translating cleanly**, i.e. flip-flop, not an error. Blocking would refuse two published plastomes, one of them an *Iris*. Gene presence genuinely varies between annotations (`ycf1` 24/34, §8.31) — auto-blocking on it would refuse real references |

**The error message is `file — reason`**, specific to the feature or coordinate, **and nothing
else** — no suggested fix. `NC_024936.1.gb: CDS rps12 — partial coordinates (<153124)`.

**The confirmation layers** for the classes no rule can settle alone: MAFFT `--adjustdirection`
over the region names which references came out reverse-oriented; a dotplot (the Phase-5 figure
system) lets the user see it; the user confirms **once**, and the confirmed exception is pinned
next to the reference set **bound to the file's sha256** — change the file and the exception is
void, and the audit starts over.

**Batch.** One run takes 1..N samples × 1..N reference records; every sample gets its own artifact
and the run gets one report. A failed sample is marked and named **without killing the run**, but
the exit code reflects that something failed. Provenance records every reference file with its
sha256 — the same sample annotated against 1 vs 7 references measurably differs (§8.15).

> **Built 09/10/2026 (P3B-09, §8.49), with two amendments this text did not have.** (1) **A batch is
> a batch of runs**: each sample goes through the same run body every other command uses and gets
> its own directory (`<out>/<accession>/<run-id>/`), its own manifest, report and artifact — so
> "the run gets one report" holds per sample, and "a failed sample is marked" has a durable place
> to be marked. (2) **The gate stops the run, and dropping a file is a declared act**:
> `plastid_annotate._reference_set` audits exactly the records it is about to search, before the
> first BLAST, and a BLOCK raises with `file — reason` and exit 2 — nothing annotated.
> `--exclude <accession>` is the mechanical form of this document's own loop, *"fix or drop the bad
> file → run again"*, recorded in the manifest and named in the run's notes; a name that matches no
> record stops the run. `--audit-only` is the flag that re-checks a set without annotating.
> Measured need for it: the pinned 34-record set blocks on `NC_024936.1`, so without the drop the
> annotator could not run against the project's own reference set at all.

**The report.** English `.docx` (the D2-04 pattern: optional extra, byte-reproducible, and every
number traced to a table cell).

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
   target**; diff each result against the published record (the answer key) and against PGA's
   reading of the same sequence, which the corpus already carries (§8.31); classify every
   difference as vocabulary or conclusion, the way §8.33 did.
3. **No silent failure**: a search that finds nothing is a failed run; every below-floor hit is
   reported with its reason.
4. **Provenance**: every tool call recorded (name · version · path · argv · exit code); the
   artifact re-reads; `/organelle` present.

**The bar is the published record, not another tool**: (1) exact, and (2) every difference named
and classified — with fewer conclusion-level differences than PGA's reading of the same sequences.

## Open decisions (a call is needed before coding)

| # | Question | Options (recommendation in bold) |
|---|---|---|
| D-a | tRNA/rRNA calling | `blastn` transfer from the reference set's tRNA/rRNA sequences (the CGAS/PGA way; what "blastn for non-coding" implies) vs **ARAGORN de novo** for tRNA (already declared, `aragorn=1.2.41`; matches the mito path; no reference needed) + `blastn` for rRNA |
| D-b | The plastid reference set | built and pinned by the taxonomy rule (same family first, like `mito-reference`); for validation, leave-one-out on the corpus. Which records — a dev call with the rule applied |
| D-e | Shared core | **extract the calling core from `core/mito_caller.py`** (the second real copy is arriving — DRY's own rule) with the mito tests as the regression guard |

D-c (phase slot) was settled on 2026-10-07: **a separate phase 3B after Phase 5, opening next**,
with the naming step (**D5-09**) as its **M3B.0** in front of the annotator — the reference set
must be name-normalized before this engine uses it (`ROADMAP.md`).

D-d (replacement policy) was settled on 2026-10-05: **Chloe is dropped**; this engine becomes the
annotate step when it passes the Done-when, and `D3-01` gains an update line then.

## Non-goals (deliberately kept out)

No de-novo gene discovery (this is reference **transfer**), no intron model beyond alignment
evidence, no RNA-editing prediction, no repeat discovery beyond the gate's need, no assembly
(Phase 7).

## Evidence this rests on

§8.15 · §8.17 · §8.20–§8.25 · §8.27 · §8.31 · §8.33 · **§8.34** (the three measurements of
2026-10-05: the golden's four overlapping pairs · 32/2 SSC orientation with clean ORFs · the
corpus reference-quality counts); the 2026-10-05 source review of PGA (commit
`0683a0ec2b2061c49c60ec6c5a2293cd08623253`, `PGA.pl` 16,320 lines, sha256
`0d55b70e53de5d677522cb8fed646898c2ea28bc3f56275a9190128bfd65c3c`; `PGA.pl:1337–1393` query
sets: blastn for tRNA/rRNA and intron-bearing genes, tblastn for CDS and per-exon; `:1518–1560`
the four OS branches; `grep -i overlap` = 0; README: "One genome should be chosen as the
reference") and of CGAS (`cgas_module2.py` passes the whole reference directory to PGA and adds no
arbitration — measured on the golden: 147 genes from 7 references, §8.15; `cgas_module3.py`
renames to one reference; `cgas_module5.py` keys everything on post-normalization names, so a
spelling its synonym table does not know splits one locus into two genes).
