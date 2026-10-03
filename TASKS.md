# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 3 — Annotation  *(PROJECT.md §13.3)*

**Phase check (§13.3):** the run's output re-reads with **our own GenBank reader** (no format
damage) · ORF checks pass · gene recount uses the fixed counting **rule**
(**CDS + tRNA + rRNA, both IR copies** — *132 is the golden genome's value, not a constant*: P2-10
measured **131** for *S. diffusum* and **133** for *C. amboinicus*, both correct).
*De-novo assemblies get the extra "verify vs close relatives / reads across LSC–IR · SSC–IR" rule —
Phase 7 only, skipped here.*

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D3-01 | Which plastome annotator | ✅ **DONE 03/10 (dev)** — **the offline tool is the pipeline; GeSeq is a one-off cross-check whose output gets pinned.** Not a preference: the toolkit must run offline and ship in Docker, and every run must record tool name · version · argv — a web service can record none of that. PGA primary (GPL-3.0, batch, emits GenBank), plann.pl fallback. GeSeq's result is **pinned** via `panel pin`, which turns one manual web step into an offline reproducible input. Verified 03/10: GeSeq has no standalone distribution. PGA/plann **ABSENT** here; Perl present, BLAST+ needed | Dev | `docs/D3-01_annotator.md` + §3 |
| D3-02 | **Semantic gene-name aliases** — `clpP1` · `pafI` · `pafII` · `pbf1` (in *C. amboinicus* `PX923121.1`) | ✅ **DONE 03/10, dev, with evidence**: all four are **the same locus** as a classic gene, proven three ways in our own 9 records — identical syntenic neighbours (`pafI` between `psaA`/`trnS-GGA` where 7 records have `ycf3`), identical length (`pbf1` and `psbN` both exactly 132 bp), and matching `/product`. **Canonical = the team's own dialect** (`pafI`·`pafII`·`pbf1`·`clpP1`): those records are **the team's** (Vu,H.-T. = cô Trang), one is the **manuscript species**, and GeSeq — their tool — emits those names. Normalising to the classic names would make the team's paper print a name that disagrees with the team's own submitted record. Classic names travel as `source_name`, and output prints **both** (`pafI (ycf3)`). `gene_groups.yaml` gains the four names in their functional groups | Dev | PROJECT.md §8.10 + §3 |
| D3-06 | The name **printed in the manuscript** and **submitted to NCBI** | recommended: **print both** — `pafI (ycf3)` — which is what published gene tables already do and what removes the reader's need to know which dialect a record uses. The PI confirms the wording | PI (cô Trang) | PROJECT.md §3 |
| D3-03 | Which translation table for reading frames | the record's own `/transl_table` — **all 1,099 CDS in the dataset declare 11** — falling back to `defaults.yaml` (1) when absent. Tables 1 and 11 share amino-acid assignments and differ only in initiation codons, so P2-03/P2-04 numbers are unaffected | Dev | `data/defaults.yaml` + PROJECT.md §3 |
| D3-04 | `.sqn`: build it, or ship the modern path | emit **`.fsa` + `.tbl`** in pure Python (always); produce `.sqn` only when `table2asn` is present — measured **ABSENT** here. NCBI's web submission takes `.fsa` + `.tbl` directly | Dev | PROJECT.md §3 |

### Pitfalls for this phase *(measured on this dataset — PROJECT.md §8, §9)*

| Trap | Instead |
|---|---|
| "A CDS must start with ATG" | **5 of the golden's 86 do not** — rps19 GTG, rpl2 ACG ×2, ndhD ATC. Whole set: ATG 746 · GTG 11 · ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1. A rule that fails 3.7 % of real plastid CDS is a false alarm, not a check |
| Assuming translation table 1 because D2-03 says so | **every CDS here declares 11**; use the record's, and fail loud on a table whose amino acids we do not know |
| Ignoring `/codon_start` | all 1,099 declare `1`; a `2` shifts the frame — refuse rather than translate a frame the record did not ask for |
| "Fixing" a truncated gene | the golden's truncated ycf1 has `gene` + `misc_feature` and **no CDS** (§9.B-9). `count_genes()` is right; a repair is wrong |
| Renaming a gene so it looks canonical | a rename changes a submitted record. Unmapped names are **listed**, with their `/product`, never guessed |
| Reading the naming difference as "four odd genes" | it is a **whole dialect**: tRNA anticodon hyphen (8 records) · brackets (1) · **bare `trnK`/`trnI`** (2); rRNA `rrn16` (6) vs `rrn16S` (3); protein classic (7+golden) vs functional (2). Normalise the syntax *and* the four semantic names, or the comparison splits one gene into two |
| Reading a tRNA anticodon at a fixed offset | measured: offset 34 matches only **106/293** named tRNAs — the anticodon is not at a constant position here. Identify a bare name by **sequence identity** against the copies that do carry one (`trnK` ≡ `trnK-UUU`, 100.0 %; `trnI` ≡ `trnI-GAU`, 100.0 %) |
| Trusting an annotator's output blindly | read it back with our own reader and require the model to round-trip; a file that does not re-read is not a result |
| Assuming an annotator is a pure function of its input | PGA deletes the contents of its output folder, breaks on filenames with spaces, and mis-annotates a multi-FASTA target. Our fresh-per-run layout, `slugify()` and one-record-per-`--input` rule already make all three impossible — keep them that way |
| Pasting a web tool's answer into the pipeline | pin it (`plastid panel pin`: accession · date · sha256, re-read with our own reader) so the manual step happens once and every later run is offline and byte-identical |
| Writing a TBL from memory | the 5-column NCBI table has exact rules (partial flags `<` `>`, `join(...)`, `codon_start`, `transl_table` per feature) — check against a real submitted record |
| Quoting "132 genes" for a new genome | 132 is *Iris domestica*'s. Write the breakdown: **132 = 86 + 38 + 8** for the golden, whatever the record actually has otherwise |

### P3-01 · ORF checks — start · stop · internal codons · M3.3 · [ ]

- [ ] read the record's own `/transl_table` and `/codon_start` into the model (`Gene` gains both);
      refuse a table we cannot translate and a `codon_start` that is not 1
- [ ] classify every CDS: `ok` · `alternative_start` · `no_start` · `no_stop` · `internal_stop` ·
      `not_multiple_of_three`; pseudogenes marked, not silently counted as failures
- [ ] report genes **annotated without a CDS** as their own class — the truncated-ycf1 shape
- [ ] never rewrite a CDS to make it pass
**Done when:** on the golden genome all 86 CDS give 0 internal stops, 0 missing stops, every length
divisible by three, and **exactly 5** classified `alternative_start` (rps19 GTG · rpl2 ACG ×2 ·
ndhD ATC); the truncated ycf1 is reported as *annotated, not coding* rather than as a failure.
**Pitfall:** the 3.7 % alternative-start rate above. **Needs:** D3-03.

### P3-02 · Gene-name check & normalization · M3.2 · [ ]

- [ ] **one** normalizer in `core/gene_names.py` — and nothing else in the codebase may rename a gene.
      It covers the measured dialects: tRNA `trnL(UAA)` → `trnL-UAA`, rRNA `rrn16S` → `rrn16`, case,
      whitespace — **plus a bare tRNA name resolved by sequence identity**, not by guessing
      (`trnK` ≡ `trnK-UUU`, `trnI` ≡ `trnI-GAU`, both 100.0 % on this dataset)
- [ ] `data/gene_aliases.yaml` for **semantic** aliases — seeded with `pafI`·`pafII`·`pbf1`·`clpP1`
      per D3-02, each entry carrying the evidence it rests on; it ships **empty of guesses**
- [ ] `data/gene_groups.yaml` gains those four names in their functional groups (PSI assembly, PSII
      biogenesis, protease) — today they land in *unclassified*, which is how this was found
- [ ] a name with no mapping is **kept** and listed in the Result with its `/product`, never guessed
- [ ] a mapped name keeps its **`source_name`**, so a reader can always trace back to the record
**Done when:** the normalizer is idempotent over all **156** distinct names, leaves the golden's
**112** unchanged, and on the panel resolves `pafI`·`pafII`·`pbf1`·`clpP1` **and** `trnK`·`trnI` —
each while still reporting the name it arrived under; `gene_groups` reports **0 unclassified** on all
11 records; any name with no mapping is reported rather than renamed.
**Needs:** nothing — D3-02 is decided (PROJECT.md §8.10). Only the *printed wording* (D3-06) waits.

### P3-03 · Annotation — plastome (PGA) · M3.1 · [ ]

- [ ] `external/pga.py` through `run_tool()`; path, version and argv land in the run info
- [ ] the annotator's GenBank output is read back by **our own reader** and must round-trip
- [ ] `modules/annotation.py`: the resulting gene set, its ORF status (P3-01) and its name mapping (P3-02)
      go through the same pipeline as a downloaded record — no second code path
**Done when:** on a machine with PGA, annotating a de-annotated golden FASTA reproduces the golden's
gene counts by the fixed rule (86 CDS + 38 tRNA + 8 rRNA) and the output re-reads clean. Here:
implemented and `PTS_REQUIRE_TOOLS=1`-gated, **no number claimed from a run that did not happen**.
**Needs:** D3-01, and PGA on the machine.

### P3-04 · Annotation — mitochondrion (MITOS2) · M3.1 · [ ]

- [ ] the mito path: no canonical IR, so nothing may call `ctx.require_ir()` on it — measured on both
      stand-in mitogenomes, which raise `IRDetectionError` and are **right** to
- [ ] `external/mitos.py` through `run_tool()`, output re-read by our reader
**Done when:** on a mitogenome the gene-set table reports plant mito gene names and the run does not
touch the IR gate; on a machine with MITOS2 the count reproduces the record's own annotation.
**Needs:** MITOS2 on the machine.

### P3-05 · NCBI submission files — FASTA · TBL · SQN · M3.4 · [ ]

- [ ] `.fsa` and the 5-column `.tbl`, in pure Python, from the *normalized* gene set
- [ ] `.sqn` through `table2asn` when it exists; a named error and exit 3 when it does not
- [ ] a first-pass check on the record: no duplicate locus tags, no feature past the sequence end,
      `transl_table` on every CDS
**Done when:** a golden record produces `.fsa` + `.tbl` whose features re-read into the same set as
the source, and the TBL's partial/join syntax matches a real submitted record field for field.
**Needs:** D3-04. **Blocked for the `.sqn` half by:** `table2asn` (absent here).

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Which gene name the manuscript and the NCBI submission print (D3-06) | cô Trang | the *wording* of the gene table, not the pipeline |
| A machine with PGA and MITOS2 (or a run from anh Điềm's toolbox) | anh Điềm | the *measured* halves of P3-03, P3-04 |
| Item (c) text for the Reference Panel notes | cô Trang | D-03 only (documentation) |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| H-01 | MAFFT 7.526 on a machine that can run it | BLOCKED | P2-07/P2-08 real-data numbers (PROJECT.md §7.8 pins 7.526) |
| H-03 | P2-10's *C. amboinicus* label | **DONE 03/10** | re-run as `own-record-coleus-amboinicus`; both draft labels now state what they are (`standin-` for someone else's record, `own-record-` for the group's) — `docs/P2-10_manuscript_draft.md` |
| H-02 | Three merged Phase-1 branches (`feat/phase-1`, `feat/p1-02-core-model`, `feat/p1-04-result-schema`) are labels on commits already in `main` | OPEN | deleted or deliberately kept — commits are reachable from `main` either way |
