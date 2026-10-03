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
| D3-01 | Which plastome annotator | **PGA** (GPL-3.0, Java, runs local) primary; **plann.pl** the offline fallback. Measured 03/10 on the dev machine: PGA **ABSENT**, plann **ABSENT**, `java` present | Dev | PROJECT.md §7.8 + §3 |
| D3-02 | **Semantic gene-name aliases** — do `clpP1` · `pafI` · `pafII` · `pbf1` get renamed? | **Not a dev call.** All four are in *C. amboinicus* `PX923121.1` — a manuscript species — with `/product` reading "ATP-dependent Clp protease proteolytic subunit 1", "photosystem I assembly factor I", "…II", "photosystem biogenesis factor 1". A rename changes what a **submitted record** says ⇒ PI | PI (cô Trang) | PROJECT.md §3 |
| D3-03 | Which translation table for reading frames | the record's own `/transl_table` — **all 1,099 CDS in the dataset declare 11** — falling back to `defaults.yaml` (1) when absent. Tables 1 and 11 share amino-acid assignments and differ only in initiation codons, so P2-03/P2-04 numbers are unaffected | Dev | `data/defaults.yaml` + PROJECT.md §3 |
| D3-04 | `.sqn`: build it, or ship the modern path | emit **`.fsa` + `.tbl`** in pure Python (always); produce `.sqn` only when `table2asn` is present — measured **ABSENT** here. NCBI's web submission takes `.fsa` + `.tbl` directly | Dev | PROJECT.md §3 |

### Pitfalls for this phase *(measured on this dataset — PROJECT.md §8, §9)*

| Trap | Instead |
|---|---|
| "A CDS must start with ATG" | **5 of the golden's 86 do not** — rps19 GTG, rpl2 ACG ×2, ndhD ATC. Whole set: ATG 746 · GTG 11 · ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1. A rule that fails 3.7 % of real plastid CDS is a false alarm, not a check |
| Assuming translation table 1 because D2-03 says so | **every CDS here declares 11**; use the record's, and fail loud on a table whose amino acids we do not know |
| Ignoring `/codon_start` | all 1,099 declare `1`; a `2` shifts the frame — refuse rather than translate a frame the record did not ask for |
| "Fixing" a truncated gene | the golden's truncated ycf1 has `gene` + `misc_feature` and **no CDS** (§9.B-9). `count_genes()` is right; a repair is wrong |
| Renaming a gene so it looks canonical | a rename changes a submitted record. Unmapped names are **listed for the PI**, with their `/product`, never guessed |
| Trusting an annotator's output blindly | read it back with our own reader and require the model to round-trip; a file that does not re-read is not a result |
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

- [ ] **one** normalizer in `core/gene_names.py` — syntax only (case, tRNA anticodon `trnL(UAA)` →
      `trnL-UAA`, `rrn16S` → `rrn16`, whitespace) — and nothing else in the codebase may rename a gene
- [ ] `data/gene_aliases.yaml` for **semantic** aliases, each entry carrying the evidence it rests on;
      it ships **empty of guesses**
- [ ] a name with no mapping is **kept** and listed in the Result with its `/product`, so the PI can rule
**Done when:** the normalizer is idempotent over all **156** distinct names in the dataset and leaves
the golden's **112** names unchanged; on the panel exactly **4** names are reported as needing a
decision (`clpP1` · `pafI` · `pafII` · `pbf1`, in the two *Coleus/Oldenlandia* records).
**Needs:** D3-02 for the semantic half only — the syntax half does not wait.

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
| Nomenclature ruling on `clpP1` · `pafI` · `pafII` · `pbf1` (D3-02) | cô Trang | the semantic half of P3-02, and every submission file for *C. amboinicus* |
| A machine with PGA and MITOS2 (or a run from anh Điềm's toolbox) | anh Điềm | the *measured* halves of P3-03, P3-04 |
| Item (c) text for the Reference Panel notes | cô Trang | D-03 only (documentation) |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| H-01 | MAFFT 7.526 on a machine that can run it | BLOCKED | P2-07/P2-08 real-data numbers (PROJECT.md §7.8 pins 7.526) |
| H-02 | Three merged Phase-1 branches (`feat/phase-1`, `feat/p1-02-core-model`, `feat/p1-04-result-schema`) are labels on commits already in `main` | OPEN | deleted or deliberately kept — commits are reachable from `main` either way |
