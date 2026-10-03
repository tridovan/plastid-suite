# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 2 — Genome reading & statistics  *(PROJECT.md §13.2)*

**Phase check (§13.2):** the golden numbers must come out of the module runs — 153,724 bp ·
**132 = 86 CDS + 38 tRNA + 8 rRNA** · IR 26,214 bp (≥ 99.9 %) · LSC 83,127 / SSC 18,169 ·
37 SSR (mono 31 · di 4 · tri 2). Any mismatch ⇒ stop and fix inside this phase.
Applies to **both plastid and mitochondrion**.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D2-01 | **Multi-genome contract** — π, SNP and later Ka/Ks take more than one genome; `run(genome, params, ctx)` is frozen in API v0.1 | extend to `run(genomes: Sequence[Genome], params, ctx) -> Result` (API v0.2) and migrate `genome_stats`/`ssr` to `genomes[0]` with an arity check | Dev | `docs/API_v0.2.md` + PROJECT.md §3 |
| D2-02 | π window size / step | 600 bp window · 200 bp step (the mVISTA-style convention in plastome papers); recorded in params + report either way | Dev | `data/defaults.yaml` |
| D2-03 | Translation table for codon/amino-acid work | NCBI table 1 (standard) for plastids; the table used is recorded in run info | Dev | `data/defaults.yaml` |
| D2-04 | Word/Excel writers | add `python-docx` + `openpyxl` as an optional extra so the minimum install stays pure-Python | Dev | `pyproject.toml` |

### Pitfalls for this phase *(every line has bitten somewhere — PROJECT.md §9, §7.9)*

| Trap | Instead |
|---|---|
| Writing "132 genes" without the breakdown | always **132 = 86 CDS + 38 tRNA + 8 rRNA**; `gene` features = 133 is a diagnostic only |
| "Fixing" the gene count by counting `gene` features | the truncated ycf1 has `gene` + `misc_feature`, no CDS — `count_genes()` is right |
| Assuming one contiguous span per CDS | rps12 is trans-spliced (3 parts, mixed strands) — build from `parts` in biological order |
| Trusting a `repeat_region` note as truth | IR comes from the gate-passing detection; the annotation is only a cross-check |
| Positional work without the gate | call `ctx.require_ir()` first; a failing gate must stop the run, not warn |
| Typing numbers into report templates | templates reference `{{table.column}}` only — the scan test enforces it |
| Skipping a file and carrying on | record the failure in the Result/manifest and name the file; exit non-zero when required |
| Hand-editing an output table | every table is regenerated from the module; if it cannot be, the module is wrong |
| Changing a window/step/threshold "just to see" | params go into the manifest and the report — always |

### P2-01 · Gene content & summary statistics · M2.1 · [ ]

- [ ] per-kind counts (CDS/tRNA/rRNA), GC, length summary — reuse `count_genes` and `gc_percent`
- [ ] run on both golden genomes and on the two diverged-IR genomes (the second must exit 2)
**Done when:** one `plastid run` reproduces 153,724 bp · 132 = 86+38+8 · GC 37.85 % for NC_050833.1
in the report and the CSV, with scope labels.

### P2-02 · Genome structure — boundaries · %GC per region · GC-skew · M2.3 · [ ]

- [ ] LSC/IRa/SSC/IRb boundaries and lengths (from the IR pair, normalized frame)
- [ ] %GC per region; cumulative GC-skew per region, plotted values tabulated
**Done when:** LSC 83,127 · SSC 18,169 · IRa = IRb = 26,214 · regions tile the genome exactly
(sum = 153,724) — the golden run asserts it.
**Pitfall:** GC-skew sign convention ((G−C)/(G+C)) and the cumulative form must be stated in the
column headers; a silent sign flip is this phase's classic quiet error.

### P2-03 · Codon usage — RSCU · ENC · CAI · CBI · Fop · M2.4 · [ ]

- [ ] CDS extraction: exons joined in biological order, strand-aware, pseudogenes excluded
- [ ] RSCU per codon + the family sums check (Σ RSCU per amino acid = family size)
**Done when:** the RSCU family sums check passes for every amino acid and the CDS set equals the
86 counted CDS. **Needs:** D2-03 (translation table).

### P2-04 · Amino acid composition (+ GRAVY · aromaticity) · M2.5 · [ ]

- [ ] translate the same CDS set used by P2-03 (one source of truth for the CDS set)
**Done when:** recomputable from the codon table in P2-03 — same CDS set, same totals.

### P2-05 · Introns — count · length · position · M2.7 · [ ]

- [ ] per gene: number of introns, lengths, positions, phase; trans-spliced genes handled
**Done when:** the intron list cross-checks with the GenBank annotation part boundaries.

### P2-06 · SSR by region and motif · M2.8 · [ ]

- [ ] summary tables by region (LSC/IRa/SSC/IRb) and by motif class, on top of the existing
      per-SSR table from `modules/ssr.py`
**Done when:** tier-B run with MISA reports 37 SSRs split by region and motif; the counts reconcile
with the per-SSR table. **Needs:** nothing — the tier-B harness already exists.

### P2-07 · SNP / substitutions by position · M2.6 · [ ]

- [ ] pairwise substitution table against a chosen reference, positions scope-labelled
**Done when:** the table states its scope (whole plastome vs CDS-only, IR in/out) and the numbers
reconcile with the alignment used. **Needs:** D2-01.

### P2-08 · Nucleotide diversity π (sliding window) · M2.9 · [ ]

- [ ] π per gene / intron / IGS + sliding window with region shading data
**Done when:** windows are reproducible from the alignments and the window/step are in the report.
**Needs:** D2-01, D2-02.

### P2-09 · Publication gene table (Word) · M2.2 · [ ]

- [ ] gene table export (Word) generated from the gene table, never hand-edited
**Done when:** regenerating the document from the module produces an identical file.
**Needs:** D2-04.

### P2-10 · Manuscript deliverable — Oldenlandia / Coleus · M2.10 · [ ]

- [ ] tables + figures + draft Results for the two species (numbers traced to table cells)
**Done when:** every number in the draft resolves to a cell in a run's tables (the queryable-number
rule, §13.5). **Needs:** the real datasets from anh Điềm/cô Trang.

### P2-11 · Sample-group labels on comparison output · M2.11 · [ ]

- [ ] every comparison table/figure carries a shared group label from `core/samples.py`
**Done when:** at least one comparison output uses a group label end-to-end.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Oldenlandia / Coleus GenBank records (real data) | anh Điềm / cô Trang | P2-10 |
| Item (c) text for the Reference Panel notes | cô Trang | D-03 only (documentation) |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
