# TASKS.md — active tasks

> **Phase 5 is active.** Phases 4 and 4a are closed; their record is `ROADMAP.md`'s ticked
> milestones, `docs/log/`, and the commits — nothing finished is kept here.
> Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 5 — Figures & reports  *(ROADMAP.md → Phase 5 · PROJECT.md §13.5)*

**Goal:** publication-grade output, regenerated from intermediates — **never by hand**, and never a
number that exists only in a picture.

**In-phase checks:** every table/figure regenerates from intermediate data · **queryable-number
test**: every number in a report traces to **exactly one cell** in a source table · the tool report
carries name, version, path, threshold.

**What already exists (do not rebuild it).** The report engine resolves `{{table.cell}}` or raises,
and a scan test already fails on a bare digit in a shipped template (`find_bare_digits`) · CSV,
JSON, Markdown and **byte-reproducible Word** exports · every module already emits tables with
scope labels and a run manifest with tool provenance. What does **not** exist: any plotting code at
all (no matplotlib in `pyproject.toml`), multi-sheet Excel, PDF, and the automated check over
*generated* reports rather than templates.

**The figure list is 11 types; 7 have a source table and 4 do not.** Measured 04/10/2026 against a
full Phase-4 run directory:

| Figure (M5.1) | Source table | Exists? |
|---|---|---|
| circular plastome map | `structure.regions` + `gene_set` | ✓ |
| RSCU heatmap | `codon_usage` | ✓ |
| π sliding window, LSC/SSC/IR shaded | `diversity` windows + `structure.regions` | ✓ |
| IR junction | `structure.regions` + the sequence | ✓ |
| mVISTA-style %IDY + synteny | `synteny_identity` · `synteny_blocks` | ✓ |
| SSR / repeat distribution | `ssr` | ✓ |
| barcode assessment | `barcodes_candidates` · `barcodes_windows` | ✓ |
| **haplotype network** | `popgen` carries a haplotype **count**, not a haplotype table | **✗ needs the table first** |
| **annotated tree** | nothing in the toolkit builds a tree | **✗ no module** |
| **NUPT/MTPT by gene** | — | **✗ M6.1, BLOCKED on a nuclear genome** |
| **mito configuration** | — | **✗ M6.2 (Phase 6)** |

The four are named here so the phase's own milestone text is not read as "eleven figures are
waiting". None of them is scheduled in this phase.

### Decisions — made where the evidence is ours, surfaced where it is not

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D5-01 | The plotting library | **matplotlib, as an optional extra** `[figures]` — the D2-04 pattern (`python-docx`): imported inside the function, missing ⇒ `ExportDependencyError` + exit 3, never a traceback. It is the only offline, pip-installable library that writes **SVG · PDF · EPS** and draws circular maps, heatmaps and statistical plots; the minimum install stays biopython + PyYAML (CLAUDE.md) | Dev | PROJECT.md §3 |
| D5-02 | The statistics library | **scipy, as an optional extra** `[stats]` — and **every statistic pinned by a hand-computed or reference vector**. D4-03 wrote `popgen` by hand because the *tools* were a Windows GUI and a VCF-only CLI; that reason does not carry to scipy. Mann–Whitney with tie correction written by hand is exactly how a plausible wrong p-value gets in | Dev | PROJECT.md §3 |
| D5-03 | Where a figure gets its numbers | **A figure is a pure function of a `Result`** and reads table cells — it never recomputes. A figure that recomputes is a second answer to a question the table already answers, which is the §8.25 failure in a new costume | Dev | PROJECT.md §3 |
| D5-04 | Scaling rule for N genomes | **Measured, not assumed**: render the same figure at N = 2, 8 and 34 against the corpus and state the rule the measurements support (series cap, or a summary shape above it) | Dev | §8, `docs/` |
| D5-05 | What a generated sentence may claim | **Only what the table's own columns carry**, with n and the test named; no causal or adaptive wording. This is the **scope of a scientific claim**, so it is the PI's call, not ours — the default above is what the code will do until told otherwise | **PI** | PROJECT.md §3 |

### Pitfalls for this phase

| Trap | Instead |
|---|---|
| a figure that recomputes instead of reading the table | two numbers for one quantity — the §8.25 lesson; a figure reads cells, and a test asserts the cells it read |
| "it rendered, so it is right" | LibreOffice rendered a Word file whose styles the pruner had destroyed (§9.C-20) — **figures are looked at**, not only asserted on |
| a PNG/SVG that differs between two runs | matplotlib stamps a timestamp and a creator string into SVG/PDF — the Word export already solved this (fixed timestamp, sorted entries); the same rule applies, or "regenerates from intermediates" is false |
| a figure carrying the only copy of a number | tables are the record; a figure is a view. Every number stays queryable in a table cell |
| 34 series on one axis | the corpus makes this cheap to check — render it and look, then write the rule from what you saw |
| an EPS/PDF that no journal accepts | state the target: SVG for editing, PDF for submission, EPS only if a journal demands it — and check what actually opens |

### P5-01 · The figure system, and the first figure · M5.1a · [ ]

- [ ] `[figures]` optional extra + `core/figure.py`: a figure is **a function of a `Result`** and the
      theme (fonts, colours, sizes) is one place, not per-figure choices
- [ ] vector export — **SVG · PDF · EPS** — byte-identical between two runs of the same input
- [ ] figures land under the run directory beside the tables, named by the layout (`ctx.out`) — a
      module still never names a path
- [ ] the first figure end-to-end: **π sliding window with LSC/SSC/IR shading** (its tables exist)

**Done when:** `plastid run` on a pinned panel writes the figure into the run directory, a second
run of the same input produces a **byte-identical** SVG, and the figure's own data is asserted equal
to the cells of the table it read. **Needs:** D5-01, D5-03.

### P5-02 · The six figures whose tables exist, and the scaling rule · M5.1b · [ ]

- [ ] circular plastome map · RSCU heatmap · IR junction · mVISTA-style %IDY + synteny ·
      SSR/repeat distribution · barcode assessment
- [ ] each one reads its table and names it; none recomputes a quantity a module already produced
- [ ] **the scaling rule (D5-04)**, measured: the same figure at N = 2, 8 and 34 on the corpus

**Done when:** all seven figures (this task's six plus P5-01's) render from a **pinned run directory
alone** with no re-analysis, and the N-scaling rule is written down with the three measurements that
produced it. **Needs:** P5-01. **Pitfall:** at N = 34 an unreadable figure is a failed figure, not a
dense one.

### P5-03 · Statistics · M5.2 · [ ]

- [ ] descriptive · non-parametric tests · correlation · **FDR/Bonferroni** — each function naming
      its test, its n and its assumption, each pinned by a hand-computed or reference vector
- [ ] a table shape that carries the statistic **with its n and its corrected p**, never a bare p
- [ ] sentence-level interpretation: generated from table cells only, under D5-05's rule
- [ ] a statistic that needs more data than exists is reported *not computable*, never a number
      (the `popgen` rule, §8.30)

**Done when:** every statistic matches its pinned vector to a stated tolerance, an FDR table exists
over a real comparison, and no generated sentence carries a claim its source table does not.
**Needs:** D5-02; D5-05 for the wording only — the mechanism does not wait on the PI.

### P5-04 · The publication exports · M5.3 · [ ]

- [ ] multi-sheet Excel (openpyxl, into the same optional extra pattern)
- [ ] a **Supplementary package** — tables + figures + a manifest, regenerated from the run
      directory, never assembled by hand
- [ ] **one-page PDF per genome** — and it must survive the corpus: 34 records, not one
- [ ] a technical QC report: what ran, with which versions, on which inputs, what was skipped
- [ ] the existing Word export keeps working (P2-09) — a regression test, not a rewrite

**Done when:** the package is regenerated from a run directory alone, the one-page PDF is produced
for **all 34 corpus records without failure**, and the QC report names every tool · version · input
hash. **Needs:** P5-01 (figures), P5-03 (statistics) for the parts that carry them.

### P5-05 · The queryable-number check, automated · M5.4 · [ ]

- [ ] extend the existing template scan to **generated reports**: every number in a report resolves
      to **exactly one** cell of a shipped table
- [ ] fail the build on an orphan number, and on a number that resolves to two rows
- [ ] a negative test: a report with a planted orphan number must fail (a guard nobody has watched
      fail is not a guard)

**Done when:** `make check` fails when a number is planted in a report that no table cell carries,
passes on every report the toolkit produces, and the check runs over the **whole** run directory of
a Phase-4 deliverable run. **Needs:** nothing — it can be written first, and should be: it is the
phase's own check.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| A nuclear genome of a manuscript species | cô Trang | the **NUPT/MTPT figure** (M6.1) — not scheduled in Phase 5 |
| A haplotype table (`popgen` counts haplotypes today, it does not list them) | **us** — a task, once a figure needs it | the **haplotype network** figure |
| A phylogeny | us (M6.4) or an external tree | the **annotated tree** figure |
| A ruling on what a generated sentence may claim (D5-05) | **cô Trang** | the wording of P5-03's interpretation sentences only; the statistics themselves do not wait |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED on cô Trang | the text exists |
| H-04 | `trnG-GCC` vs `trnG-UCC` (§10.10) | OPEN | cô Trang decides the name printed in the manuscript table |
| H-05 | PMGA web service re-test (re-tested 03/10: still refuses connections) | OPEN | the `curl` in ROADMAP answers; then one web run + `panel pin --file` — nothing waits on it |
| H-06 | `data/mito_gene_names.yaml` holds **31** loci; the references carry **10 standard genes more** (`rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` · `rps14` · `rps19` · `sdh3`) | OPEN, measured 03/10 | each new locus carries the reference `/product` string that justifies it, and a test proves the vocabulary covers the set — it changes naming for **every** mito record, so it is its own task |
| H-07 | π per **intron/IGS** — the last piece of M2.9 | OPEN | intron spans are derived inside `modules/introns.py` and IGS spans nowhere: extract that derivation to the shared layer first (DRY's second copy) and feed `diversity`'s region machinery; a run whose regions cover LSC/IR/SSC **plus** every intron and IGS, cross-checked against `introns`' own count (17 on the golden genome) |
| H-08 | The run cache keys on `module_version` but not on the module's source: a behaviour change without a version bump serves a stale payload (measured 04/10, §9.C-21) | OPEN | either the cache key carries a source hash, or CI fails a changed module whose version did not move |
| H-09 | No module declares **which organelle it serves**, and no run checks that its **inputs** are single-organelle: the registry's only signal is `requires_ir`, which answers "does this need the IR gate", not "is this plastid or mito" — running `mito_genes`/`mito_annotate` on a plastome, or a comparison run over a **mixed** plastome + mitogenome input set, is blocked by nothing, while the other mismatch direction already stops at the gate | OPEN, decided 05/10 (D5-06, PROJECT.md §13 "Hai bào quan") | the registry carries each module's group (plastid · mito · both); a run whose inputs **disagree on organelle** exits with a named error; and a run naming a module outside the record's own organelle exits with a named error — negative tests for the mixed input set and for **both** mismatch directions. The Phase-6 transfer analysis (nucleus ↔ plastid ↔ mito) is the named exemption, and the check must not forbid it. The run manifest also records the run's organelle — today it carries `ir` only, so a run's own report does not say which organelle it described |
| H-10 | Seven shared-by-nature analyses are still welded to plastid regions, so they stop at the IR gate on a mitogenome: `ssr` · `introns` · `codon_usage` · `amino_acid` · `substitutions` · `diversity` · `synteny` | OPEN, measured 03/10 (§8.6: both stand-in mitogenomes raise `IRDetectionError`) | each module either runs on a pinned mitogenome with its region column explicitly *not applicable*, proven by a real run on the project's own `AP045469.1` annotation — or is recorded as staying plastid-only, with the measured reason |
| H-11 | An artifact this toolkit writes does not declare its organelle: `plastid mito annotate` copies the input's (empty) `source`, so its GenBank carries **no `source` feature and no `/organelle`** (measured 05/10: 0 of each in `.panel/mito-prototype/AP045469.1.gb`) — `aliases_for` then falls back to the plastid vocabulary, and `submission` refuses a record with no source organism (`core/submission.py:377`) | OPEN, measured 05/10 (D5-06) | every path that produces a record writes a `source` carrying at least `/organelle` (mitochondrion on the mito path), proven by reading the artifact back in a test; the `mito-prototype` artifact is rebuilt and **re-pinned**, with the new sha256 recorded |
