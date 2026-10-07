# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phases 1 – 5 and 4a — COMPLETE and merged to `main`.** **Phase 5 (figures & reports) closed
05/10/2026**: five tasks, each Done-when measured, recorded in `docs/log/phase-5.md` (§8.35–§8.39)
and ticked in `ROADMAP.md`. `TASKS.md` holds no open phase — only the Housekeeping table.

**Immediate next action:** the next board is **Phase 3B — the self-built plastid annotator**
(D5-08), whose slot is an open decision — and with it the slot for the **D5-09 naming step**
(reference-DB audit + `gene_normalize`, which stands in front of 3B: the annotator's reference set
must be name-normalized first). Designs and open questions: `docs/D5-08_plastid_annotator.md` ·
`docs/D5-09_reference_naming.md`; ROADMAP.md → Open decisions. Nothing in Phase 5 waits on either,
and nothing else is scheduled: Phase 6 waits on mito/nuclear data, Phase 7 on real reads.

**What Phase 5 left in the toolkit.** Three optional extras, none of them on the minimum install:
`[figures]` (matplotlib), `[stats]` (scipy), `[excel]` (openpyxl) — each imported inside the
function, each missing one a named `ExportDependencyError` + exit 3 · **eight figure types** in
`result/figures/`, rendered at run level from the run's own `result.json` · `plastid run --figures`,
`plastid figures`, `plastid package`, `plastid qc` · `core/stats.py` + the `stats` module ·
`report/audit.py` (the queryable-number check) · `read_result`, which makes a run directory
readable without re-running anything.

**The phase's own numbers:** SVG · PDF · EPS byte-identical between two runs (a cached second run
included) · the N-scaling rule **measured** at 2 · 8 · 34 (a series cap of 8, then the span the
genomes cover) · 320 statistical tests over the stand-in panel with 0 significant after
Benjamini–Hochberg, and 3 of the panel's 11 records carrying no annotation at all · **34/34**
corpus records producing a one-page PDF, none failing · the audit clean over the whole Phase-4
deliverable run (12 reports, 44 numbers).

**Three layout defects were found by looking at rendered figures**, not by tests — the §9.C-34
lesson, paid three times : overlapping region bands (§8.36) · a barcode band shading the whole axis
plus a GC-skew ring with 50 % gaps (§8.37) · summary lines drawn inside the GC-skew panel (§8.39).

**New housekeeping item, measured 05/10:** **H-14** — `--docx` is silently lost on a cache hit
(two runs of `gene_table`, the second `cached` with 0 `.docx` files), and the manifest's `command`
does not carry `--docx` either. Same family as H-08.

**On hold by D5-07 (plastid first):** H-06 · H-10 · H-12. **Open:** H-05 · H-07 · H-08 · H-09 ·
H-11 · H-13 · H-14. **Waiting on the PI:** H-04 (`trnG`) and D5-05 (what a generated sentence may
claim — the mechanism is built and follows the conservative default).

## Log

### 2026-10-07 · Claude · dev WSL · session 22 (doc alignment with D5-09 — and the cache defect the review exposed)

**Did.** Read `PROJECT.md` and `CLAUDE.md` in full against the new D5-09 design record (the
reference-driven naming step) and fixed every contradiction the review found.

* **CLAUDE.md** — the pipeline sentence read `assembly → annotate → gene-name unification →
  analyses`; D5-09 reverses it (names are unified **on the reference DB, before annotation**; the
  step after annotation only **verifies**). Rewritten, plus a D5-09 paragraph in the briefing (the
  naming reference is a sixth kind, chosen by the end user from the audit report; D3-02 is
  reversed; `gene_normalize` is a new module) and a "who decides" row for the ★ choice. Also fixed
  **"Current state: Phase 5 in progress"** — stale since 05/10; the file now says phases 1–5 + 4a
  complete, next board 3B, with its slot (and D5-09's) open.
* **PROJECT.md §3** — D3-02 marked `[ĐÃ THAY 07/10/2026]`; the D4a-05 alias row points at the
  reversal; the **D5-09 decision row** added (it did not exist). §13.3's naming bullet marked as
  moved. §13.0's `[5] RNA editing` said **PREPACT** while the table beneath it said "chưa chọn" and
  §7.4/§7.5 measured PREPACT unreachable (`http=000`) — the diagram now says "chưa chọn"; the D5-09
  copy synced and re-diffed **byte-identical**.
* **D5-09** — §5's status column used "the reference" for two different things; `measured` /
  `not measured` now say **in A**, `outside vocabulary` says **in B**. §9 gained the three gaps the
  review found: the **exit code** for a conflict stop (2 is already the IR gate's — and D5-08's
  audit's), the **organelle group** `gene_normalize` must declare in §13's table, and the fact that
  **§7.2's cache promise cannot hold today**.
* **ROADMAP** — a second open-decision row: the naming step's slot, which interacts with 3B's.
* **The cache defect, measured and named H-15.** `Cache.key` sorts input hashes and `_input_hash`
  hashes the sequence only, so the corpus's own "one sequence, two annotations" pair shares a key in
  either order — while six modules take `genomes[0]` as their reference. Measured on
  `tests/data/NC_050833.1.gb` (86 CDS) vs the corpus PGA copy (85 CDS): equal `_input_hash`, equal
  key with the inputs swapped. The run directory's fingerprint **is** order-sensitive (`cli.py` —
  "ordered: order is meaning"); only the cache is not. → TASKS.md **H-15**.
* **Docs corrected as D5-09 §8 asked**: `gene_name_variants.md` §6 (direction note; the two CGAS
  citations fixed to the measured lines — `return candidates[0]` at **211 and 214**, not 208; the
  prefix match at **149**, not 150–153) · `D5-08` line 54 (pointer: the reference set must be
  name-normalized before the annotator uses it) · `gene_aliases.yaml` header (D3-02's direction
  marked as reversing when the naming step lands — values unchanged until then, so code and tests
  stay consistent).

* **Two more false statements the stale sweep found, corrected in place**: `README.md` still said
  "Phase 5 in progress" and "Nineteen analysis modules" (it is 20 — `stats` arrived in Phase 5);
  `core/gene_names.py`'s docstring still stated D3-02's direction as current.
* **Three more places still stated the old pipeline order** (CLAUDE.md had been fixed first):
  ROADMAP's Phase-3B goal, D5-08's pipeline diagram and its "next steps" bullet, and the D5-08
  decision row in PROJECT.md §3.
* **The four dead items removed**, one commit: `docs/PHASE3_BLOCKERS.md` (its "still blocked — a
  decision on the mito annotator" item was false since 03/10; its one live ask — who `Vu,H.-T.` is
  / whether `PX923121.1` is the manuscript's specimen — moved to `TASKS.md` → "Waiting on data /
  people", and `environment.yml`'s pointer to the file was re-pointed) · `TASKS.md` **H-05** (PMGA
  web re-test — nothing waits on it; the mito path is settled by D3A-02 and frozen by D5-07) ·
  PROJECT.md **§13.10's** 26/9 snapshot table (heading kept — a log entry cites it; every row of it
  now lives in §8 or §7.1) · PROJECT.md **§12's** unpinned install block (the §9.B-14 hazard;
  `environment.yml` is the pinned source, and §7.8's pointer was re-pointed). `docs/tools/pmga.md`
  **stays** — it is the marked record of a measured dead end and keeps a later agent from
  re-investigating PMGA.

**Validation.** `make check` green (lint · 668 fast tests · golden with named skips). The
§13.0 ↔ D5-09 §1 diagram diff is empty. Committed as two commits: the design record as received,
then the alignment.

**Residual risks.** (1) The D3-02 reversal is **documented, not built** — the alias table and the
three modules that read the canonical key still run the old direction until `gene_normalize` lands.
(2) `gene_normalize`'s exit code and its organelle group are named open in D5-09 §9, not decided.
(3) H-15 is measured but not fixed.

**Blocked.** Nothing.

**Next.** The two slots (3B · D5-09) are the one open decision; then the Phase 3B board.

### 2026-10-05 · Claude · dev WSL · session 21 (Phase 5 — figures, statistics, publication exports; and the queryable-number check)

**Did.** The whole Phase 5 board, task by task, each committed and pushed on `feat/phase-5`, then
merged to `main`. Order followed the board's own advice: the phase's check first (P5-05), then the
figure system, the rest of the figures, the statistics, the exports.

* **P5-05 · the queryable-number audit.** `report/audit.py` requires every number in a generated
  report to be a bounded number inside the rendered value of a cell of that module's own tables.
  The strict rule was run against the Phase-4 deliverable first and failed with **10 findings over
  its 12 reports, in 6 classes** — every one a quantity colliding with a number printed inside a
  *quoted text* column (a block id `…#10`, `mean + 1 sd …`, `Coleus: 6, Oldenlandia: 5`, a note).
  A `str` value containing a letter is prose: set aside by name, never counted as a cell carrying
  a number. After that rule the deliverable run is clean, and planting an orphan in a real report
  is caught by name. `read_result` came with it — the JSON round-trips into a `Result`.
* **P5-01 · the figure system.** A figure is a pure function of a `Result`: `series()` reads cells
  into plain data, `build()` draws through one theme. The plan said `core/figure.py`; the layer map
  lets `core/` import only `plastid.core`, and two figures read tables from **two different
  modules**, so the system lives in `result/` and renders at **run level from the stored
  `result.json`** — which is why a fully cached run still gets its figures and `plastid run
  --figures` and `plastid figures <run-dir>` are the same function. Determinism was **measured, not
  assumed**: SVG carries a per-process salt *and* a `dc:date` (both pinned), PDF has a metadata
  hook, the PostScript backend has none so its `%%CreationDate` is rewritten.
* **P5-02 · six more figures and the N rule.** `identity_synteny` was rendered at N = 2 · 8 · 34:
  at 8 the series read, at 34 the legend covered the panel and 33 lines merged into one mass. The
  rule is a **cap of 8**, above which the panel draws the span the genomes cover — both edges still
  cell values — and the canvas grows with N.
* **P5-03 · statistics.** `core/stats.py` (scipy, `[stats]`) with every function pinned by a vector
  computed **by hand**; writing those vectors caught two of my own errors. The `stats` module runs
  a real comparison over the stand-in panel: 80 loci × 4 measures = 320 tests, 17 not computable,
  **0 significant after Benjamini–Hochberg** — the honest answer at 4 vs 4 observations per locus.
* **P5-04 · the exports.** Byte-identical multi-sheet Excel, `plastid package` (the Supplementary
  package, regenerated from a run directory), the technical QC report, and a PDF-only one-page
  sheet per genome: **34/34** corpus records, none failing, each one page, measured at 114 s for
  the whole corpus.

**Validation.** `make check` green after every task. Final: **668 fast tests**, golden **60 passed,
4 skipped** (PGA · chloe · julia · table2asn absent on this machine) with the tool environment on
PATH. Every figure was rendered from a real run and **looked at**.

**Residual risks.** (1) The figure and export extras are not installed by a bare `make env` —
`uv sync --all-extras` is what this tree uses. (2) The one-page PDF is A4 portrait with a fixed
layout; a genome with very many regions would crowd it, and no such record is in the corpus.
(3) The N cap of 8 is measured on *this* figure at N = 34; another figure may support more series,
and the rule is stated per figure rather than globally.

**Blocked.** Nothing.

**Next.** Phase 3B (the self-built plastid annotator, D5-08) once its slot is confirmed — or the
open housekeeping items, of which **H-14** was found and measured this session.

### 2026-10-05 · Claude · dev WSL · session 20 (doc ↔ code sweep: false statements fixed in place, §9 renumbered, the log split)

**Did.** A consistency review of `PROJECT.md` · `CLAUDE.md` · `ROADMAP.md` · `TASKS.md` · `docs/`
against each other and against the code, then the fixes it found. One rule clarified in CLAUDE.md:
**a statement that has become false is corrected in place** — history stays in §8/§9 (numbers
untouched), but a stale sentence in a spec is a wrong instruction to the next agent.

* **Superseded numbers leave the live text.** §13's Ts/Tv bullet and §13.4's in-phase check now
  carry **2.6917 · 79 genes · 7 blocks · 97.37** (was 2,686 · 76 · 21 blocks); §10.1 and §12
  follow; §7.8's MAFFT row states what was measured without quoting the old golden number; §13.10's
  two stale rows carry inline provenance like §7.1/§7.4 already did.
* **"Applies to both organelles" is true as written again** — §13's chapeau, §13.2/§13.4 and
  ROADMAP's cross-phase rule + Phase 2/4 headings now say the seven transitional analyses are
  plastid-only today (H-10, on hold). §13.3's "(132) recount" reads like ROADMAP's Phase 3 check:
  the *rule* is fixed, the *number* is the golden genome's.
* **§9 renumbered to unique, continuous numbers 1–40** (A 1–8 · B 9–26 · C 27–36 · J 37–40). The
  duplicates were real and had already bitten: `20` appeared three times, and `§9.C-20` was written
  in one place for the NCBI-filter item and in another for the style-pruner item. Every citation
  repo-wide updated (~25 sites incl. code docstrings); `test_project_docs.py` gained a guard that
  fails on a repeat or a gap.
* **Decisions that had no row:** **D2-03** (compiled table 1 for codon/AA; §3) and **D5-05** (§10.11
  as an open PI question — TASKS claimed it was "recorded in §3" and it was not).
* **New §9.J — four measurements** made while auditing the two-organelle rule in code:
  `gene_names` never loads the mito vocabulary (0 renamed / 0 aliased vs `mito_genes` naming 31 of
  51); the translation-table column implied the record's own table (5 of the golden's 86 CDS start
  GTG/ACG/ATC); our mito artifact has **0** source features and **0** `/organelle`; and
  `submission` on that artifact reports **24** first-pass problems (23 = no `/product`) while still
  writing the files. **Only the table label was fixed** (code, no number moved); the rest are named
  — **H-12** (new), plus H-09 and H-11 — and **on hold with H-10/H-06 by D5-07**.
* **Housekeeping:** `docs/handoff-archive.md` holds sessions 1–13 (a pure move, re-concatenation
  verified byte-identical); `pmga.md` gets the 05/10 update note its siblings already had; README's
  "Phase 1" status; `data/…` → `src/plastid/data/…`; the "Nhật ký & tài liệu" pointer → "The run
  log"; ROADMAP's duplicate M1.6/M1.8 lines, three decided rows and M4a.3 (ticked with its named
  exceptions); PHASE3_BLOCKERS' stale "can start immediately" list.

**Validation.** `ruff check` + `ruff format --check` clean · **567 fast** tests · golden with the
tool env on PATH: **57 passed, 4 skipped** (chloe · julia · table2asn absent on this machine).

**Residual risks.** (1) The §9 renumbering touched citations in code, tests and `docs/` — all
updated, and the new guard fails on a duplicate or a gap. (2) §9.J.38 *states* the initiator
convention (V/T/I where a published protein begins with M) rather than changing it; changing it
moves 5 residues of the golden composition — a decision, not a fix. (3) The HANDOFF archive is a
move; a citation by **content** rather than session number should name the archive file.

**Blocked.** Nothing.

**Next.** Phase 5 — **P5-05** (queryable-number check over generated reports) or **P5-01** (figure
system + first figure).

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


### Older sessions — moved 05/10/2026

Sessions 1–13 (02/10 – 03/10/2026) are in [`docs/handoff-archive.md`](docs/handoff-archive.md) —
a pure move of the text; no entry was retyped or renumbered. This file keeps the resume point and
the sessions from Phase 3A onward.
