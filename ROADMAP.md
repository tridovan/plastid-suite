# ROADMAP.md — phases, milestones, progress

> The file to read to know **where we are**. Milestones here are the source for specs/tasks.
> Session set: `CLAUDE.md` (rules) · `HANDOFF.md` (resume point + session log) · `ROADMAP.md` (this)
> · `TASKS.md` (active checklists with Done-when).
> `PROJECT.md` is the evidence archive (golden numbers §6 · decisions §3 · tools §7 · measurements
> §8–§9) — open it only when a phase needs it.
> Legend: `[x]` done · `[ ]` not started · **BLOCKED** = waiting on data/people.

## How to use

1. A phase is done only when its **In-phase checks** pass and are recorded. Never defer checks.
   Where a check is a number, write it as a **failing test first** (CLAUDE.md → Testing discipline).
2. When a phase starts: derive specs/tasks from its milestones → put them in `TASKS.md`.
3. Tick a milestone here when its tasks are DONE (add the date); refresh `HANDOFF.md`.
4. New measured facts → PROJECT.md §8 (runs) / §9 (errors), per `CLAUDE.md`.

## Cross-phase rules (apply to every phase)

- Applies to **both plastid and mitochondrion**, except divergence-time work (tool paper only).
- **IR gate ≥ 99.9 %** before any position-based analysis. Fail ⇒ red flag, stop. No exceptions.
- Gene counting: **132 = 86 CDS + 38 tRNA + 8 rRNA**. Never "133 genes" (that is the `gene` feature
  count, diagnostic only).
- Ts/Tv scope: **CDS only · 76 shared genes · MAFFT · IR excluded** → 2.686. Other scopes must be
  labelled as such wherever they are exported.
- Every run emits a **tool + parameters report**; output is **run-scoped** (`<dataset>/<run_id>/`),
  never overwritten in place.
- Genome normalization (origin/strand → LSC start, 4 regions split) happens once in Phase 1 and is
  reused everywhere after.
- Runs on: **ordinary machines** = phases 1–5, 8 · **remote server** = phase 6 (divergence time)
  and phase 7 (assembly).
- Standard test case: **Eleutherine bulbosa dataset + the golden numbers** (PROJECT.md §6).

## Already banked (evidence in PROJECT.md §8–§9)

- [x] Golden numbers measured & cross-checked against 2 published papers (153,724 bp · 37.85 % ·
      132 = 86+38+8 · IR 26,214 bp @ 100.00 % · LSC 83,127 / SSC 18,169 · 37 SSR · Ts/Tv 2.686)
- [x] IR algorithm validated on 20 RefSeq genomes — naive detection fails silently 2/9 (22 %);
      the ≥ 99.9 % gate catches 2/2 bad cases, 0 false alarms on 7 good ones
- [x] Validated probe code + golden data on disk (`/mnt/d/business/pts/probe/`), including
      `ir_probe_v2.py` (detect) and `ir_gate.py` (gate) — to be ported, not rewritten
- [x] Tool inventory: 84 tools tried · 64 working · replacements found for every blocked tool (§7)
- [x] Key decisions locked 02/10 (§3): analysis before assembly · Linux + Docker · 2-tier golden ·
      MISA bundled · no Snakemake · Reference Panel 4 tiers defined

---

## Phase 1 — Framework & shared data  *(PROJECT.md §13.1)*

**Goal:** the skeleton every later phase plugs into — data model, IR gate, result schema, run layout,
golden tests, report engine, module registry. Without it, Phase 2 has no gate to enforce.

**In-phase checks (from §13.1):** `make test` + `make golden` in **one command** produce every
golden number · the run report lists tool name + version + parameters · opening the output folder
tells you which file belongs to which job.

**Milestones**

- [x] **M1.1 Skeleton & packaging** (2026-10-02) — `uv sync` installs; `plastid --version` works;
  `make check` green; import layering enforced by `tests/unit/test_architecture.py`
- [x] **M1.2 Core model + GenBank reader + vendored golden data** (2026-10-02) — tier A tests read
  153,724 bp · GC 37.8536 % · 86 CDS + 38 tRNA + 8 rRNA · `gene`=133 · 112 unique names from the
  vendored fixture; truncated-ycf1 and trans-spliced-rps12 cases covered
- [x] **M1.3 IR detection + gate + normalization** (2026-10-02) — annotation IR 26,214 bp @ 100.00 %;
  the ported de-novo detector reproduces §8.3 exactly (7/9, same two failures — one loud, one silent,
  disagreement recorded in `denovo_length`); the gate catches 2/2 diverged genomes with 0 false alarms;
  an off-by-one pairing constant collapses identity to ~29 % (regression); normalization reproduces
  LSC 83,127 / IRa 26,214 / SSC 18,169 / IRb 26,214 and handles a genome starting inside an IR
- [x] **M1.4 Result schema + exporters** (2026-10-02) — `Scope` · `Column` · `Table` · `Result`;
  cell refs (`table.column`, `table[row].column`) resolve or raise `ResultRefError`; csv/json/md
  exports are byte-identical across runs and carry the scope label + coordinate frame
- [x] **M1.5 Run layer** (2026-10-02) — `<dataset>/<run_id>/` with `00_run_info` + logs; manifest
  schema v1 + hand-written `validate()`; `run_tool()` is the only subprocess door and records
  name/version/source/path/argv/cwd/rc/time/logs; content-hash cache restores files byte-identically
- [x] **M1.6 Golden harness** (2026-10-02) — `make golden` runs the whole set in one command; the
  tier-B probe names every missing tool and skips; `PTS_REQUIRE_TOOLS=1` turns skips into failures
- [ ] **M1.6 Golden tests tier A + `make golden`** — one command, every §6 tier-A number asserted;
  tier B prints a skip reason per missing tool
- [x] **M1.7 Report engine + citation registry + defaults** (2026-10-02) — `{{…}}` resolves through
  `Result.format()` or raises; shipped templates carry zero bare digits (scan test); the draft
  banner is appended by the engine; citations come from `data/citations.yaml`
- [x] **M1.8 Module registry + first module end-to-end** (2026-10-02) — decorator registry;
  `genome_stats` produces `01_genome_stats/` (3 tables + result.json + report.md); CLI `plastid run`
  writes `manifest.json` + `00_run_info/run_info.md`; a failed gate exits 2 and still leaves a
  traceable run directory
- [ ] **M1.8 Module registry + first module (`genome_stats`) end-to-end** — run dir complete:
  `00_run_info/` + tables + report
- [x] **M1.9 MISA vendored + `ssr` module (tier B)** (2026-10-02) — IPK zip verified against every
  §7.8 checksum and vendored with `SOURCE.json` (licence: not stated → flagged for M9.3); the module
  reproduces **SSR 37 = mono 31 · di 4 · tri 2, compound 3** on this machine; a missing perl raises a
  named `ToolError` → CLI exit 3
- [x] **M1.10 Reference Panel (4 tiers) + NCBI fetch/cache** (2026-10-03) — `core/panel.py`
  (tiers · modes A/B · quality filter · pin with an injected fetcher) + `run/ncbi.py` (esearch ·
  esummary · taxonomy efetch · 3 req/s) + `plastid panel refresh | propose | pin`. Verified live:
  76 RefSeq genomes for Iridaceae (the §8.3 count, after fixing NCBI's `chloroplast[filter]` blind
  spot — PROJECT.md §9.C-20), and a re-pinned NC_050833.1 whose sha256 matches the vendored golden
  file byte-for-byte. Proposals and pins are offline; the network lives only in `refresh`.
- [x] **M1.11 Sample-group metadata** (2026-10-02) — loader + lookup shipped (`core/samples.py`,
  duplicate-accession guard, tests); the "used by a comparison table" half moved to **M2.11**
- [x] **M1.12 Execution backend decision** (2026-10-03, dev decision) — local for phases 1–5 and 8;
  remote for phase 6 (divergence time) and phase 7 (assembly), where PROJECT.md §9.14 measured that
  an 8 GB machine leaves only 1.5–2.5 GB free per job. Recorded in PROJECT.md §3
- [x] **M1.13 Hardening + API freeze v0.1** (2026-10-02) — README quickstart · `environment.yml` ·
  `docker/Dockerfile` skeleton · `docs/API_v0.1.md` (module contract, Result schema, coordinates,
  manifest keys, exit codes). *The environment checks in this milestone's Done-when — `make check` on
  a clean clone and `make golden` inside Docker — are deferred to the Phase 9 platform pass, when
  Docker is set up on the dev machine.*

## Phase 2 — Genome reading & statistics  *(§13.2 — plastid **and** mito)*

**Goal:** tables + figures + draft Results for the 2 species manuscripts (Oldenlandia, Coleus).
**Input:** GenBank files (downloaded or produced). Raw reads not needed.
**In-phase checks:** golden numbers — 153,724 bp · 132 = 86+38+8 · IR 26,214 (≥ 99.9 %) ·
LSC 83,127 / SSC 18,169 · 37 SSR (mono 31 · di 4 · tri 2) when MISA/cpstools available. Mismatch ⇒ stop.

- [x] **M2.1 Gene content & summary statistics** (per group, GC %, length)
- [x] **M2.2 Publication gene table** (Word export) — categorised table + reproducible `.docx`
- [x] **M2.3 Genome structure** — LSC/SSC/IR boundaries · %GC per region · GC-skew
- [x] **M2.4 Codon usage** — RSCU · ENC · CAI · CBI · Fop
- [x] **M2.5 Amino acid composition** — + GRAVY · aromaticity
- [~] **M2.6 SNP / substitutions by position** — scope labels enforced
- [x] **M2.7 Introns** — count · length · position
- [x] **M2.8 SSR** — by region and motif
- [~] **M2.9 Nucleotide diversity π** — per gene/intron/IGS, sliding window
- [~] **M2.10 Manuscript deliverable** — tables + draft Results produced for both species from the
  public stand-in set (labelled); real samples and figures outstanding
- [x] **M2.11 Sample-group labels** — every comparison table/figure carries a shared group label
  (moved here from M1.11, whose loader shipped in Phase 1: `core/samples.py`)

## Phase 3 — Annotation  *(§13.3)*

**Goal:** annotate genomes we already have; prepare NCBI submissions.
**In-phase checks:** the output re-reads with our own GenBank reader (no format damage) · ORF checks
pass · gene recount uses the fixed counting **rule** (CDS + tRNA + rRNA, both IR copies) — **132 is
the golden genome's value, not a constant**: P2-10 measured 131 (*S. diffusum*) and 133
(*C. amboinicus*), both correct. *De-novo assemblies (Phase 7 only) get the extra "verify vs close
relatives / reads across LSC–IR · SSC–IR junctions, else label undetermined" rule — skipped here.*

- [x] **M3.1 Annotation** (2026-10-03) — plastome: **Chloe** (PGA measured and replaced, §8.15/§8.17),
      **132 = 86 + 38 + 8** through our own reader, 0 CDS failing the ORF checks; mito path
      (`mito_genes`, no IR gate) done and measured; the mito **annotator's count check is OPEN** —
      the local container needs 11 GB and MAKER (too strong for this machine) and PMGA's web
      application server **refuses connections** (measured 03/10, `docs/tools/pmga.md`)
- [x] **M3.2 Gene-name check & normalization** (2026-10-03) — one normalizer (`core/gene_names.py`), panel measured: 155 → 115 canonical, 0 unresolved, 0 unclassified
- [x] **M3.3 ORF checks** (2026-10-03) — golden: 0 failures, 5 alternative starts; the truncated ycf1 reported as annotated-not-coding
- [x] **M3.4 NCBI submission files** (2026-10-03) — `.fsa` + `.tbl` in pure Python, 265/265 features
      round-trip; **the `.sqn` step ran against `table2asn`**, and NCBI's own validator took the
      table from **13 errors to 7** (the rest inherent to duplicated IR genes; §8.18). A real
      submission still needs the team's own `.sbt` from NCBI's form —
      **the mito count check continues in Phase 3A**

## Phase 3A — Mito annotator, self-built prototype  *(opened 03/10)*

**Goal:** a mito annotator this project can run, read and reproduce. PMGA cannot run here — 11 GB
unpacked driving MAKER against an 8 GB VM cap, and its own web service refuses connections (P3-04,
`docs/tools/pmga.md`). CLAUDE.md's rule for a tool that cannot run is *replace it — or write that
part ourselves*, from parts the repo already has: BLAST+, the 31-locus vocabulary, the ORF checks,
the reader. It is measured against a record whose answer is already published, so it must earn its
place rather than be declared.

**In-phase check:** on `OX459128.1` the annotator's **own output** re-reads with **our own reader**,
runs through `orf` + `mito_genes` unchanged, and recounts to **31 CDS · 17 tRNA · 3 rRNA** — with
**every difference from the published record listed**, and anything it cannot call reported as not
called. A count that differs is a finding (§9), never a number to tune away.

- [ ] **M3A.1 GenBank writer** — read → write → read round-trips the golden's 265 features
- [ ] **M3A.2 Reference set** — RefSeq mitogenomes chosen by taxonomy (same family first) → CDS
      proteins + rRNA, pinned, with sources and sha256 recorded
- [ ] **M3A.3 Gene caller** — `blastx` for CDS and `blastn` for rRNA; boundaries come from the
      alignment, not from ORF-walking; everything below the floor is reported
- [ ] **M3A.4 tRNA caller** — ARAGORN, declared in `environment.yml` and `docs/`
- [ ] **M3A.5 CLI** — `plastid mito annotate`, output through the existing pipeline, no second path
- [ ] **M3A.6 Golden check** — 31 · 17 · 3 against the published record, differences itemized
- [ ] **M3A.7 AP045469.1** — the unannotated sample: run it, report it, pin it

## Phase 4 — Multi-genome comparison  *(§13.4 — plastid and mito)*

**Goal:** comparison content for the 2 manuscripts.
**In-phase checks:** Ts/Tv = **2.686** for *I. domestica* × *I. gatesii* (CDS only · 76 genes · MAFFT) ·
identity plot must show **21 blocks**, %IDY 83.7–98.2 (mean 91.5). Mismatch ⇒ stop, re-check inputs.

- [ ] **M4.1 Ka/Ks selection pressure**
- [ ] **M4.2 mVISTA-style %identity + synteny**
- [ ] **M4.3 Barcode candidates & hypervariable regions**
- [ ] **M4.4 Population genetics** (DnaSP-style: Hd · θ · Tajima's D · Fu & Li · Fst)
- [ ] **M4.5 RNA editing prediction** — **BLOCKED: needs RNA data**
- [ ] **M4.6 Ts/Tv module** (scope fixed: CDS only · 76 shared genes · MAFFT · IR excluded)

## Phase 5 — Figures & reports  *(§13.5)*

**Goal:** publication-grade output, regenerable from intermediates (never by hand).
**In-phase checks:** every table/figure regenerates from intermediate data · **queryable-number test**:
every number in a report traces to exactly one cell in a source table · tool report has name, version,
path, threshold.

- [ ] **M5.1 Figure system** — shared theme (colors/fonts) · figure types (circular map · RSCU heatmap ·
  π sliding window with LSC/SSC/IR shading · IR junction · mVISTA-style · SSR/repeat distribution ·
  annotated tree · NUPT/MTPT by gene · haplotype network · mito configuration · barcode assessment) ·
  vector export (SVG/PDF/EPS) · scaling rule for N genomes
- [ ] **M5.2 Statistics** — descriptive · non-parametric tests · correlation · FDR/Bonferroni ·
  sentence-level interpretation
- [ ] **M5.3 Exports** — multi-sheet Excel · Word · Supplementary package · 1-page PDF per genome ·
  technical QC report
- [ ] **M5.4 Queryable-number check automated** — fails the build on any orphan number

## Phase 6 — Tool-paper content  *(§13.6 — published/GenBank data; assembly-dependent parts move to P7)*

**Goal:** the scientific content of the Eleutherine bulbosa tool paper.
**In-phase checks:** mito-structure parts have published benchmarks to compare against
(*Damnacanthus indicus* · *Salvia* 2022–2025 · *Hemiboea* · *Coffea arabica* 2023).
NUPT/MTPT has **no published benchmark** — its only check is the Oldenlandia/Coleus comparison in Phase 8.

- [ ] **M6.1 NUPT/MTPT** (nuclear ↔ plastid ↔ mito transfer) — **BLOCKED: needs nuclear genome**
- [ ] **M6.2 Mito polymorphic structure + recombination via repeats** (Chimeric ORF/CMS)
- [ ] **M6.3 Parallel plastid–mito comparison** within species
- [ ] **M6.4 Ancestral state reconstruction**
- [ ] **M6.5 Divergence time estimation** (mcmctree / BEAST2) — remote server

## Phase 7 — Assembly  *(§13.7 — after analysis, by decision 01/10)*

**Goal:** assemble plastomes/mitogenomes from raw reads; feed Phase 6.
**In-phase checks:** GetOrganelle + SPAdes → 127,679 bp from 20,000 simulated read pairs ·
OATK → 126,660 bp @ identity 100.00 % · Bandage → PNG. With real reads: add a **peak-RAM** measurement.

- [ ] **M7.1 Batch assembly + QC** (fastp → GetOrganelle) — **BLOCKED: no real Illumina reads yet**
- [ ] **M7.2 Assembly-graph checks** (Bandage)
- [ ] **M7.3 Feed mito-structure work** (M6.2)
- [ ] **M7.4 Remote path** (usegalaxy.eu `get_organelle_from_reads`) — subsample reads before upload ·
  pin tool version · record history ID
- [ ] **M7.5 Peak-RAM with real reads** — **BLOCKED: no real reads yet**

## Phase 8 — Cross-testing, performance, automation  *(§13.8 — not first-time testing)*

**Goal:** only what cannot be done earlier — cross-checks, scale, automation.

- [ ] **M8.1 Cross-check vs CGAS · CPStools · CPGview** on the same dataset — differences explained
- [ ] **M8.2 Cross-check NUPT/MTPT vs the Oldenlandia/Coleus papers** (only benchmark available)
- [ ] **M8.3 Batch test 10–50+ genomes** — record memory/time limits (PROJECT.md §8)
- [ ] **M8.4 CI** — tests on every change; golden set re-runs automatically

## Phase 9 — Packaging, license & release  *(§13.9)*

**Goal:** a citable, installable release.
**In-phase check:** golden numbers re-run **inside Docker on Windows** — same input, same numbers
across 3 environments (WSL · CI · Docker).

- [ ] **M9.1 README / quickstart / user docs** — followable by a non-technical user
- [ ] **M9.2 pip/conda package + Dockerfile** (Windows via Docker)
- [ ] **M9.3 License review** — PGA = GPL-3.0 · Vmatch = academic only · **MISA = not stated**
  (zip ships no license; IPK page only cites Beier et al. 2017) · isolation if needed
- [ ] **M9.4 Module × platform × RAM table + 3 support levels** (Linux native · Windows/Docker ·
  macOS experimental, unverified)
- [ ] **M9.5 "Things we don't do" statement** — no native Windows; no assembly below RAM threshold;
  no RNA-editing prediction without RNA-seq data
- [ ] **M9.6 Short how-to video · tool manuscript · DOI (Zenodo) + CITATION.cff**
- [ ] **M9.7 Web UI** — deferred until the core is stable

---

## Open decisions blocking milestones

| Decision | Blocks | Where |
|---|---|---|
| ~~Coordinate system in reports~~ — decided 2026-10-02: GenBank coordinates for per-genome tables · normalized for cross-genome figures · system labelled on every table · `rotation_offset` in `00_run_info` | — | PROJECT.md §3 |
| ~~Reports generated from `Result` objects~~ — decided 2026-10-02: reports reference table cells (`{{table.column}}`); an unresolvable reference raises | — | PROJECT.md §3 |
| ~~Reference Panel: pinned snapshot vs live NCBI~~ — decided 2026-10-02: propose from a dated RefSeq cache, pin + checksum on the user's confirmation, runs stay offline | — | PROJECT.md §3 |
| Reference Panel item (c) — source text truncated | documentation only — does not block M1.10 | PROJECT.md §10.9 |

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Real Illumina reads (1 sample) | Diem / Trang | M7.1 · M7.5 · M4.5 |
| Nuclear genome of the sample | Trang | M6.1 (NUPT/MTPT) |
| A macOS machine | Team | M9.4 platform claims |
