# ROADMAP.md — the whole project

> The macro view. Rewritten when a milestone closes or a phase opens or closes — never appended to.
> Closed work is one line: the detail is in git (commit messages), the evidence in `docs/log/` —
> `PROJECT.md` §8 indexes it by `§8.N`. Legend: `[x]` done · `[~]` partly · `[ ]` open.

## Now — 17/10/2026

- **Closed and on `main`:** Phases 1–5, 4a, 3B, 2B, 8, 8B and **9** — every analysis module, the
  in-house annotator, the cross-checks against CGAS / CPStools and the Crataegus answer key, CI, and
  the release itself: MIT, the wheel, the Docker image and the user docs. The analysis round is done,
  except RNA editing (needs RNA-seq) and the parked mito/nuclear round.
- **Open:** **Phase 10** on `feat/phase-10` — the paper's figures. anh Điềm's bundle of R scripts and
  patched CPJSdraw drew the Iridaceae manuscript's Figures 2–9. **M10.1 is done**: the junction
  figure is a command, and it reproduces the bundle's own output byte for byte. `TASKS.md` holds the
  next task — the R engine.
- **Waiting on people:** anh Điềm's P9-08 report — the code is on `main` now, so that is the pull —
  the release accounts, and a macOS machine.
- **After that:** Phase 7, assembly — judged on a published paper's public reads against its assembled
  plastome; it waits on no one.

## Scope of this round

**Plastid only** (D5-07): one organelle per run, never a mixed input set. The mito/nuclear milestones are
parked by name (below) and never picked up; the mito code that exists stays as it is and stays tested.
Divergence time and ancestral state are the PI's analyses, not modules. The tool paper (NUPT/MTPT, mito
structure) waits with the parked round.

## Open milestones

**Phase 9 — packaging, licence, release** (`feat/phase-9`; M9.1 – M9.5 are done, below)

- [~] **M9.6 Release artifacts** — `CITATION.cff` written and validated (CFF 1.2.0) and `docs/release.md`
  is the procedure. Left: the DOI (Zenodo), the PyPI upload, the how-to video and the authors to name.
  Each is outward-facing and permanent, so each waits for the dev's explicit go.
- [ ] **M9.7 Web UI** — deferred until the core is stable

**Phase 10 — the paper's figures** (`feat/phase-10`; the manuscript's Figures 2–9, drawn outside the
toolkit until now — seven R scripts and a patched CPJSdraw, in `feedback/`)

- [x] **M10.1 The IR-junction figure (Figure 2)** — `core/junctions.py` reads the record's own four
  regions in the GenBank frame, `external/cpjsdraw.py` draws them, `plastid ir-junction` drives both
  into a run directory. Measured 17/10/2026 on the bundle's eight Iridaceae records: the
  configuration is reproduced **8/8**, the SVG is **byte-identical** to the bundle's own verified
  output, and the table's eight region and junction columns agree on every row.
- [x] **M10.2 The R figure engine** — `external/rscript.py`, the seven scripts under
  `data/r_scripts/`, `docs/tools/r.md` and the `environment.yml` pins. Measured 17/10/2026: all
  seven run through the adapter in 35 s and each writes its PNG, on 4.3.3 here and 4.1.2 for the
  scripts' author. The one change to the scripts is in `00_ham_chung.R`: input paths are made
  absolute before `setwd()`, without which a relative `--inputdir` finds nothing
- [ ] **M10.3 The figures PTS lacks** — amino-acid heatmap (paper 4), substitution types (7), Ka/Ks (8)
- [ ] **M10.4 The pivots** — RSCU (3), SSR (5) and π (6) need the long→wide shape their R scripts read
- [ ] **M10.5 The tree (9 + S1–S3)** — stays the user's script over the user's own tree; PTS builds none

**Phase 7 — assembly** (not built; judged on a published paper — its public reads in, its assembled plastome
as the answer; so far run only on simulated pairs)

Check: GetOrganelle + SPAdes → 127,679 bp from 20,000 simulated read pairs · OATK → 126,660 bp at
100.00 % identity · Bandage → PNG; with real reads, add a peak-RAM measurement. Assembly runs locally
when the user's machine has the RAM, otherwise on the Galaxy remote path (confirmed 08/10/2026).

- [ ] **M7.1 Batch assembly + QC** (fastp → GetOrganelle · BWA/minimap2 map-back) — needs a published
  paper with public Illumina reads and its assembled plastome as the answer key. Neither `fastp` nor `BWA`
  is declared in the repo yet: declare each in the session it is first used.
- [ ] **M7.2 Assembly-graph checks** (Bandage) — after M7.1
- [ ] **M7.4 Remote path** (usegalaxy.eu `get_organelle_from_reads`) — subsample reads before upload ·
  pin the tool version · record the history ID · data-sharing policy for unpublished reads is **open**
  (PROJECT.md §10.12)
- [ ] **M7.5 Peak-RAM with real reads** — measured on the M7.1 run
- [ ] **M7.6 Long-read path** (OATK · raven-assembler, ONT/HiFi) — needs a published paper with public
  long reads and its assembly; both tools are declared in `environment.yml` and were run in the survey
  (§7.1)

**Earlier phases, still open**

- [~] **M3B.0 Naming (D5-09)** — built and measured (`gene_normalize`; Ts/Tv 2.6917 unchanged, §8.45).
  **Open: the ★ pick** — the end user chooses the project's naming record from the audit report
  (§8.43); a recorded override stands until then.
- [ ] **M4.5 RNA editing prediction** — needs a published paper with public RNA-seq and its list of
  editing sites
- [~] **M8B.4 RSCU vs Table S5** — measured, not reproduced: 1–2 % away under every one of 11 counting
  rules tested; the sequences are proven identical (Table S7 reproduces exactly), so the residual is the
  authors' counting rule (§8.70) — unexplained, reported as measured.

**Parked by scope — never picked up this round** (D5-07; the mito/nuclear round, with its follow-ups H-06 · H-10 · H-12)

- [ ] **M6.1 NUPT/MTPT** (nuclear ↔ plastid ↔ mito transfer) — needs a nuclear genome when it reopens
- [ ] **M6.2 Mito polymorphic structure + recombination via repeats** (chimeric ORF / CMS)
- [ ] **M6.3 Parallel plastid–mito comparison** within a species
- [ ] **M7.3 Feed mito-structure work** (M6.2)
- [ ] **M8.2 Cross-check NUPT/MTPT against the Oldenlandia/Coleus papers** — checks a module (M6.1) this
  round does not build, so it is parked, not waiting

## Known defects — named, not fixed

- **H-25** An IR pair whose two annotated copies differ in length by *d* measures as unrelated.
  `measure_ir` takes the mirror constant from the two annotated boundaries (`first.start +
  second.end - 1`), so a copy *d* bases longer shifts every base *d* positions away from its true
  partner: measured 17/10/2026 on eight Iridaceae records, *d* = 1 (Chasmanthe NC_072562), 2 (Iris
  hermona NC_081931) and 4 (Trimezia NC_072662) give identities of **0.2924 / 0.2636 / 0.2778** —
  all under the 99.9 % gate — while the constant that accounts for *d* gives **1.0000 / 0.9985 /
  0.9998**. The gate therefore refuses three genomes whose inverted repeats are exact, and
  `plastid` cannot process them at all. Records whose difference exceeds
  `IR_COPY_LENGTH_TOLERANCE` (Freesia 37, Olsynium 40) fail the annotation path and are rescued by
  de-novo detection; those inside the tolerance are not. Note the interaction: a copy written
  `join(...)` across the origin arrives as two intervals and the small piece is dropped, which is
  one way *d* becomes non-zero. A task of its own with a regression test, never a silent fix.
- **H-24** `trnV-UAC` is never called: ARAGORN reads that locus as `tRNA-Stop(tta)` in the Crataegus
  corpus and as `tRNA-Cys(aca)` on the golden, so an intron-containing tRNA is missing from every
  `plastid annotate` run. Three fixes are possible — transfer tRNAs from the reference, as CDS and rRNA
  already are · a second engine · report and keep — and none is chosen (§8.71;
  `docs/corpus_crataegus_2026.md`).
- **H-23** `introns` publishes CDS introns only, so the golden's 8 tRNA introns appear in no module's
  output. Extending it touches `core/spans.py` and the π-per-intron column of `diversity`: a task, not a
  patch (§8.62; the limit is stated in `docs/platforms.md`).
- **H-22** ARAGORN's `-i` searches introns up to 3,000 bp and says nothing when it skips a longer one; the
  longest measured is 2,671 bp, a 12 % margin (§8.52, §8.55).
- **H-21** The audit's partial-coordinate BLOCK refuses `NC_024936.1`, whose `rps19` is complete and which
  produced the published call — the rule's stated reason does not hold for this shape. Changing a BLOCK
  rule decided in P3B-03 means re-measuring all of it (§8.49).

## Waiting on people

| What | Who | Blocks |
|---|---|---|
| **P9-08** — pull the code from git (branch `feat/phase-9` until Phase 9 is merged into `main`) and test it by the process that suits the tool: input in, output out | anh Điềm | Phase 9's close. One machine's green is not a release's green; anything found becomes a regression test before it is fixed |
| Authors to name in `CITATION.cff` · a Zenodo account · the repository's visibility (private today; Zenodo archives public repositories) · the PyPI name `plastid-suite` · the how-to video | the team · the dev | M9.6 — and a publish needs the dev's explicit go, each time |
| A macOS machine | the team | macOS stays *unverified* in `docs/platforms.md` |

Nothing else waits on people: the tools are judged on a published paper's public input against the
paper's own results (CLAUDE.md → How the tools are judged). The open design questions for the PI and the end
user — the naming ★, which `trnG` name prints, what a generated sentence may assert — are in PROJECT.md §3
and §10, and none of them blocks evaluation.

## The pipeline — every step, and who owns it

*The diagram below is a **copy** of the canonical one in `PROJECT.md` §13.0 — edit it there
first, then copy it back (the rule `docs/D5-09_reference_naming.md` §1 also follows). Tools are
named per step; a step whose tool is not declared in the repo yet says so in its milestone.*

```
   NGUỒN THAM CHIẾU
   NCBI RefSeq  ·  bản ghi nhóm tự nộp
                │
                │  plastid panel pin        ← accession · ngày · sha256
                ▼
   ┌──────────────────────────────────────┐
   │  DB THAM CHIẾU      (nhiều file .gb) │
   └──────────────────┬───────────────────┘
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐÁNH GIÁ CHẤT LƯỢNG                 │
   │   · file nào KHÁC so với phần còn lại│
   │   · file nào không dùng được (BLOCK) │
   └──────────────────┬───────────────────┘
                      │  báo cáo
                      ▼
   ┌──────────────────────────────────────┐
   │  NGƯỜI DÙNG ĐỌC VÀ CHỌN              │
   │  1 file .gb tốt nhất                 │
   └──────────────────┬───────────────────┘
                      ▼
                ★ CHUẨN TÊN    →    từ vựng A
                      │
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐỒNG NHẤT TÊN CHO CẢ DB THAM CHIẾU  │   gene_normalize · D5-09
   └──────────────────┬───────────────────┘
                      ▼
        DB THAM CHIẾU ĐÃ CHUẨN HOÁ   (A)
                      │
                      │              reads (Illumina · ONT/HiFi)
                      │                │   fastp · BWA/minimap2 · Bandage
                      │                ▼
                      │       ┌─────────────────────────┐
                      │       │ 1  LẮP RÁP              │  GetOrganelle + SPAdes (short read)
                      │       │                         │  OATK · raven (long read)
                      │       └────────────┬────────────┘
                      │                    │  FASTA
                      │                    ▼
                      └───────────►┌─────────────────────────┐
                                   │ 2  CHÚ GIẢI             │  công cụ tự viết (D5-08)
                                   │                         │  blastx · blastn  [BLAST+]
                                   └────────────┬────────────┘
                                                │  GenBank — tên theo A
                                                ▼
                                   ┌─────────────────────────┐
                                   │ 3  KIỂM TRA TÊN         │  gene_normalize (chế độ kiểm)
                                   └────────────┬────────────┘
                                                ▼
                        ╔═══════════════════════════════════════╗
                        ║   BỘ INPUT CHUẨN HOÁ   (từ vựng A)    ║
                        ╚═══════════════════┬═══════════════════╝
                                            │
                             ┌──────────────┴───────────────┐
                             ▼                              ▼
                    [4]  PHÂN TÍCH                [5]  RNA EDITING (C.3)
                    cấu trúc · gen · codon         ◄── reference RIÊNG (B)
                    SSR (MISA) · repeat (vmatch)        ghép cặp có bằng chứng
                    synteny (MUMmer4)                   → alignment → chưa chọn
                    Ka/Ks (PAML) · Ts/Tv (MAFFT)
                    barcode · popgen · pycirclize
                             │
                             ▼
                    [6]  XUẤT BẢN
                    table2asn (.sbt) · Excel · Supplementary
```

**Coverage** — every box, and every tool inside a box, names its owner (the rule this section
exists for; check it whenever the diagram changes):

| Box in the diagram | Owner | State |
|---|---|---|
| Reference supply · `plastid panel pin` | M1.10 | done 03/10 |
| Quality assessment (the audit) | M3B.1 | **done** — P3B-01 · P3B-03 · `--audit-only` (P3B-09), §8.40 · §8.43 · §8.49 |
| The end user picks the naming reference ★ | **M3B.0** | **waiting on the end user** — the audit report is the evidence (§8.43) |
| Unify the reference DB (`gene_normalize`) | **M3B.0** | **done** — P3B-04, §8.44 (the ★ override stands until it is picked) |
| reads · `fastp` · `BWA`/`minimap2` (trim · map back) | M7.1 | **not built** — judged on a published paper's public reads; neither tool declared yet |
| reads · `Bandage` (assembly graph) | M7.2 | **not built** — after M7.1 |
| 1 Assembly — short read (`GetOrganelle` + `SPAdes`) | M7.1 | **not built** — judged on a published paper's public Illumina reads |
| 1 Assembly — long read (`OATK` · `raven`) | **M7.6** | **not built** — judged on a published paper's public long reads; declared · run in the survey (§7.1) |
| 2 Annotation (self-built · BLAST+) | M3B.2–M3B.5 | **done 10/10/2026** — P3B-07 … P3B-10, §8.46–§8.55; close = P3B-11 · P3B-12 |
| 3 Verify names (`gene_normalize`, verify mode) | **M3B.0** | **done** — P3B-05, §8.45 |
| [4] Analyses — structure · gene content · codon · SSR (`MISA`) · **repeats** (`vmatch`) · synteny (`MUMmer4`) · Ka/Ks (`PAML`) · Ts/Tv (`MAFFT`) · barcode · popgen · pycirclize | M2.3 · M2.1/M2.2 · M2.4 · M2.8 · **M2.12** · M4.2 · M4.1 · M2.6/M4.6 · M4.3 · M4.4 · M5.1 | **done** — M2.9 · M2.10 · M2.12 closed in Phase 2B (§8.57–§8.59) |
| [5] RNA editing (a reference of its own) | M4.5 | **not built** — judged on a published paper's public RNA-seq and its editing sites |
| [6] Submission (`table2asn` · Excel · Supplementary) | M3.4 · M5.3 | done — except the team's own `.sbt` (see Waiting) |

## Done — one line each

### Phase 1 — Framework & shared data
- [x] **M1.1 Skeleton & packaging**
- [x] **M1.2 Core model + GenBank reader + vendored golden data**
- [x] **M1.3 IR detection + gate + normalization**
- [x] **M1.4 Result schema + exporters**
- [x] **M1.5 Run layer**
- [x] **M1.6 Golden harness**
- [x] **M1.7 Report engine + citation registry + defaults**
- [x] **M1.8 Module registry + first module end-to-end**
- [x] **M1.9 `ssr` module (tier B)** — MISA is obtained, no longer vendored (M9.3)
- [x] **M1.10 Reference Panel (4 tiers) + NCBI fetch/cache**
- [x] **M1.11 Sample-group metadata**
- [x] **M1.12 Execution backend decision**
- [x] **M1.13 Hardening + API freeze v0.1**

### Phase 2 — Genome reading & statistics
- [x] **M2.1 Gene content & summary statistics**
- [x] **M2.2 Publication gene table**
- [x] **M2.3 Genome structure**
- [x] **M2.4 Codon usage**
- [x] **M2.5 Amino acid composition**
- [x] **M2.6 SNP / substitutions by position**
- [x] **M2.7 Introns**
- [x] **M2.8 SSR**
- [x] **M2.9 Nucleotide diversity π**
- [x] **M2.10 Manuscript deliverable** — the recipe is `docs/P2-10_manuscript_draft.md`
- [x] **M2.11 Sample-group labels**
- [x] **M2.12 Repeats (dispersed + palindromic)**

### Phase 3 — Annotation
- [x] **M3.1 Annotation** — plastid: superseded by Phase 3B (Chloe and PGA retired); mito: Phase 3A
- [x] **M3.2 Gene-name check & normalization**
- [x] **M3.3 ORF checks**
- [x] **M3.4 NCBI submission files**

### Phase 3A — Mito annotator, self-built prototype
- [x] **M3A.1 GenBank writer**
- [x] **M3A.2 Reference set**
- [x] **M3A.3 Gene caller**
- [x] **M3A.4 tRNA caller**
- [x] **M3A.5 CLI**
- [x] **M3A.6 Golden check**
- [x] **M3A.7 AP045469.1**

### Phase 3B — Plastid annotator, self-built
- [x] **M3B.1 Reference-set audit**
- [x] **M3B.2 Calling engine**
- [x] **M3B.3 Arbitration & export**
- [x] **M3B.4 CLI & batch**
- [x] **M3B.5 Golden + corpus validation**

### Phase 4 — Multi-genome comparison
- [x] **M4.1 Ka/Ks selection pressure**
- [x] **M4.2 mVISTA-style %identity + synteny**
- [x] **M4.3 Barcode candidates & hypervariable regions**
- [x] **M4.4 Population genetics**
- [x] **M4.6 Ts/Tv module**

### Phase 4a — The Điềm corpus as the test dataset
- [x] **M4a.0 Corpus in the working tree**
- [x] **M4a.1 Pinned**
- [x] **M4a.2 Vendored fixtures + Phase-1 checks**
- [x] **M4a.3 Scale (M8.3)**
- [x] **M4a.4 Annotator at scale**

### Phase 5 — Figures & reports
- [x] **M5.1 Figure system**
- [x] **M5.2 Statistics**
- [x] **M5.3 Exports**
- [x] **M5.4 Queryable-number check automated**

### Phase 8 — Cross-testing, performance, automation
- [x] **M8.1 Cross-check vs CGAS · CPStools · CPGview**
- [x] **M8.3 Batch test 10–50+ genomes**
- [x] **M8.4 CI**

### Phase 8B — The Crataegus answer-key corpus
- [x] **M8B.1 The corpus on disk, verified** — provenance, and how to prove a copy is the same: `docs/corpus_crataegus_2026.md`
- [x] **M8B.2 Structure & GC vs Table 1**
- [x] **M8B.3 SSR vs Table S6**
- [x] **M8B.5 Ka/Ks vs Table S7**
- [x] **M8B.6 π and hotspots vs the paper**
- [x] **M8B.7 Annotation vs the published GeSeq annotation**
- [x] **M8B.8 The defects the corpus forced**

### Phase 9 — Packaging, license & release
- [x] **M9.1 README / quickstart / user docs**
- [x] **M9.2 pip/conda package + Dockerfile (Windows via Docker)**
- [x] **M9.3 License review**
- [x] **M9.4 Module × platform × RAM table + 3 support levels**
- [x] **M9.5 "Things we don't do" statement**

*Open questions: PROJECT.md §10. Decisions: PROJECT.md §3. Original per-phase scope text: PROJECT.md §13.*
