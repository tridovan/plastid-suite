# ROADMAP.md — the whole project

> The macro view. Rewritten when a milestone closes or a phase opens or closes — never appended to.
> Closed work is one line: the detail is in git (commit messages), the evidence in `docs/log/` —
> `PROJECT.md` §8 indexes it by `§8.N`. Legend: `[x]` done · `[~]` partly · `[ ]` open · **BLOCKED** = waits on data or people.

## Now — 17/10/2026

- **Closed and on `main`:** Phases 1–5, 4a, 3B, 2B, 8 and 8B — every analysis module, the in-house
  annotator, the cross-checks against CGAS / CPStools and the Crataegus answer key, CI. The analysis
  round is done, except RNA editing (needs RNA-seq) and the parked mito/nuclear round.
- **Open:** Phase 9 on `feat/phase-9`. The packaging, the licence review (MIT), the Docker image and the
  user docs are done; the 1.0.0 candidate is with the team (P9-08), and the DOI, PyPI and the video wait
  on them.
- **Next, unblocked:** fixes to the figure tool from anh Điềm's package (`feedback/`, untracked) — then
  close Phase 9 on the team's report.
- **After that:** Phase 7, assembly — blocked on real Illumina/ONT reads.

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

**Phase 7 — assembly** (waits on real reads; so far run only on simulated pairs)

Check: GetOrganelle + SPAdes → 127,679 bp from 20,000 simulated read pairs · OATK → 126,660 bp at
100.00 % identity · Bandage → PNG; with real reads, add a peak-RAM measurement. Assembly runs locally
when the user's machine has the RAM, otherwise on the Galaxy remote path (confirmed 08/10/2026).

- [ ] **M7.1 Batch assembly + QC** (fastp → GetOrganelle · BWA/minimap2 map-back) — **BLOCKED: no real
  Illumina reads.** Neither `fastp` nor `BWA` is declared in the repo yet: declare each in the session
  it is first used.
- [ ] **M7.2 Assembly-graph checks** (Bandage) — **BLOCKED: no real reads**
- [ ] **M7.4 Remote path** (usegalaxy.eu `get_organelle_from_reads`) — subsample reads before upload ·
  pin the tool version · record the history ID · data-sharing policy for unpublished reads is **open**
  (PROJECT.md §10.12)
- [ ] **M7.5 Peak-RAM with real reads** — **BLOCKED: no real reads**
- [ ] **M7.6 Long-read path** (OATK · raven-assembler, ONT/HiFi) — **BLOCKED: no ONT/HiFi reads**; both
  tools are declared in `environment.yml` and were run in the survey (§7.1)

**Earlier phases, still open**

- [~] **M3B.0 Naming (D5-09)** — built and measured (`gene_normalize`; Ts/Tv 2.6917 unchanged, §8.45).
  **Open: the ★ pick** — the end user chooses the project's naming record from the audit report
  (§8.43); a recorded override stands until then.
- [ ] **M4.5 RNA editing prediction** — **BLOCKED: needs RNA-seq data**
- [~] **M8B.4 RSCU vs Table S5** — measured, not reproduced: 1–2 % away under every one of 11 counting
  rules tested; the sequences are proven identical (Table S7 reproduces exactly), so the residual is the
  authors' counting rule (§8.70) — a question for anh Điềm.

**Parked by scope — never picked up this round** (D5-07; the mito/nuclear round, with its follow-ups H-06 · H-10 · H-12)

- [ ] **M6.1 NUPT/MTPT** (nuclear ↔ plastid ↔ mito transfer) — needs a nuclear genome when it reopens
- [ ] **M6.2 Mito polymorphic structure + recombination via repeats** (chimeric ORF / CMS)
- [ ] **M6.3 Parallel plastid–mito comparison** within a species
- [ ] **M7.3 Feed mito-structure work** (M6.2)
- [ ] **M8.2 Cross-check NUPT/MTPT against the Oldenlandia/Coleus papers** — checks a module (M6.1) this
  round does not build, so it is parked, not waiting

## Known defects — named, not fixed

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

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| **P9-08** — run the 1.0.0 candidate (`docs/handoff_test_1.0.0.md`, levels 1–3) and read a `report.md` | anh Điềm (levels 2–3) · cô Trang (the report) | Phase 9's close. One machine's green is not a release's green; anything they find becomes a regression test before it is fixed |
| Authors to name in `CITATION.cff` · a Zenodo account · the repository's visibility (private today; Zenodo archives public repositories) · the PyPI name `plastid-suite` · the how-to video | the team · the dev | M9.6 — and a publish needs the dev's explicit go, each time |
| The ★ naming record, picked from the audit report | the end user (the PI signs off the names a manuscript prints) | M3B.0 |
| Two decisions that are the PI's: which `trnG` name prints in tables and submissions (H-04 — the golden calls two different tRNAs `trnG-UCC`; PROJECT.md §10.10) · what a machine-written sentence may assert (D5-05; PROJECT.md §10.11, the cautious default is in force) | cô Trang | manuscript tables, submission files, report wording |
| Four open questions from the Crataegus corpus (§8.69 … §8.73) | anh Điềm | M8B.4 — none of them blocking |
| Written permission from IPK Gatersleben to redistribute `misa.pl` | cô Trang / anh Điềm (one mail) | nothing — the release does not wait on it |
| Real Illumina reads (1 sample) | anh Điềm / cô Trang | M7.1 · M7.5 · M4.5 |
| ONT/HiFi reads (1 sample) | anh Điềm / cô Trang | M7.6 |
| The team's own `.sbt` submitter template (NCBI's form; a stand-in is vendored for tests) | the team | M3.4 · M5.3 — a real `.sqn` submission only |
| A macOS machine | the team | macOS stays *unverified* in `docs/platforms.md` |

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
| reads · `fastp` · `BWA`/`minimap2` (trim · map back) | M7.1 | **blocked — no real reads yet**; neither tool declared yet |
| reads · `Bandage` (assembly graph) | M7.2 | **blocked — no real reads yet** |
| 1 Assembly — short read (`GetOrganelle` + `SPAdes`) | M7.1 | **blocked — no real Illumina reads yet** |
| 1 Assembly — long read (`OATK` · `raven`) | **M7.6** | **blocked — no ONT/HiFi reads**; declared · run in the survey (§7.1) |
| 2 Annotation (self-built · BLAST+) | M3B.2–M3B.5 | **done 10/10/2026** — P3B-07 … P3B-10, §8.46–§8.55; close = P3B-11 · P3B-12 |
| 3 Verify names (`gene_normalize`, verify mode) | **M3B.0** | **done** — P3B-05, §8.45 |
| [4] Analyses — structure · gene content · codon · SSR (`MISA`) · **repeats** (`vmatch`) · synteny (`MUMmer4`) · Ka/Ks (`PAML`) · Ts/Tv (`MAFFT`) · barcode · popgen · pycirclize | M2.3 · M2.1/M2.2 · M2.4 · M2.8 · **M2.12** · M4.2 · M4.1 · M2.6/M4.6 · M4.3 · M4.4 · M5.1 | **done** — M2.9 · M2.10 · M2.12 closed in Phase 2B (§8.57–§8.59) |
| [5] RNA editing (a reference of its own) | M4.5 | **blocked — needs RNA-seq** |
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
