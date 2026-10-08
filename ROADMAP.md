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
5. **Coverage**: every box of the pipeline diagram (canonical copy: PROJECT.md §13.0 — repeated in
   "The pipeline" below) — **and every tool named inside a box** — names the milestone that owns
   it, or is marked *surveyed — not adopted* with the measurement that dropped it. A box or tool
   with no owner is a **gap**: name it in "Open decisions", never leave it silent. The check runs
   at tool level because the box-level pass of 07/10/2026 missed three tools and one team artifact.

## Cross-phase rules (apply to every phase)

- Serves **both plastid and mitochondrion** — but the seven *transitional* analyses (see PROJECT.md
  §13, "Hai bào quan") run on plastid only today: on a mitogenome they stop at the IR gate. Opening
  them is **H-10, on hold** by D5-07. Divergence-time work is tool-paper only.
- **Plastid first** (team decision 05/10/2026, D5-07): plant mito data is still scarce, so new module
  work serves the plastome. The existing mito path stays as it is — not extended — and the
  two-organelle structure is kept so mito analyses can be integrated later without rework
  (PROJECT.md §13 → "Hai bào quan"; H-06 · H-10 on hold by D5-07).
- **IR gate ≥ 99.9 %** before any position-based analysis. Fail ⇒ red flag, stop. No exceptions.
- Gene counting: **132 = 86 CDS + 38 tRNA + 8 rRNA**. Never "133 genes" (that is the `gene` feature
  count, diagnostic only).
- Ts/Tv scope: **CDS only · 79 shared genes · MAFFT · IR excluded** → 2.6917. Other scopes must be
  labelled as such wherever they are exported. *(Was 2.686 over 76 genes; the scope definition is
  unchanged — the shared set was re-measured by locus on 04/10/2026 after three differently-spelled
  loci turned out to be shared, §8.27.)*
- Every run emits a **tool + parameters report**; output is **run-scoped** (`<dataset>/<run_id>/`),
  never overwritten in place.
- Genome normalization (origin/strand → LSC start, 4 regions split) happens once in Phase 1 and is
  reused everywhere after.
- Runs on: **ordinary machines** = phases 1–5, 8 · **remote server** = phase 6 (divergence time).
  Assembly (phase 7) runs **locally when the user's machine has the RAM**, otherwise the Galaxy
  remote path — confirmed 08/10/2026; data-sharing policy for unpublished reads: **open**
  (PROJECT.md §10.12).
- Standard test case: **Eleutherine bulbosa dataset + the golden numbers** (PROJECT.md §6).

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
| Quality assessment (the audit) | M3B.1 | planned — Phase 3B |
| The end user picks the naming reference ★ | **M3B.0** | planned — Phase 3B |
| Unify the reference DB (`gene_normalize`) | **M3B.0** | planned — Phase 3B |
| reads · `fastp` · `BWA`/`minimap2` (trim · map back) | M7.1 | **blocked — no real reads yet**; neither tool declared yet |
| reads · `Bandage` (assembly graph) | M7.2 | **blocked — no real reads yet** |
| 1 Assembly — short read (`GetOrganelle` + `SPAdes`) | M7.1 | **blocked — no real Illumina reads yet** |
| 1 Assembly — long read (`OATK` · `raven`) | **M7.6** | **blocked — no ONT/HiFi reads**; declared · run in the survey (§7.1) |
| 2 Annotation (self-built · BLAST+) | M3B.2–M3B.5 | planned — Phase 3B |
| 3 Verify names (`gene_normalize`, verify mode) | **M3B.0** | planned — Phase 3B |
| [4] Analyses — structure · gene content · codon · SSR (`MISA`) · **repeats** (`vmatch`) · synteny (`MUMmer4`) · Ka/Ks (`PAML`) · Ts/Tv (`MAFFT`) · barcode · popgen · pycirclize | M2.3 · M2.1/M2.2 · M2.4 · M2.8 · **M2.12** · M4.2 · M4.1 · M2.6/M4.6 · M4.3 · M4.4 · M5.1 | done, except **M2.12** (open) and M2.9 · M2.10 `[~]` |
| [5] RNA editing (a reference of its own) | M4.5 | **blocked — needs RNA-seq** |
| [6] Submission (`table2asn` · Excel · Supplementary) | M3.4 · M5.3 | done — except the team's `.sbt`, a TASKS.md waiting row |

*Checked 07/10/2026: two boxes had no owner — the naming step (now **M3B.0**) and **repeats**
(now **M2.12**). Making that check a rule ("How to use" 5) is what keeps this table true.*

*Checked again 07/10/2026, now **at tool level**: three tools had no owner — the long-read
assemblers (`OATK` · `raven` → **M7.6**), the read map-back (`BWA`/`minimap2` → named in **M7.1**),
and `Sibelia` (surveyed and run in the 26/9 survey, never adopted: the project's synteny route is
`nucmer`, D4-02/M4.2 — so it left the diagram); plus one team artifact with no row at all (the
submitter template `.sbt` → `TASKS.md` → "Waiting on data / people"). `Flye` also left the diagram
— measured to die on RAM here and replaced by the two long-read assemblers (§7.1). PROJECT.md
§13.0's dated note is the canonical record.*

## Already banked (evidence in PROJECT.md §8–§9)

- [x] Golden numbers measured & cross-checked against 2 published papers (153,724 bp · 37.85 % ·
      132 = 86+38+8 · IR 26,214 bp @ 100.00 % · LSC 83,127 / SSC 18,169 · 37 SSR ·
      Ts/Tv 2.6917 — re-measured 04/10 over the 79 shared loci, §8.27)
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
- [x] **M1.7 Report engine + citation registry + defaults** (2026-10-02) — `{{…}}` resolves through
  `Result.format()` or raises; shipped templates carry zero bare digits (scan test); the draft
  banner is appended by the engine; citations come from `src/plastid/data/citations.yaml`
- [x] **M1.8 Module registry + first module end-to-end** (2026-10-02) — decorator registry;
  `genome_stats` produces `01_genome_stats/` (3 tables + result.json + report.md); CLI `plastid run`
  writes `manifest.json` + `00_run_info/run_info.md`; a failed gate exits 2 and still leaves a
  traceable run directory
- [x] **M1.9 MISA vendored + `ssr` module (tier B)** (2026-10-02) — IPK zip verified against every
  §7.8 checksum and vendored with `SOURCE.json` (licence: not stated → flagged for M9.3); the module
  reproduces **SSR 37 = mono 31 · di 4 · tri 2, compound 3** on this machine; a missing perl raises a
  named `ToolError` → CLI exit 3
- [x] **M1.10 Reference Panel (4 tiers) + NCBI fetch/cache** (2026-10-03) — `core/panel.py`
  (tiers · modes A/B · quality filter · pin with an injected fetcher) + `run/ncbi.py` (esearch ·
  esummary · taxonomy efetch · 3 req/s) + `plastid panel refresh | propose | pin`. Verified live:
  76 RefSeq genomes for Iridaceae (the §8.3 count, after fixing NCBI's `chloroplast[filter]` blind
  spot — PROJECT.md §9.C-27), and a re-pinned NC_050833.1 whose sha256 matches the vendored golden
  file byte-for-byte. Proposals and pins are offline; the network lives only in `refresh`.
- [x] **M1.11 Sample-group metadata** (2026-10-02) — loader + lookup shipped (`core/samples.py`,
  duplicate-accession guard, tests); the "used by a comparison table" half moved to **M2.11**
- [x] **M1.12 Execution backend decision** (2026-10-03, dev decision) — local for phases 1–5 and 8;
  remote for phase 6 (divergence time) and phase 7 (assembly), where PROJECT.md §9.28 measured that
  an 8 GB machine leaves only 1.5–2.5 GB free per job. Recorded in PROJECT.md §3.
  *(Updated 08/10/2026: assembly runs locally when the user's machine has the RAM — the Galaxy
  remote path is for machines that do not. Divergence time stays remote.)*
- [x] **M1.13 Hardening + API freeze v0.1** (2026-10-02) — README quickstart · `environment.yml` ·
  `docker/Dockerfile` skeleton · `docs/API_v0.1.md` (module contract, Result schema, coordinates,
  manifest keys, exit codes). *The environment checks in this milestone's Done-when — `make check` on
  a clean clone and `make golden` inside Docker — are deferred to the Phase 9 platform pass, when
  Docker is set up on the dev machine.*

## Phase 2 — Genome reading & statistics  *(§13.2 — plastid; the mito half waits on H-10)*

**Goal:** tables + figures + draft Results for the 2 species manuscripts (Oldenlandia, Coleus).
**Input:** GenBank files (downloaded or produced). Raw reads not needed.
**In-phase checks:** golden numbers — 153,724 bp · 132 = 86+38+8 · IR 26,214 (≥ 99.9 %) ·
LSC 83,127 / SSC 18,169 · 37 SSR (mono 31 · di 4 · tri 2) when MISA/cpstools available. Mismatch ⇒ stop.

- [x] **M2.1 Gene content & summary statistics** (per group, GC %, length)
- [x] **M2.2 Publication gene table** (Word export) — categorised table + reproducible `.docx`
- [x] **M2.3 Genome structure** — LSC/SSC/IR boundaries · %GC per region · GC-skew
- [x] **M2.4 Codon usage** — RSCU · ENC · CAI · CBI · Fop
- [x] **M2.5 Amino acid composition** — + GRAVY · aromaticity
- [x] **M2.6 SNP / substitutions by position** (2026-10-03) — scope labels enforced; the **CDS
      scope reproduces §6's golden Ts/Tv exactly**. Originally 317 Ts / 118 Tv = 2.686 over 76
      shared genes and 66,066 columns (§8.25); **re-measured 04/10/2026** as **323 Ts / 120 Tv =
      2.6917 over 79 shared genes and 68,448 columns** after three differently-spelled loci were
      recognized as shared (§8.27). The whole-alignment default is a *different* number (1.8416)
      and every table says which scope it is
- [x] **M2.7 Introns** — count · length · position
- [x] **M2.8 SSR** — by region and motif
- [~] **M2.9 Nucleotide diversity π** — **real-data run done 03/10** (§8.25): π 0.007181 whole,
      LSC 0.009185 · IRa 0.001682 · IRb 0.001681 · SSC 0.014066, per gene and 600/200 sliding
      window. **Still missing: per intron/IGS** — those spans need the intron derivation extracted
      from `modules/introns.py` into the shared layer first (DRY: it would be its second copy), so
      it is named rather than half-built here (H-07)
- [~] **M2.10 Manuscript deliverable** — tables + draft Results produced for both species from the
  public stand-in set (labelled); **the draft predates the π/Ts/Tv sections** (MAFFT has been
  present since 03/10, §8.25 — the draft simply has not been re-run with them) and figures; real
  samples outstanding
- [x] **M2.11 Sample-group labels** — every comparison table/figure carries a shared group label
  (moved here from M1.11, whose loader shipped in Phase 1: `core/samples.py`)

- [ ] **M2.12 Repeats (dispersed + palindromic)** — `repeat-match` + Vmatch on the golden fixture;
      the module reports every match ≥30 bp with its coordinates and region, the IR pair included;
      the count is **measured, not inherited** — the 26/9 survey's "9 matches" was measured on other
      data (§9.H's trap). *In scope since §1 and §7.1 and named in the canonical diagram (§13.0,
      `repeat (vmatch)`); it had no milestone until 07/10/2026.*

## Phase 3 — Annotation  *(§13.3)*

**Goal:** annotate genomes we already have; prepare NCBI submissions.
**In-phase checks:** the output re-reads with our own GenBank reader (no format damage) · ORF checks
pass · gene recount uses the fixed counting **rule** (CDS + tRNA + rRNA, both IR copies) — **132 is
the golden genome's value, not a constant**: P2-10 measured 131 (*S. diffusum*) and 133
(*C. amboinicus*), both correct. *De-novo assemblies (Phase 7 only) get the extra "verify vs close
relatives / reads across LSC–IR · SSC–IR junctions, else label undetermined" rule — skipped here.*

- [x] **M3.1 Annotation** (2026-10-03) — plastome: **Chloe** (PGA measured and replaced, §8.15/§8.17),
      **132 = 86 + 38 + 8** through our own reader, 0 CDS failing the ORF checks; mito path
      (`mito_genes`, no IR gate) done and measured; the mito **count check closed in Phase 3A** by
      this project's own annotator — PMGA never ran here (11 GB + MAKER against an 8 GB cap; its web
      server still refuses connections, re-measured 03/10) and nothing depends on it now (§8.22).
      *The annotate step itself is being replaced by this project's own tool — **Phase 3B** (D5-08):
      Chloe is dropped; this milestone records what Phase 3 shipped.*
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

**In-phase check — MET 03/10/2026.** On `OX459128.1` the annotator's **own output** re-reads with
**our own reader**, runs through `orf` + `mito_genes` unchanged, and recounts to **31 CDS · 17 tRNA
· 3 rRNA** — with **every difference from the published record listed** (13 CDS spans within 75 bp,
1 missing and 1 extra tRNA, 1 anticodon disagreement), and anything it cannot call reported as not
called. The counts agree; the *members* of two of the three sets do not, and that is what the
golden file prints. On the unannotated `AP045469.1` the same command yields 23 CDS of 31, pinned —
with the eight it could not call reported and measured (§8.24).

- [x] **M3A.1 GenBank writer** — read → write → read round-trips the golden's 265 features
- [x] **M3A.2 Reference set** — RefSeq mitogenomes chosen by taxonomy (same family first) → CDS
      proteins + rRNA, pinned, with sources and sha256 recorded
- [x] **M3A.3 Gene caller** (2026-10-03) — `blastx` for CDS and `blastn` for rRNA; boundaries come from the
      alignment, not from ORF-walking; everything below the floor is reported
- [x] **M3A.4 tRNA caller** (2026-10-03) — ARAGORN, declared in `environment.yml` and `docs/`
- [x] **M3A.5 CLI** (2026-10-03) — `plastid mito annotate`, output through the existing pipeline, no second path
- [x] **M3A.6 Golden check** (2026-10-03) — 31 · 17 · 3 against the published record, differences itemized
- [x] **M3A.7 AP045469.1** (2026-10-03) — the unannotated sample: run it, report it, pin it

## Phase 3B — Plastid annotator, self-built  *(slot decided 07/10/2026: this phase, with the naming step M3B.0 in front of the annotator; **board open 07/10/2026** — `TASKS.md` P3B-01…P3B-10, branch `feat/phase-3b`)*

**Goal:** the annotate step of the pipeline (assembly → annotate → analyses, with gene names
unified **on the reference set before** annotation and verified after it — D5-09) as this
project's own tool — **CGAS's flow with our search-and-decide core**: a pinned
reference set, per-sample annotation, an audit gate in front of every run. **Chloe is dropped**
(D5-08). The mito counterpart is Phase 3A; this is its plastid sibling. Design and the other open
decisions (D-a tRNA caller · D-b reference set · D-e shared core): `docs/D5-08_plastid_annotator.md`.

**In-phase checks (proposed):** the golden *sequence* annotates to exactly **132 = 86 + 38 + 8**
with **every difference from the published record listed** · the 34-record corpus runs
leave-one-out, each result diffed against the published record and every difference classified ·
a failed reference audit stops the run with `file — reason`, nothing annotated · no silent
failures: a search that finds nothing is a failed run, every below-floor hit reported.

- [~] **M3B.0 Naming (D5-09)** — the step that stands in front of the annotator: the audit's
      report is what the end user reads to pick the naming reference ★ (D5-09 §7.1); the three
      numbers of D5-09 §9.3 measured **before** coding (the margin · whether neighbours separate
      `rpl2`/`rpl23` and the two `trnL` · the `identity_floor`); `gene_normalize` — normalize the
      reference DB, verify the annotated samples — with every table in vocabulary **A**; the
      `gene_aliases.yaml` direction flips with its guard test; **Ts/Tv re-measured** (§6's 2.6917
      is a prediction until this lands)
      *P3B-02 (the three numbers, §8.41) and P3B-04 (`gene_normalize`, §8.44) are done 08/10/2026,
      and the alias direction is flipped. Left: **P3B-06** — Ts/Tv re-measured under the new
      direction — and the **★ pick**, which is the end user's (a recorded override stands until
      then).*
- [~] **M3B.1 Reference-set audit** — the gate (structural errors block; flip-flop and
      majority-deviation report), `--audit-only`, errors printed as `file — reason`, confirmed
      exceptions pinned by sha256; MAFFT `--adjustdirection` per region + a dotplot for the
      user's confirmation; **its report is the evidence the end user picks the naming reference ★
      from** (D5-09 §7.1)
      *P3B-01 (the set, §8.40) and P3B-03 (the gate, §8.43) are done 08/10/2026. What is left is
      the flag itself: `--audit-only` belongs to `plastid annotate`, the command it is a flag of,
      which arrives in **P3B-09** — the audit it will call (`core/reference_audit.audit_files`) is
      what P3B-03 built.*
- [ ] **M3B.2 Calling engine** — orthogonal searches (`blastx` on published `/translation` ·
      `blastn` for rRNA · IR self-blast behind the structural gate), candidates carrying their
      evidence, the six measured rules of §8.22
- [ ] **M3B.3 Arbitration & export** — structure first and score within one program; overlap kept
      and reported (the golden itself carries four overlapping pairs); many references → one call
      per locus with the losing candidates listed; ties → one call when identical, both with a
      note otherwise; GenBank that re-reads, `/organelle` on the source, a reason per rejected hit
- [ ] **M3B.4 CLI & batch** — 1..N samples × 1..N reference records; organelle declared at input
      (H-09); a failed sample is marked without killing the run; docx report (D2-04 pattern)
- [ ] **M3B.5 Golden + corpus validation** — the phase's Done-when at both scales, every
      difference itemized

## Phase 4 — Multi-genome comparison  *(§13.4 — plastid; the mito half waits on H-10)*

**Goal:** comparison content for the 2 manuscripts. **Deliverable: the three runs and their table
sets — `docs/P4-05_comparison_tables.md`** (two species over the pinned stand-in panel + the mito
sketch), every table carrying the stand-in label.
**In-phase checks:** Ts/Tv = **2.6917** for *I. domestica* × *I. gatesii* (CDS only · 79 genes · MAFFT, §8.27) ·
identity plot **7 blocks**, %IDY 89.66–99.11 (mean 97.37), coverage 99.09 % — measured 04/10
(§8.26). *The earlier "21 blocks · 83.7–98.2 · mean 91.5" could not be reproduced from the probe's
inputs and is recorded as a finding, §9.H.* Mismatch ⇒ stop, re-check inputs.

- [x] **M4.1 Ka/Ks selection pressure** (2026-10-04) — PAML yn00 (Yang & Nielsen 2000) per shared
      CDS on a codon alignment; 79 shared genes, 78 scored, pooled Ka/Ks 0.2167; unscoreable genes
      listed with a reason (§8.27)
- [x] **M4.2 mVISTA-style %identity + synteny** (2026-10-04) — `nucmer --maxmatch → delta-filter -1
      → show-coords -rcl`, MUMmer 4.0.1 pinned; block table + 600/200 % identity series + raw
      coords/delta as run artifacts
- [x] **M4.3 Barcode candidates & hypervariable regions** (2026-10-04) — π reused from the
      diversity module (773/773 windows identical in a cross-check), 26 regions and 10 candidate
      loci against the record's own annotation (§8.29)
- [x] **M4.4 Population genetics** (2026-10-04) — written here in pure Python (D4-03): Hd · θ(S) ·
      θ(π) · Tajima's D · Fu & Li's D*/F* · Nei's Gst, per gene and group; cross-checked against
      scikit-allel and hand-computed vectors (§8.28)
- [ ] **M4.5 RNA editing prediction** — **BLOCKED: needs RNA data**
- [x] **M4.6 Ts/Tv module** (2026-10-04) — the CDS-scope module (M2.6) is this phase's Ts/Tv
      deliverable; its shared-gene set is measured by locus, and the re-measurement that followed
      the three alias loci (§8.27) is what the manuscript tables quote: **2.6917 = 323 Ts / 120 Tv
      over 79 shared genes**, with the whole-alignment scope a different, labelled number

## Phase 4a — The Điềm corpus as the test dataset  *(opened 04/10/2026)*

**Goal:** the phases so far were tested on data this project picked itself. anh Điềm's corpus — **34
Iridaceae plastomes × 2 annotations** (published NCBI · PGA) **+ the bare sequences**, 28 MB in
`.data/thanhdiem_data/` (gitignored, in the working tree) — is the first dataset nobody here chose,
and it carries its own answer key three ways. It is a **test corpus**: Iridaceae, not the manuscript
species, so nothing from it becomes a manuscript table.

**In-phase checks (all measured 04/10/2026):** the golden genome's own file is **in the corpus**
(set 01 `Iris_domestica_NC_050833.1.gb`, sha256 `a946f23e…` — byte-identical to the vendored
fixture) · set 02 is the *same sequence* annotated by PGA, **85 CDS where the contract says 86**,
with **3/34 records failing the IR gate at 26–29 %** while de-novo finds 100.00 % on every one ·
**34/34** records spell at least one locus differently between the two annotations — §8.27's lesson
at 34× scale.

- [x] **M4a.0 Corpus in the working tree** (2026-10-04) — `.data/thanhdiem_data/`, 136 files, every
      sha256 identical to `D:\business\pts\thanhdiem_data`; gitignored; ~100× faster than `/mnt/d`
      on the small-file writes a run is made of (measured)
- [x] **M4a.1 Pinned** (2026-10-04) — `.panel/iridaceae-34/` (0 gate failures) and
      `.panel/iridaceae-34-pga/` (3, named) with `PINNED.json` as the verification table, every
      sha256 · md5 · byte count and measured value re-read; the `.fasta` half as `fasta.json`
      (68 entries) because the pin route is GenBank-only (measured), with all 68 sequences proven
      identical to the pinned records. §8.31
- [x] **M4a.2 Vendored fixtures + Phase-1 checks** (2026-10-04) — five files · 1.3 MB in
      `tests/data/corpus/` (three pairs; the golden's published half is not duplicated), and
      `tests/golden/test_corpus.py` asserts each one: the golden re-annotated is one CDS short,
      *Iris hermona* passes the gate as published and fails it at 26.36 % under PGA on the same
      sequence (exit 0 / exit 2 through the CLI), the group's own record 87 → 85, and the
      `ndhI`/`ndhL` alias resolving both spellings to one key. A `SHA256SUMS` completeness check
      now guards every vendored file. §8.31
- [x] **M4a.3 Scale (M8.3)** (2026-10-04) — wall clock and peak RSS for every module run at N = 34;
  `diversity` and `barcodes` are **named as not measured** (they align whole genomes like
  `substitutions`, whose curve *is* measured). `docs/M8.3_scale.md` + PROJECT.md §8.32
- [x] **M4a.4 Annotator at scale** (2026-10-04) — 34/34 bare sequences annotated by Chloe, 0 CDS
      failing the ORF checks, and three annotations of one sequence side by side: **no record has
      all three agreeing (0/34)**, and the differences cluster on `ycf1` (24/34) and the
      trans-spliced `rps12` (22/34). Needed a FASTA door in `plastid run`, which exposed a third
      bug of the same family. §8.33

## Phase 5 — Figures & reports  *(§13.5)*

**Goal:** publication-grade output, regenerable from intermediates (never by hand).
**In-phase checks:** every table/figure regenerates from intermediate data · **queryable-number test**:
every number in a report traces to exactly one cell in a source table · tool report has name, version,
path, threshold.

- [x] **M5.1 Figure system** (2026-10-05) — eight figure types in `result/figures/`, one theme, SVG ·
  PDF · EPS byte-identical between runs (§8.36), and the N-scaling rule **measured** at 2 · 8 · 34:
  a series cap of 8, above which the identity panel draws the span the genomes cover (§8.37)
- [x] **M5.2 Statistics** (2026-10-05) — `core/stats.py` (scipy, `[stats]`): descriptive ·
  Mann–Whitney U · Kruskal–Wallis H · Spearman · Pearson · Benjamini–Hochberg · Bonferroni, every
  one pinned by a hand-computed vector; the `stats` module runs a real two-group comparison and
  reports a test the data cannot support as *not computable* (§8.38)
- [x] **M5.3 Exports** (2026-10-05) — multi-sheet Excel (byte-identical) · `plastid package` (the
  Supplementary package, regenerated from a run directory) · one-page PDF per genome, **34/34 corpus
  records** · `plastid qc` (the technical QC report) · the Word export unchanged (§8.39)
- [x] **M5.4 Queryable-number check automated** (2026-10-05) — `report/audit.py` fails on an orphan
  number and on a number two cells disagree about; clean over the whole Phase-4 deliverable run (§8.35)

*Of M5.1's eleven figure types, **seven have a source table today and four do not** (measured
04/10/2026): the haplotype network needs a haplotype table (`popgen` counts haplotypes,
it does not list them), the annotated tree needs a phylogeny nothing builds yet, NUPT/MTPT is M6.1
(blocked on a nuclear genome) and the mito configuration is M6.2. Phase 5 shipped the seven, plus a
PDF-only one-page sheet per genome (P5-04).*

**Phase 5 — CLOSED 05/10/2026**, merged to `main` (fast-forward). Five tasks, each with its
Done-when measured: the figures (§8.36–§8.37), the statistics (§8.38), the exports (§8.39) and the
queryable-number audit (§8.35). Three layout defects were found by *looking at* the rendered
figures rather than by a test — the lesson of §9.C-34, paid three times in one phase.

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

*The canonical diagram (PROJECT.md §13.0) names `fastp` · `BWA/minimap2` · `Bandage` for read QC;
neither `fastp` nor `BWA` is declared in the repo yet — declare each in the session it is first
used (the §7.8b lesson); the map-back step itself is M7.1's. The long-read assemblers (`OATK` ·
`raven-assembler`, both declared) belong to **M7.6**; `Flye` left the diagram — measured to die on
RAM here and replaced (§7.1, PROJECT.md §13.0's dated note).*

- [ ] **M7.1 Batch assembly + QC** (fastp → GetOrganelle · BWA/minimap2 map-back) — **BLOCKED: no real Illumina reads yet**
- [ ] **M7.2 Assembly-graph checks** (Bandage)
- [ ] **M7.3 Feed mito-structure work** (M6.2)
- [ ] **M7.4 Remote path** (usegalaxy.eu `get_organelle_from_reads`) — subsample reads before upload ·
  pin tool version · record history ID · data-sharing policy for unpublished reads: **open**
  (PROJECT.md §10.12)
- [ ] **M7.5 Peak-RAM with real reads** — **BLOCKED: no real reads yet**
- [ ] **M7.6 Long-read path** (OATK · raven-assembler — ONT/HiFi) — **BLOCKED: no ONT/HiFi reads**;
  both tools are declared (`environment.yml`) and were run in the survey (§7.1); feeds M7.3

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
| Reference Panel item (c) — source text truncated | documentation only — does not block M1.10 | PROJECT.md §10.9 |

*Decided 07/10/2026 (dev): the self-built annotator's slot is **Phase 3B, opened next** (after
Phase 5, mirroring 3A), and the naming step (D5-09) sits **inside it as M3B.0** — the reference set
must be name-normalized before the annotator uses it (D5-08 §7.3). `docs/D5-08`'s D-c and
`docs/D5-09`'s §9.5 are answered there; D-a · D-b · D-e stay open.*

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Real Illumina reads (1 sample) | Diem / Trang | M7.1 · M7.5 · M4.5 |
| ONT/HiFi reads (1 sample) | Diem / Trang | M7.6 (long-read path) |
| Plant mito data — a real mitogenome dataset | anh Điềm / cô Trang | M6.2 · M6.3 · H-06 · H-10 · H-12 (on hold by D5-07) |
| Nuclear genome of the sample | Trang | M6.1 (NUPT/MTPT) |
| A macOS machine | Team | M9.4 platform claims |
