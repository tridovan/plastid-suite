# D5-09 — reference-driven gene naming (`gene_normalize`)

**Status: design record, 07/10/2026.** Sections 1–8 are decided; §9 lists what is still open and
what must be measured before coding. No code exists for this step yet.

## 1. The pipeline, and the tools

*The canonical copy of this diagram is `PROJECT.md` §13.0. It is repeated here so this record reads
standalone — change `PROJECT.md` first, then copy it back.*

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

**How to read it.** The left column is the **reference supply chain** — not a sequential step of the
sample flow but a source feeding into step 2. The long arrow on the left *is* "the normalized
reference DB is what annotates the samples".

Three things the picture is meant to make hard to forget:

1. **Name unification happens before annotation, and runs on the reference DB.** The samples
   inherit vocabulary A through annotation rather than being renamed afterwards.
2. **Step 3 is the same module run again in verify mode.** If the annotator takes its names from an
   already-normalized reference, its output must already be in A — so **any change at step 3 is a
   finding**, not routine.
3. **RNA editing is the only branch with a reference of its own**, and it sits *after* the
   normalized input — so its reference B has to be paired back onto A (§3).

### Tools, step by step

| # | Step | Tools | In → Out | Reference needed |
|---|---|---|---|---|
| 0 | Reference supply | `plastid panel refresh` · `propose` · `pin` | NCBI / a team file → a pinned DB (`PINNED.json`, sha256) | — |
| 0b | Quality assessment | the audit (`docs/D5-08_plastid_annotator.md`) | DB → report: BLOCK · REPORT | — |
| — | **Unify the reference DB** | **`gene_normalize`** (this file) | DB + ★ the chosen file → DB in vocabulary A | **★ the naming reference** |
| 1 | Assembly | GetOrganelle + SPAdes (short read) · OATK · raven-assembler (long read); fastp · BWA/minimap2 · Bandage for QC | reads → FASTA | seed DB `embplant_pt` · HMM DB `embryophyta_pltd.fam` |
| 2 | Annotation | **self-built** (D5-08) · BLAST+ 2.17.0 | FASTA → GenBank | **the reference set, 1..N records** |
| 3 | Verify names | **`gene_normalize`, verify mode** | GenBank → checked GenBank | **★ the naming reference** |
| 4 | Analyses | MAFFT · MISA · MUMmer4 · PAML · vmatch · pycirclize | A → tables · figures | kind 1, chosen by taxonomy (Appendix B) |
| 5 | RNA editing | *not yet chosen* (C.3) | A → editing sites | **a separate reference B** (Appendix B) |
| 6 | Submission | table2asn | tables → `.sqn` | submitter template `.sbt` |

## 2. What the step does

The step takes **any GenBank that enters the pipeline** — the reference DB before annotation, the
team's samples after it, and the separate reference a downstream tool needs. It writes each back out
**named after the project's naming reference**, plus the tables that say what it changed and why.

**Two runs, one module, two expectations:**

| Run | On | Mode | Expected outcome |
|---|---|---|---|
| 1 | the reference DB | `normalize` | renames happen; the DB comes out in A |
| 2 | the annotated samples | `verify` | **nothing to change** — a change means the annotator named something itself, or called a gene the reference does not have |

Two principles make this more than a rename pass:

- **Nothing is guessed.** Every rename carries the evidence that decided it (§4); anything that
  cannot be decided honestly is reported rather than resolved (§5).
- **It is cross-checking, not merging.** The goal is not to make two sources agree; it is to
  **surface the places where they disagree**. A disagreement is a finding, not a mess to clean up.

## 3. Vocabulary: A and B

| | What it is | Where it may appear |
|---|---|---|
| **A** | The project's naming reference — one GenBank record, chosen by the end user (§7.1) | **Every result table, every manuscript table** |
| **B** | Any other source that enters: a tool's reference set, a foreign record, a comparison genome | **Only in the `pairing` table**, in a column of its own |

**B lives at the input end and nowhere else.** It exists to *find* the counterpart of a gene; once
the pair is found its job is done. It stays in `pairing` so a reader can see what was matched
against what — and so it is provable that nothing was matched by guesswork.

> **No table may mix two vocabularies.** Either the whole table is A, or a column says which
> vocabulary each name belongs to. A reader must never have to guess.

This is the rule the whole step exists to serve: A is the join key. If one analysis printed B, that
table would not join to the gene table, to Ka/Ks, or to gene content — and every later join is a
place a silent error can enter.

### Three edge cases

| Case | Handling |
|---|---|
| A gene exists in B but not in A | print **B's name** plus an explicit **`outside vocabulary` flag column** — never silently re-labelled to the nearest A name. This is the only place B may appear in a result table, and it must be a column, not a font or a colour |
| A gene exists in A but not in B | `not measured` — see §5 |
| B spells a name "better" than A | **change A** — the project's naming reference — and re-run everything. Never a local exception for one analysis |

That last row is a rule, not advice: **the vocabulary is a project-level decision, not a
per-analysis one.** If each analysis could pick its own vocabulary from whatever reference it used,
three analyses would produce three name sets — which is what this step was built to prevent.

## 4. Arbitration: the name proposes, the sequence disposes

Four tests, in order, stopping at the first that decides:

| # | Test | Rules out |
|---|---|---|
| 1 | **Name** — syntax rules, then the alias table, → a canonical name | nothing; it produces the hypothesis |
| 2 | **Sequence ranking** — compare this gene against every candidate; the **best match must be the gene carrying that name** | a name that points at the wrong locus |
| 3 | **Margin** — `best − runner-up` must be wide enough | a ranking too close to trust |
| 4 | **Neighbours** — the genes flanking it, *within the same region* (LSC / IRa / SSC / IRb) | paralogues sequence alone cannot separate |

| Name matches? | Sequence matches? | Verdict |
|---|---|---|
| yes | yes | **confirmed** |
| yes | **no** | **conflict** — same name, different locus → stop (§5) |
| no | yes, uniquely | **proposed rename**, printed with its identity → needs a decision |
| no | no | `outside vocabulary` / `not measured` |

**Why not sequence alone.** Matching purely on "best sequence hit" would silently *correct* a
mis-annotated record and never say so. This project has met exactly that: NC_024936.1 spells
`atpI`·`psaI`·`psbB` as `atpL`·`psaL`·`bspB` (PROJECT.md §8.27). Sequence-only matching would have
fixed it and reported nothing. **What is wanted is not the correction — it is knowing where the two
sources disagree.**

**Why not name alone.** That is CGAS module 3, and it is the defect this design exists to avoid
(Appendix A).

**Why sequence is not enough either.** **[ĐÃ THAY 08/10/2026 — the two examples this paragraph
was written around do not survive the measurement; the rule does. Measured (PROJECT.md §8.41,
`tests/golden/test_naming_numbers.py`): `rpl2` vs `rpl23` are **22.03–44.44 %** identical across
the 34 corpus records and `trnL-UAA` vs `trnL-CAA` **61.18 %** — neither pair is near-identical,
and the ranking separates both. The measured cases that need test 4 are: the **serine tRNAs**
(26 of 1,182 cross-record tRNA comparisons rank a different locus first, by up to 15.05 identity
points), **two genes under one name** — the golden itself annotates two trnG features as
`trnG-UCC`, 73.24 % identical — and names the vocabulary cannot resolve at all (72 tRNA name
occurrences and 1 CDS name across 33 records). Original text kept below.]**
`rpl2`/`rpl23` and the two `trnL` genes are near-identical;
protein identity is high against both and ranking cannot separate them. That is what test 4 is for —
the project's existing "structure before score" principle (D5-08).

The column `rule` records which test decided, so a later reader can tell "confirmed by name and
sequence" from "confirmed only by gene order".

**What is compared.** CDS → **protein**, using the record's published `/translation` where it has
one, never a self-translation (§8.21). tRNA and rRNA → nucleotide (too short to translate).

## 5. Conflicts, and the gate

| Tier | What falls in it | Behaviour |
|---|---|---|
| **Automatic** | case folding, `trnL(UAA)`, `rrn16S`, and renames confirmed by test 2 or 3 | applied, recorded in the report |
| **Stop — the user decides** | same name / different locus · a rename that would create a duplicate · a name that cannot be resolved · a reference that fails its audit | the run **stops**, naming every case with its coordinates |
| **Reported, never stopping** | extra genes · missing genes · `/product` differs · CDS length and intron differences | written to the tables; the run continues |

Extra and missing genes sit in the third tier deliberately: a lost gene is a **biological result**,
not a tool failure. Stopping on it would make a real finding indistinguishable from a broken input.

**How "the user decides" works, without an interactive prompt.** The module is a pure function and
must run under Docker and in scripts, so the decision is two invocations, not a question:

- **default** — stop with a named error and a non-zero exit code, listing every conflict with its
  coordinates (the project already uses exit 2 for the IR gate and exit 3 for a missing tool);
- `--set on_conflict=continue` — proceed, writing the normalized GenBank artifacts **and** a
  `conflicts` table carrying contig, coordinates and strand.

A yes/no prompt in `cli.py` (TTY only, degrading to "stop") may be added as a convenience, but it is
not the mechanism.

**Every gene carries a status, and three values must stay distinguishable:**

| Status | Meaning |
|---|---|
| `measured` | a counterpart was found in A; there may or may not be a result |
| `not measured` | **no counterpart was found in A** ← must never be read as "none found" |
| `outside vocabulary` | present in B, absent from A |

Without this column `0` and `not measured` look identical, and every per-gene aggregate is wrong.

## 6. The module

A **new module**, `gene_normalize` (arity 2, `requires_ir=False`, modes `normalize` and `verify`),
registered in `modules/__init__.py`. **`modules/gene_names.py` is not touched.**

| | extend `gene_names` (rejected) | new `gene_normalize` (chosen) |
|---|---|---|
| CLAUDE.md — *"Extend by adding, not editing … Never grow an existing module to cover a new job"* | violated | satisfied |
| Input contract | one module carrying two: a **symmetric** panel today, a **privileged** record tomorrow | one contract each |
| The gate's home | a gate inside a *reporting* module | a gate inside a *transforming* step |
| Existing tests, P3-02 goldens | module identity changes under them | untouched |

The cost is zero for the user: `-m` is repeatable, so `plastid run -m gene_names -m gene_normalize …`
is still one command and one run directory.

- **All renaming logic lives in `core/gene_names.py`** — it stays the only place a gene is renamed
  (P3-02's rule).
- **The normalized `.gb` is a `Result.artifacts` entry**; `runner.py` already writes artifacts, so
  the export layer is untouched.
- **The comparison tables** (extra / missing / length / intron) live here rather than in a third
  module: once every gene has been matched, those columns are a **by-product**, and a separate
  module would compute the same matching twice.
- **Verify mode reports and does not write.** It produces the same tables; a non-empty change set is
  a finding, and whether it stops the run follows §5's tiers.

## 7. The project's naming reference

### 7.1 The user chooses it, from the audit report

**Which file becomes the standard is not a team decision.** It is the **end user's** — a researcher
in biotechnology / bioinformatics, the people this toolkit is built for (defined 08/10/2026) — made
after reading the report produced by the reference-set quality assessment — the audit in
`docs/D5-08_plastid_annotator.md` §"The reference-set audit", whose REPORT class already compares
gene sets as **name sets** and names the records that deviate from the majority. That is exactly the
evidence a choice needs.

So the toolkit's job is not to pick a file. It is to make the choice **possible, informed and
recorded**:

- the audit produces the report and the verdicts, and blocks only what cannot be used at all;
- the user names the chosen record;
- the run pins it, records it, and refuses to proceed without it (§7.2).

### 7.2 What is recorded

| | |
|---|---|
| **Identity is content, not a path** | record `{path, accession, sha256}` — a path string alone cannot tell whether the reference changed |
| **The digest enters the cache key** | otherwise a run can be served a result computed against a different reference (Appendix C) |
| **A run may override the standard; the override is recorded** | so two runs are comparable, and a run that used something else says so |
| **Pinning reuses existing conventions** | `tests/data/PROVENANCE.md` + `SHA256SUMS` for anything committed; `.panel/*/PINNED.json` for anything fetched |
| **The parameter must not be called `reference`** | that key already carries two incompatible meanings (Appendix C); proposed name in §9 |

### 7.3 This reference is a sixth kind, and must not be merged with the others

Surveyed 07/10/2026 (Appendix B). The plastid path has **three adjacent stages with three different
answers** to "what is the reference":

| Stage | Reference | Chosen by |
|---|---|---|
| annotate — today | Chloe's own curated set; a different one is **refused** | the tool |
| annotate — after D5-08 | a reference **set** of 1..N records, same family first | the team, by taxonomy |
| **name unification — this file** | **one record, fixed for the project** | **the end user, from the audit report** |

The naming reference is a kind of its own: **it is not chosen by taxonomy per run**, and it is not
the tool's. Writing that down is what stops a later change from quietly collapsing the three into
one.

**One consequence for D5-08.** Its annotator now has a step standing in front of it: the reference
set must be name-normalized before it is used. D5-08 does not say this; it needs a line pointing
here.

## 8. Consequences of reversing D3-02

Canonical direction stops being the group's dialect and becomes the naming reference's. Measured
against the golden `NC_050833.1`, **five of the nine measured loci flip**:

| Locus | Golden spells | Canonical today (D3-02) | Under this decision |
|---|---|---|---|
| ycf3 | `ycf3` | `pafI` | `ycf3` |
| ycf4 | `ycf4` | `pafII` | `ycf4` |
| psbN | `psbN` | `pbf1` | `psbN` |
| clpP | `clpP` | `clpP1` | `clpP` |
| ndhL | `ndhL` | `ndhI` | `ndhL` |
| atpI · psaI · psbB · trnI-GAU | match | match | unchanged |

Carry these:

- **Three analysis modules read the canonical key** — `kaks`, `substitutions` (this is Ts/Tv) and
  `popgen`, all through `cds_by_canonical_name`. Each locus still maps one-to-one, so the *shared
  gene set* should not change and **Ts/Tv 2.6917 should hold — but that is a prediction, and it
  must be re-measured, not assumed.**
- **Printed names change** in the gene table, and therefore in manuscript tables. D3-02 also
  recorded that the names *printed in a manuscript* and *submitted to NCBI* are the PI's call
  (D3-06) — this decision changes what the default print is.
- **A bad reference becomes a bad standard.** Were the chosen reference NC_024936.1, then
  `atpL`·`psaL`·`bspB` — three names this project measured as wrong — would become canonical for
  every query. This is why §7.1's audit is load-bearing rather than decorative.
- **Documents corrected in place** 07/10/2026 with the project's `[ĐÃ THAY — …]` convention:
  `PROJECT.md` §3 (D3-02, plus the D4a-05 alias row and a D5-09 decision row),
  `docs/gene_name_variants.md` §6 (the direction note, and two CGAS citations that were off:
  `return candidates[0]` sits at **211 and 214**, not 208, and the two-way prefix at **149**, not
  150–153), and `docs/D5-08_plastid_annotator.md` line 54 (a pointer that the reference set must be
  name-normalized before the annotator uses it, and that the unification step now writes the
  normalized GenBank).

## 9. Open, and what must be measured

1. **The parameter's name** — proposed `naming_reference`, to avoid a third meaning for the key
   `reference`.
2. **Whether `MANIFEST_SCHEMA_VERSION` is bumped** for the new per-module reference entry.
3. **Three numbers, measured before coding — MEASURED 08/10/2026, PROJECT.md §8.41.** They are
   now thresholds rather than questions, and the answers are not the ones this file guessed:
   - the **margin** of test 3 is a **score** margin, in bits — a margin in identity points is
     meaningless when the runner-up's alignment is partial (`psaB` against `psaA` reads a
     *negative* margin in points). Measured over 2,838 cross-record CDS queries: smallest margin
     **19.2 bits**, median 298.0, none under 5. Recommendation: **margin ≥ 10 bits**.
   - **test 4 is not needed for `rpl2`/`rpl23` or the two `trnL`** (44 % and 61 % identical —
     see §4's correction). It is needed for the serine tRNAs, for two genes sharing one name,
     and for names the vocabulary cannot resolve;
   - the **`identity_floor`**: with **coverage ≥ 0.7** the highest identity a genuinely different
     locus reaches is **44.93 %** and the lowest a true ortholog reaches is **67.75 %**.
     Recommendation: **identity ≥ 50 % AND coverage ≥ 0.7** — the mito path's identity floor
     holds, its coverage floor (0.5) does not: plastid CDS are short (`psbJ` is 36 aa) and a
     21-aa alignment of it reaches 52.38 % against `matK`. **For tRNA no identity floor works**
     (53.76–100 % true vs up to 73.97 % wrong) — the decision must come from the name, the
     anticodon or the neighbours.
4. **PI sign-off** on the printed-name change, since it reaches manuscript tables.
5. **Phase slot** — **decided 07/10/2026: inside Phase 3B, as M3B.0**, in front of the annotator
   (`ROADMAP.md`); D5-08's D-c is answered there too.
6. **The cache defect of Appendix C is not part of this design and should not wait for it.** It
   returns silently wrong tables today. Whether it is fixed inside this work or as its own task is
   the dev's call, but the fix belongs to `run/cache.py` and `runner.py`, not to `gene_normalize`.
7. **The exit code a conflict stop uses.** The project's ladder today is 1 error · 2 usage (the IR
   gate — and D5-08's audit already claims 2) · 3 missing tool (`src/plastid/cli.py`). A conflict
   stop needs a number that stays distinguishable from the IR gate's; pick it here, not in the code.
8. **The organelle group `gene_normalize` declares.** The group table in `PROJECT.md` §13 is the
   single source and a new module must be declared there — "shared" looks right (like `gene_names`),
   but the declaration is what counts.
9. **§7.2's cache promise needs Appendix C's fix to be true.** Measured 07/10/2026: `Cache.key`
   sorts the input hashes and `_input_hash` hashes the sequence only, so the ★ record's *annotation*
   cannot reach the key — on the corpus pair (one sequence, two annotations: 86 vs 85 CDS) the
   `_input_hash` values are equal and the key is identical with the two inputs swapped. Either the
   content hash lands with this work, or §7.2 is recorded as unmet until item 6's fix does.

---

# Appendix A — CGAS module 3: what is taken, what is not

The copy under review is `module3_source/cgas_module3.py`, md5 `2851c57a6ca31aa568c57caec5cf195b` —
byte-identical to CGAS commit `52b82e60…` and to this repo's own review checkout
(`.data/cgas-review/cgas/cgas_module3.py`, same md5).

**Taken (the workflow):** one reference GenBank plus N query records · the reference is the
vocabulary · each query is written out as a normalized GenBank · the run produces a comparison ·
one command.

**Not taken (measured 07/10/2026, on the four *Erigeron* records shipped with it — `BK070967`
*compositus* · `BK070968` *denalii* · `BK070970` *elongatus* · `BK070975` *acris*):**

| Defect | Measurement |
|---|---|
| The bare-tRNA fallback is **file order** | `match_trna` ends `return candidates[0]` (line 214) over a dict filled in reference *feature order*. The three query records spell `trnS` bare; the reference lists `trnS-GCU` first (line 123). Reordering the reference's trnS features in a faithful re-implementation returns `trnS-UGA` or `trnS-GGA` instead. The answer was right **by accident** |
| The anticodon is **destroyed before it is read** | `match_trna` strips `\([^)]*\)` at line 159, then looks for an anticodon after a hyphen — so the bracket form `trnL(UAA)` loses the only information that could disambiguate it and falls through to the same guess. This project has measured **30 names in that form** in one record (NC_058328.1) |
| **Two-way prefix match** | line 149: `raw.startswith(ref) or ref.startswith(raw)` → first hit in reference order. `ycf1-like` → `ycf1`, silently |
| Renames leave **no trace in the file** | the output contains **0** occurrences of the bare `trnH` it renamed; the rename exists only in the Excel |
| The reference is **never audited** | whatever the reference spells becomes the standard for every query |

Its own reported numbers were reproduced exactly: 8 features changed per target (4 tRNA + 4 `gene`),
24 in total, no duplicate created, no gene gained or lost.

**Checked and found sound** — the README's claim that the output keeps everything else unchanged
holds: sequence, every feature location and every non-`/gene` qualifier are identical after the
Biopython round-trip; only the file's *text* is re-serialized (PUBMED key spacing, `/translation`
line wrapping, location notation). Measured feature by feature, not assumed.

**Two observations about the *Erigeron* records, recorded but not explained** — do not put them in a
manuscript table before they are checked. All four carry `129 = 85 CDS + 36 tRNA + 8 rRNA` against
the golden's `132 = 86 + 38 + 8`, with the same 112 distinct names. All four lack **`trnT-GGU`** (a
72 bp probe taken from the golden finds its best window at 83.3 %, against a 100 % control in both
*Iris* records — so it is absent rather than unannotated) and all four carry a **`trnL-UAG`** the
golden does not have. Whether that is lineage biology or an annotation artefact is **not measured**.

# Appendix B — reference mechanisms of the plastid tools

Surveyed 07/10/2026 from the source and from `PROJECT.md` §7. Only tools the project actually calls
are listed; tools declared in `environment.yml` with no code path (cpstools · mreps · vmatch ·
pycirclize · sibelia · bandage · oatk · getorganelle · spades · flye · novoplasty · minimap2 ·
samtools · wgsim · beast2 · treePL) are not.

| Kind | Tools | Reference is |
|---|---|---|
| 1 — a species / a record set | Chloe · PGA · the self-built annotator (D5-08) · GeSeq · mVISTA (→ nucmer + matplotlib) · PAML mcmctree · BEAST2 | a genome, chosen by taxonomy |
| 2 — a seed or database | GetOrganelle (`embplant_pt`) · OATK (`embryophyta_pltd.fam`) | a downloaded data file |
| 3 — the user's own genome | minimap2 + samtools (RNA editing measured from RNA) · BWA | the record itself |
| 4 — an artifact of the team | table2asn (`.sbt` template) | not a genome |
| 5 — none | MAFFT · MISA · repeat-match · Vmatch · nucmer · Sibelia · pycirclize · PyPopART · Bandage · SPAdes · raven-assembler | — |
| **6 — the naming reference (this file)** | `gene_normalize` | **one record, fixed for the project, chosen by the end user (§7.1)** |

Kind 6 is a kind of its own. In particular it is **not** kind 1: the naming reference is not
re-selected per run by taxonomy.

**Where a tool's reference cannot be substituted.** Chloe reads its reference set from its clone's
`Artifacts.toml` and **refuses** a different one — `external/chloe.py:73` pins the tree hash
`9163915d99fdfbc7d8510f0479b3d08905682ca5`, checked at `:200-207` with the message *"a different
reference set is a different annotation"*. The module records what was used as
`its own curated set (Chloe_refs 9163915d99fd)` (`modules/annotation.py:273`) — the pattern to
follow: when a tool ignores the user's reference, the run has to say so.

**The RNA-editing case, and why the pairing is ours.** The two published plastid RNA-editing
predictors (PREP-Cp, PREPACT) are both web services and both measured unreachable (`http=000`,
PROJECT.md §7). PREPACT's *alignment* mode accepts a **user-defined reference**, but it matches by
**alignment position, not by gene name** — so supplying our own reference moves the whole naming
problem into the step that builds the alignment, which is ours, and it fails silently: a gene with
no counterpart yields **0 predicted sites**, indistinguishable from a gene with no editing sites.
That is the concrete reason §4's tests and §5's status column are not optional.

# Appendix C — the repo's own reference mechanisms, and one measured defect

The three mechanisms, and who consumes them:

- **Reference Panel** (`core/panel.py`, `plastid panel`): a data source **and** an integrity ledger
  (`PINNED.json` with sha256/md5/read-back measurements). **No production code reads it** — in
  `src/` only `cli.py:31` imports `core.panel`, and `load_pinned()` has exactly one caller,
  `tests/unit/test_panel.py:170`. Humans pass the files to `plastid run -i`.
- **Annotation reference** — a directory of GenBank records. `annotation.py:284` and
  `mito_annotate.py:102` read the **same run-global key `params["reference"]`** for two incompatible
  purposes; `--set` is one flat dict per run, so a run holding both modules cannot give them
  different values.
- **The naming pool** (`core/gene_names.py` `references_from()`): `(canonical name, tRNA sequence)`
  pairs drawn from **the run's own records** — not a genome.

And an unenforced positional convention: six modules take `genomes[0]` as the reference with no
parameter, no manifest field and no validation — `kaks.py:63` · `substitutions.py:51,198` ·
`popgen.py:66` · `synteny.py:49` · `diversity.py:51` · `barcodes.py:68`.

**Measured defect — the cache is blind to input order.** `Cache.key` **sorts** the input hashes
(`run/cache.py:39`); `_input_hash` hashes **only the sequence** (`run/runner.py:51-52`); the IR
fingerprint is also a property of the sequence (`runner.py:55-58`); and the cache lives at
`out_root/.cache`, shared across runs (`run/layout.py:58-59`). So two records carrying the same
sequence produce the same key in either order — and this repository has exactly such pairs on
purpose: the corpus fixtures are *"one sequence, two annotations"*. `plastid run -m kaks -i
published.gb -i pga.gb` and the same command with the inputs swapped therefore share a cache key,
and the second is served the first's tables (`runner.py:143-146`, status `cached`) — while `kaks`
treats `genomes[0]` as its reference and the two annotations have **different gene sets** (85 vs 86
CDS on the golden pair). The run directory's own fingerprint *is* order-sensitive
(`cli.py:476-478`); only the cache is not.

This is not this design's bug, but it is this design's problem: a step whose whole output is "named
after the reference" must not be cache-blind to which reference that was. Hence §7.2.
