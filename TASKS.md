# TASKS.md — active tasks

> **Phase 8B is OPEN** on `feat/phase-8b` — the **Crataegus answer-key corpus**: the toolkit's
> numbers beside a **published paper's own result tables**, on an independent Rosaceae corpus the
> team just supplied (anh Điềm's `db_test_new`, 10/10/2026). Phase 8 (tool cross-checks, scale, CI)
> is closed and merged to `main`.
> Milestones: `ROADMAP.md` → Phase 8B · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

> **PHASE 8B IS CLOSED (10/10/2026)** — tagged `phase-8b`, squashed to `main` as one commit. All
> nine tasks passed their Done-when. What remains are **questions for anh Điềm**, not work on this
> board (see "Waiting on data / people"). **Next: Phase 9** (packaging · licence review · release),
> which inherits the cross-check tools' licence questions.
> The headlines, so nobody re-derives them: **S7 reproduces to the last printed digit**
> (*Amelanchier cusickii* + MLWL, §8.69) · **Table 1 and the 50 barcodes reproduce exactly**
> (§8.67 · §8.72) · **S6 reproduces at MISA's defaults** (14/49 sections whole, all 49 compared at
> both settings, §8.73) · **S5 does not reproduce** — 1–2 % under all 11 counting rules tested,
> with the sequences proven identical (§8.70) · **our annotation differs from the published GeSeq
> one in 341 conclusion-level places, and two thirds of them are the deposits' own gaps or one
> reference record's outlier calls** (§8.71) · and **H-24**: `trnV-UAC` is never called (§8.71).

## Phase 8B — The Crataegus answer-key corpus (M8B.1 · … · M8B.8)

**Goal:** every number this toolkit publishes for a plastome has now been checked against another
*tool* (Phase 8). This phase checks them against another **laboratory**: Sun et al. 2026
(*Front. Plant Sci.* 17:1767012) printed their own tables for the same 49 genomes, and the Methods
name every tool and parameter they used — so a difference can be attributed instead of guessed.
**In-phase checks:** every published table quoted with its own scope · every difference with a named
cause (parameter · counting rule · class definition · coordinate convention) · every output labelled
as third-party published data · the corpus's provenance recorded **in the repo**, not only on `D:`.

- [x] **P8B-01 · The MISA compound flag** *(the defect the corpus found first)* — **DONE
      10/10/2026 (§8.65)**: MISA marks a cluster's components two ways and the parser read only
      one, so a component on its own line said "no" — 0 of PX413282's 16 components flagged while
      MISA's own statistics counted 8 clusters. Test written first, module 0.3.0; the golden still
      reads 37 = 31 · 4 · 2, and its table now flags its **6** components (was 2).

- [x] **P8B-02 · The corpus in the working tree, verified** — **DONE 10/10/2026**: 117 files,
      **byte-identical** to `/mnt/d/business/pts/thanhdiem_data/db_test_new` (sha256 per file),
      gitignored under `.data/thanhdiem_data/db_test_new/`. 100/100 GenBank records read by our own
      reader; 50 plastomes all declare `/organelle` and pass the IR gate at 100.00 %.
      **Done when:** the copy is byte-identical (sha256 per file) and 100/100 records read —
      **met 10/10/2026**. The repo-side provenance note (accessions · retrieval date · DOI ·
      manifest sha256) is carried by **P8B-09**, because it is a record, not a copy.

- [x] **P8B-03 · Structure & GC vs Table 1** *(M8B.2)* — **DONE 10/10/2026 (§8.67)**: 44
      plastomes run; **18/18 Table-1 rows match our four sizes exactly**; per-region GC matches
      to the paper's printed precision except **one cell** — *C. jozana* IR GC 45.62 (paper) vs
      42.64 (ours), with that row's other three GC values and its IR length matching exactly.
  - Run `structure` + `genome_stats` over the 18 new and the 26 published plastomes; compare
    region lengths and per-region GC against Table 1 row by row. The 18 new ones already reproduce
    **exactly** (checked by hand 10/10/2026); the published ones have no table to compare against,
    so they get our own range check only.
  - Name what the paper's **text** says and its **table** says where they differ: the text's IR
    range starts at 26,311 bp, Table 1's own minimum is 26,331 (measured), and 26,311 does not
    appear in any of the 44 genomes we measured.
  **Done when:** a §8 entry holds the side-by-side table for the 18 (with the counting rule each
  column uses) and names every row that differs.

- [x] **P8B-04 · SSR vs Table S6** *(M8B.3)* — **DONE 10/10/2026 (§8.68 + §8.73)**: at MISA's
      defaults our table reproduces **every row** of 14 of the 49 sections (the extras are
      compound components). **The full measurement is now in** (§8.73): MISA run at **both**
      settings over all 50 plastomes, one genome assigned per section by their own SSR
      coordinates and read for both settings — **defaults 2,581/3,546 rows exact (72.8 %), 14/49
      sections complete · their stated thresholds 3,070/3,546 (86.6 %), 4/49 complete**; all 14
      complete sections are in the 18-new half. The mechanism is measured row by row: in the ZWSZ
      section 9 of the 9 non-exact rows are a compound cluster their table writes as **one
      spanning row** (`(242,267)` theirs = `(242,252)+(256,267)` ours). **Which run produced each
      section of S6** stays a question for anh Điềm, but it no longer blocks the comparison.
  - Measured 10/10/2026: our `ssr` (MISA defaults 10,6,5,5,5,5) reproduces **57/57** of their
    *C. altaica* rows once the compound components are set aside, while their stated Method
    thresholds (10,5,4,3,3,3) give 97 and reproduce only 55/57 — i.e. **their Methods and their own
    table disagree**. Run the comparison at both settings over all 49 sections of S6 and report both.
  **Done when:** the §8 entry states which parameter setting reproduces S6, at what agreement over
  all genomes, and names the compounds (our table lists components; theirs appears to list none).

- [x] **P8B-05 · RSCU vs Table S5** *(M8B.4)* — **DONE 10/10/2026 (§8.70)**, with an open
      question for anh Điềm. Measured: **max |ΔRSCU| 0.1039 · mean 0.0148 · 1.98 % relative**
      (our rule), **0.1018 · 0.0106 · 1.40 %** counting one copy per gene name. The difference is
      **systematic** (same sign across all 18 samples) and **classed by record**: 74–75-CDS
      records differ 1.55 %, 83–85-CDS records 2.6 % (→ 1.10 % under one-copy counting) — so the
      IR-duplicate counting rule is *part* of it. Eleven candidate rules were tested and each was
      rejected **with its number** (listed in §8.70): no gene-set variant, no added single gene,
      no tRNA/rRNA, no intron span, no whole genome, no averaging rule reproduces it.
  - **Why it is not a data difference, measured**: the same deposited CDS reproduce Table S7
    (Ka/Ks) **to the last printed digit** (§8.69), and our own RSCU side is validated against
    cpstools on the golden (§8.61: 0 codons differ). So the residual is their **counting rule**,
    which S5 does not state (Methods say only "PhyloSuite … screen the PCGs" + an RSCU-Plot web
    tool) — now a row in "Waiting on data / people".

- [x] **P8B-06 · π, hotspots and the barcode region** *(M8B.6)* — **DONE 10/10/2026 (§8.72)**.
      `diversity` (600/200 — their DnaSP setting and our default) over the **44 plastomes**: π
      0.001633 whole, LSC 0.002044 · SSC 0.003533 · IR ≈0.00027. **The hotspot ranking does not
      agree at the top**: their three named regions rank **4th, 5th and 9th of 91 spacers** in our
      measurement (all inside our top 10), while our first is **`petG-trnW-CCA` at 0.02042** — more
      than twice the runner-up. Stated with the scope difference: they ranked by sliding windows
      over whole plastomes, we by per-spacer π over the same 44.
  - Barcode half: **50/50 records reproduce Table S2 exactly** (length, and GC to the printed
    digit — max difference 0.0000); the missing row is `seq45` *C. scabrifolia* YNSZ3, the 51st
    sample, never deposited — named as missing, which is the whole of what is available. π inside
    the locus: **0.018759** over 455 usable columns, against **0.00742** for the same spacer across
    the 44 plastomes — the finer sampling is ~2.5× more diverse.
  - Measured on the way: the Pitfalls row was wrong — `diversity` **and** `barcodes` both require
    the IR gate, so neither accepts the amplicons; π went through the same core function and MAFFT
    adapter the module uses, and the row is corrected in place.

- [x] **P8B-07 · Annotation vs the published GeSeq annotation** *(M8B.7)* — **DONE 10/10/2026
      (§8.71)**: 18/18 leave-one-out runs through `plastid annotate`, classified with
      `tools/corpus_bar.py`'s own rule (imported). **341 conclusion-level differences** (missing
      21 · extra 136 · extent 184) over 4,448 published features against 2,339 called; 288
      boundary, 305 vocabulary, strict 629; 18.9 per record against the Iridaceae bar's 12.3
      (418/34). **Two of the three components are not our reading being wrong**: 110 of the 136
      extras are the IR genes the deposited records lack, and 105 of the 184 extents are six genes
      of **one reference record** (*M. germanica*'s `psbM` is 270 bp where the other 17 say 105) —
      the transfer is faithful to its reference, as designed. The third is ours and is named:
      **H-24**, `trnV-UAC` is never called (18/18 and on the golden too — ARAGORN reads that locus
      as `tRNA-Stop(tta)`, and as `tRNA-Cys(aca)` on the golden), and the run does say so.
  - Counting rules as required: their Table 1 writes "total (unique)"; §6's rule is CDS + tRNA +
    rRNA with both IR copies — both printed where the two meet, in §8.71 and §8.69.

- [x] **P8B-08 · Ka/Ks vs Table S7** *(M8B.5)* — **DONE 10/10/2026 (§8.69)**: the table
      **reproduces**, and the design it actually used is not the one the paper's text states.
      Measured: **KaKs_Calculator 2.0.1 `-m MLWL -c 11`** on MAFFT-aligned CDS pairs, reference
      **`MN068257` *Amelanchier cusickii*** (a record **outside** the study group) — **629/714
      comparable cells match with zero error** to the last printed digit. The text's "relative to
      *C. kansuensis*" is ruled out by measurement, twice: against `MF784433` the values are wholly
      different (rpoC2 *C. monogyna* 0.6667 vs the table's 2.10854), and the table's own structure
      refutes it (rpoB of altaica · GSSZ · ZWSZ is **byte-identical** to `MF784433`'s yet carries a
      value, where that comparison prints `NA`). 77 cells are `0` for an undefined ratio (Ks = 0 →
      tool prints `NA`), 347 are `0` because Ka = 0, 2 differ (clpP *M. germanica* 0.4 %, ycf1
      *M. germanica* 10.8 %) — and 42 are not comparable because **12/18 deposited records lack the
      whole IR gene set**, so **the analysis annotation is not the deposited annotation** (input to
      P8B-07). Version note: their v3.0 vs our 2.0.1 does not show at printed precision.
  - Done when clause met: **reproduced**, with the design question answered by measurement rather
    than left open — and the paper's own text/table disagreement recorded as a finding (D-l).

- [x] **P8B-09 · The record: provenance note + the defects the corpus forced** *(M8B.8)* — **DONE
      10/10/2026**: `docs/corpus_crataegus_2026.md` names the paper (DOI `10.3389/fpls.2026.1767012`,
      CC-BY), the delivery (10/10/2026, `D:\business\pts\thanhdiem_data\db_test_new`), the four
      accession groups, **how to prove a copy is the same** (`sha256(gb_data/manifest.tsv)` =
      `0a34e679…`, and the whole tree in one line = `4ee47625…`), the labels this corpus gets wrong
      (`PX413286` is **chlorosarca in GenBank and jozana in Table S1**) and what it is not (not
      manuscript content; not a benchmark to tune against).
  - **The defects it forced, each fixed test-first or named open with evidence:**
    | defect | where | status |
    |---|---|---|
    | MISA's compound flag missed the per-component encoding (0 of 16 flagged on `PX413282`) | `ssr` | **fixed**, module 0.3.0, §8.65 |
    | `ssr`'s `Context` column read a trans-spliced gene's outer span | `ssr` | **fixed** in Phase 8, §9.M-1 |
    | **`trnV-UAC` is never called** — ARAGORN misreads that locus (`tRNA-Stop(tta)` here, `tRNA-Cys(aca)` on the golden), so an intron-containing tRNA is missing from every run, 18/18 and the golden too | `plastid_annotate` + ARAGORN | **open — H-24** (ROADMAP M3B.5), evidence in §8.71 |
  **Done when clause met:** two fixed with regression tests, one named open with its measurement.

## What the paper states — read 10/10/2026, so nobody opens the PDF again for these

anh Điềm's answer to the six asks was "read the paper", and he was right: four of the six are
answered there, just never where a Methods section would put them. **Read this table before asking
him anything again.**

| Ask | Where it is | What it says |
|---|---|---|
| **Ka/Ks design** | **Results**, *not* Methods — and **the table disagrees with it** *(measured 10/10/2026, §8.69)* | the text says "relative to those of *C. kansuensis*", but S7 **reproduces exactly** (629/714 cells, zero error) against **`MN068257` *Amelanchier cusickii*** with **MLWL** — and *C. kansuensis* is ruled out by the table's own structure. The text's rpoC2 names (2.11 · 1.83 · 1.83) are S7's own cell values ✓ |
| **π** | Methods + Results | window **600 / step 200** (our `diversity` default); the values live in Figure 6 only, **but the text names its three most diverse regions**: `ndhC~trnV-UAC` · `psbZ~trnG-UCC` · `trnR-UCU~atpA` ⇒ the **ranking** is comparable without their numbers |
| **IRscope junctions** | Figure 3 only | the numeric junction comparison goes through **Table 1's region sizes** instead — already done (§8.67) |
| **MISA thresholds** | Methods | **stated**: (10,5,4,3,3,3) with compound 100 — *and their own Table S6 disagrees with them* (defaults reproduce 57/57 rows, their stated thresholds 55/57). An **inconsistency, not a gap** — D-l covers how to report it |
| **51st barcode sample** | Data availability statement | covers `PX676094–PX676143` = **50 records** — *C. scabrifolia* was never deposited ⇒ the barcode half compares **50 of 51**, and that is the whole of what is available |
| **raw reads** | not deposited | anh Điềm's answer: they need the assembly step. **Phase 7's blocker is unchanged by this corpus** |

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Recorded in |
|---|---|---|---|
| D-k | What the Crataegus corpus is to the project | **A test corpus with a published answer key — never manuscript content.** Rosaceae, while the manuscripts are Oldenlandia and Coleus; the M4a rule. Every output is labelled third-party published data | ROADMAP Phase 8B · §8 entries |
| D-l | How to compare when their Methods and their tables disagree | **Report both settings, and say which one reproduces their table.** Measured: their SSR table follows MISA's defaults, not the thresholds their Methods state — matching the *stated* parameters would look like disagreement with a table that is actually reproducible | §8 entry for P8B-04 |
| D-m | Where the corpus lives | **`.data/` (gitignored) + a provenance note in the repo** — **executed 10/10/2026**: `docs/corpus_crataegus_2026.md` (accessions · retrieval date · DOI · `sha256(gb_data/manifest.tsv)` and the whole-tree digest, so a fresh clone can prove a copy is the same). The paper is CC-BY but 22 MB of third-party PDF belongs in the note, not in git — the same shape as `docs/tools/*.md` for tools | `docs/corpus_crataegus_2026.md` |
| D-n | Whether the toolkit should be able to run their *stated* SSR parameters | **Answered by measurement 10/10/2026 — not needed. `ssr` keeps pinning MISA's defaults, and no `--set` is added**: at defaults 14/49 sections reproduce whole (72.8 % of rows), at their stated thresholds 4/49 do (86.6 % of rows) — the extra rows found there are matched by extra rows produced, so the parameter would buy coverage, not reproduction. If it is ever wanted it stays a module change that must keep §6's golden at 37 | §8.73 · §3 if it is decided |
| D-o | How to use the paper's tables when their design is not stated | **Run ours first with every parameter printed and justified; use their table as an answer key to *infer* a design (labelled as inference); never tune a parameter until a table matches** — a setting chosen because it makes a number agree is *fitting, not reproduction*, and the entry must say which one it is. anh Điềm's own suggestion ("chạy output thô ra rồi so") with these two guards | §8 entries · this board |

## Pitfalls — this phase's own traps

| Trap | What to do instead |
|---|---|
| A published table taken as a benchmark | Their Methods and their own Table S6 **disagree** (measured 10/10/2026). A table is data to explain, exactly like CGAS's labels |
| Their counting rule mixed with ours | Table 1 writes "119 (105)" = total (unique); §6's rule is **132 = 86 CDS + 38 tRNA + 8 rRNA**, both IR copies. Print both rules wherever they meet |
| The 50 barcode records fed to a plastome module | They are ~600 bp fragments with **no IR**: the gate stops them by design. That is the gate working — but **measured 10/10/2026 (§8.72): `barcodes` and `diversity` both declare `requires_ir` and every genome in a run must pass, so they refuse the amplicons too**. Measure π through `core.align.nucleotide_diversity` on a MAFFT alignment instead, and say so |
| Their SSR count read as a count | It lists **no compound components** (57 vs our 73 on *C. altaica*); ours lists every component. The difference is a class definition, not a miscount |
| A number in the paper's text quoted as measured | The text's IR minimum (26,311) appears in **no** genome we measured, and not in their own Table 1 (26,331) — quote the table, or the measurement |
| Running the drop's `fetch_gb.py` | It hits NCBI; the files are already here, and this project's tests never touch the network. It stays as provenance, not as a step |
| A corpus number reaching a manuscript table | It is Rosaceae and it is third-party: label every output (D-k), the same rule as a stand-in |

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| The **51st barcode sample** — Table S2's last row is *Crataegus scabrifolia* with GenBank ID `----` | anh Điềm | P8B-06's barcode half |
| *(closed 10/10/2026 — by measurement, not by asking)* the Ka/Ks design: S7 is **KaKs_Calculator MLWL against *Amelanchier cusickii* `MN068257`** (§8.69). What *is* worth telling anh Điềm: their own deposited records **lack the IR genes** while S7 has values for them — **which annotation produced the tables?** | anh Điềm (optional) | P8B-07's reading of the annotation |
| **How their RSCU was counted** — Table S5 is 1–2 % away from ours under every one of 11 rules we tested (§8.70), and S5 states no gene set; with the sequences proven identical (S7 reproduces exactly), only the counting rule can explain it | anh Điềm | P8B-05's residual (reported as measured, not as agreement) |
| **Which annotation produced the tables** — 12/18 deposited records lack the whole IR gene set (`ycf1` · `ycf2` · `ndhB` · `rpl2` · `rpl23` · `rps7`), yet S7 carries `ndhB` values for 11 of them (§8.69) | anh Điềm | P8B-07's reading |
| **How each section of Table S6 was produced** — its two halves encode a compound cluster differently (the published half writes one spanning row, the new half omits clusters). **No longer blocking**: all 49 sections are compared at both settings (§8.73); an answer would let the comparison be *stricter*, and would explain why the published half's rows line up better at the thresholds their own Methods state | anh Điềm | — (P8B-04 closed) |
| *(no longer blocking)* their **π sliding-window numbers** — figures only, and the **ranking** comparison is done without them (§8.72: their three regions rank 4th · 5th · 9th of 91 for us). Their numbers would let the *magnitudes* be compared too, which nothing now waits on | anh Điềm | — (P8B-06 closed) |
| *(closed 10/10/2026)* the IRscope junction table — the paper prints Figure 3 only; the numeric junction comparison is done through **Table 1's region sizes** instead (§8.67) | — | — |
| *(separate ask)* raw Illumina reads for ≥1 of the 18 samples — Table S1 shows 2,786× coverage, so the reads exist | anh Điềm / cô Trang | **Phase 7** (M7.1 · M7.5 · M7.6 — still the project's oldest blocker) |
