# TASKS.md — active tasks

> **Phase 8B is OPEN** on `feat/phase-8b` — the **Crataegus answer-key corpus**: the toolkit's
> numbers beside a **published paper's own result tables**, on an independent Rosaceae corpus the
> team just supplied (anh Điềm's `db_test_new`, 10/10/2026). Phase 8 (tool cross-checks, scale, CI)
> is closed and merged to `main`.
> Milestones: `ROADMAP.md` → Phase 8B · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

> **START HERE → P8B-05 (RSCU vs Table S5).** P8B-08 closed 10/10/2026 (§8.69) by
> reproducing S7 exactly — *Amelanchier cusickii* + MLWL, not *C. kansuensis* as the text
> says. Then **P8B-06** (π hotspot + the 50 barcode records) → **P8B-07** (annotation vs
> GeSeq — the biggest) → **P8B-04**'s remaining 35 sections, which waits on anh Điềm. Read
> "What the paper states" below before asking him anything — and note P8B-08 measured that
> **12/18 deposited records lack the IR genes** while the paper's tables analyse them, which
> P8B-07 must read as its first fact, not as a surprise.

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

- [~] **P8B-04 · SSR vs Table S6** *(M8B.3)* — **partly done 10/10/2026 (§8.68)**: at MISA's
      defaults our table reproduces **every row** of 14 of the 49 sections (the extras are
      compound components); the other 35 cannot be compared yet because **their table encodes
      a compound cluster as one spanning row** in the published half and **omits clusters** in
      the new half. **Blocked on:** which run produced each section of S6 (see Waiting).
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

- [ ] **P8B-06 · π, hotspots and the barcode region** *(M8B.6)*
  - Their sliding-window analysis is DnaSP v6 at **600/200** — our `diversity` default — over whole
    plastomes, PCGs and intergenic regions; their hotspot is *ndhC~trnV-UAC*, amplified for 51
    individuals. Run `diversity` on the 44 plastomes and `barcodes` on the same set, then compare
    **the hotspot choice** (do we rank the same region highest?) and the 50 deposited barcode
    records (length/GC vs Table S2, then π within the locus).
  - **Needs:** the 51st sample (*Crataegus scabrifolia*, no GenBank ID in Table S2) before the
    barcode half can be complete.
  **Done when:** the §8 entry names whether our hotspot ranking agrees, and reports the barcode
  comparison for the 50 records we have, with the 51st named as missing.

- [ ] **P8B-07 · Annotation vs the published GeSeq annotation** *(M8B.7)*
  - Their 18 records were annotated by **GeSeq** (with Prodigal/HMMER/aragorn behind it) and are
    deposited in GenBank — so the corpus carries an answer key for the annotator, in a family ours
    has never seen (phase 3B's bar was Iridaceae). Run `plastid annotate` leave-one-out over the
    corpus and classify every difference the way `tools/corpus_bar.py` does.
  - Their counting rule is **not** ours: Table 1 writes "119 (105)" = total (unique), while §6's
    rule is CDS + tRNA + rRNA with both IR copies. Both must be printed wherever they meet.
  **Done when:** a §8 entry holds the bar (per-record differences vs the published annotation,
  classified) on the Crataegus corpus, beside the Iridaceae bar's 418-vs-PGA-588 result — same
  measure, different family, stated as such.

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

- [ ] **P8B-09 · The record: provenance note + the defects the corpus forced** *(M8B.8)* — a
      `docs/` note naming the accessions, the retrieval date, the paper's DOI and the manifest's
      sha256, so a fresh clone can re-fetch and prove it got the same 100 files (D-m); and P8B-01
      is the first defect on the list. Keep a running
      list in the §8 entries and fix each one test-first; a defect found here is worth more than a
      comparison that agrees, because the corpus is the first one nobody here chose.
  **Done when:** every defect found by the corpus is either fixed with a regression test or named
  as open on this board with its evidence.

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
| D-m | Where the corpus lives | **`.data/` (gitignored) + a provenance note in the repo** (accessions · retrieval date · DOI · manifest sha256). The paper is CC-BY but 22 MB of third-party PDF belongs in the note, not in git — the same shape as `docs/tools/*.md` for tools | README of the drop · `docs/` |
| D-n | Whether the toolkit should be able to run their *stated* SSR parameters | **Yes, but as a parameter, not a code change**: `ssr` pins MISA's defaults today. A `--set` for thresholds is a module change that must keep §6's golden at 37 — decide it when P8B-04 shows how many of the 49 sections need it | §3 · P8B-04 |
| D-o | How to use the paper's tables when their design is not stated | **Run ours first with every parameter printed and justified; use their table as an answer key to *infer* a design (labelled as inference); never tune a parameter until a table matches** — a setting chosen because it makes a number agree is *fitting, not reproduction*, and the entry must say which one it is. anh Điềm's own suggestion ("chạy output thô ra rồi so") with these two guards | §8 entries · this board |

## Pitfalls — this phase's own traps

| Trap | What to do instead |
|---|---|
| A published table taken as a benchmark | Their Methods and their own Table S6 **disagree** (measured 10/10/2026). A table is data to explain, exactly like CGAS's labels |
| Their counting rule mixed with ours | Table 1 writes "119 (105)" = total (unique); §6's rule is **132 = 86 CDS + 38 tRNA + 8 rRNA**, both IR copies. Print both rules wherever they meet |
| The 50 barcode records fed to a plastome module | They are ~600 bp fragments with **no IR**: the gate stops them by design. That is the gate working — use `barcodes`/`diversity`, never `structure` |
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
| **How each section of Table S6 was produced** — its two halves encode a compound cluster differently (the published half writes one spanning row, the new half omits clusters), and the published half carries SSRs our MISA does not find; the measured suspicion is that those sections were copied from each record's own publication | anh Điềm | P8B-04's remaining 35 sections |
| *(reduced 10/10/2026)* their **π sliding-window numbers** — still figures only, **but the paper names its three most diverse regions** (ndhC~trnV-UAC · psbZ~trnG-UCC · trnR-UCU~atpA), so the hotspot **ranking** can be compared without them | anh Điềm | P8B-06 (values) |
| *(closed 10/10/2026)* the IRscope junction table — the paper prints Figure 3 only; the numeric junction comparison is done through **Table 1's region sizes** instead (§8.67) | — | — |
| *(separate ask)* raw Illumina reads for ≥1 of the 18 samples — Table S1 shows 2,786× coverage, so the reads exist | anh Điềm / cô Trang | **Phase 7** (M7.1 · M7.5 · M7.6 — still the project's oldest blocker) |
