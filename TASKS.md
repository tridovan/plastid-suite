# TASKS.md — active tasks

> **Phase 8B is OPEN** on `feat/phase-8b` — the **Crataegus answer-key corpus**: the toolkit's
> numbers beside a **published paper's own result tables**, on an independent Rosaceae corpus the
> team just supplied (anh Điềm's `db_test_new`, 10/10/2026). Phase 8 (tool cross-checks, scale, CI)
> is closed and merged to `main`.
> Milestones: `ROADMAP.md` → Phase 8B · log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

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

- [ ] **P8B-03 · Structure & GC vs Table 1** *(M8B.2)*
  - Run `structure` + `genome_stats` over the 18 new and the 26 published plastomes; compare
    region lengths and per-region GC against Table 1 row by row. The 18 new ones already reproduce
    **exactly** (checked by hand 10/10/2026); the published ones have no table to compare against,
    so they get our own range check only.
  - Name what the paper's **text** says and its **table** says where they differ: the text's IR
    range starts at 26,311 bp, Table 1's own minimum is 26,331 (measured), and 26,311 does not
    appear in any of the 44 genomes we measured.
  **Done when:** a §8 entry holds the side-by-side table for the 18 (with the counting rule each
  column uses) and names every row that differs.

- [ ] **P8B-04 · SSR vs Table S6** *(M8B.3)*
  - Measured 10/10/2026: our `ssr` (MISA defaults 10,6,5,5,5,5) reproduces **57/57** of their
    *C. altaica* rows once the compound components are set aside, while their stated Method
    thresholds (10,5,4,3,3,3) give 97 and reproduce only 55/57 — i.e. **their Methods and their own
    table disagree**. Run the comparison at both settings over all 49 sections of S6 and report both.
  **Done when:** the §8 entry states which parameter setting reproduces S6, at what agreement over
  all genomes, and names the compounds (our table lists components; theirs appears to list none).

- [ ] **P8B-05 · RSCU vs Table S5** *(M8B.4)*
  - Their RSCU came from PhyloSuite's PCG screening; ours from `codon_usage` (all CDS, stops
    excluded, IR copies counted). Expect the §8.61 story: identical counts once both sides count
    one gene set — measure it rather than assume it.
  **Done when:** the §8 entry reports per-codon agreement (or the named cause of each difference)
  for the 18 genomes, with both sides' gene-set rules printed.

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

- [ ] **P8B-08 · Ka/Ks vs Table S7** *(M8B.5)* — **BLOCKED: the comparison design is not in the
      paper or the supplement** (each column is a species; the partner is unstated). Ask anh Điềm
      for the design, or their `.axt`/codon alignments, or the exact `KaKs_Calculator` invocation.
  **Done when:** either the §8 entry reproduces S7 within a stated tolerance, or it records what was
  tried and what the design question is — never a silent skip.

- [ ] **P8B-09 · The record: provenance note + the defects the corpus forced** *(M8B.8)* — a
      `docs/` note naming the accessions, the retrieval date, the paper's DOI and the manifest's
      sha256, so a fresh clone can re-fetch and prove it got the same 100 files (D-m); and P8B-01
      is the first defect on the list. Keep a running
      list in the §8 entries and fix each one test-first; a defect found here is worth more than a
      comparison that agrees, because the corpus is the first one nobody here chose.
  **Done when:** every defect found by the corpus is either fixed with a regression test or named
  as open on this board with its evidence.

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Recorded in |
|---|---|---|---|
| D-k | What the Crataegus corpus is to the project | **A test corpus with a published answer key — never manuscript content.** Rosaceae, while the manuscripts are Oldenlandia and Coleus; the M4a rule. Every output is labelled third-party published data | ROADMAP Phase 8B · §8 entries |
| D-l | How to compare when their Methods and their tables disagree | **Report both settings, and say which one reproduces their table.** Measured: their SSR table follows MISA's defaults, not the thresholds their Methods state — matching the *stated* parameters would look like disagreement with a table that is actually reproducible | §8 entry for P8B-04 |
| D-m | Where the corpus lives | **`.data/` (gitignored) + a provenance note in the repo** (accessions · retrieval date · DOI · manifest sha256). The paper is CC-BY but 22 MB of third-party PDF belongs in the note, not in git — the same shape as `docs/tools/*.md` for tools | README of the drop · `docs/` |
| D-n | Whether the toolkit should be able to run their *stated* SSR parameters | **Yes, but as a parameter, not a code change**: `ssr` pins MISA's defaults today. A `--set` for thresholds is a module change that must keep §6's golden at 37 — decide it when P8B-04 shows how many of the 49 sections need it | §3 · P8B-04 |

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
| The **Ka/Ks comparison design** — which sequence each species column of Table S7 is compared against (or their `.axt` files / the `KaKs_Calculator` invocation) | anh Điềm | P8B-08 |
| *(if available)* their **π sliding-window numbers** and the **IRscope junction table** — both are figures today | anh Điềm | P8B-06 · P8B-03 (numeric rather than figure comparison) |
| *(separate ask)* raw Illumina reads for ≥1 of the 18 samples — Table S1 shows 2,786× coverage, so the reads exist | anh Điềm / cô Trang | **Phase 7** (M7.1 · M7.5 · M7.6 — still the project's oldest blocker) |
