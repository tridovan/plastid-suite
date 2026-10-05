# HANDOFF — archive: sessions 1–13 (02/10 – 03/10/2026)

**What this is.** The oldest part of `HANDOFF.md`'s session log, moved here 05/10/2026 — the same
reason §8 moved to `docs/log/` on 04/10: `HANDOFF.md` is the first file every session reads, and
it had reached 61 KB / 813 lines, most of it history no session needs to start.

**This is a pure move.** The kept part of `HANDOFF.md` plus this file is byte-identical to the
file it was split from — no entry retyped, renumbered or trimmed. Sessions keep their numbers
("session 5" is still "session 5"), so a citation like "HANDOFF session 5" resolves here.

---

### 2026-10-03 · Claude · dev WSL · session 13 (P3A-03 — the mito gene caller, measured on the answer key)

**Did.** Built the caller that transfers an annotation from the pinned mito reference set:
`core/mito_caller.py` (pure: chaining, exon merging, floors, refinement), `external/blast.py`
(BLAST+ through `run_tool`, with the outfmt-6 parser beside the argv that produces it),
`modules/mito_annotate.py` (wiring + six tables + the GenBank artifact), `defaults.yaml`'s `mito:`
section, the report template, `docs/tools/blast.md`, and `environment.yml` repinned to
`blast=2.17.0`. 47 unit + 4 tier-B tests, written test-first.

**Validation.** `make check` green (424 fast + 30 golden, 8 tier-B skips with named reasons). Real
run on `OX459128.1` with BLAST+ 2.17.0: **31 CDS + 3 rRNA called, nothing not called**; the golden
test pairs every call with the published CDS its exons overlap (bijection over 31) and prints
**18/31 spans byte-identical, 13 differing by ≤ 75 bp, strand 31/31**, ORF profile
**21 ok · 3 no_start · 7 no_stop**.

**Five measurements changed the design** (all in PROJECT.md §8.22 / §9.D, each with a test):
BLAST truncates `sseqid` at the first whitespace; the committed reference FASTA holds 8 duplicate
headers — 6 identical repeat copies and **2 symbols carrying two different sequences**
(Arabidopsis `atp6` 385/349 aa), so an HSP belongs to the smallest reference length containing it
and two sequences under one key are never chained; `-seg no` (SEG halves atp9's alignment,
137 → 67 bits); `-parse_seqids` **fails** on our reference set and upper-cases ids; and ranking
candidates by bitscore **builds chimeras** (`atp6` pulled 381 bp upstream onto a 75.6 % alignment)
⇒ the winner is the highest-identity candidate that passes the floors.

**A bug outside this module:** `run_tool`'s version probe tried `--version`/`-v`/`-V`, but BLAST+
accepts **only** `-version` — so the tool the module depends on was recorded as "unknown" while
every probe "succeeded" in failing. `-version` added, with a regression test.

**Residual risks.** (1) Two published exons of nad1/nad5 (19 and 7 codons) are not recoverable from
any alignment; the caller reports 4 exons and says so rather than inventing a fifth — P3A-06 owns
listing that difference. (2) The caller cannot know a gene is a pseudogene (the answer key marks
4), so its ORF profile reads 21 ok · 3 no_start · 7 no_stop against the record's 21 · 4 · 3 · 3.
(3) The floors are calibrated on one record; a divergent sample is exactly what P3A-07
(`AP045469.1`) will test. (4) BLAST version mismatch is a hard refusal (the Chloe precedent) — if
that proves too strict for anh Điềm's machine it must be decided in PROJECT.md §3, not loosened
quietly.

**Blocked.** Nothing.

**Next.** P3A-04 — ARAGORN 1.2.41 for tRNA (bioconda; `environment.yml` + `docs/tools/aragorn.md`
in the same session), then P3A-05's CLI (`plastid mito annotate`), P3A-06's itemized golden check
(31 · 17 · 3) and P3A-07 (`AP045469.1`, then pin it).

### 2026-10-03 · Claude · dev WSL · session 12 (Phase 3A opened — a mito annotator we build ourselves)

**Did.** Opened **Phase 3A** (`ROADMAP.md`, `TASKS.md`, branch **`feat/phase-3a`** off
`feat/phase-3`) after measuring what a self-built mito annotator would cost. Then **P3A-01** and
**P3A-02**.

* **Why the phase exists.** P3-04's count check cannot close by either PMGA route: the container is
  11 GB unpacked driving MAKER against an 8 GB WSL cap, and PMGA's web service refuses connections.
  CLAUDE.md's rule for a tool that cannot run is *write that part ourselves* — and the estimate is
  **under 10 MB on disk and seconds of runtime**, against PMGA's 14 GB.
* **P3A-01 · GenBank writer — DONE (§8.20).** `to_genbank()` beside the reader. Round trip on the
  golden: **265/265** features, parts in transcript order, sequence and source intact; 324 unit
  tests green, so the `.tbl` path verified against NCBI's own tables is unchanged. Two findings:
  `/pseudo` and `/trans_splicing` are **bare** qualifiers, and BioPython writes `=""` for an empty
  string while reading its own output back happily (§9.B-26 — a round trip through the same library
  is not evidence); and the qualifier builder moved to **`model.gene_qualifiers()`** so the table
  and GenBank writers cannot disagree.
* **P3A-02 · Reference set — DONE (§8.21).** Eight RefSeq mitogenomes pinned by the panel's own
  taxonomy rule (three Lamiaceae for *Coleus*, three Rubiaceae for *Oldenlandia*, plus Arabidopsis
  · Nicotiana · Oryza), 5.2 MB **not** committed; **287 proteins + 25 rRNAs (149 KB)** are, and a
  test rebuilds them from the pinned panel and compares. Coverage **31/31**.
  Two measurements changed the design: **self-translating is impossible** — all 530 reference CDS
  publish a `/translation`, none declares a `/transl_table`, and plant mito CDS begin with `ACG`
  which no table accepts as an initiator, so CDS-mode translation raises; and **the references are
  polluted** — 236 anonymous ORFs annotated as CDS (121 in *Nicotiana* alone), plastid-derived
  `rbcL`/`petG`/`atpB` inside the mitogenome, and a *Psychotria* "5S rRNA" that is 722 bp.

**Validation.** `make check` green after each task (362 unit + 30 golden; tier-B skips name their
missing tool). Every claim above is a measurement on the pinned records, not a reading of a paper.

**Residual risks.** (1) **The vocabulary is narrower than the gene set** — 23 reference genes have
no entry, ten of them standard (`rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` ·
`rps14` · `rps19` · `sdh3`). They are kept in the reference with `-` so a run can say it found one;
named as **H-06** and deliberately *not* fixed here. (2) The writer is **not** a lossless GenBank
editor — `locus_tag`, `note`, `protein_id`, `translation` and `misc_feature` are outside the model
and stated as such in its docstring. (3) The golden round trip has **no partial coordinates** to
exercise; the partial path is covered by a synthetic test only.

**Blocked.** Nothing in Phase 3A. P3-04's count check stays open (H-05) and is now *also* the thing
the phase is trying to make unnecessary.

**Next.** **P3A-03 — the gene caller** (`blastx` → best hit per gene → codon boundaries → `Gene`;
`blastn` for rRNA), then P3A-04 (ARAGORN), P3A-05 (CLI), P3A-06 (the 31·17·3 golden check),
P3A-07 (`AP045469.1`, then pin it). All measured facts and the pitfalls are in `TASKS.md`; the
design decision that matters most is **D3A-02 — boundaries come from the reference alignment, never
from ORF-walking**, because the answer key's own 31 CDS score 21 ok · 4 pseudogene · 3 no_stop ·
3 no_start.

### 2026-10-03 · Claude · dev WSL · session 11 (Phase 3 closed)

- **Did:** closed P3-04 by decision — the mito annotator is PMGA on its web service, pinned
  (`plastid panel pin --file`, §8.19), not a local adapter, because the machine cannot run the
  annotator (11 GB unpacked + MAKER vs an 8 GB WSL cap and 42 GB free on C:). Deleted the 14 GB
  container after reading its real CLI out of it. Board, ROADMAP M3.1, PROJECT.md §3 and
  `docs/tools/pmga.md` all say so, with the trigger that would reopen the adapter.
- **Validation:** `make check` green — 345 fast + 30 golden. Phase 3's own check (§13.3): every
  run's output re-reads with our reader, the ORF checks pass, the gene recount uses the fixed rule
  — all three hold for both the plastome and the mito path.
- **Risks / notes:** (1) the mito **count check** against a real PMGA run has still never been
  made — it needs one web-service run by the team; until then `mito_genes` is verified against a
  downloaded record, not against an annotator's output; (2) the disk: deleting inside WSL returns
  nothing to Windows (`ext4.vhdx` does not shrink; WSL refuses its own sparse mode) — a manual
  `diskpart compact vdisk` from Windows is what gives the 14 GB back.
- **Blocked:** nothing. Open for the PI: D3-06 (which gene name is printed in the manuscript) and
  the `trnG` naming question (§10.10).
- **Next:** Phase 4 — M4.1 Ka/Ks · M4.2 %identity + synteny · M4.3 barcodes · M4.4 population
  genetics · M4.6 Ts/Tv (**scope fixed: CDS only · 76 shared genes · MAFFT · IR excluded**).
  M4.5 (RNA editing) needs RNA data nobody has yet. Note MAFFT is not on this machine — Phase 4's
  golden checks need it installed and pinned (`docs`/§7.8, measured earlier at 7.526).

### 2026-10-03 · Claude · dev WSL · session 10 (PMGA downloaded, and the backend rule completed)

- **Did:** downloaded and unpacked the PMGA container (Figshare 27201798), read its real command
  line out of its own POD, and closed the half of the plan's backend rule that was missing:
  `plastid panel pin --file/--note` (§8.19), tested, so a tool that ran on a web service or
  another machine can be pinned and consumed offline. `docs/tools/pmga.md`, TASKS.md P3-04,
  PROJECT.md §3 + §8.19.
- **Validation:** `make check` green — **345 fast + 30 golden**. `panel pin --file` measured on the
  stand-in mitogenome: 381,931 bytes · sha256 `e31da1d3…` · 258,274 bp · 31 CDS + 17 tRNA +
  3 rRNA · circular (the §8.14 numbers).
- **Risks / notes:** (1) the container is 11 GB unpacked and made the WSL VM balloon toward its
  8 GB cap — the user saw `vmmemWSL` at 4 GB after the extraction; dropping the page cache inside
  the guest does **not** shrink the VM (that takes `wsl --shutdown` or idle time with
  `autoMemoryReclaim=gradual`). (2) PMGA is still unrun: peak RAM and wall clock are unmeasured,
  and that is the number that decides whether it can ever be a requirement rather than an option.
- **Blocked:** nothing. The container is on this machine and can be deleted whenever the disk or
  the VM footprint is unwanted — the web-service route does not need it.
- **Next:** P3-04 via the pin route (pin a PMGA/GeSeq mito annotation, run `mito_genes`), or
  Phase 4's comparison modules.

### 2026-10-03 · Claude · dev WSL · session 9 (P3-05's `.sqn` half — run against NCBI's own converter)

- **Did:** ran `table2asn` for real (a stand-in `.sbt` is vendored at
  `tests/data/standin_template.sbt`), which found that our table was **not submittable**: fixed the
  missing `source` feature, the dropped `/exception`, and the partial-start rule for a CDS whose
  codon 1 is not an initiator. `Gene` gained `/exception`, `Genome` gained the record's source
  qualifiers, the module reads NCBI's `.val` into its tables, and a tier-B test runs the whole
  `.sqn` path on a machine that has the tool.
- **Validation:** `make check` green — **343 fast + 30 golden**, plus the real-tool tier-B tests
  (`PTS_TABLE2ASN`, `PTS_CHLOE`, `PTS_JULIA`). Measured: NCBI's validator on the golden table
  **13 errors → 7**, the remainder being the class a real accepted plastome submission also has
  (control run: *Dioscorea alata*'s published table → 9 errors, 7 of them the same class).
- **Risks / notes:** (1) the partial-start rule is applied by the writer and named in the report —
  it adds a marker, never moves a coordinate, but it *is* a claim, so a genome whose CDS start is
  wrong for another reason will look "partial" rather than "broken" until someone reads the note;
  (2) the 7 remaining errors are inherent to duplicated IR genes in this format and cannot be
  removed from our side; (3) the stand-in template must never be used for a real submission.
- **Blocked:** nothing. The **team's own `.sbt`** is still needed *to submit* (NCBI's form, cô
  Trang's details) — that is now a submission prerequisite, not a code one.
- **Next:** P3-04's mito annotator (PMGA container — record its real entry point), or Phase 4's
  comparison modules, which can key on the canonical names.

### 2026-10-03 · Claude · dev WSL · session 8 (Chloe adapter — P3-03 closes)

- **Did:** `core/gff.py` (GFF3 reader), `external/chloe.py`, `modules/annotation.py` (annotator
  selection: `chloe` default · `pga` for non-angiosperms), `tests/unit/test_gff.py`,
  `tests/unit/test_chloe.py`, two tier-B tests, a real GFF fixture + provenance, `docs/tools/chloe.md`,
  D3-01 updated, PROJECT.md §8.17 · §9.B.21-22 · §10.10 · §3.
- **Validation:** `make check` green — **333 fast + 30 golden** + the two Chloe tier-B tests run for
  real (`PTS_CHLOE`/`PTS_JULIA`). Measured through our reader on Chloe's own file:
  **132 = 86 + 38 + 8**, ORF **80 `ok` + 6 `alternative_start`, 0 failures**; end-to-end CLI 22 s,
  manifest carries `chloe 0.1.18` · `julia version 1.13.0` · the full argv with `--no-transform`.
- **Risks / notes:** (1) the transform guard is a FASTA comparison — a future Chloe that transforms
  without writing that file would slip past it (the release pin and the reference-set pin are the
  second line); (2) `--no-transform` is measured on one genome; (3) `core/gff.py` is ours because
  BioPython has no GFF3 parser — strict, but it has met one producer's output shape.
- **Blocked:** nothing. Open for the PI: the `trnG` name (§10.10). Licence review (M9.3): the
  checkout has no LICENSE file, the Zenodo deposit declares MIT.
- **Next:** P3-04's annotator (PMGA container — its entry point must be recorded, not invented) ·
  P3-05's `.sqn` (needs the team's `.sbt`) · Phase 4 comparison modules, which can now key on the
  canonical names.

### 2026-10-03 · Claude · dev WSL · session 7 (Phase 3: P3-01 · P3-02 · P3-05 · P3-03)

- **Did:** P3-01 (ORF checks), P3-02 (gene-name normalization + the whole printed form),
  P3-05's Python half (`.fsa` + 5-column `.tbl`, artifacts through the Result/export/cache path),
  P3-03 (PGA adapter + annotation module, gated). P3-04's tool decided (PMGA) but not implemented.
- **Validation:** `make check` green — **276 fast + 30 golden + 1 skipped** (the PGA test, with a
  named reason). Real CLI runs produced, for the golden genome: the ORF report (0 failures, 5
  alternative starts), the 11-record panel name report (155 → 115, 0 unresolved, both bare tRNA
  names resolved), and `submission.fsa` + `submission.tbl` (265/265 features round-trip, cache
  restores them byte-identically).
- **Six bugs found by the work, each now a regression test:** (1) a `gene` feature named `trnK`
  spans the intron, so resolving bare names per *feature* compared an intron against a tRNA and
  reported a false "unresolved" — resolution now runs from the tRNA feature only; (2) the CLI
  refused more than `arity` genomes, contradicting API v0.2 and making an 11-record panel
  impossible; (3) `annotated-not-coding` matched by position alone and `ndhF`'s CDS "covered" the
  truncated ycf1; (4) the TBL writer put `>` after the number (`1050>`) where NCBI writes it as a
  prefix (`>1050`); (5) a table key not present in its own columns (caught by the template fixture);
  (6) `name_mappings` keyed on a column that did not exist.
- **Residual risks:** PGA has never run here, so `external/pga.py`'s argv is from PGA's README and
  unverified; the `.tbl`'s partial/qualifier layout is per NCBI's documentation and their verbatim
  examples, not against a real submitted `.tbl`; the `.sqn` path has never run (no `table2asn`, no
  `.sbt`); template prose cannot contain digits at all, so gene names with digits must come from
  table cells.
- **Blocked:** P3-04 on the PMGA container; P3-03's measured half on a machine with PGA; the `.sqn`
  half on a submitter template (`.sbt`) — `table2asn` itself is installed and declared.
- **Late in the session, three more things landed:** the **mito path** (`modules/mito_genes.py`:
  a plant mitogenome described with no IR gate, 51 = 31 + 17 + 3 on the stand-in, every CDS named
  by a symbol); **`table2asn` 1.29.324 installed and declared** (`docs/tools/table2asn.md`, md5
  recorded) with both of its error paths tested; and **three real NCBI-produced `.tbl` files
  vendored** (`tests/data/real_*.tbl`), which closed the "field for field" clause for everything
  but joins — and which corrected the writer (three fields per interval line, not five) and the
  parser (a table file may hold several records; one VADR file carried eight and they had been
  merged silently).
- **Next:** see the Resume point.

### 2026-10-03 · Claude · dev WSL · session 6 (blocker clearance)

- **Did:** cleared what could be cleared before Phase 3 starts. Installed **MAFFT 7.526** from
  conda-forge and verified it byte-for-byte against §7.8's pin, plus BLAST+ 2.17.0. Found that two
  named blockers were not blockers, and that one tool in the plan is the wrong kingdom. Wrote
  `docs/PHASE3_BLOCKERS.md`.
- **Validation:** the MAFFT tarball is 2,610,460 B with md5 `cf3b2acf649bc31e576ed43a2fa4aaa3` —
  both numbers are §7.8's, so this is the pinned build and not merely a 7.526. `environment.yml`
  dry-runs clean and the solver picks `mafft 7.526 h4bc722e_0` from conda-forge rather than
  bioconda's 7.525.
- **A mistake I made and fixed:** I installed tools on this machine and changed nothing in the repo.
  Every number measured here would have been reproducible by nobody. Worse, the repo's own
  `environment.yml` listed a **bare `mafft`** while §9.B-14 had already measured that the aligner
  changes Ts/Tv (2.667 vs 2.686) — so following our own install file produced a different golden
  number. Both fixed; the rule is now an invariant in CLAUDE.md.
- **Residual risks:** `table2asn` is only partly downloaded (NCBI's FTP timed out twice). PGA is
  cloned nowhere useful — running freshly-cloned third-party code was refused by the permission
  classifier, correctly, since no one had named the source in-session; one manual `perl PGA.pl`
  clears it. PMGA's benchmark numbers come from a web summary, **not the paper**.
- **Blocked:** `docs/PHASE3_BLOCKERS.md` — four asks, each with a named owner.
- **Next:** Phase 3 A→Z. Start with P3-01 and P3-02 (both pure Python, nothing blocks them).

### 2026-10-03 · Claude · dev WSL · session 5 (Phase 2 closed, Phase 3 board)

- **Did:** pushed `feat/phase-2`; opened branch `feat/phase-3` (stacked — Phase 2 is not merged to
  `main` yet); cleared `TASKS.md` of Phase 2 and built the Phase 3 board from the milestones
  (M3.1–M3.4) with the `make-tasks` skill: 5 tasks, each with a measurable Done-when.
- **Validation:** the board's numbers are measured, not guessed — the whole dataset was swept for
  start codons (`/transl_table=11` and `/codon_start=1` on all 1,099 CDS; starts ATG 746 · GTG 11 ·
  ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1) and for gene names the grouping file does not claim
  (exactly four, all in the two manuscript-species records). `TASKS.md` is 108 lines.
- **Corrected a wrong phase check:** ROADMAP's Phase 3 check read "gene recount uses the fixed
  counting rule (**132**)". 132 is *Iris domestica*'s value, not a constant — P2-10 measured 131 and
  133 for two other species. The *rule* is fixed; the *number* is not. Left uncorrected this is the
  "133 genes" trap in a new costume.
- **Did (cont.):** **D3-02 decided** — the four unusual gene names are the same loci as
  `ycf3` · `ycf4` · `psbN` · `clpP`, proven three ways (syntenic neighbours, identical length — `pbf1`
  and `psbN` both exactly 132 bp — and `/product`). PROJECT.md §8.9.
- **Did (cont.):** the user identified `Vu,H.-T.` as **cô Trang**, which made §8.9's conclusion
  incomplete and flipped its main recommendation. Measured afterwards: **two of the eleven records are
  the group's own** (`PX923121.1` *C. amboinicus* — the manuscript species — and `PX977793.1`); the
  naming difference is a **whole dialect** (tRNA hyphen/bracket/**bare**, rRNA `rrn16`/`rrn16S`, protein
  classic/functional), not four genes; and a bare tRNA name resolves by **sequence identity**, not by
  reading an anticodon at a fixed offset (offset 34 matched only 106/293). The comparison key moved to
  the **team's own names**, because normalising to the classic ones would make the team's paper print a
  name that disagrees with the team's own submitted record. PROJECT.md §8.10.
- **Residual risks:** P3-03/P3-04 (PGA, MITOS2) have no tool on this machine, so they can only be
  implemented and `PTS_REQUIRE_TOOLS=1`-gated — the same honest position P2-07/P2-08 are in.
- **Process note worth keeping:** I first routed D3-02 to the PI *because I had no evidence*, then
  kept it there after measuring — which is the wrong half to keep. The mapping is a dev call with
  evidence; only the name **printed in the manuscript / submitted to NCBI** (D3-06) is the PI's.
  A dev-decidable item parked on the PI blocks silently until someone asks.
- **Blocked:** PGA/MITOS2 (anh Điềm's toolbox) block the measured halves of P3-03/P3-04. A
  confirmation of whether `Vu,H.-T.` — submitter of the two records that use the new names — is this
  group would settle whether those names are ours to explain.
- **Next:** P3-01.

### 2026-10-03 · Claude · dev WSL · session 4 (P2-09, D2-04, P2-10)

- **Did:** finished Phase 2. **P2-09** — a categorised publication gene table
  (`data/gene_groups.yaml` is the grouping, `modules/gene_table.py` builds it) plus a Word export
  (`result/docx_export.py`, `plastid run … --docx`) that is **byte-reproducible** and pruned of the
  template's dead weight. **D2-04** — `python-docx` as an optional extra, imported inside the
  function, so the minimum install stays biopython + PyYAML. **P2-10** — `-m` is repeatable, so a
  run holds several modules and writes one combined `report/report.docx`; draft produced for both
  manuscript species from the stand-in set.
- **Validation:** `make check` green — **208 fast + 28 golden**. Every document was rendered
  (LibreOffice → PDF → PNG) and **looked at**, not just asserted on. Measured: two builds of one run
  are byte-identical, and so is a fully cached re-run (including the merged document).
- **Four real bugs found by this work, each now a regression test:**
  1. the gene table's own reconciliation reported **0** on its Total row — `by_kind.get("Total")`
     on a Counter holding only CDS/tRNA/rRNA returns the default. The column caught it on the
     first run, which is the whole reason it exists (PROJECT.md §9.B-20).
  2. the OOXML style pruner read styles as **attributes** (`element.get("w:pStyle")`) when they are
     **child elements** (`<w:pStyle w:val="…"/>`), so it removed nearly every style the document
     used — including `Title` and `Heading 1`. **LibreOffice still rendered something plausible**;
     only a structural test ("every style the document references must still be defined") caught it
     (PROJECT.md §9.C-34).
  3. a table's own NUL-looking allocation bug: scaling every column by one factor starved the short
     ones, so a header broke mid-word ("Gene / s"); columns now get a floor before the rest is
     shared out.
  4. `Table.emphasis` (a bold totals row) accepted any string — a typo would have been invisible in
     every export except the Word document's font weight. It now must name an existing key.
- **Residual risks:** the Word layout is verified by eye against LibreOffice; Word itself is not
  installed here, so a Word-only rendering difference would not have been seen. The stand-in drafts
  are labelled but are **public records of the same two species** — same species is not the same
  specimen.
- **Blocked:** P2-07/P2-08 numbers still need MAFFT 7.526 (absent).
- **Next:** Phase 3, or collect the three open sub-items in the resume point.

### 2026-10-03 · Claude · dev WSL · session 3 (P2-03 … P2-08, P2-11, D2-01/02/03)

- **Did:** completed seven Phase-2 tasks. `core/codon.py` + `modules/codon_usage.py` (RSCU · ENC ·
  CAI · CBI · Fop), `core/protein.py` + `modules/amino_acid.py` (GRAVY · aromaticity),
  `modules/introns.py`, SSR by region + motif, `modules/diversity.py` (π) and
  `modules/substitutions.py` (Ts/Tv), plus group labels on comparison output.
  D2-01 (multi-genome contract, `docs/API_v0.2.md`), D2-02 (π 600/200), D2-03 (NCBI table 1).
- **Validation:** `make check` green — **169 fast + 22 golden**. Every new number is either a
  §6 golden number or an analytic check (ENC = 61 unbiased / 20 at the extreme; CAI = Fop = CBI = 1
  for a self-referenced sequence; π of a known alignment; the two IR copies' skew mirrors exactly).
- **Six real bugs found by the new work, each now a regression test:**
  1. `_parse_gff` read only `monomeric_repeat` — MISA writes every longer unit as `microsatellite`
     and a compound as one `repeat_region` line, so **7 of 37 SSRs were silently dropped**. The old
     test fixture invented a type MISA never writes, which is why it passed. P2-06's reconciliation
     column caught it.
  2. `coding_sequence` — re-sorting a trans-spliced gene's parts or reading them all forward gives
     a wrong ORF; only per-part strand + annotation order yields rps12's `ATG…TAA`.
  3. Two independent copies of a `Table` key were allowed (`rpl2.1` twice) — references resolved
     silently to one row. `Table` now rejects a repeated key.
  4. `load_groups` let `yaml.YAMLError` through (not a `ValueError`) — a malformed file crashed
     with a traceback instead of a message.
  5. `NcbiClient.resolve` could leak a `KeyError` where `NcbiError` belongs.
  6. `entry_from_summary` claimed `annotated=True` for every record; GenBank submissions are not
     annotated by default — 3 of the 11 stand-in plastomes have zero gene features.
- **Residual risks:** the comparison modules are unit-verified but **no real-data number was
  measured** — MAFFT is absent and installing a different version than the pinned 7.526 would
  measure the wrong thing. π per intron/IGS is still missing from P2-08.
- **Blocked:** nothing. P2-10 still needs the team's real samples.
- **Next:** P2-09 (see the binary-export note in the resume point), then P2-10.

### 2026-10-03 · Claude · dev WSL · session 2 (stand-in data, P2-01, P2-02)

- **Did:** **D-04** — a public stand-in dataset for Oldenlandia/Coleus: `data/*.yaml` definitions
  (13 records), `panel pin --groups` (new) + `--index` now optional + live accession resolution via
  the new `NcbiClient.resolve`; `PinnedGenome` measures gene counts · topology · IR outcome at pin
  time, so `PINNED.json` *is* the verification table. Hardened what the first cut got wrong
  (all found by a review pass): malformed YAML raised a traceback (`yaml.YAMLError` is not a
  `ValueError`), `resolve` could leak a `KeyError`, `annotated=True` was assumed for GenBank
  records. **P2-01** done (e2e tests: golden numbers in report *and* CSV · 2nd golden genome ·
  both diverged-IR genomes exit 2). **P2-02** done — `modules/structure.py` + `gc_skew()` +
  `structure.en.md` + `tests/unit/test_structure.py`.
- **Validation:** `make check` green — **124 fast + 18 golden**. All 13 sha256 re-verified
  byte-identical on a second pin. Structure on the golden genome reproduces LSC 83,127 ·
  IRa = IRb = 26,214 · SSC 18,169 · tiling 153,724, and the two IR copies skew exactly opposite
  (reverse complements) — a real correctness signal for the sign convention.
- **Residual risks:** `PZ428230.1` is a bad record (15,756 bp IR, `linear`) — excluded from
  position-based analysis, and nothing in the tool would have caught it automatically (§9.A-8).
  The API doc said `table[rowkey].column` but the implementation parses `table.rowkey.column`;
  no template had used a row key until now — docs corrected to the tested syntax, code unchanged.
  3 of 11 stand-in plastomes have no annotation (0 CDS) — sequence-level analyses only.
- **Blocked:** nothing. P2-10 still needs the team's real samples.
- **Next:** P2-03 (decide D2-03 first), then P2-04 … P2-11.

### 2026-10-03 · Claude · dev WSL

- **Did:** finished Phase 1. P1-10 complete: `core/panel.py` (4 tiers, modes A/B, quality filter,
  `pin()` with an injected fetcher that re-reads the file with our own reader) · `run/ncbi.py`
  (esearch/esummary/taxonomy efetch, 3 req/s, injectable transport) · `plastid panel refresh |
  propose | pin`. Closed P1-12 as a **dev decision** (local for phases 1–5/8, remote for 6–7) and
  added the "Who decides what" rule to CLAUDE.md so engineering calls stop being queued for the PI.
- **Validation:** `make check` (94 fast + 13 golden). Live NCBI run: 76 RefSeq genomes for
  Iridaceae; a re-pinned NC_050833.1 whose **sha256 matches the vendored golden file byte-for-byte**
  (`a946f23e…`), NC_024936.1 → `fd441727…`. Measured and recorded in PROJECT.md §9.C-27: NCBI's
  `chloroplast[filter]` misses "plastid"-worded records (71 vs 76) — including our own comparison
  genome.
- **Residual risks:** the panel's tier model relies on NCBI lineage strings (genus = last component);
  a record with an unusual lineage would land in the wrong tier — the user confirms the final list,
  and `propose` prints the tier for every candidate. `pin()` records unparseable downloads as
  `parse_error` rather than failing the batch — a caller must check it (the CLI prints UNVERIFIED).
- **Blocked:** nothing. Waiting on people only: real reads + nuclear genome (Phase 6/7), item (c).
- **Next:** Phase 2, starting with M2.1 (gene content & summary statistics).

### 2026-10-02 · Claude · dev WSL

- **Did:** created the session set (CLAUDE.md · ROADMAP.md · TASKS.md · HANDOFF.md); squashed git
  history into one `init` commit and force-pushed (amended to remove Claude attribution — repo rule
  + memory); built **P1-01** scaffold, **P1-02** core model + BioPython GenBank reader + vendored
  golden data (sha256-pinned), **P1-03** IR detection + mandatory gate + coordinate normalization.
  §10.6 decided and recorded (GenBank coords per genome · normalized for cross-genome figures).
- **Validation:** `make check` (40 fast + 7 golden, no skips); the 20-genome sweep reproduces
  PROJECT.md §8.3 — 7/9 vs GenBank, gate catches 2/2 diverged genomes with 0 false alarms; golden
  tier-A numbers asserted in tests. Recorded in PROJECT.md §8.5.
- **Did (cont.):** §10.7 decided (reports resolve table cells; an unresolvable reference raises) →
  PROJECT.md §3/§10. **P1-04 green** — `result/result.py` (`Scope` · `Column` · `Table` · `Result`,
  refs `table.column` / `table[row].column`) and `result/export.py` (csv with provenance comments,
  json with raw values, markdown); exports byte-identical across runs.
- **Residual risks:** de-novo detection is still silently short on 1/9 known cases (NC_063920 →
  22,998 bp vs the true 26,161); the annotation cross-check *records* the 3,163 bp disagreement but
  the run does not fail on it — decide the policy later. The de-novo detector is unvalidated on
  mitochondria. `csv` exports carry provenance as `#` comments — readers must strip them (stdlib
  `csv` cannot).
- **Blocked:** §10.8/§10.9 (P1-10) · D-01. Nothing blocks P1-05.
- **Next:** P1-05 test-first — run layout, manifest, `run_tool()` (the single subprocess door), cache.
