# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 2 — Genome reading & statistics  *(PROJECT.md §13.2)*

**Phase check (§13.2):** the golden numbers must come out of the module runs — 153,724 bp ·
**132 = 86 CDS + 38 tRNA + 8 rRNA** · IR 26,214 bp (≥ 99.9 %) · LSC 83,127 / SSC 18,169 ·
37 SSR (mono 31 · di 4 · tri 2). Any mismatch ⇒ stop and fix inside this phase.
Applies to **both plastid and mitochondrion**.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D2-01 | **Multi-genome contract** — π, SNP and later Ka/Ks take more than one genome; `run(genome, params, ctx)` is frozen in API v0.1 | ✅ **DONE 03/10**: contract thứ hai `run(genomes, params, ctx)` chọn bằng `arity`, giữ nguyên v0.1 (không migrate module nào) | Dev | `docs/API_v0.2.md` + PROJECT.md §3 |
| D2-02 | π window size / step | ✅ **DONE 03/10**: 600 bp window · 200 bp step, vào params + report | Dev | `data/defaults.yaml` |
| D2-03 | Translation table for codon/amino-acid work | NCBI table 1 (standard) for plastids; the table used is recorded in run info | Dev | `data/defaults.yaml` |
| D2-04 | Word/Excel writers | ✅ **DONE 03/10**: extra tuỳ chọn `docx = ["python-docx>=1.2"]`, import bên trong hàm ⇒ **bản cài tối thiểu vẫn thuần Python**. `openpyxl` **hoãn** tới khi có bảng Excel thật (Giai đoạn 5) — chưa có gì dùng thì khai báo chỉ là chỗ giữ chỗ | Dev | `pyproject.toml` + PROJECT.md §3 |

### Pitfalls for this phase *(every line has bitten somewhere — PROJECT.md §9, §7.9)*

| Trap | Instead |
|---|---|
| Writing "132 genes" without the breakdown | always **132 = 86 CDS + 38 tRNA + 8 rRNA**; `gene` features = 133 is a diagnostic only |
| "Fixing" the gene count by counting `gene` features | the truncated ycf1 has `gene` + `misc_feature`, no CDS — `count_genes()` is right |
| Assuming one contiguous span per CDS | rps12 is trans-spliced (3 parts, mixed strands) — build from `parts` in biological order |
| Trusting a `repeat_region` note as truth | IR comes from the gate-passing detection; the annotation is only a cross-check |
| Positional work without the gate | call `ctx.require_ir()` first; a failing gate must stop the run, not warn |
| Typing numbers into report templates | templates reference `{{table.column}}` only — the scan test enforces it |
| Skipping a file and carrying on | record the failure in the Result/manifest and name the file; exit non-zero when required |
| Hand-editing an output table | every table is regenerated from the module; if it cannot be, the module is wrong |
| Changing a window/step/threshold "just to see" | params go into the manifest and the report — always |

### P2-01 · Gene content & summary statistics · M2.1 · [x] DONE

- [x] per-kind counts (CDS/tRNA/rRNA), GC, length summary — reuse `count_genes` and `gc_percent`
- [x] run on both golden genomes and on the two diverged-IR genomes (the second must exit 2)
**Done when:** one `plastid run` reproduces 153,724 bp · 132 = 86+38+8 · GC 37.85 % for NC_050833.1
in the report and the CSV, with scope labels.

### P2-02 · Genome structure — boundaries · %GC per region · GC-skew · M2.3 · [x] DONE

- [x] LSC/IRa/SSC/IRb boundaries and lengths (from the IR pair, normalized frame)
- [x] %GC per region; cumulative GC-skew per region, plotted values tabulated
      — `modules/structure.py`, tables `regions` + `gc_skew_windows` (window/step params);
      measured on the golden genome: LSC 35,98 % / IR 43,05 % / SSC 31,45 % GC, skew
      LSC −0,0274 · IRa +0,0346 · SSC −0,0618 · IRb −0,0346 (the two IR copies mirror exactly,
      as reverse complements must); `PZ428230.1` investigated and excluded — PROJECT.md §8.6, §9.A-8
**Done when:** LSC 83,127 · SSC 18,169 · IRa = IRb = 26,214 · regions tile the genome exactly
(sum = 153,724) — the golden run asserts it.
**Pitfall:** GC-skew sign convention ((G−C)/(G+C)) and the cumulative form must be stated in the
column headers; a silent sign flip is this phase's classic quiet error.

### P2-03 · Codon usage — RSCU · ENC · CAI · CBI · Fop · M2.4 · [x] DONE

- [x] CDS extraction: exons joined in biological order, strand-aware, pseudogenes excluded
      — `coding_sequence()` reads each part on its **own** strand in annotation order; measured on
      both golden rps12 copies: only this reading gives `ATG…TAA` (re-sorting → `AAT…TAT`,
      all-forward → `ATA…ATT`). `cds_features()` is the one source of truth for the CDS set
- [x] RSCU per codon + the family sums check (Σ RSCU per amino acid = family size)
      — `core/codon.py` + `modules/codon_usage.py`; golden run: 86 CDS, 26.263 codons,
      ENC 49,69 · CAI 0,5375 · CBI 0,2705 · Fop 0,5118 (reference = the CDS set itself, named in
      the output), 0 internal stops, all 20 families `ok`; **D2-03 decided**: NCBI table 1
**Done when:** the RSCU family sums check passes for every amino acid and the CDS set equals the
86 counted CDS. **Needs:** D2-03 (translation table).

### P2-04 · Amino acid composition (+ GRAVY · aromaticity) · M2.5 · [x] DONE

- [x] translate the same CDS set used by P2-03 (one source of truth for the CDS set)
      — `cds_features()` in `core/model.py` is that one source, used by both modules;
      golden run: 26.263 residues = the P2-03 codon total exactly, 86 terminal stops excluded,
      GRAVY −0,1034 · aromaticity 0,1118 (Kyte–Doolittle · Lobry & Gautier, named in the output)
**Done when:** recomputable from the codon table in P2-03 — same CDS set, same totals.

### P2-05 · Introns — count · length · position · M2.7 · [x] DONE

- [x] per gene: number of introns, lengths, positions, phase; trans-spliced genes handled
      — `modules/introns.py`: golden 17 introns in 12 genes · 12.547 bp · longest 1.055 bp ·
      phase = preceding-exon bases mod 3. **Trans-splicing junctions are separated from introns**
      by the region test (parts in different plastid regions) and cross-checked against the
      record's own `/trans_splicing` qualifier: 2 junctions, both rps12, `agrees`.
      Also hardened: `Table` now rejects a repeated key column (the IR gives several genes a
      second copy — `rpl2.1.1` vs `rpl2.2.1` — and a repeated key resolved silently)
**Done when:** the intron list cross-checks with the GenBank annotation part boundaries.

### P2-06 · SSR by region and motif · M2.8 · [x] DONE

- [x] summary tables by region (LSC/IRa/SSC/IRb) and by motif class, on top of the existing
      per-SSR table from `modules/ssr.py` — `ssr_by_region` is computed **from the per-SSR rows**,
      and the summary carries a `GFF vs statistics` reconciliation column.
      **The reconciliation immediately caught a real silent bug**: MISA writes monomeric SSRs as
      `monomeric_repeat` but every longer unit as `microsatellite` (class only inside the Note),
      and a compound cluster as one `repeat_region` line — the parser read only the first
      vocabulary, so **7 of 37 SSRs were being dropped silently** (30 rows). Fixed; the golden
      run now reports 37 = mono 31 · di 4 · tri 2 and `agrees`; LSC 27 · SSC 10 · IR 0.
      Regression tests in `tests/unit/test_misa.py` use **real MISA line shapes** — the old
      fixture invented a `dinucleotide_repeat` type MISA never writes, which is why it passed
**Done when:** tier-B run with MISA reports 37 SSRs split by region and motif; the counts reconcile
with the per-SSR table. **Needs:** nothing — the tier-B harness already exists.

### P2-07 · SNP / substitutions by position · M2.6 · [~] IMPLEMENTED — real-data run pending

- [x] pairwise substitution table against a chosen reference, positions scope-labelled
      — `modules/substitutions.py` (arity 2+): per-genome and per-region Ts/Tv against row 0,
      gap columns skipped and counted, Ts/Tv blank when undefined. Maths in `core/align.py`,
      unit-tested analytically.
- [ ] **verify on real data — needs MAFFT 7.526, absent from this machine.** `apt` would install
      a different version, and §7.8 pins 7.526, so no number was measured here rather than a
      number measured with the wrong aligner. `PTS_REQUIRE_TOOLS=1` on a machine that has it.
**Done when:** the table states its scope (whole plastome vs CDS-only, IR in/out) and the numbers
reconcile with the alignment used. **Needs:** D2-01.

### P2-08 · Nucleotide diversity π (sliding window) · M2.9 · [~] IMPLEMENTED — real-data run pending

- [x] π per gene / intron / IGS + sliding window with region shading data
      — `modules/diversity.py` (arity 2+): whole-alignment π, sliding windows (600/200, D2-02),
      per normalized region and **per gene built from its exons** (an intron is not its gene).
      Every π carries its usable-column count; gaps dropped, never scored.
- [ ] **verify on real data — needs MAFFT 7.526** (see P2-07). Pending: π per **intron** and per
      **IGS**, which need the inter-feature intervals mapped to columns.
**Done when:** windows are reproducible from the alignments and the window/step are in the report.
**Needs:** D2-01, D2-02.

### P2-09 · Publication gene table (Word) · M2.2 · [x] DONE

- [x] gene table export (Word) generated from the gene table, never hand-edited
      — `modules/gene_table.py` + `data/gene_groups.yaml` + `result/docx_export.py` + `--docx`.
      The table is the **categorised** one a paper prints, not a per-feature dump; the grouping
      is data (gene → functional group), and it carries its own reconciliation against
      `count_genes()` — CDS 86/86 · tRNA 38/38 · rRNA 8/8 · Total 132/132, all `agrees`
      (a first run reported **0** on the Total row; the reconciliation column is what caught it).
- [x] reproducible bytes: the `.docx` is rebuilt with fixed zip timestamps and pruned of the
      template's 164-style catalogue (~800 kB → **13 kB**), so two builds of one run are the
      same file — **verified through the cache too** (same sha256 on a cached re-run)
- [x] rendered and **looked at** (LibreOffice → PDF → PNG): title, provenance grid, table
      captions above, scope line, header row repeating per page, rows unsplit, numbers
      right-aligned, bold totals row, `Page X of Y` footer
**Done when:** regenerating the document from the module produces an identical file. ✅
**Needs:** D2-04 (DONE).

### P2-10 · Manuscript deliverable — Oldenlandia / Coleus · M2.10 · [~] DRAFT PRODUCED — real data pending

- [x] **one command per species → one Word draft** — `-m` is now repeatable, so a run holds
      several modules and writes `report/report.docx` (prose + tables + provenance, joined from
      the modules' own documents so a cache hit cannot produce an empty one);
      recipe: `docs/P2-10_manuscript_draft.md`
- [x] draft produced for both species from the public stand-in set, labelled `standin-…`
      throughout (header, footer, manifest): **131 = 86 + 37 + 8** (*S. diffusum* `PP078674.1`)
      and **133 = 88 + 37 + 8** (*C. amboinicus* `PX923121.1`), IR gate PASS 100.00 % on both,
      gene table reconciles `agrees` on both — PROJECT.md §8.8
- [x] **the stand-in set turns out to contain both manuscript species**: `PP078674.1` is
      *Scleromitrion diffusum* (the current name of *O. diffusa*) and `PX923121.1` is
      *C. amboinicus*. Same species, **not** the team's specimen — the label stays
- [ ] real datasets substituted, then re-run — numbers are never patched in
- [ ] π and Ts/Tv sections (need MAFFT 7.526, absent here) · figures (Phase 5, M5.1)
**Done when:** every number in the draft resolves to a cell in a run's tables (the queryable-number
rule, §13.5). ✅ for the structure — **the real samples are still what makes it the team's draft**.
**Needs:** the real datasets from anh Điềm/cô Trang.

### P2-11 · Sample-group labels on comparison output · M2.11 · [x] DONE

- [x] every comparison table/figure carries a shared group label from `core/samples.py`
      — `plastid run --groups FILE` labels every input; the label travels in params (a module
      never reads a file) and lands in the comparison table's `Group` column. An input with no
      label stops the run (exit 1) rather than producing an unlabelled row.
      End-to-end test with a stubbed aligner: `tests/test_smoke.py`
**Done when:** at least one comparison output uses a group label end-to-end.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Oldenlandia / Coleus GenBank records (real data) | anh Điềm / cô Trang | P2-10 only — **bản công khai thay thế đã có** (D-04, PROJECT.md §8.6) nên P2-07/08/11 dev được ngay |
| Item (c) text for the Reference Panel notes | cô Trang | D-03 only (documentation) |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| D-04 | Public stand-in dataset — Oldenlandia / Coleus (CLAUDE.md §Data) | DONE | 13 records pinned with sha256; measured table recorded in PROJECT.md §8.6 |
| D2-04 | Word writer as an optional extra | DONE | `pip install 'plastid-suite[docx]'`; the base install imports no third-party beyond biopython + PyYAML |
