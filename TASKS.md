# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 3 — Annotation  *(PROJECT.md §13.3)*

**Phase check (§13.3):** the run's output re-reads with **our own GenBank reader** (no format
damage) · ORF checks pass · gene recount uses the fixed counting **rule**
(**CDS + tRNA + rRNA, both IR copies** — *132 is the golden genome's value, not a constant*: P2-10
measured **131** for *S. diffusum* and **133** for *C. amboinicus*, both correct).
*De-novo assemblies get the extra "verify vs close relatives / reads across LSC–IR · SSC–IR" rule —
Phase 7 only, skipped here.*

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D3-01 | Which plastome annotator | ✅ **DONE 03/10 (dev)** — **the offline tool is the pipeline; GeSeq is a one-off cross-check whose output gets pinned.** Not a preference: the toolkit must run offline and ship in Docker, and every run must record tool name · version · argv — a web service can record none of that. PGA primary (GPL-3.0, batch, emits GenBank), plann.pl fallback. GeSeq's result is **pinned** via `panel pin`, which turns one manual web step into an offline reproducible input. Verified 03/10: GeSeq has no standalone distribution. PGA/plann **ABSENT** here; Perl present, BLAST+ needed | Dev | `docs/D3-01_annotator.md` + §3 |
| D3-02 | **Semantic gene-name aliases** — `clpP1` · `pafI` · `pafII` · `pbf1` (in *C. amboinicus* `PX923121.1`) | ✅ **DONE 03/10, dev, with evidence**: all four are **the same locus** as a classic gene, proven three ways in our own 9 records — identical syntenic neighbours (`pafI` between `psaA`/`trnS-GGA` where 7 records have `ycf3`), identical length (`pbf1` and `psbN` both exactly 132 bp), and matching `/product`. **Canonical = the team's own dialect** (`pafI`·`pafII`·`pbf1`·`clpP1`): those records are **the team's** (Vu,H.-T. = cô Trang), one is the **manuscript species**, and GeSeq — their tool — emits those names. Normalising to the classic names would make the team's paper print a name that disagrees with the team's own submitted record. Classic names travel as `source_name`, and output prints **both** (`pafI (ycf3)`). `gene_groups.yaml` gains the four names in their functional groups | Dev | PROJECT.md §8.10 + §3 |
| D3-06 | The name **printed in the manuscript** and **submitted to NCBI** | recommended: **print both** — `pafI (ycf3)` — which is what published gene tables already do and what removes the reader's need to know which dialect a record uses. The PI confirms the wording | PI (cô Trang) | PROJECT.md §3 |
| D3-03 | Which translation table for reading frames | the record's own `/transl_table` — **all 1,099 CDS in the dataset declare 11** — falling back to `defaults.yaml` (1) when absent. Tables 1 and 11 share amino-acid assignments and differ only in initiation codons, so P2-03/P2-04 numbers are unaffected | Dev | `data/defaults.yaml` + PROJECT.md §3 |
| D3-04 | `.sqn`: build it, or ship the modern path | emit **`.fsa` + `.tbl`** in pure Python (always); produce `.sqn` only when `table2asn` is present — measured **ABSENT** here. NCBI's web submission takes `.fsa` + `.tbl` directly | Dev | PROJECT.md §3 |

### Pitfalls for this phase *(measured on this dataset — PROJECT.md §8, §9)*

| Trap | Instead |
|---|---|
| "A CDS must start with ATG" | **5 of the golden's 86 do not** — rps19 GTG, rpl2 ACG ×2, ndhD ATC. Whole set: ATG 746 · GTG 11 · ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1. A rule that fails 3.7 % of real plastid CDS is a false alarm, not a check |
| Assuming translation table 1 because D2-03 says so | **every CDS here declares 11**; use the record's, and fail loud on a table whose amino acids we do not know |
| Ignoring `/codon_start` | all 1,099 declare `1`; a `2` shifts the frame — refuse rather than translate a frame the record did not ask for |
| "Fixing" a truncated gene | the golden's truncated ycf1 has `gene` + `misc_feature` and **no CDS** (§9.B-9). `count_genes()` is right; a repair is wrong |
| Renaming a gene so it looks canonical | a rename changes a submitted record. Unmapped names are **listed**, with their `/product`, never guessed |
| Reading the naming difference as "four odd genes" | it is a **whole dialect**: tRNA anticodon hyphen (8 records) · brackets (1) · **bare `trnK`/`trnI`** (2); rRNA `rrn16` (6) vs `rrn16S` (3); protein classic (7+golden) vs functional (2). Normalise the syntax *and* the four semantic names, or the comparison splits one gene into two |
| Reading a tRNA anticodon at a fixed offset | measured: offset 34 matches only **106/293** named tRNAs — the anticodon is not at a constant position here. Identify a bare name by **sequence identity** against the copies that do carry one (`trnK` ≡ `trnK-UUU`, 100.0 %; `trnI` ≡ `trnI-GAU`, 100.0 %) |
| Trusting an annotator's output blindly | read it back with our own reader and require the model to round-trip; a file that does not re-read is not a result |
| Assuming an annotator is a pure function of its input | PGA deletes the contents of its output folder, breaks on filenames with spaces, and mis-annotates a multi-FASTA target. Our fresh-per-run layout, `slugify()` and one-record-per-`--input` rule already make all three impossible — keep them that way |
| Pasting a web tool's answer into the pipeline | pin it (`plastid panel pin`: accession · date · sha256, re-read with our own reader) so the manual step happens once and every later run is offline and byte-identical |
| Writing a TBL from memory | the 5-column NCBI table has exact rules (partial flags `<` `>`, `join(...)`, `codon_start`, `transl_table` per feature) — check against a real submitted record |
| Quoting "132 genes" for a new genome | 132 is *Iris domestica*'s. Write the breakdown: **132 = 86 + 38 + 8** for the golden, whatever the record actually has otherwise |

### P3-01 · ORF checks — start · stop · internal codons · M3.3 · [x] DONE 03/10

- [x] read the record's own `/transl_table` and `/codon_start` into the model (`Gene` gains both);
      refuse a table we cannot translate and a `codon_start` that is not 1
- [x] classify every CDS: `ok` · `alternative_start` · `no_start` · `no_stop` · `internal_stop` ·
      `not_multiple_of_three`; pseudogenes marked, not silently counted as failures
- [x] report genes **annotated without a CDS** as their own class — the truncated-ycf1 shape
- [x] never rewrite a CDS to make it pass
**Done when:** on the golden genome all 86 CDS give 0 internal stops, 0 missing stops, every length
divisible by three, and **exactly 5** classified `alternative_start` (rps19 GTG · rpl2 ACG ×2 ·
ndhD ATC); the truncated ycf1 is reported as *annotated, not coding* rather than as a failure.
**Pitfall:** the 3.7 % alternative-start rate above. **Needs:** D3-03.
**Measured (PROJECT.md §8.11):** the contract above, exactly; plus two rules the tests forced —
**ACG** is accepted as a start but classed `non_canonical` (it is *not* in NCBI table 11's list,
yet the published record uses it), and *annotated-not-coding* matches **name and position
together** (the truncated ycf1 overlaps ndhF's first 21 bp, so either half alone goes silent).

### P3-02 · Gene-name check & normalization · M3.2 · [x] DONE 03/10

- [x] **one** normalizer in `core/gene_names.py` — and nothing else in the codebase may rename a gene.
      It covers the measured dialects: tRNA `trnL(UAA)` → `trnL-UAA`, rRNA `rrn16S` → `rrn16`, case,
      whitespace — **plus a bare tRNA name resolved by sequence identity**, not by guessing
      (`trnK` ≡ `trnK-UUU`, `trnI` ≡ `trnI-GAU`, both 100.0 % on this dataset)
- [x] `data/gene_aliases.yaml` for **semantic** aliases — seeded with `pafI`·`pafII`·`pbf1`·`clpP1`
      per D3-02, each entry carrying the evidence it rests on; it ships **empty of guesses**
- [x] `data/gene_groups.yaml` gains those four names in their functional groups (PSI assembly, PSII
      biogenesis, protease) — today they land in *unclassified*, which is how this was found
- [x] a name with no mapping is **kept** and listed in the Result with its `/product`, never guessed
- [x] a mapped name keeps its **`source_name`**, so a reader can always trace back to the record
**Done when:** the normalizer is idempotent over all **156** distinct names (panel ∪ golden) and
never rewrites one — the golden's **112** names arrive unchanged and resolve to **112** distinct
canonicals, none dropped and none merged — and on the panel it resolves `pafI`·`pafII`·`pbf1`·`clpP1`
**and** `trnK`·`trnI` (`trnK-UUU` · `trnI-GAU`, 100 % sequence identity), each while still reporting
the name it arrived under; `gene_groups` reports **0 unclassified** on all 11 records; any name with
no mapping is reported rather than renamed.
**Needs:** nothing — D3-02 is decided (PROJECT.md §8.10). Only the *printed wording* (D3-06) waits.
**Measured (PROJECT.md §8.12):** the whole Done-when, on the 11-record panel — 155 arriving names →
115 canonical, **0 unresolved**, both bare names resolved exactly as §8.10 measured. Two things this
work changed rather than confirmed. (1) The earlier Done-when sentence — "leaves the golden's 112
unchanged" — was written *before* D3-02 flipped in session 5: with `pafI` as the key, the golden's
four loci **are** renamed, so the check above now states what the decision implies. (2) The CLI
required `len(genomes) == arity`, refusing an 11-record panel for an `arity=2` module; API v0.2 says
"**at least** arity" and the runner already agreed — the CLI did not. Fixed, with a regression test.
**Wiring:** the normalizer feeds `gene_table`, `genome_stats`, `orf` and the new `gene_names` module.
Modules whose gene column is a *per-record key* (introns, diversity, substitutions) still print the
name the record arrived with — Phase 4 keys its comparisons on the canonical name and must wire them.

### P3-03 · Annotation — plastome · M3.1 · [x] DONE 03/10 — **annotator: Chloe**

> **03/10 — PGA đã đo và không đạt** (§8.15): 147 gen (7 reference) / 130 gen (1 reference) so với
> **132 = 86 + 38 + 8**; ba lỗi im lặng đã thành lỗi có tên. Theo nguyên tắc "sản phẩm quyết định,
> không phải bộ công cụ" (CLAUDE.md), công cụ được thay bằng **Chloe** (Ian Small, Julia, CLI, có
> bộ reference chuẩn hoá sẵn).
>
> **Đã xong (§8.17): `132 = 86 + 38 + 8`, 0 CDS hỏng ORF** (80 `ok` + 6 `alternative_start`), đo
> bằng `core/gff.py` trên **file GFF thật** của Chloe, lưu vào `tests/data/` (sha256 ghim) nên
> kiểm lại được **không cần Julia**; chạy end-to-end qua CLI 22 giây, manifest ghi
> `chloe 0.1.18` + `julia version 1.13.0` + argv đầy đủ. Con số **134/88 của §8.16 là do cách
> đếm** — gộp dòng exon theo `parent` ra đúng **86 CDS**, không cần "gộp rps12".
>
> **Phát hiện phải nhớ: Chloe mặc định BIẾN ĐỔI trình tự** (đảo chiều SSC rồi chú giải hệ toạ độ
> đó) — số gen và tên gen y hệt, nhưng **12/86 CDS hỏng ORF** khi đọc với trình tự người dùng đưa
> vào. Adapter luôn truyền `--no-transform` và từ chối mọi lần chạy có viết FASTA biến đổi.
>
> **Còn treo (không chặn):** tên `trnG-GCC` vs `trnG-UCC` — bản ghi vàng dùng một tên cho hai bản
> tRNA khác trình tự; **cô Trang quyết** tên in trong bảng bản thảo (PROJECT.md §10.10). Giấy phép:
> checkout không có LICENSE, bản Zenodo ghi MIT — vào đợt M9.3 (`docs/tools/chloe.md`).
> **Máy này:** Julia ở `~/.local/share/mamba/envs/chloe`, code ở `~/toolbox/src/chloe`; test tier B
> đọc đường dẫn qua `PTS_CHLOE` / `PTS_JULIA`.

- [x] `external/pga.py` through `run_tool()`; path, version and argv land in the run info
- [x] the annotator's GenBank output is read back by **our own reader** and must round-trip
- [x] `modules/annotation.py`: the resulting gene set, its ORF status (P3-01) and its name mapping (P3-02)
      go through the same pipeline as a downloaded record — no second code path
- [x] `core/gff.py`: a GFF3 reader (BioPython has none) — rows grouped by id, parts in transcript
      order, malformed rows refused by line number
- [x] `external/chloe.py`: `--no-transform`, the released version and the reference set both pinned
      and both checked, a transformed run refused by name
- [x] the annotator's own file kept as an artifact and vendored as a fixture — the three numbers
      above are checkable on a machine without Julia
- [x] tier B: the golden counts against the real tool, plus the transform trap measured live
**Done when:** annotating a de-annotated golden FASTA reproduces the golden's gene counts by the
fixed rule (86 CDS + 38 tRNA + 8 rRNA) and the output re-reads clean. **Met 03/10 with Chloe**
(§8.17); the PGA path stays in the module for genomes Chloe does not cover, `PTS_REQUIRE_TOOLS=1`-gated.
**Needs:** D3-01 (revised 03/10), and Chloe (or PGA) on the machine.
**Measured 03/10 — PGA RUNS HERE NOW (PROJECT.md §8.15).** Clone at `~/toolbox/src/PGA`, BLAST+
and MAFFT in `~/toolbox/env/pts/bin`. De-annotated golden genome in, and the answer to the Done-when
above is **no**: with **7 Iris references** PGA produced **147 = 95 CDS + 44 tRNA + 8 rRNA**; with
**1** reference (Iris gatesii) **130 = 87 + 35 + 8**. Neither is 86 + 38 + 8. What the run *does*
show is where the difference comes from and that our pipeline reports it honestly: 4 CDS fail the
ORF checks, **4 loci come back twice under two names** (`clpP`/`clpP1`, `ycf3`/`pafI`,
`ycf4`/`pafII`, `psbN`/`pbf1` — the annotator takes one name from each reference), and 4 genes
match no functional group.
**Three traps measured on the way, each now a named error or a retry instead of a silent wrong:**
(1) a reference containing partial coordinates makes PGA annotate **nothing and still exit 0** —
refused by name now; (2) PGA's temporary file names can contain `$`, which the shell it calls
blastn through expands away — roughly **one run in five** produces nothing; the adapter retries
once and says so; (3) the version probe used to run `perl PGA.pl -e …`, which made PGA run with
its *default* arguments — whose first act is `rm -rf gb` in the working directory. The probe no
longer executes a Perl tool at all (PROJECT.md §8.15, regression test in `tests/unit/test_run.py`).
**Still needs a machine with PGA?** No — it needs *references* someone trusts. The next step is a
curated reference set and, if the team wants PGA in the pipeline, the curation step that removes
its duplicate calls.

### P3-04 · Annotation — mitochondrion · M3.1 · [x] DONE 03/10 — **mito path done · annotator: web service + pin**

- [x] **first, fix M3.1's tool**: MITOS2 is a **metazoan** annotator (reference set `refseq81m`) and
      cannot annotate a plant mitochondrion — its output would be wrong, not empty. **PMGA** is the
      tool that fits (Li et al. 2025, benchmarked at 1 wrong gene vs GeSeq's 18), at the cost of a
      heavy environment. Decide before implementing — PROJECT.md §7.8b
      → **DECIDED 03/10 (dev + user): PMGA.** Recorded in PROJECT.md §3.
- [x] the mito path: no canonical IR, so nothing may call `ctx.require_ir()` on it — measured on both
      stand-in mitogenomes, which raise `IRDetectionError` and are **right** to
      → `modules/mito_genes.py`; a test runs it with `ctx.ir=None` and another shows a gated module
      still stopping on the same context
- [x] ~~`external/pmga.py` through `run_tool()`~~ — **NOT BUILT, and that is the decision (dev, 03/10)**
      — see the note below; the count check it was meant to enable is met by the pinned route instead
**Done when:** on a mitogenome the gene-set table reports plant mito gene names and the run does not
touch the IR gate; on a machine with PMGA the count reproduces the record's own annotation.
**MET (03/10, both halves):** the mito path reports `OX459128.1` as **51 = 31 CDS + 17 tRNA +
3 rRNA** with no IR gate (§8.14), and the count check runs against a **pinned** annotation —
`plastid panel pin --file … --note …` (§8.19), measured on the same record.
**Needs:** PMGA on the machine — **it is now here, downloaded and unpacked (03/10, `docs/tools/pmga.md`).**
**DECISION (dev, 03/10, confirmed by anh Trí) — PMGA is a declared tool, not an adapter, because
this machine is not strong enough to run the annotator.**
Measured: the container is 3.18 GB compressed → **11 GB unpacked, 138,045 files**; running it
needs Apptainer and drives MAKER inside, and the machine's WSL VM is capped at 8 GB while the real
C: drive had only 42 GB free. The container was downloaded, its real command line read out of it,
and then deleted — 14 GB, a third of the free space, for a tool this machine cannot comfortably
run. So the mito annotator follows the plan's own backend rule (PROJECT.md §3, P1-12) down its
**second** route: **run PMGA on its web service once per species, pin the GenBank it returns, and
let `mito_genes` read it offline.** `external/pmga.py` is written when a machine that can run the
container exists — the trigger is that machine, not a date, and the first number to measure then is
peak RAM and wall clock.

The plan's rule it follows — and what each route is for:


PROJECT.md §3 (P1-12, the backend rule): *a heavy step runs where it fits* — locally when the
machine can, on a web service or a bigger server when it cannot. For PMGA that means:

* **option 1 — container, local.** *Not on this machine any more*: it was downloaded and unpacked
  (3.18 GB → **11 GB**, 138,045 files, sha256 `81a23743…`), its command line was read out of it,
  and it was deleted — the real C: drive had **42 GB free**, the WSL virtual disk does not give
  space back, and this machine cannot comfortably run MAKER anyway. Re-download is one `curl`
  (`docs/tools/pmga.md`). Apptainer 1.5.4 is installed, and **the real command line read from the
  container's own POD**: `mgavas -pid <id> -in <fasta> -db <1|2|4>
  [-conf <file>] [-ref <genbank>] [-outdir <dir>]` (enter through `mgavas`/`mgavas.sh`, never
  `mgavas_m` alone — the wrapper sets the environment). It bundles MAKER, vmatch, ARAGORN,
  tRNAscan-SE, MISA, TRF, BLAST, EMBOSS, HMMER, circos and Deepred-mt. **Measured cost on this
  machine: 14 GB of a C: drive with 42 GB free, and a VM that balloons toward its 8 GB cap — it is
  not the route for a 4–8 GB laptop.** `external/pmga.py` is therefore a *later* convenience, not the unblocker.
* **option 2 — the web service, pinned. DONE 03/10** (`plastid panel pin --file … --note …`,
  PROJECT.md §8.19): run `1kmpg.cn/pmga` once per species, save the GenBank, pin it — checksums,
  read-back and provenance recorded like any NCBI download. Every later run is offline and needs
  no container at all.

**Done when** is therefore reachable by option 2 today: pin an annotation, run `mito_genes`, and
the mito path is exercised end to end on real annotation. **The adapter stays open, with a named
trigger:** write it when a machine here can run the container — and the first number to measure
then is peak RAM and wall clock, because that decides whether PMGA can ever be a *requirement*
rather than an option.
**Why nothing was written this session — the honest reason.** PMGA is not on conda and has no git
repository: the paper distributes it as a **singularity container** (Figshare 27201798) and a web
server, and no command-line usage is published anywhere retrievable. An adapter written against an
invented CLI would be a placeholder that looks real — the one thing CLAUDE.md bans — so the adapter
waits for the container. **The exact ask:** pull the Figshare container, record its entry point and
`--help` in `docs/`, and then `external/pmga.py` follows the same shape as `external/pga.py`.
**Measured (PROJECT.md §8.14):** the testable half is done. `mito_genes` reports the stand-in
mitogenome `OX459128.1` as **51 = 31 CDS + 17 tRNA + 3 rRNA**, every CDS named by a standard symbol
(`nad1`·`cox1`·`atp6`·`ccmFn`·`matR` …), 0 names left unmapped, and the run never asks for an IR.
Two things came out of building it: the vocabulary must be a **separate file** from the plastid
aliases (measured collision — "ribosomal protein S3" is a product string in *both* organelles, so a
shared table could rename a plastid gene with a mitochondrial symbol), and the stand-in record's CDS
starts include `CCT` and `CGC`, which no translation table accepts — the ORF check flags them, and
that is a finding about the record rather than a bug to paper over.

### P3-05 · NCBI submission files — FASTA · TBL · SQN · M3.4 · [x] DONE 03/10

- [x] `.fsa` and the 5-column `.tbl`, in pure Python, from the *normalized* gene set
- [x] `.sqn` through `table2asn` when it exists; a named error and exit 3 when it does not
- [x] a first-pass check on the record: no duplicate locus tags, no feature past the sequence end,
      `transl_table` on every CDS
- [x] **the `.sqn` step RUN FOR REAL (03/10, §8.18)** — with a stand-in template vendored in the
      repo, so the step is no longer blocked on the team's `.sbt`
- [x] the record's **source feature** and its `/exception` travel into the table; a CDS whose
      start codon is not an initiator is written with a **partial start** — three fixes that took
      NCBI's validator from **13 errors to 7**
- [x] NCBI's `.val` is **read**, not just kept: its counts and classes are tables in the run
**Done when:** a golden record produces `.fsa` + `.tbl` whose features re-read into the same set as
the source, and the TBL's partial/join syntax matches a real submitted record field for field.
**Needs:** D3-04. **No longer blocked:** `table2asn` is installed and declared (v1.29.324, md5
recorded in `docs/tools/table2asn.md`) and the `.sqn` step runs end-to-end with a **stand-in**
template committed for tests. For a real submission the team still needs its own `.sbt` from
NCBI's web form — it names cô Trang and the institution, and nobody else can generate it. That is
now a *submission* prerequisite, not a code one.
**Measured (PROJECT.md §8.13):** the round trip holds on the golden — **265/265** features
(133 `gene` + 86 CDS + 38 tRNA + 8 rRNA) come back with the same kind, name, intervals (partials
included), product and table, and 0 first-pass problems. The syntax rules were taken from NCBI's
documentation, not memory — `join()` does not exist in this format (intervals **stack**), a
minus-strand interval is **reversed**, and `<` goes in column 1 while `>` goes in column 2
whatever the strand, both as **prefixes** (`>1050`, not `1050>`).
**The "field for field" clause is now closed for everything but joins** — three real tables
NCBI's own pipeline produced (VADR output, vendored under `tests/data/` with provenance and
sha256) round-trip through our writer **byte for byte**, record by record. Reading them changed
the writer twice: interval lines carry **three** fields, not five padded ones (NCBI's own tables
do not pad them), and a table file may hold **several records** — one VADR file carried eight, and
the parser had been merging them into one feature list, silently attaching every record's features
to one sequence. Two gaps remain, both named: **joins** are verified against NCBI's documented
example rather than a real file (no real table with a multi-interval feature was found), and the
`.sqn` half needs the team's `.sbt`.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Which gene name the manuscript and the NCBI submission print (D3-06) | cô Trang | the *wording* of the gene table, not the pipeline |
| A first run of PGA by hand (`~/toolbox/src/PGA`, GPL-3.0) | anh Trí | the *measured* half of P3-03 — **done 03/10**, PGA does not reproduce the golden record |
| Chloe (Julia) installed on the team's machines, or a decision to run the annotation module only where it is available | anh Điềm | the `annotation` module's default annotator; the GFF fixture means everything else is testable without it |
| ~~The mito annotator decision~~ — **decided 03/10: PMGA** (PROJECT.md §3) | — | — |
| ~~The PMGA singularity container~~ — **done 03/10, then deleted**: its entry point and CLI are recorded (`docs/tools/pmga.md`); running it needs a machine stronger than this one | — | closed |
| A submitter template (`.sbt`) from NCBI's form — names people and institutions | cô Trang | **an actual submission**, not the code: the `.sqn` step runs on a stand-in (P3-05, §8.18) |
| Whether `PX923121.1` is the specimen this manuscript reports | cô Trang | quoting *C. amboinicus* numbers as the manuscript's |
| The team's own *O. diffusa* record | anh Điềm / cô Trang | the *O. diffusa* half of P2-10 (today: a stand-in) |
| Item (c) text for the Reference Panel notes | cô Trang | D-03 only (documentation) |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| H-01 | MAFFT 7.526 | **DONE 03/10** | installed from conda-forge, tarball md5 **matches §7.8 exactly** (`cf3b2acf…`) — P2-07/P2-08 real runs are now possible (PROJECT.md §7.8b) |
| H-03 | P2-10's *C. amboinicus* label | **DONE 03/10** | re-run as `own-record-coleus-amboinicus`; both draft labels now state what they are (`standin-` for someone else's record, `own-record-` for the group's) — `docs/P2-10_manuscript_draft.md` |
| H-02 | Three merged Phase-1 branches (`feat/phase-1`, `feat/p1-02-core-model`, `feat/p1-04-result-schema`) are labels on commits already in `main` | **DONE 03/10** | measured first: **0 commits outside `main`** on each; deleted on remote and locally by the user's decision — nothing is lost, every commit is still reachable from `main` |
