# TASKS.md — active tasks

> **Current phase only.** Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 3 — closed 03/10 (record only; details in ROADMAP + PROJECT.md)

| Task | Outcome | Record |
|---|---|---|
| P3-01 ORF checks | golden: 0 failures, 5 alternative starts; the truncated ycf1 reported, not repaired | §8.11 |
| P3-02 Gene-name check | 155 arriving names → 115 canonical, 0 unresolved | §8.12 |
| P3-03 Annotation, plastome | **132 = 86 + 38 + 8** with **Chloe**; PGA measured and replaced | §8.15, §8.17 |
| P3-04 Annotation, mito | mito path done and measured; **the count check is OPEN — Phase 3A is how it closes** | §8.14 |
| P3-05 Submission files | 265/265 features round-trip; `.sqn` ran for real; NCBI's validator 13 → 7 errors | §8.13, §8.18 |

**P3-04's open half, kept here verbatim so it is not lost.** Both routes were tried; both fail, and
only one of them is ours. **Local:** this machine is not strong enough (11 GB unpacked + MAKER vs an
8 GB WSL cap). **Web:** PMGA's application server **refuses connections** — `www.1kmpg.cn/pmga/`
answers 200 but only redirects to `47.96.249.172:16084`, which answers `Connection refused` on
16084/16081/16086. It stays a **re-test**, not a task:
`curl -sS -o /dev/null -w '%{http_code}\n' http://47.96.249.172:16084/home`
— the moment that answers, one web run + `panel pin --file` closes it.

## Phase 3A — Mito annotator, self-built prototype  *(ROADMAP.md → Phase 3A)*

**Why this phase exists.** P3-04 chose PMGA and then could not run it. CLAUDE.md's rule for exactly
that case is *replace it — or write that part ourselves*, and the parts are already here: BLAST+
2.17.0, the 31-locus vocabulary (`data/mito_gene_names.yaml`, measured **1:1** with the answer
key's 31 CDS), the ORF checks, the reader. What was missing is a **public reference set** and a
**tRNA caller** — both located, both small (below). The prototype must **earn** its place: it is
measured against a record whose answer is published, not declared correct.

**Phase check:** on `OX459128.1` the annotator's own output re-reads with our reader, runs through
`orf` + `mito_genes` unchanged, and recounts to **31 CDS · 17 tRNA · 3 rRNA**, every difference from
the published record listed.

**Measured cost, so nobody has to guess:** ~3–5 MB of reference GenBank (not committed) → a
~150 KB committed protein FASTA + a 147 KB ARAGORN package + a <1 MB BLAST DB. **Under 10 MB on
disk**, seconds of runtime, no MAKER. PMGA was 14 GB.

### Decisions to make before the tasks that need them

| ID | Decision | Recommended default (evidence) | Who | Recorded in |
|---|---|---|---|---|
| D3A-01 | The tRNA caller | **ARAGORN 1.2.41** — 147 KB on bioconda, and it is the tool PMGA itself runs inside. Upstream `ansikte.se` answers **HTTP 455** (blocked), so bioconda is the route. A hand-written structural caller is **rejected**: no independent validation, and a caller tuned until it hits 17 is overfitting | Dev | `docs/tools/aragorn.md` + `environment.yml` |
| D3A-02 | How a CDS's boundaries are decided | **From the alignment to reference proteins — never from ORF-walking.** Evidence: the answer-key record's own 31 CDS score **21 ok · 4 pseudogene · 3 no_stop · 3 no_start**. A "find ATG…stop" caller cannot reproduce it, and one that "fixes" it would be inventing a record NCBI did not publish | Dev | PROJECT.md §3 |
| D3A-03 | The reference set | RefSeq mitogenomes by the **same taxonomy rule as the plastid panel** — same family first. Measured available: Lamiaceae `NC_023209.1` *Salvia* · `NC_023103.1` *Ajuga* · `NC_049064.1` *Rotheca*; Rubiaceae `NC_066984.1` *Psychotria* · `NC_057654.1` *Scyphiphora*; models `NC_037304` · `NC_006581` · `NC_011033` | Dev | `data/mito_reference.md` |
| D3A-04 | Trans-spliced genes (`nad1`, `nad5`) | **Exons called, joins refused.** Measured spans **28,320–240,089** and **39,406–104,374** — the prototype reports the exons it can align and never invents a join it cannot show | Dev | PROJECT.md §9 |
| D3A-05 | What this prototype is **not** | Annotation + submission only: **no** de-novo gene discovery, **no** intron model, **no** RNA-editing prediction, **no** repeats. A gene absent from the reference set is reported **absent**, never guessed — those stay PMGA's job if a machine ever runs it | Dev | PROJECT.md §3 |

### Pitfalls for this phase *(measured on this project's own records)*

| Trap | Instead |
|---|---|
| "A CDS must start ATG and end in a stop" | measured: **10 of the answer key's 31 do not** — 3 `no_stop`, 3 `no_start` (ACG · CGC · ATA), 4 pseudogenes. These are RNA-editing-dependent ends, not errors. Report the ORF profile; a "repair" is a fabricated record |
| Treating a BLAST hit as a gene | a plant mitogenome is stuffed with plastid-derived (MTPT) sequence, so a plastid copy can outscore the mitochondrial one. Require the hit to be in the mito vocabulary, keep a stated coverage floor, and **list what was dropped** |
| Taking a reference protein from a `/pseudo` or partial CDS | it poisons the set silently — this is §8.15's measured PGA trap. Exclude pseudogenes and partial coordinates when building the reference |
| Counting tRNAs by anticodon | measured: **17 = 14 unique + `trnH-GUG` · `trnK-UUU` · `trnM-CAU` twice**. Deduplicate and the count is wrong |
| Hand-writing GenBank | the writer must round-trip through **our own reader** before anything is pinned; a file that does not re-read is not a result |
| Quoting a count as "the answer" | **31 · 17 · 3 is *this record's*** — another mitogenome gets its own numbers, and the module already says so |
| A tool that exists only on this machine | ARAGORN ships in the **same session** as the `environment.yml` line and `docs/tools/aragorn.md` — otherwise no number here is reproducible by anh Điềm |

### P3A-01 · GenBank writer · M3A.1 · [x] DONE 03/10

- [x] `to_genbank(genome) -> str` beside the reader in `core/genbank.py` — the source feature, every
      gene kind, qualifiers, partial flags in the location, wrapped sequence
- [x] partials (`<` / `>`) and `transl_table` written the way NCBI writes them — checked against the
      golden record's own text, not memory

**Done when:** `NC_050833.1` read → written → read back gives the same `Genome`: **265/265**
features with the same kind, name, span, product, translation table and partial flags, and the
sequence byte-identical. **Met 03/10** (§8.20), 10 tests in `tests/unit/test_genbank.py`, 324 unit
tests green — the `.tbl` path is untouched. **Needs:** nothing.
**Measured, and it changed the design twice:** `/pseudo` and `/trans_splicing` are **bare**
qualifiers in every real record, and BioPython writes `=""` for an empty string but a bare
qualifier for `None` (§9.B-25); and the qualifier builder became **one shared function**
(`model.gene_qualifiers`) because the table writer already had it — two writers, one decision about
what a feature says. Also decided there: **`/gene` is written only when it says something
`/product` does not**, because `OX459128.1` names its features by a product *sentence*.

### P3A-02 · Reference set · M3A.2 · [x] DONE 03/10

- [x] pin the RefSeq mitogenomes of D3A-03 with `plastid panel pin` — 8 records, **5.2 MB**,
      sha256 in `PINNED.json`
- [x] `core/mito_reference.py` extracts each record's **own `/translation`** and its rRNA into
      `data/mito_reference_proteins.faa` (**287 entries, 100 KB**) and `data/mito_reference_rrna.fa`
      (**25 entries, 49 KB**), headers `accession|gene|vocabulary|organism`
- [x] every exclusion is **returned and printed**, never dropped in silence
- [x] `data/mito_reference.md` records the sources, the checksums and the rebuild line

**Done when:** every one of the 31 vocabulary loci has at least one reference protein behind it.
**Met 03/10** — coverage **31/31**, and the committed FASTA is held against the pinned panel by
`test_the_committed_set_is_what_the_pinned_panel_builds`. **Needs:** D3A-03.
**Measured, and two of them changed the design:** (1) **self-translating is impossible, not just
worse** — all 530 reference CDS publish a `/translation` and none declares a `/transl_table`, and
plant mito CDS begin with `ACG`, which no table accepts as an initiator, so a CDS-mode translation
raises; the reference is therefore the curator's own protein. (2) The records are **polluted** —
236 anonymous ORFs (*Nicotiana* alone annotates 121), plus plastid-derived `rbcL`/`petG`/`atpB` and
horizontally transferred `citB`/`tuf`, plus a *Psychotria* "5S rRNA" that is 722 bp. A naive build
would have swallowed all of it.
**Left open on purpose (H-06):** 23 reference genes have **no vocabulary entry** — ten of them
standard (`rpl2` · `rpl16` · `rps7` · `rps12` · `rps14` · `rps19` · `rpl10` · `rpl23` · `rps2` ·
`sdh3`). They stay in the set marked `-` so a run can *say* it found one.

### P3A-03 · Gene caller — CDS and rRNA · M3A.3 · [x] DONE 03/10

- [x] `blastx` (target genome × reference proteins) through `run_tool()`; best hit per gene above a
      **stated** coverage/identity floor — floors 50 % identity / 50 % coverage, in `defaults.yaml`,
      printed in every run's summary
- [x] refine each hit to codon boundaries on the nucleotide sequence; a first codon that is not an
      initiator becomes a **partial start**, never a shifted one (the submission module already
      writes those). The stop codon is taken only where the sequence confirms it, one codon past
      the alignment — never searched for (D3A-02)
- [x] `blastn` against the reference rRNA for `rrn18` · `rrn26` · `rrn5` — strand from the subject's
      direction, which is the inverse of the blastx convention
- [x] hits below the floor, vocabulary genes with **no** hit, and genes whose symbol the vocabulary
      lacks (H-06) are each reported in their own table

**Done when:** on `OX459128.1` the caller returns the **same 31 symbols** the record carries, each
with coordinates and strand, plus the explicit list of anything not called.
**Met 03/10** — **31 CDS + 3 rRNA called, `not_called` empty**; `tests/golden/test_mito_annotate.py`
pairs every call with the published CDS its **exons** overlap (a bijection over all 31 — spans
cannot be used, a trans-spliced gene's span covers 100 kb of its neighbours) and prints the
differences: **18/31 spans identical, 13 differ by ≤ 75 bp, strand 31/31**, ORF profile
**21 ok · 3 no_start · 7 no_stop**. 47 unit tests + 4 tier-B. **Needs:** D3A-02, D3A-04.
**Measured, and five of them changed the design** (§8.22): BLAST truncates `sseqid` at the first
whitespace; the committed reference FASTA holds 8 duplicate headers, two of them **different
sequences under one symbol** (Arabidopsis `atp6` 385/349 aa) ⇒ an HSP belongs to the smallest
reference length containing it; **`-seg no`** (the mask halves atp9's alignment, 137 → 67 bits);
**no `-parse_seqids`** (it fails outright on our reference set and upper-cases the ids);
**ranking candidates by bitscore builds chimeras** — on `atp6` it pulled the call 381 bp upstream
onto a 75.6 %-identity alignment, so the winner is the highest-**identity** candidate that passes
the floors (total coordinate deviation 310 bp vs 655 bp). Also fixed in `run_tool`: BLAST+ accepts
only `-version`, so the tool this module depends on was being recorded as "unknown".
**Left open on purpose:** the two 19- and 7-codon exons of the published nad1/nad5 are not
recoverable from any alignment, and the caller reports 4 exons rather than inventing a fifth.

### P3A-04 · tRNA caller · M3A.4 · [x] DONE 03/10

- [x] `external/aragorn.py` through `run_tool()`; its tabular output parsed; names through the
      existing normalizer (`core/gene_names.py`), never a private rename
- [x] `environment.yml` gains `aragorn=1.2.41`; `docs/tools/aragorn.md` carries source, licence,
      install line and why it is not vendored

**Done when:** on `OX459128.1` the count is reported against the published **17**, with every
difference — extra or missing — listed by name, position and anticodon. A count of 16 or 18 is
**reported as measured**, not adjusted. **Needs:** D3A-01.
**Met 03/10** — ARAGORN 1.2.41 calls **17**, the record carries **17**, and they are **not the
same 17**: `trnA-UGC 13244-13348` is in the record and not called (105 bp, longer than a tRNA);
`trnL-AAG 36355-36427` is called and not in the record; the tRNA-His at 11741-11819 is
**`trnH-AUG` to ARAGORN and `trnH-GUG` to the record** — an anticodon disagreement, printed rather
than resolved. 9 spans differ by 1–6 bp, the largest being the one tRNA that **crosses the
origin** (`[258274,73]`, written as the two parts a `join()` needs). 12 unit tests + the golden
check. **Needs:** D3A-01.
**Measured, and two of them changed the code:** ARAGORN's `-V` **is not a version flag** — it
exits 0 and prints "No sequence file specified", which `run_tool` recorded as the version of the
tool calling tRNA genes; the probe now requires a digit in the answer, and the adapter verifies the
release from `-h`'s banner itself. And a name is not an identity: the record's 17 tRNAs carry only
**14 distinct names**, so a comparison that deduplicates by name reports 14.

### P3A-05 · CLI + pipeline wiring · M3A.5 · [ ]

- [ ] `plastid mito annotate --input <fasta> --reference <dir> --out <dir>` — FASTA in, GenBank out,
      plus the standard run report (tool · version · argv · time)
- [ ] the output goes through `orf` + `mito_genes` **unchanged** — no second code path, per code seam 3

**Done when:** the de-annotated `OX459128.1` sequence in, a GenBank out that re-reads with our
reader and drives both module reports end to end through the normal runner.

### P3A-06 · Golden check · M3A.6 · [ ]

- [ ] `tests/golden/test_mito_annotate.py`: the 31 CDS symbols · tRNA against 17 · 3 rRNA, against
      the published record, with the differences itemized rather than summarized

**Done when:** the check runs and prints its counts **next to** the published ones. A count that
differs is recorded in PROJECT.md §9 as a finding and the task says so — it is **never** tuned away,
and a red check that reports honestly is a finished task, not a failed one. **Needs:** P3A-03…05.

### P3A-07 · The sample nobody has annotated · M3A.7 · [ ]

- [ ] run on `AP045469.1` (*Coleus hadiensis*, 313,585 bp, **0 features** — measured)
- [ ] report the gene set found **and** the vocabulary genes not found; ORF profile through the
      existing module; scope label stating this is a prototype run, not a submission
- [ ] `plastid panel pin --file <out.gb> --note "prototype mito annotate, <date>"`

**Done when:** `PINNED.json` holds the file with sha256 and a **clean read-back**, and the gene-set
table runs offline from the pin — which is the same shape every later mito run will take.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| Which gene name the manuscript and the NCBI submission print (D3-06) | cô Trang | the *wording* of the gene table, not the pipeline |
| Chloe (Julia) on the team's machines, or a decision to run the annotation module only where it exists | anh Điềm | the `annotation` module's default annotator; the vendored GFF keeps everything else testable |
| A submitter template (`.sbt`) from NCBI's form — it names people and institutions | cô Trang | **an actual submission**, not the code: the `.sqn` step runs on a stand-in (§8.18) |
| Whether `PX923121.1` is the specimen this manuscript reports | cô Trang | quoting *C. amboinicus* numbers as the manuscript's |
| The team's own *O. diffusa* record | anh Điềm / cô Trang | the *O. diffusa* half of P2-10 (today: a stand-in) |
| Whether a third party may receive the team's **unpublished** mitogenomes (PMGA's web service, if it ever answers) | cô Trang | nothing in Phase 3A — the prototype runs offline, which is why it is the route |

## Housekeeping

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED | cô Trang provides the text |
| H-04 | `trnG-GCC` vs `trnG-UCC` — the golden uses one name for two different tRNA sequences | OPEN | cô Trang decides the name printed in the manuscript table (PROJECT.md §10.10) |
| H-05 | PMGA web service re-test | OPEN | `curl` above answers; then one web run + `panel pin --file` closes P3-04 |
| H-06 | `data/mito_gene_names.yaml` holds **31** loci; the references carry **10 standard genes more** (`rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` · `rps14` · `rps19` · `sdh3`) | OPEN, measured 03/10 | each new locus carries the reference `/product` string that justifies it; a test proves the vocabulary covers the set — it changes naming for **every** mito record, so it is its own task, not a side effect |
