# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.

## Resume point

**Phase 1 — COMPLETE.** **D-04 (public stand-in dataset) — COMPLETE** (13 records pinned, sha256
stable, PROJECT.md §8.6). **D2-04 (Word writer as an optional extra) — COMPLETE.**

**Phase 2 status:** P2-01 … P2-06 · P2-09 · P2-11 **DONE**. P2-10 **draft produced and labelled**
(PROJECT.md §8.8, `docs/P2-10_manuscript_draft.md`) — the real samples are still what makes it the
team's draft. P2-07 · P2-08 **implemented and unit-verified, real-data run pending MAFFT 7.526**.

**Phase 3 — code COMPLETE, one verification OPEN (4 full + P3-04's count check), branch `feat/phase-3` (stacked on `feat/phase-2`, which is not yet
merged to `main`).** The board is `TASKS.md`; the one task that closed by *decision* rather than by
code is P3-04's annotator (see below). **Next phase: Phase 4 — multi-genome comparison**
(`ROADMAP.md`; M4.1–M4.6, of which M4.5 needs RNA data the team does not have yet).

- **P3-01 ORF checks — DONE.** `core/orf.py` + `modules/orf.py`; `Gene` now carries the record's own
  `/transl_table` and `/codon_start`. Golden: 86 CDS, 81 `ok`, **5 `alternative_start`**
  (rps19 GTG ×2 · rpl2 ACG ×2 · ndhD ATC), 0 failures; the truncated ycf1 is reported as
  *annotated, not coding*. Two rules the tests forced: **ACG is a start but `non_canonical`** (not in
  NCBI table 11 yet used by published records), and *annotated-not-coding* matches **name and
  position together** (the truncated ycf1 overlaps ndhF's first 21 bp). PROJECT.md §8.11.
- **P3-02 gene names — DONE.** One normalizer (`core/gene_names.py`) + `data/gene_aliases.yaml`
  (4 loci, each with its evidence) + `modules/gene_names.py`. On the 11-record panel: 155 arriving
  names → 115 canonical, **0 unresolved**, `trnK` → `trnK-UUU` and `trnI` → `trnI-GAU` by 100 %
  sequence identity, 0 unclassified. Output prints **`canonical (synonym)`** everywhere. Also fixed
  a real CLI bug: `arity ≥ 2` is "**at least**" (API v0.2) but the CLI demanded exactly. §8.12.
- **P3-05 submission files — DONE, and the `.sqn` half ran for real (03/10, late).** `.fsa` + 5-column `.tbl`, round-tripping
  **265/265** golden features; `.sqn` is opt-in and needs both `table2asn` and a `.sbt` template.
  The `.tbl` rules came from NCBI's docs, not memory (`join()` does not exist — intervals stack;
  minus-strand intervals reverse; `<` in column 1, `>` in column 2, both as prefixes). §8.13.
  **Then the `.sqn` step was run against `table2asn` with a stand-in template (no longer blocked on
  the team's `.sbt`): our table had 13 ERROR-level findings and would have been refused at upload**
  — no `source` feature, an RNA-edited start without its `/exception`, and that start without a
  partial marker (plus two source qualifiers the format has no place for). All three fixed; the
  validator now reports 7, every one of them the `GeneXrefStrandProblem` class that a **published,
  accepted** plastome submission also has — proved by running that real table through the same
  validator. NCBI's `.val` is read into the run's tables rather than kept unopened. §8.18.
- **P3-03 plastome annotation — DONE (03/10, late).** Annotator replaced: PGA → **Chloe**.
  Written: `core/gff.py` (a GFF3 reader — BioPython has none), `external/chloe.py`,
  `modules/annotation.py` grown an `annotator` parameter (`chloe` default, `pga` for genomes
  Chloe does not cover), `docs/tools/chloe.md`, a real GFF fixture in `tests/data/`, tier-B tests
  for the real tool. **Measured through our own reader: 132 = 86 + 38 + 8, 0 CDS failing the ORF
  checks** (80 `ok` + 6 `alternative_start`; the golden record itself is 81 + 5) — §8.17. The
  earlier "134 / 88 CDS" was a *counting* artefact: group the exon rows by the `parent` they
  share and the file carries 86 CDS, so no rps12 merge was ever needed.
  **The finding to remember: Chloe transforms the sequence by default.** With the transform on,
  the single-copy region is reverse-complemented and the run returns the *same 132 genes under the
  same names* — while 12 of 86 CDS stop translating and ndhH sits where ndhH is not. The adapter
  always passes `--no-transform` and refuses any run that writes a transformed FASTA. Also
  resolved: `ndhI` (Chloe) = `ndhL` (golden record) — same position, identical protein, so a fifth
  alias entered `gene_aliases.yaml`. Open for the PI: `trnG-GCC` vs `trnG-UCC`, where the golden
  record uses one name for two different tRNA sequences (PROJECT.md §10.10).
  Julia: `~/.local/share/mamba/envs/chloe`; Chloe: `~/toolbox/src/chloe`. Tier-B tests find both
  through `PTS_CHLOE` / `PTS_JULIA`.
- **P3-03 PGA annotation — RUN FOR REAL, and the Done-when does not hold.** The "waiting on anh Trí"
  item was a *permission* boundary, not a capability one: same machine, and once anh Trí granted
  `Bash(perl *)` in `.claude/settings.local.json` (gitignored), PGA ran. Measured: 7 Iris references
  → **147 = 95 CDS + 44 tRNA + 8 rRNA**; 1 reference → **130 = 87 + 35 + 8**; truth **132 = 86 + 38
  + 8**. So PGA does not reproduce the golden counts, and the module says so rather than hiding it.
  Three traps measured and turned into named errors or a retry: a reference with **partial
  coordinates** makes PGA annotate nothing and exit 0; PGA's temp names can contain `$` which the
  shell eats (≈1 run in 5 is empty, retried once); and the version probe used to run `perl PGA.pl
  -e …`, which made PGA run with defaults — whose first act is `rm -rf gb` in the working directory
  (fixed: the probe never executes a Perl tool now). Also measured: 3 of the RefSeq Iris references
  use the **functional dialect** too, which strengthens D3-02. §8.15.
- **P3-04 mito — the mito path is DONE; the annotator's count check is OPEN, and both routes were
  tried and measured (03/10).** Local: the container needs 11 GB + MAKER, this machine cannot.
  Web: `www.1kmpg.cn/pmga/` is a redirect to `47.96.249.172:16084`, and that host **refuses
  connections** on 16084/16081/16086 (`curl -v`, measured 03/10) — the authors' service is down.
  Closing it later is a re-test, not a task: when that port answers, one web run +
  `plastid panel pin --file` finishes it. (I first closed this task on the *decision* alone, before
  running the web route; that was wrong and is corrected in the board.) The mito path works
  and is measured (`mito_genes`: `OX459128.1` → 51 = 31 CDS + 17 tRNA + 3 rRNA, no IR gate). The
  annotator itself is **PMGA on its web service, pinned** — because **this machine is not strong
  enough to run the annotator**: the container is 3.18 GB → 11 GB unpacked and drives MAKER, the
  WSL VM is capped at 8 GB, and the real C: drive had 42 GB free (14 GB of container was a third
  of it, now deleted). `external/pmga.py` is deliberately **not written**; its trigger is the first
  machine that can run the container, and the first number to measure then is peak RAM + wall
  clock. `docs/tools/pmga.md`, PROJECT.md §3.
- **P3-04 mito — the tool's two options are now in place (03/10, late).** The plan's backend rule
  (PROJECT.md §3 P1-12: a heavy step runs where it fits) was half-implemented — `panel pin` took
  only NCBI accessions, so a tool that ran elsewhere could not be pinned. `plastid panel pin
  --file … --note …` closes that (§8.19), and PMGA now has both routes: the container is
  downloaded here (3.18 GB → **11 GB unpacked**, 138,045 files, Apptainer 1.5.4 installed, real
  CLI read from its own POD: `mgavas -pid <id> -in <fasta> -db <1|2|4> [-conf] [-ref] [-outdir]`)
  but needs more RAM than a 4–8 GB laptop should give it; the web service + pin route works today
  and needs no container. `external/pmga.py` is a later convenience, not a blocker. `docs/tools/pmga.md`.
- **P3-04 mito — the testable half is DONE.** `modules/mito_genes.py` describes a plant mitogenome
  **without asking for the IR gate** (a test runs it with `ctx.ir=None`; a companion test shows a
  gated module still stopping on the same context). On the stand-in `OX459128.1`:
  **51 = 31 CDS + 17 tRNA + 3 rRNA**, every CDS named by a standard symbol, 0 unmapped. The
  vocabulary is a **separate file** (`data/mito_gene_names.yaml`) because the same product string
  occurs in both organelles — measured: "ribosomal protein S3" is a plastid product *and* a
  mitochondrial one. **PMGA itself is still not implemented**: it ships as a singularity container
  with no published CLI, and an invented CLI would be a placeholder that looks real.
- **P3-05 submission — `table2asn` is now installed and declared** (`docs/tools/table2asn.md`,
  v1.29.324, md5 recorded) and both of its error paths are tested (no tool → exit 3; no template →
  exit 3, named). What the `.sqn` still needs is the **team's `.sbt`**: a hand-written one is
  rejected (`Error loading template file`, measured twice), so it comes from NCBI's web form.

**Immediate next action:** write `external/chloe.py` + the GFF read-back and finish P3-03 with Chloe
(the measurement is done; the adapter is not). Then: the two remaining asks are a **curated PGA
reference set** (the annotator's
over-calling is a property of its references — with 7 uncurated Iris references it returns 147
genes for a 132-gene genome) and the **team's `.sbt`** (`cô Trang`) for the `.sqn` half; PMGA's
container is now documented (`singularity exec … mgavas_m -pid … --in … --db 1 --outdir …`, from the
Figshare record) so `external/pmga.py` can be written from the publisher's own instructions — it
needs singularity, which is installable here, and the container (3.2 GB). Everything else Phase 3 can do here has run.
The next session should therefore either chase those two, or start Phase 4 by wiring the remaining
name-bearing modules (introns · diversity · substitutions) to the canonical name, since Phase 4
keys its comparisons on it. **Blockers: `docs/PHASE3_BLOCKERS.md`** (the map of who owes what).

Blocker state, short version:

- **Cleared:** MAFFT **7.526** is installed and verified against §7.8's pin (tarball md5
  `cf3b2acf649bc31e576ed43a2fa4aaa3`, byte-identical) — so **P2-07/P2-08 can finally be measured**
  instead of staying "implemented, run pending". BLAST+ 2.17.0 and Perl 5.32.1 present.
- **Two "blockers" were wrong blockers:** MITOS2 **cannot** annotate a plant mitochondrion (it is a
  metazoan tool, reference set `refseq81m`) — PMGA is the tool that fits; and the gene-naming
  question was a dev call, not the PI's.
- **Still open, each with a named owner and a one-line ask:** `docs/PHASE3_BLOCKERS.md`.
- **New invariant in CLAUDE.md:** a tool that exists only on one machine is not part of the toolkit.
  It must be declared in the repo — `environment.yml` (version pinned), `vendor/` + `SOURCE.json`
  (small scripts), or `docs/` — **in the same session as the install**. `environment.yml` was
  repinned because it listed a bare `mafft`, and §9.B-13 had already measured that the aligner
  changes Ts/Tv.

**Next action: P3-01 (ORF checks).** It is pure Python, fully verifiable on this machine, and it is
the phase's core value. The measurement that shapes it is already made: **5 of the golden genome's
86 CDS do not start with ATG** (rps19 GTG · rpl2 ACG ×2 · ndhD ATC), and across the whole dataset the
starts are ATG 746 · GTG 11 · ACG 6 · CTG 5 · ATA 3 · TTG 2 · ATC 1 · ATT 1 — so a "must be ATG"
rule would fail 3.7 % of real plastid CDS. Every CDS declares `/transl_table=11` and
`/codon_start=1`; the reader currently ignores both.

Then P3-02 (gene-name normalization — **decided**, no longer waiting on anyone), then P3-03/04 (PGA,
MITOS2 — neither tool is on this machine), then P3-05 (NCBI submission files).

**Two facts from session 5 that change earlier work — read these before touching data:**

1. **`PX923121.1` and `PX977793.1` are the group's OWN records** (`Vu,H.-T.` = cô Trang, Jan–Feb
   2026), and `PX923121.1` **is the manuscript species** *C. amboinicus*. So P2-10's
   `standin-coleus-amboinicus` label is **wrong** — see H-03. The *O. diffusa* record `PP078674.1`
   is by Shen,C. and **is** a genuine stand-in.
2. **The naming difference is a whole dialect, not four genes** — and it comes from GeSeq, the
   team's own annotator (GetOrganelle → GeSeq). The comparison key is now the **team's** names
   (`pafI`·`pafII`·`pbf1`·`clpP1`), with the classic names kept as `source_name` and both printed.
   Full evidence: PROJECT.md §8.9–8.10.

Open items carried over from Phase 2, in the order they matter:

1. **P2-07/P2-08 real-data run.** Needs MAFFT **7.526** (§7.8 pins it; `apt` installs a different
   version and §9.B-13 records that the aligner changes Ts/Tv). `PTS_REQUIRE_TOOLS=1` on a machine
   that has it. π per **intron/IGS** is still missing from P2-08.
2. **Nothing cross-checks the de-novo IR *length*** on an unannotated record — `PZ428230.1` passed
   the gate on a 15,756 bp repeat while every congener is ~25.7 kb (PROJECT.md §9.A-8).
3. **`--docx` on a machine without the extra** — verified by test (exit 3, named error), never yet
   run on a machine that genuinely lacks python-docx.
4. **Three Phase-1 branches are labels on commits already in `main`** (`feat/phase-1`,
   `feat/p1-02-core-model`, `feat/p1-04-result-schema`) — deletable, nothing is lost; the user has
   not decided.

**Before any P2-07/08 number is quoted:** the alignment must come from MAFFT 7.526, and §6's golden
Ts/Tv (2.686) has the narrower scope *CDS only · 76 shared genes · IR excluded* — the module's
default whole-alignment scope is a different number and must never be quoted as that one.

**Waiting on people:** a nomenclature ruling on four gene names → P3-02 · real Oldenlandia/Coleus
samples → P2-10 (the stand-in set already contains both species, as *public records*, which is not
the same as the team's specimens) · real Illumina reads + nuclear genome → Phase 6/7 · item (c) →
docs.

**Branches (all pushed to `origin`):** `main` = **Phase 1 complete** · `feat/phase-2` (nine commits,
not yet merged) · **`feat/phase-3` = current**, stacked on `feat/phase-2`. The three Phase-1 label
branches (`feat/phase-1`, `feat/p1-02-core-model`, `feat/p1-04-result-schema`) were **deleted
03/10** — measured first: each had 0 commits outside `main`, so nothing was lost.

Every commit on `feat/phase-2` and `feat/phase-3` was `make check`-verified on its own, so
`git bisect` is usable across both branches.

**Data note:** `.panel/` is gitignored; re-fetch the stand-in set with the two `panel pin --groups`
commands in PROJECT.md §8.6.

## Log

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
  D3-01 updated, PROJECT.md §8.17 · §9.B.20-21 · §10.10 · §3.
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
  `environment.yml` listed a **bare `mafft`** while §9.B-13 had already measured that the aligner
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
     first run, which is the whole reason it exists (PROJECT.md §9.B-19).
  2. the OOXML style pruner read styles as **attributes** (`element.get("w:pStyle")`) when they are
     **child elements** (`<w:pStyle w:val="…"/>`), so it removed nearly every style the document
     used — including `Title` and `Heading 1`. **LibreOffice still rendered something plausible**;
     only a structural test ("every style the document references must still be defined") caught it
     (PROJECT.md §9.C-20).
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
  (`a946f23e…`), NC_024936.1 → `fd441727…`. Measured and recorded in PROJECT.md §9.C-20: NCBI's
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
