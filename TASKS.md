# TASKS.md — active tasks

> **Phase 3B is OPEN** on `feat/phase-3b` — the self-built plastid annotator, with the naming step
> (**D5-09**) as **M3B.0** in front of it. Design: `docs/D5-08_plastid_annotator.md` ·
> `docs/D5-09_reference_naming.md` · milestones: `ROADMAP.md` → Phase 3B.
> Progress: `ROADMAP.md` · session log: `HANDOFF.md` · rules: `CLAUDE.md`.
> Markers: `[ ]` todo · `[~]` doing · `[x]` done · **BLOCKED** = a named decision or dataset missing.
> **Every task ends with a "Done when" — a runnable, measurable check. No check, no done.**
> How this file is produced: the `make-tasks` skill (`.claude/skills/make-tasks/SKILL.md`).

## Phase 3B — the plastid annotator (M3B.0–M3B.5)

**Goal:** the annotate step as this project's own tool — CGAS's flow with our search-and-decide core:
a pinned reference set, per-sample annotation, an audit gate in front of every run, and the naming
step in front of the annotator. **Chloe is dropped** (D5-08). **In-phase checks:** the golden
*sequence* annotates to exactly **132 = 86 + 38 + 8** with every difference from the published
record listed · the 34-record corpus runs leave-one-out, every difference classified · a failed
audit stops the run with `file — reason`, nothing annotated · no silent failures.
*Order: the three naming numbers are measured **before** any `gene_normalize` code (D5-09 §9.3), and
the ★ pick needs the audit's report — so tasks run in the sequence below, not by milestone number.*

- [x] **P3B-01 · The reference set, built and pinned** *(M3B.1 · D-b)* — DONE 08/10 (§8.40)
  - Choose by the taxonomy rule (same family first — the `mito-reference` pattern); the pinned
    Iridaceae corpus (`.panel/iridaceae-34/`, 0 gate failures) is the worked example.
  - Pin accession · sha256 · date in a `PINNED.json`; write the **leave-one-out** design down.
  **Done when:** the pin exists, a test re-verifies every sha256, and a LOO manifest shows the
  target is excluded for all 34 corpus records.

- [x] **P3B-02 · The three numbers, before coding** *(M3B.0 · D5-09 §9.3)* — DONE 08/10 (§8.41);
  margin **≥ 10 bits** · floor **50 % identity + 0.7 coverage** (not the mito path's 0.5) ·
  test 4's measured cases are the **serine tRNAs**, **two genes under one name**, and names the
  vocabulary cannot resolve — *not* `rpl2`/`rpl23` or the two `trnL`
  - The **margin** of test 3, including the genus-distance case (Appendix A: one base in 72).
  - Whether **test 4 (neighbours)** separates `rpl2`/`rpl23` and the two `trnL`; the **`identity_floor`**.
  **Done when:** all three are measured on data already in the repo + the corpus and recorded as a
  new §8 entry — **before** any `gene_normalize` code lands.

- [ ] **P3B-03 · The reference-set audit** *(M3B.1)*
  - The BLOCK / REPORT classes exactly as D5-08's table; `--audit-only`; every message
    `file — reason`, nothing else; confirmed exceptions pinned by sha256 (change the file, void).
  - MAFFT `--adjustdirection` per region + a dotplot for the user's confirmation; name-set report (N ≥ 3).
  **Done when:** negative tests — a CDS missing `/translation` blocks (exit 2, nothing annotated,
  the reason in the manifest); on the corpus the audit reproduces §8.34 — 32/34 one SSC orientation,
  **2/34 inverted and REPORTED, not blocked**.

- [ ] **P3B-04 · `gene_normalize` — the module** *(M3B.0)*
  - New module, arity 2, `requires_ir=False`, modes `normalize` / `verify`; renaming logic stays in
    `core/gene_names.py`; `gene_names` untouched.
  - The four tests in order + the `rule` column; the **status** column (`measured` · `not measured` ·
    `outside vocabulary`); B appears only in the `pairing` table, in its own column.
  - Conflict tiers + the chosen exit code; `--set on_conflict=continue` writes the `conflicts` table;
    the normalized `.gb` is a `Result.artifacts` entry; verify mode reports and writes nothing.
  - Flip `gene_aliases.yaml`'s direction (canonical = the naming reference's spelling) with its guard test.
  **Done when:** over the corpus pair (one sequence, two annotations) the normalize run writes the
  `.gb` in A plus a pairing table; verify on an unchanged file reports **zero** changes; a planted
  conflict stops with the chosen code, naming coordinates; each of the four tests has a negative case.
  **Needs:** the ★ pick (Waiting table) for the real-DB run — a recorded override is allowed meanwhile.

- [x] **P3B-05 · The cache stops silently changing the answer** — DONE 08/10 (§8.42):
  H-15 · H-08 · H-14 fixed, plus §7.2's promise (a file-valued parameter is content-addressed) *(H-15 + H-08 · H-14 · D5-09 §7.2, §9.9)*
  - The cache key carries what changes the output: the **annotation** (names, not only the sequence),
    the input **order**, the **module source**, and the output options (`--docx`).
  - The ★ record's sha256 reaches the key — §7.2's promise, impossible today.
  **Done when:** the three regression tests (the corpus pair with inputs swapped · a changed module
  source · a second run asking `--docx`) fail on today's code and pass after.

- [ ] **P3B-06 · Ts/Tv re-measured under the new direction** *(M3B.0)*
  **Done when:** `substitutions` re-run over the CDS scope (**CDS only · 79 shared genes · MAFFT ·
  IR excluded**) reproduces **2.6917 = 323 Ts / 120 Tv** — or the difference is explained and §6
  updated with fresh measured evidence (it is a prediction until then); the `kaks` ·
  `substitutions` · `popgen` goldens still pass.

- [ ] **P3B-07 · The calling engine** *(M3B.2 · D-e)*
  - Extract the calling core from `core/mito_caller.py` (DRY's second copy) — the mito tests are the guard.
  - `blastx` against the reference's own published `/translation` (**never self-translate**), `blastn`
    for rRNA, IR self-blast behind the structural gate.
  - The six measured rules of §8.22 as code: `-seg no` · never stitch two sequences' HSPs · boundaries
    from the alignment · the minimum-length subject on key collisions · strand from the subject for
    rRNA · identity + coverage floors (**measured on the answer key**).
  **Done when:** the mito tests pass unchanged against the extracted core; each of the six rules has a
  test naming it; a below-floor hit is reported with its reason.

- [ ] **P3B-08 · Arbitration & export** *(M3B.3)*
  - Structure first, score only **within one program**; **overlap kept and reported**; one call per
    locus with the losing candidates listed; ties → one call when identical, both with a note otherwise.
  - GenBank that re-reads, `/organelle` on the source, a reason per rejected hit.
  **Done when:** the golden output carries the four measured overlapping pairs (psbC×psbD 53 bp ·
  ndhC×ndhK 10 · atpB×atpE 4 · rps12×rps12 114); two references disagreeing at one locus give one
  call + a listed loser; the artifact re-reads and carries `/organelle`.

- [ ] **P3B-09 · CLI & batch** *(M3B.4 · H-09)*
  - 1..N samples × 1..N reference records; one artifact per sample, one report per run; a failed
    sample is marked **without killing the run** (the exit code reflects it).
  - The organelle declared at input; a mito input refused **by name**; docx report (the D2-04 pattern).
  - The engine registered as the annotate step, `annotation`'s chloe · pga set retired with D3-01's
    update line.
  **Done when:** a batch with one broken sample completes, names it, exits non-zero; a mito FASTA is
  refused by name; the run directory holds one artifact per sample.

- [ ] **P3B-10 · Golden + corpus validation** *(M3B.5 — the phase's Done-when)*
  **Done when:** the golden sequence annotates to exactly **132 = 86 + 38 + 8** with every difference
  from the published record listed (kind + span, nothing summarized away); all 34 records ran
  leave-one-out, each result diffed against the published record **and** against PGA's reading, every
  difference classified (vocabulary vs conclusion) with **fewer conclusion-level differences than
  PGA's**; a search that finds nothing is a failed run.

## Decisions to make

The recommendation is the default unless the dev says otherwise (CLAUDE.md → who decides).

| # | Question | Recommended default, with its evidence | Recorded in |
|---|---|---|---|
| D-a | tRNA/rRNA calling | **ARAGORN de novo for tRNA** (declared, `aragorn=1.2.41`, matches the mito path, no reference needed) + `blastn` for rRNA (P3A-03's pair) | PROJECT.md §3 · docs/D5-08 |
| D-b | Which records form the reference set | the **taxonomy rule** (same family first) over pinned records; the corpus is validation, never production | the pin · §8 |
| D-e | Shared calling core | **extract from `core/mito_caller.py`** — second real copy arriving (DRY), its tests the guard | code · PROJECT.md §3 |
| — | `gene_normalize`'s organelle group | **plastid** for now — the ★ is one plastid record, and normalizing a mito record against it would print plastid spellings for mito genes (H-12's lesson) | §13 group table |
| — | Conflict-stop exit code | **4** — the ladder is 1 error · 2 usage/IR gate (and the audit) · 3 missing tool | `cli.py` · API doc |
| — | Parameter name | **`naming_reference`** — `reference` already carries two meanings (§7.2) | `defaults.yaml` · API doc |
| — | Manifest schema | **bump `MANIFEST_SCHEMA_VERSION`** — a new per-module reference entry | `run/manifest.py` |
| — | Replacing Chloe | a **new module mirroring `mito_annotate`**; chloe · pga retires when P3B-10 passes; D3-01 gains its update line | PROJECT.md §3 |

## Pitfalls — the project's own measured traps

- **PGA counted `gene` features and reported 147 genes where the golden has 132, with three silent
  failure modes (§8.15)** → the convention is **132 = 86 + 38 + 8**; never print an annotation rate
  that congratulates itself.
- **Never self-translate a CDS** (§8.21) → always `blastx` against the record's own `/translation`.
- **Blast's measured traps** (§8.22) → `-seg no` (masking cost `atp9` half its bits) · never
  `-parse_seqids` (it errors on our reference set) · boundaries from the alignment, never ORF-walking.
- **The golden carries four overlapping pairs in the contract** → a rule that drops the weaker
  overlap deletes real genes; keep both, report.
- **2/34 corpus records carry an inverted SSC with clean ORFs (§8.34)** → report, never block.
- **34/34 corpus records spell at least one locus differently (§8.31)** → A and B never mix in one
  table; `outside vocabulary` is a column, not a colour.
- **The cache is blind to the annotation and the input order (H-15)** → a `cached` run can be the
  wrong answer; P3B-05 lands before any normalize run is trusted.
- **The mito artifact writes no `/organelle` (H-11)** → the plastid annotator's `source` carries it,
  and P3B-08 checks the artifact, not the intention.

## Waiting on data / people

| What | Who | Blocks |
|---|---|---|
| The **★** pick — one record, the project's naming reference | **the end user**, from the audit report (D5-09 §7.1); the PI signs off on printed names | the real-DB normalize run (P3B-04) and every annotation run after it — a recorded override is allowed until then |
| PI sign-off on the names a manuscript prints (D5-09 §9.4) | **cô Trang** | only the printed-name change; the naming step itself does not wait |
| Plant mito data — a real mitogenome dataset | anh Điềm / cô Trang | **H-06 · H-10 · H-12** (on hold by D5-07) and Phase 6's mito content (M6.2 · M6.3) |
| One RNA-seq sample (any tissue) | team | M4.5 RNA editing |
| The team's own Oldenlandia/Coleus samples | anh Điềm / cô Trang | the manuscript tables being *the team's* — the stand-in tables exist today |
| Confirmation of who `Vu,H.-T.` is, and whether `PX923121.1` is the *C. amboinicus* manuscript's specimen | cô Trang | the *C. amboinicus* numbers being quoted as the manuscript's |
| A nuclear genome · a haplotype table · a phylogeny | cô Trang · **us**, once a figure needs them | the NUPT/MTPT (M6.1) · haplotype-network · annotated-tree figures |
| Real Illumina reads (1 sample) · ONT/HiFi reads (1 sample) | anh Điềm / cô Trang | M7.1 · M7.5 · M7.6 |
| The team's own `.sbt` submitter template (NCBI's form) | cô Trang | only the `.sqn` half of M3.4; `.fsa` + `.tbl` submits without it |
| A ruling on what a generated sentence may claim (D5-05) | **cô Trang** | the wording of the `stats` interpretation sentences only |

## Housekeeping

*(H-15 · H-08 · H-14 moved into P3B-05. H-16 is the repeats milestone M2.12, on its own branch when picked up.)*

| ID | Task | Status | Done when |
|---|---|---|---|
| D-03 | Fill Reference Panel item (c) in PROJECT.md §13.1 | BLOCKED on cô Trang | the text exists |
| H-04 | `trnG-GCC` vs `trnG-UCC` (§10.10) | OPEN | cô Trang decides the name printed in the manuscript table |
| H-06 | `mito_gene_names.yaml` holds **31** loci; the references carry **10 standard genes more** | ON HOLD — D5-07; measured 03/10 | each new locus carries the reference `/product` string that justifies it, and a test proves coverage; reopen with plant mito data |
| H-07 | π per **intron/IGS** — the last piece of M2.9 | OPEN | intron spans derived in the shared layer first (DRY's second copy), then `diversity`'s region machinery; cross-checked against `introns`' count (17 on the golden) |
| H-09 | No module declares **which organelle it serves**; no run checks its **inputs** are single-organelle. **Partly in P3B-09** | OPEN, decided 05/10 (D5-06) | the registry carries each module's group; a **mixed** input set and **both** mismatch directions exit with a named error; the Phase-6 transfer analysis stays the exemption |
| H-11 | An artifact this toolkit writes declares no organelle — `mito annotate` copies the input's empty `source`; the same artifact fails `submission`'s first pass 24× | OPEN, measured 05/10 | every path that writes a record writes a `source` with `/organelle`, proven by reading it back; `mito-prototype` rebuilt and re-pinned |
| H-12 | `gene_names` is declared **shared** but loads only the plastid alias table (0 renamed on a mito record, measured 05/10) | ON HOLD — reopen with plant mito data | the module picks the vocabulary per record — or the §13 row says the *layer* is organelle-aware |
| H-13 | **An unannotated record's IR is only half-checked** — identity, never length (`PZ428230.1`: 15,756 bp vs ~25.7 kb, §9.A-8) | OPEN, named 05/10 | IR length checked against ≥ 3 same-genus records; a length outlier stops the run by name; the corpus record is the negative test |
| H-16 | **Repeats (`repeat-match` + Vmatch) has no module** — milestone **M2.12** | OPEN — own branch when picked up (Phase 2 is closed + merged) | a `repeats` module reports every match ≥30 bp with coordinates and region (the IR pair included), its count a new §8 entry |
