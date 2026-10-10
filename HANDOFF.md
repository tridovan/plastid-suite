# HANDOFF.md — session log

> Purpose: a fresh agent (or a human) picks up here when a session ends or runs out of context.
> **Newest first. Terse. English.** One entry per session: Did · Validation · Residual risks ·
> Blocked · Next. The **Resume point** is the single answer to "where do I start?" — keep it current.
> **An open phase lives on its branch.** `main` is the closed-phase state, one commit per phase;
> the live session docs — this file, `TASKS.md`, `PROJECT.md` — are on `feat/phase-N`
> (`git branch --list 'feat/phase-*'`). See `CLAUDE.md` → Commits.

## Resume point

**Phases 1 – 5, 4a, 3B, 2B, 8 and 8B — COMPLETE and merged to `main`.** **Phase 9 is OPEN on
`feat/phase-9`** (board: `TASKS.md`), and six of its seven tasks are done:

| Task | What it produced |
|---|---|
| **P9-01** (§8.74) | **MIT** for the toolkit (D-p) · **MISA left the package** (D-q) — its zip states no licence — and is obtained now, with the download line and both hashes in `docs/tools/misa.md` · `docs/licensing.md`: 14 tools, licences **measured** from conda metadata |
| **P9-02** (§8.75) | version **1.0.0**, wheel/sdist built and **installed into a clean venv**, golden numbers reproduced from the installed console script |
| **P9-03** (§8.76) | `docker/Dockerfile` real: pinned tools + MISA fetched with its md5, **golden suite run at build time**; the first container run found a real defect — `gene_normalize` could not run from an installed package (**§9.N-1**, fixed with two regression tests) |
| **P9-04 · P9-05 · P9-06** (§8.77) | README **executed verbatim** in a clean venv (five exit codes measured) · `docs/platforms.md` (22 modules × tool × RAM × 3 support levels) · `docs/things_we_do_not_do.md` (12 refusals with their reasons) |
| **P9-07** | **half**: `CITATION.cff` written and validated (CFF 1.2.0) · `docs/release.md` the procedure |

**What Phase 9 still needs is the team, not code** — the four items are on the board's "Waiting on
data/people": **who is named as an author** in `CITATION.cff` (it says "PlastidToolSuite
contributors" on purpose) · a **Zenodo** account, and the **repository is private** while Zenodo's
GitHub integration archives public repos · the **PyPI** name `plastid-suite` · the **how-to video**.
Publishing is outward-facing and permanent (a DOI is minted per version, a PyPI filename is taken on
upload), so it needs the dev's explicit go — `docs/release.md` has the order and the check between
each step. **IPK Gatersleben** may also be asked (one mail) whether MISA may be redistributed; until
an answer, the release does not wait on it.

**Still open from earlier phases, unchanged:** Phase 7's assembly waits on **real Illumina/ONT
reads** (the project's oldest blocker) · the whole **mito round is parked by scope** (Phase 6, M7.3,
M8.2) · **macOS is unverified** — no machine · and the four questions for anh Điềm from the
Crataegus corpus (§8.69–§8.73) are still open, none of them blocking.

## Log

### 2026-10-11 · Claude · dev WSL · session 39 (Phase 9 opened and run: licence → packaging → Docker → release-facing docs)

**Did.** Phase 8 and 8B were already closed on `main`, so the session opened **Phase 9** on
`feat/phase-9` — board first (the `make-tasks` skill, P9-01 … P9-07), then the tasks in dependency
order, each committed on the branch.

* **P9-01 (§8.74) — the licence review, and MISA left the package.** The trigger was a measurement:
  `uv build --wheel` shipped **146 files including `plastid/vendor/misa/misa.pl`**, and MISA's own
  `SOURCE.json` says its zip states **no licence**. No licence means all rights reserved, so
  bundling it was the one act not permitted. `LICENSE` is **MIT** now (D-p), the wheel is **145
  files** with `License-Expression: MIT` (PEP 639), and MISA is **obtained** (D-q): `docs/tools/
  misa.md` carries the URL, the zip md5 and both file hashes with a verified download line, the
  module takes `--set misa=/path/to/misa.pl`, reads the release from the script's own header and
  refuses anything but the pinned 2.1, and `misa.ini` is the toolkit's own file now.
  `docs/licensing.md` inventories **14 tools** with licences **measured from conda metadata**, and
  answers the obvious question (PAML and ARAGORN are GPL — and neither is redistributed).
  `tests/unit/test_licensing.py` was **watched failing on the old tree** before it passed. Golden
  unchanged: 37 = mono 31 · di 4 · tri 2.
* **P9-02 (§8.75) — the package, and a clean-venv install.** Version **1.0.0** (D-r) + classifiers,
  keywords, URLs — and **no Windows classifier**, deliberately. Verified where it matters: a fresh
  venv (`python3 -m venv` has no ensurepip on this machine, so `uv venv --seed`), `pip install` of
  the built wheel, then a real run — **153,724 bp · 37.85 % · 132 = 86 + 38 + 8**, exit 0, on a
  minimum install of biopython + numpy + PyYAML.
* **P9-03 (§8.76) — the Docker image, and the defect only the container could see.** The image now
  carries the pinned tools + MISA (md5-checked) and runs the **golden suite at build time** under
  `PTS_REQUIRE_TOOLS=1`. The first successful container run failed three checks: two for want of
  Vmatch/table2asn (REQUIRE_TOOLS doing its job — Vmatch is installed now), and one **real defect**:
  `gene_normalize` could not run **from an installed package** (`repo_root()` walked up from
  `site-packages` looking for `pyproject.toml`). Fixed with two regression tests — **§9.N-1** — and
  reproduced independently on the wheel from a user's own directory. table2asn is deliberately
  *not* in the image: it is optional and its 56 MB download crawled through the Docker VM, and an
  optional tool must never hang a user's build; its one check is deselected **by name, reason
  printed**.
* **P9-04 · P9-05 · P9-06 (§8.77) — the documents a user actually reads.** `README.md` rewritten
  for a non-expert and then **executed verbatim in a clean venv** (every command, the run-directory
  picture, and all five exit codes — `2` verified on a real amplicon refusal); its status line is
  now a version a test holds to `pyproject.toml`. `docs/platforms.md` puts **all 22 modules** beside
  their tool, their measured cost (scope named) and the three support levels, with the seven
  unmeasured modules named rather than implied — and it found that the registry's `citations` field
  is stored but never read, and has drifted. `docs/things_we_do_not_do.md` is twelve refusals, each
  with its measured reason.
* **P9-07 (partly).** `CITATION.cff` written and **validated with `cffconvert`** (CFF 1.2.0);
  `docs/release.md` is the procedure. **The DOI, the PyPI upload and the how-to video are not
  done** — they need the team's decisions and accounts, see below.

**Validation.** `make check` green — **861 fast + 85 golden** with the pinned tools present
(`PTS_MISA_DIR` now required: MISA is no longer in the tree). The wheel installs and runs in a fresh
venv. The container build runs the golden suite and only builds if it reproduces — **54 passed,
30 skipped (all `.panel/`, each naming its reason), 1 deselected (table2asn, named)**.

**CI found what the local run could not** (§9.N-2): the fix commit cites `§8.76` in its own
comments, and §8.76 was written one commit later — so that commit does not stand on its own, and
per-commit CI said so while `make check` on the working tree was green. The guard
(`test_every_cited_log_entry_exists`) was working; the order was wrong. Rule for next time: an
entry a commit cites goes in **that** commit, and the cheap check is the docs test run **after
staging**.

**Residual risks.** (1) The registry's `citations` drift is **named, not fixed** — the fix is either
to feed it from each module's `Result` or to drop the field. (2) A user who installs the wheel and
runs `gene_normalize` without a reference now gets a message (D-t) — correct per D5-09, but it is a
first-run experience worth watching. (3) The Docker image is built and verified **on this machine**;
no second machine has built it yet. (4) macOS is still untested — no machine.

**Blocked / waiting.** The **release itself**: the authors to name in `CITATION.cff` (the team), a
**Zenodo** account and the repository's visibility (it is **private**, and Zenodo's GitHub
integration archives public repos), the PyPI name check, and the **how-to video** — plus an explicit
"go" before anything is published, because a DOI and a PyPI filename are permanent. Phase 7's
assembly still waits on the raw reads.

### 2026-10-10 · Claude · dev WSL · session 38 (Phase 8B run whole — six tasks, and every answer table put beside ours)

**Did.** The remaining Phase 8B board, task by task, each committed and pushed on `feat/phase-8b`.
The phase's question was "does our reading of a plastome reproduce another laboratory's published
tables?" — and the answer is now per table, with the reason for every difference.

* **P8B-08 (§8.69) — the tightest comparison, and it found the design.** Table S7 (42 PCGs × 18
  samples) **reproduces with zero error on 629 of 714 comparable cells** — under **KaKs_Calculator
  2.0.1 `-m MLWL -c 11`** against **`MN068257` *Amelanchier cusickii***, a record from outside the
  study group. The paper's Results say the comparison is "relative to *C. kansuensis*"; that is
  ruled out twice — against `MF784433` the numbers are wholly different (rpoC2 *C. monogyna*
  0.6667 vs the table's 2.10854), and the table's own structure refutes it (41/42 rows are constant
  within each group of byte-identical CDS, and for rpoB and matK the samples whose CDS *is*
  `MF784433`'s carry a value where that comparison prints `NA`). Getting there took five
  measurements, each recorded: the structural test · the byte-identical CDS test · a 99-record brute
  force · five *Amelanchier* records × 8 groups × 5 methods · then the full 42 genes. **And the
  paper's own named outliers reproduce with our estimator**: the four species it calls rpoC2 > 1
  are exactly the four our yn00 calls > 1 against Amelanchier (1.5263 · 1.3684 · 1.3684 · 1.3684),
  while against *C. kansuensis* **none** exceeds 1 (pooled 0.15–0.27). A finding fell out of the
  last 42 cells: **12 of the 18 deposited records lack the whole IR gene set**, yet S7 carries
  `ndhB` values for 11 of them — the annotation analysed is not the annotation deposited.
* **P8B-05 (§8.70) — the one table that does not reproduce, and it is not the data.** RSCU is 1–2 %
  away (max |ΔRSCU| 0.1039 · mean 0.0148) under our counting rule, and 0.1018 · 0.0106 counting one
  copy per gene name. The difference is systematic — the big codons differ **in the same direction
  in all 18 samples** — and it is classed by record (74–75-CDS records 1.55 %, the complete 83–85
  CDS records 2.6 %, falling to 1.10 % under one-copy counting), so the IR-duplicate rule is part
  of it and not all of it. **Eleven candidate rules were tested and every one rejected with its
  number** (gene-set variants, adding any single gene from the corpus, tRNA/rRNA, intron spans,
  whole genome, averaging rules). That it is *not* a data difference is measured, not assumed: the
  same deposited CDS reproduce S7 to the last digit, and our RSCU side matches cpstools on the
  golden codon for codon (§8.61). Their counting rule is now a question for anh Điềm.
* **P8B-06 (§8.72) — the barcode half is exact, the hotspot ranking is not.** `diversity` (600/200,
  their DnaSP setting) over the 44 plastomes: π 0.001633 whole. Their three named hotspot regions
  rank **4th · 5th · 9th of 91** spacers for us — inside our top ten, but **our first is
  `petG-trnW-CCA` at 0.02042**, more than twice the runner-up; the scope difference (their ranking
  comes from whole-genome sliding windows) is stated with it. The **50 barcode records reproduce
  Table S2 exactly** — length 50/50, GC to the printed digit, max difference 0.0000 — and the
  missing row is `seq45`, the 51st sample, never deposited. π inside the locus: 0.018759 against
  0.00742 for the same spacer across the plastomes.
* **P8B-07 (§8.71) — the biggest, and its composition is the finding.** 18 leave-one-out runs
  through `plastid annotate`, classified with `tools/corpus_bar.py`'s own rule (imported, not
  re-implemented): **341 conclusion-level differences** (missing 21 · extra 136 · extent 184) over
  4,448 published features against 2,339 called — 18.9 per record against the Iridaceae bar's 12.3
  (418/34). **Two of the three components are not our reading being wrong**: 110 of the 136 extras
  are the IR genes the deposited records lack, and 105 of the 184 extents are **six genes of one
  reference record** (`PX413299` annotates `psbM` as 270 bp where the other 17 say 105 — the
  transfer is faithful to its reference, which is the limit P3B-10 named as a residual risk). The
  third is ours, and it is a defect of the annotator: **`trnV-UAC` is never called** — ARAGORN
  *does* find a tRNA at that locus but reads its identity wrong (`tRNA-Stop(tta)` here,
  `tRNA-Cys(aca)` on the golden), the vocabulary cannot name a Stop tRNA, so the call is dropped and
  an intron-containing tRNA is missing from every run. Recorded as **H-24** (ROADMAP M3B.5).
* **P8B-04 (§8.73) — the comparison finished rather than left half-done.** MISA run at **both**
  settings over all 50 plastomes, one genome assigned per section by their own SSR coordinates and
  read for both: **defaults 2,581/3,546 rows exact (72.8 %) and 14/49 sections complete; their own
  Methods' thresholds 3,070/3,546 (86.6 %) but only 4 sections complete** — the two numbers answer
  two different questions and both are stated. The mechanism is measured row by row: in the ZWSZ
  section **9 of 9 non-exact rows are a compound cluster their table writes as one spanning row**
  (`(242,267)` theirs = `(242,252)` + `(256,267)` ours).
* **P8B-09 (M8B.8) — the record.** `docs/corpus_crataegus_2026.md`: the paper and its DOI, the four
  accession groups, **how to prove a copy is the same** (`sha256(gb_data/manifest.tsv)` =
  `0a34e679…`; the whole tree in one line = `4ee47625…`), the labels this corpus gets wrong
  (`PX413286` is *C. chlorosarca* in GenBank and *C. jozana* in Table S1), what it is not (not
  manuscript content; not a benchmark to tune against), and the running defect list.

**Validation.** `make check` green on the branch at each commit; the doc-edits run was **853 fast
passed, 1 skipped**, plus the golden suite (83 passed, 2 skipped — `vmatch` and `table2asn` absent
in that shell) on the tree before the last doc edits. **The phase's closing run then failed once**,
and the failure was real rather than flaky noise: `test_the_workbook_holds_one_sheet_per_table_and_
is_byte_identical` was red in the full run and green in isolation. Measured: **3 exports in 300
were not byte-identical** — openpyxl's `write_workbook` assigns `properties.modified =
datetime.now()` *after* the caller pinned it, so a second-boundary race leaked the wall clock into
`docProps/core.xml` while the module's docstring claimed both stamps were pinned. Fixed in the
written entry (not the workbook), **test written first**
(`test_the_workbook_pins_both_timestamps`), 0/300 after — recorded as **§9.C-38**, the same family
as §9.C-37: a clock assumption inside something that must not depend on the clock.
**The final run is green with every tool present**: ruff clean · **846 fast passed, 1 skipped** ·
**85 golden passed, 0 skipped** (mafft · blast · mummer · paml · aragorn · vmatch · table2asn).

**And CI then caught what that run could not**: the negative-case test was written *while* the
check was already past its lint step, so `ruff format --check` never saw it — red on `main` and on
the branch at the close, for a wrapped line. Fixed on both refs (`style:`, cherry-picked to `main`,
then `main` merged back into the branch so the two trees are identical), and **both are green**:
`main` [38044476678](https://github.com/tridovan/plastid-suite/actions/runs/38044476678) ·
`feat/phase-8b` [38044509941](https://github.com/tridovan/plastid-suite/actions/runs/38044509941).
The lesson is mechanical: re-run the lint as part of committing, not as part of starting a check. `tests/unit/test_project_docs.py` moved 68 →
**73** entries in the same commits, as it requires. Every comparison quotes its own parameters: MISA
2.1 at two settings · KaKs_Calculator 2.0.1 with its argv and method · MAFFT 7.526 · BLAST+ 2.17.0 ·
ARAGORN 1.2.41 · yn00 4.10.10. The corpus's copy is proven by per-file sha256, and the working
manifest (`src/plastid/data/plastid_reference.json`) was **restored and verified clean** after the
leave-one-out runs — the recorded set is one at a time, which is itself a measured limitation.

**Residual risks.** (1) **H-24** is open: `trnV-UAC` has no call in any run, and three fixes are
possible (transfer tRNAs from the reference as CDS/rRNA already do · a second engine · report and
keep) — none chosen, and the phase did not need to choose. (2) **P8B-05's question** — their counting
rule — is the one thing the corpus did not answer; the entry says so rather than smoothing 1–2 %
into agreement. (3) The corpus's own internal inconsistencies are now a list, not a surprise: text
vs Table 1's IR range, one GC cell, S6's two encodings, S7's reference, the RSCU rule, and the
deposited annotations that lack the IR genes. (4) The corpus is third-party published data — nothing
from it may reach a manuscript table (D-k), and every output carries the label. (5) The two species
manuscripts still rest on **stand-in** panels — a real sample means re-running, never patching.

**Blocked.** Nothing in this phase. Phase 7 (real reads) unchanged; the mito round parked by scope.

**Next.** Close **Phase 8B** — `git tag phase-8b feat/phase-8b`, squash-merge to `main` as one
commit — then **Phase 9** (packaging · licence review · release), which now inherits the licence
questions from the cross-check tools (cpstools · KaKs_Calculator · CGAS · Vmatch academic-only).
The four asks for anh Điềm are on the board's Waiting table; none of them blocks Phase 9.

### 2026-10-10 · Claude · dev WSL · session 37 (Phase 8B opened — the Crataegus answer-key corpus)

**Did.** anh Điềm delivered `D:\business\pts\thanhdiem_data\db_test_new`: the paper (Sun et al. 2026,
*Front. Plant Sci.* 17:1767012, CC-BY), **its own supplementary tables** (S1–S7 + 5 data sheets) and
**100 GenBank records** (18 plastomes sequenced for the study, 26 published *Crataegus*, 50 barcode
amplicons, 6 outgroup/reference), with `groups.json`, `manifest.tsv` and a re-fetch script.
(1) **Verified and copied** into `.data/thanhdiem_data/db_test_new/` — 117 files, sha256 per file
against the D: original, byte-identical; 100/100 read by our own reader; all 50 plastomes declare
`/organelle` and pass the IR gate at 100.00 % (§8.66). (2) **Reproduced their Table 1 row for row**
on the 18 new plastomes (e.g. *C. altaica* 159,664 · 87,665 · 19,245 · 26,377 · GC
34.37/30.42/42.65). (3) Found that the paper's **text** IR minimum (26,311) is reproduced by
nothing we measure, including its own Table 1 (26,331). (4) Found that their **SSR table follows
MISA defaults**, not the thresholds their Methods state — 57/57 rows match at defaults, 55/57 at
their stated ones. (5) The first run on the corpus found a defect of ours: **the MISA compound
flag** (§8.65) — MISA writes a cluster two ways and `_parse_gff` read one, so 0 of PX413282's 16
component rows were flagged while MISA counted 8 clusters, and the golden hid it by carrying only
the other encoding. Fixed test-first, module `ssr` 0.3.0; the golden's own numbers are unchanged
(37 = 31 · 4 · 2) and its table now flags its 6 components. (6) Opened Phase 8B with a board of nine
tasks and four decisions, added to `ROADMAP.md`.

**Validation.** `make check` green on the branch (843 fast + 85 golden, before the doc-only edits
that the same-day re-run covers). The golden contract re-measured through the corpus work: 37 ·
31/4/2 · compound 3 (MISA's cluster count) with 6 flagged component rows. The corpus's copy is
proven by per-file sha256 against the source, not by a file count.

**Residual risks.** (1) The paper's Methods and its own tables disagree on SSR parameters (measured);
every comparison must say which setting it used and why (D-l). (2) The corpus is **third-party
published data**: nothing from it may reach a manuscript table, and every output is labelled (D-k).
(3) Two of the comparisons wait on anh Điềm — the 51st barcode sample and the Ka/Ks design — and
P8B-08 must not be quietly skipped. (4) `db_test_new` is a *test* corpus in a family the
manuscripts do not use; the two species manuscripts still rest on stand-ins.

**Blocked.** P8B-08 (the Ka/Ks design). Phase 7 (real reads). **Parked, not blocked:** the mito round
— Phase 6 entirely, M7.3, M8.2 (scope decision 10/10/2026).

**Two comparisons are in, one of them only half-way, and both said something.**
**P8B-03** (§8.67): all **18** of the paper's Table-1 rows reproduce our four sizes exactly
(genome · LSC · SSC · IR), 17 machine-read and the 18th read by hand; per-region GC agrees to the
paper's printed precision except **one cell** — *C. jozana*'s IR GC is 45.62 there and 42.64 here,
while that row's other three GC values and its IR length match exactly, and no other genome in the
table (ours included) leaves 42.62–42.67. Recorded as measured. **P8B-04** (§8.68): at MISA's
defaults our table contains **every row** of 14 of S6's 49 sections — and the other 35 cannot be
compared yet, because their table **encodes a compound cluster as one spanning row** in the
published half and **omits clusters** in the new half. That half-way state is the honest one: D-l is
answered (compare at defaults — they reproduce 57/57 against 55/57 for the thresholds their own
Methods state), and the remaining question is which run produced each section of S6.

**START HERE next session → `TASKS.md` → P8B-08 (Ka/Ks vs Table S7).** The design is settled and
that is new: the paper states it in its **Results**, not its Methods — every species is compared
**relative to *C. kansuensis***, whose record (`MF784433`) is in the corpus, and the paper prints
three of its own values to check ours against (rpoC2 **2.11** in *C. monogyna*, **1.83** in
*C. laevigata*, **1.83** in *C. songarica*; S7 carries 2.10854 · 1.82925 · 1.82925). It is the
tightest comparison in the phase, so it goes first. Then **P8B-05** (RSCU) → **P8B-06** (π hotspot —
the paper names its three most diverse regions — plus the 50 barcode records) → **P8B-07**
(annotation vs GeSeq, the biggest) → **P8B-04**'s remaining 35 sections, which waits on anh Điềm.

**The rule for every comparison in this phase is on the board as D-o**, and it came from anh Điềm
himself: run our own output first and print every parameter with its reason; use their table as an
answer key to **infer** a design they did not state, labelled as inference; **never tune a parameter
until a table matches** — that is fitting, not reproduction, and the entry must say which one it is.
`TASKS.md` also carries a "What the paper states" table, so nobody opens the PDF again for the six
asks: four are answered there, and only the **51st barcode sample** (never deposited — the data
availability statement covers PX676094–PX676143, 50 records) and the **raw reads** are genuinely
absent.

### 2026-10-10 · Claude · dev WSL · session 36 (Phase 8 run whole — five tasks, and the phase closed)

**Did.** (1) **P8-01** CPGview: surveyed and **not adopted** (§8.60) — the repo is a README-only
landing page, there is no pip/conda/CLI to pin, the web app's port does not answer (controlled
against port 80 on the same host and against `portquiz.net:16085`), and the only installable form
is a 2.87 GB container for a figure tool. (2) **P8-02** the CPStools 3.0 cross-check (§8.61),
installed pip-only in its own env with `docs/tools/cpstools.md` in the same session: IR identical,
SSR 37/37, RSCU counts identical codon by codon once both sides count one gene set, π identical for
100 of 104 genes to 5 dp, Ka/Ks differences traced to small Ks denominators — and the comparison
**found a real bug of ours**: the `ssr` module's `Context` column read a trans-spliced gene's
*outer* span, so `rps12`'s 68,975–139,634 made 11 of the golden's 37 rows say "Intron" — including
SSRs sitting inside `ycf1`, `ndhA` and three spacers. Fixed by the rule `core/spans.py` already
states (both ends in one region), regression test written first (`tests/unit/test_ssr.py`), module
0.1.0 → 0.2.0 because the run cache keys on the version, H-08. Recorded as **§9.M-1**. (3) **P8-03**
the CGAS 1.0.1 cross-check (§8.62): structure, RSCU and CDS introns identical; SSR differs by one
**threshold** parameter (its defaults 54 vs MISA 37; with `-t "10,6,5,5,5,5"` it prints 37 =
27·10·0 · 31/4/2, identical). Two findings: CGAS prints `LSC/SSC/IR = 0` for `NC_024936.1` with no
warning while still counting it "analysed" (our paths find IR 26,184 @ 100.00 % there), and CGAS
prints **8 tRNA introns** the golden really has and **no module of ours publishes** → **H-23** on
ROADMAP's M2.7. R is absent here, and that measured something worth keeping: CGAS modules 8/9/12/13
compute numbers in Python and call R only for plots. (4) **P8-04** the two modules
`docs/M8.3_scale.md` had honestly left *not measured* (§8.63): `diversity` **13 min 14 s · 138 MB**
(π 0.040086 over 120,964 of 191,707 columns), `barcodes` **12 min 27 s · 138 MB** (859 windows, 13
hypervariable regions, 10 candidates, 8 found) — both finished, so the old prediction was wrong in
the useful direction: it came from a run *stopped at a cap*, and a cap is not a measurement. (5)
**P8-05** CI (§8.64): two jobs, every tool pinned to §7.8's versions, `PTS_REQUIRE_TOOLS=1`,
table2asn fetched and md5-checked. Its first run was red on **three real defects** — scipy 1.18
returning `pvalue=nan` where 1.17 returns 1.0 (and `kruskal` returning nan on *both* builds) while
the `stats` module exported that nan into its FDR table; a test hardcoding `/root/toolbox/bin/
table2asn`, which raised `PermissionError` on the runner instead of skipping; and the micromamba
env missing from `PATH` in a `run:` step. All fixed, then the build was deliberately broken to
watch it fail on the golden number (`37.863555723244254 == 37.8536 ± 0.001`) and reverted. (6)
Closed the phase: tagged `phase-8`, squashed to `main`.

**Validation.** `make check` green on the phase's tree. The golden contract re-measured through
the cross-checks: Ts/Tv 2.6917 · SSR 37 · 132 = 86 + 38 + 8, and every cross-check run recorded
both tools' names, versions and argv. CI green on `feat/phase-8`
(https://github.com/tridovan/plastid-suite/actions/runs/38025149827) after being red both by
accident and on purpose.

**Residual risks.** (1) **H-23** — the `introns` module publishes CDS introns only; extending it
to tRNA introns touches `core/spans.py` and the π-per-intron column of `diversity`, so it is a
task, not a patch. (2) **Vmatch is academic-only** and present on this machine alone (M9.3 must say
so in the release notes), and `cpstools`/`KaKs_Calculator`/CGAS are now declared cross-check tools
— three more licences for that same review. (3) The **tier-B CI job** installs from conda-forge
and bioconda at run time: a channel-side change could turn it red without anything in this repo
changing, and the job prints every pinned version so that failure is readable. (4) The species
manuscripts still rest on **stand-in** panels — a real sample means re-running, never patching.

**Blocked.** M7.1 · M7.5 · M7.6 (real reads). **Parked, not blocked:** the mito round — Phase 6
entirely, M7.3, M8.2 (scope decision 10/10/2026).

**Next.** Phase 9 (packaging, licence review, release) is the next unblocked phase; before it, the
H-23 decision and the `docs/tools/` licence list are the two things it will want.

### 2026-10-10 · Claude · dev WSL · session 35 (Phase 2B run whole — three tasks, and the phase closed)

**Did.** (1) **P2B-01**: the intron and IGS span derivation moved out of `modules/introns.py` into
`core/spans.py` — DRY's second real copy, as H-07 said it would be — and `diversity` now measures π
per intron and per intergenic spacer over exactly the spans `introns` publishes (§8.57).
(2) **P2B-02**: the **repeats** module (M2.12, H-16), the last analysis box without one: maximal
exact matches ≥30 bp, dispersed and palindromic, each with both copies' coordinates and plastid
region, running on MUMmer's `repeat-match` with **Vmatch 2.3.1** as the independent second engine
(§8.58). (3) **P2B-03**: both species' manuscript drafts re-run — now with π, Ts/Tv, repeats and
figures (§8.59). (4) Closed the phase: tagged `phase-2b`, squashed to `main`.

**Validation.** `make check` green: **833 fast tests + 84 golden** (1 skip: `table2asn` is not on
this machine). The golden contract re-measured this session: Ts/Tv **2.6917**, SSR **37**,
**132 = 86 + 38 + 8**. The two engines agree on the golden row for row (9 matches, same
coordinates, same classes) — a test, not a story. The species tables reproduce the 03/10 draft's
numbers **cell for cell** (131/133 genes, SSR 32/32, introns 16/17) and add the parts it lacked.

**Residual risks.** (1) **Vmatch is academic-only** and exists on this machine only; the
cross-check test skips with a named reason elsewhere, and M9.3's licence review has to say so in
the release notes. (2) **CPGview is still unnamed in §7.8** — it is the first thing Phase 8 needs
(above). (3) The two species' panels are **stand-ins**: a real sample means re-running, never
patching numbers. (4) The 11-genome two-genus comparison runs stay in `.runs/phase2b/` as what was
measured first; the draft's numbers come from the **genus-restricted** runs — quoting the former
without the asymmetry finding would be quoting a measurement whose IR halves disagree by 3.4×.

**Blocked.** M7.1 · M7.5 · M7.6 (real reads). **Parked, not blocked:** the mito round — Phase 6 entirely, M7.3, M8.2 (scope decision 10/10/2026; this round is plastid only). **Next.** Phase 8, board P8-01.

### 2026-10-09 · Claude · dev WSL · session 34 (Phase 2B opened: the last unblocked analysis work)

**Did.** Closed Phase 3B (P3B-12: tagged `phase-3b`, squash-merged to `main` as one commit — the
merge conflicted in `CLAUDE.md` with `main`'s own "one commit per phase" edit, resolved by keeping
main's changes and taking the branch's newer annotate paragraph), then opened **Phase 2B** on
`feat/phase-2b` with a board built by the `make-tasks` skill: three tasks (π per intron/IGS ·
the repeats module · the manuscript draft re-run), two decisions, five pitfalls, the waiting-on
table.

**Validation.** `make check` green *on main's own tree* before the push (805 fast + 79 golden).
The board's self-check is done by hand, as the skill asks: every task has a runnable Done-when,
every blocker names a dataset or a person, nothing finished is left on it, 89 lines.

**Residual risks.** (1) P2B-02 needs **Vmatch**, declared but not installed here — the install
ships in the session that first uses it (§7.8), and its academic-only licence keeps it out of the
minimum install. (2) P2B-03's output stays a **stand-in** draft until real samples arrive; a real
sample means re-running, never patching numbers. (3) **CPGview** is named in M8.1/§13.8 with no
§7.8 record — named here rather than silently inherited.

**Blocked.** Phases 6 and 7 (nuclear genome · plant mito data · real reads). **Next.** P2B-01.

### 2026-10-09 · Claude · dev WSL · session 33 (P3B-10 done — the bar is met; and P3B-11, the two annotators retire)

**Did.** (1) **P3B-10**: the exon-structure transfer (`core/exon_structure.py`), the stop codon read
past an alignment that stopped short, and H-20 — the corpus bar went 681 → **418 against PGA's
588**. (2) **P3B-11**: retired `annotation`'s chloe · pga set (9 files, 2,406 lines).

**Validation.** `make check` green (lint · 778 fast tests · golden). The corpus: 34 records
leave-one-out, **30/34 better, 4 worse, 0 tied**. The golden: 132 = 86 + 38 + 8, **35** differences
listed and pinned as a set. The mito golden tripwire moved **18 → 19** when the stop-codon rule
reached the shared core — read before touching it, and it is `rps1` gaining its stop codon to match
the published record; §8.22's 18/31 stays as what it measured then.

**Residual risks.** (1) **The verdict rests on the 1 bp tolerance** — quoting "418 against 588"
without the strict line (1070 against 699) is quoting half the measurement; it is stated in
TASKS.md, §8.55, HANDOFF, D5-08 and CLAUDE.md. (2) The transfer is faithful to the chosen
reference, so where a reference's own structure differs from the published record by 2–3 bp at an
intron boundary, our call differs too — classified, not queued. (3) The corpus is not a test: it is
`make corpus-bar`, minutes long, outside `check`.

**Blocked.** Nothing. **Next.** P3B-12 — tag, squash-merge, one commit.

### 2026-10-09 · Claude · dev WSL · session 32 (P3B-10: the corpus measured with the instrument in the repo — bar not met, and the golden's 86 is right by accident)

**Did.** (1) **Re-ran the 34-record corpus** with the code from session 31's two fixes (34/34,
0 failed runs, `--exclude NC_024936.1`); the composition changed exactly as predicted: **missing
167 → 47, extra 161 → 41** — the whole IR-copy family gone. (2) **Committed the measuring
instrument** (`tools/corpus_bar.py`, `make corpus-bar`, `19 tests` in
`tests/unit/test_corpus_bar.py`, `pythonpath` in `pyproject.toml`) with its rule stated in one
place. (3) **Declared and justified the boundary tolerance** (1 bp): 719/1,264 matched pairs move
exactly one base, 634 of them with our call containing the record's, and ARAGORN run by hand on the
golden returns the coordinates we emit — two methods ending a tRNA, not two conclusions. Applied to
both readings; the no-tolerance totals printed beside it. (4) **Ran the bar**: **811 against
PGA's 588 — NOT MET under either reading** (0/34 records better). (5) **Found three defects, none
fixed**: `lhbA` is a *duplicate call* at the same span as `psbZ` (33 of 41 extras) and is what
makes the golden's 86 — the published 86 counts `rps12` twice, ours counts `lhbA` twice, so both
read 85 loci and **the golden number is right by accident**; `rps12` is missing its 26 bp middle
exon and its whole IRa copy; the first exons of `petB`/`petD`/`rpl16`/`rps16` are below blastx's
seed (136 of the 242 moves > 30 bp).

**Validation.** `make check` green at the commit: **818 fast passed** (36 deselected) · **77 golden
passed / 4 skipped** · ruff clean. The corpus numbers regenerate with `make corpus-bar` (needs
`.panel/`, `TOOLS=/root/toolbox/env/pts/bin`); the harness's PGA pairing is by **sequence digest**,
so a broken pair stops rather than comparing two different genomes.

**Residual risks.** (a) The rule has a threshold in it, so it is reported with its sensitivity, not
just its verdict — at 3 bp the bar would pass without any fix, which is exactly why the bar is read
at 1 bp. (b) Writing the rule as code exposed two holes in it (envelope-only comparison; strand not
compared) — both fixed with tests, but it says the rule was never as settled as it looked. (c) The
stale first corpus run's numbers are still in §8.50; §8.51 says which run supersedes them.

**Blocked.** Nothing. H-04 (cô Trang) and H-21 (a rule change) stay open by their own rows.

**Next.** In order: the reference's **exon structure** in place of the alignment's span (fixes the
four first exons *and* `rps12`); **H-18** together with the `rps12` IR copy (they move the golden's
count together); **H-20**; then the golden's difference list as a test. The commit `e2c268c`'s
subject says "session 31" — this entry is session 32; the subject is off by one and was left
uncorrected rather than force-pushing a pushed commit.

**Also done in the same session (§8.54, H-18 + the `rps12` IR copy).** `lhbA` and `psbZ` are one
locus by all three §8.9 signals, measured exact this time — **the same `/product` string**, the
same 189 bp / 62 aa, the same neighbours `psbC … rps14`; 33 pinned records write `psbZ`, one
(`NC_042827.1`) writes `lhbA`. The alias takes the pinned vocabulary **80 → 79 names**, and the
caller stops writing both at the golden's 35895-36083. The `rps12` half is the mirror rule:
`_ir_copies` now mirrors the exons of a **straddling** gene and keeps the shared exon in both
features, which is how the published record writes a trans-spliced gene across the IR — the
counter-case `ycf1` (one call reaching into the repeat gets **no** second copy) is pinned in
`tests/unit/test_ir_copies.py`. **The measured cost, named rather than hidden: 22 of 34 published
records write `rps12` once**, so our second copy reads as an extra there; the corpus moved 692 →
**681** (missing 47 → 34, extra 41 → 30) and PGA stays 588. Even if that representation difference
were excluded on both sides it is 659 against 576 — the bar is not met, and §8.54 says so.

**Added later in the same session (§8.52).** The corpus measurement itself found a **silent-wrong**
and it is fixed: ARAGORN's intron notation **`i(a,b)`** has `b` as the intron's **length**, not its
end. Read as an end (what the parser did from 09/10), every spliced tRNA on the golden came out
**104–134 nt** — no tRNA is that long — and `trnL-UAA`'s exons came out 35 + 85 where the published
record's own join is **35 + 50**. Read as a length they come out 72–95 nt, and `trnL-UAA` (85) and
`trnA-UGC` (73) land **exactly** on the record. The test that asserted the wrong reading was
written from the same misreading; it is **corrected with the numbers in it**, not loosened, and a
second test pins all three cases to a plausible tRNA length. Re-measured: **811 → 692** (the 119
gone are 57 + 27 + 34 + 1, the intron tRNAs), PGA re-measured over all 34 records (**588**, and the
earlier run had missed `PX516262.1` because PGA names that file after the organism — the harness now
pairs by **sequence digest**). **The bar is still not met: 692 against 588.** Measuring the introns
also produced **H-22**: ARAGORN's `-i` searches introns **≤ 3,000 bp** and says nothing when it
skips one, and the longest measured here is **2,671 bp** — a 12 % margin, named rather than
silently patched.

### 2026-10-09 · Claude · dev WSL · session 31 (P3B-10 opened: two silent-wrong fixes, and the corpus bar measured — not met)

**Did.** P3B-10 is **in progress**, not done. Three measurements, two of them fixes:

* **ARAGORN skips intron-containing tRNAs in silence** unless `-i` is passed — its own help says
  so and it reports nothing when it skips them. The golden read **30 tRNA where the published
  record has 38**, and the eight missing were exactly the intron ones (`trnK-UUU` · `trnL-UAA` ·
  `trnV-UAC` · one `trnG` · both copies of `trnI-GAU` and `trnA-UGC`). With `-i` and the intron
  parsed into the two exons a `join` needs, the artifact is **132 = 86 CDS + 38 tRNA + 8 rRNA** —
  the contract. `-i` is a *caller's* flag, not a default: on the mito answer-key record it takes
  17 tRNAs to 21, four of which the published record does not have.
* **The IR copies were written one IR length from where they belong** (26,214 bp — outside the
  repeat), because the inversion was computed from the first repeat's *start* instead of its
  *end*. The gene **count** stayed 86, so every count-based check passed. Fixed; `rpl23` · `ndhB` ·
  `ycf2` · `rps7` · `rrn16` · `rrn5` now match the published record base for base, and a new golden
  test checks the mirrored copies **against the record**, not against their own arithmetic.
* **The corpus ran leave-one-out, 34/34, 0 failed runs** (`--exclude NC_024936.1`, the audit gate's
  condition). Classified against the published record: **ours 1,579 conclusion-level differences ·
  PGA 637**. The bar is *not* met. The split says why the two mechanisms differ: ours has far more
  presence errors (167 missing · 161 extra) and PGA far fewer (39 · 4), while PGA's differences are
  almost all extent (594) — PGA copies coordinates from its reference, this engine derives them
  from the alignment.

**Validation.** `make check` green: **799 fast · 77 golden passed / 4 skipped**. New tests: the
intron split (both strands, plus the wrapped case), a `tRNA-Stop` line reported rather than fatal,
and the IR-copy guard.

**Residual risks.** The corpus composition is unmeasured: what the 161 extra and 167 missing calls
*are* is the next question, and until it is answered the phase's bar cannot be argued. Four genes
lose a 6–40 bp first exon (`rps16` · `petB` · `petD` · `rpl16`) — blastx cannot align two codons.
`rps19` −54 bp is **H-20**; `lhbA` is **H-18**; eight tRNA names follow ARAGORN's anticodon reading
on intron genes and are **H-04**'s family.

**Blocked.** Nothing.

**Next.** P3B-10: classify the extra/missing calls, decide what to do about the tiny first exons,
commit the classification harness (`/tmp/p3b10_diff.py` — `compare` · `summarise`), then the
phase's own bar can be argued.

### 2026-10-09 · Claude · dev WSL · session 30 (P3B-09: `plastid annotate`, and the gate in front of the run)

**Did.** **P3B-09** (M3B.4 · H-09) — the annotate command, the batch, and the organelle. Five
pieces, each with a measured reason:

* **`plastid annotate`** — 1..N samples × 1..N reference records. A batch is a **batch of runs**:
  every sample goes through `_execute_run` (seam 3) into `<out>/<accession>/<run-id>/` with its own
  manifest, report and artifact. A failed sample is named, its own run directory records
  `status: failed` + the reason, the batch continues, and the exit code is non-zero. The manifest
  gained an optional `dataset` key (schema 2 unchanged) so a run that failed before writing
  anything still says which sample it was.
* **The audit gate, in the module, before the first BLAST.** `_reference_set` audits exactly the
  records the run will search (H-19 named this hole, §8.48). Measured on the real pinned set:
  exit 2 in **13 s**, two `file — reason` lines (`NC_024936.1.gb — CDS rps19: partial coordinates
  (>83020)` / `(<153124)`), no artifact. `--exclude <accession>` is D5-08's "drop the bad file" as
  a mechanism — declared, recorded in the manifest, named in the reference label and in the run's
  notes; a name matching nothing **stops the run**. `--audit-only` (P3B-03's flag) audits and stops.
* **H-09 done.** `ModuleSpec.group` (`shared` · `plastid` · `mito`, **required**) on all 22 modules
  per §13's table; `core/organelle.py` reads `/organelle` and refuses **both** mismatch directions
  **and** a mixed set, by name, pointing at the command that serves it. A bare sequence declares
  nothing and is not guessed (D5-06: the path declares).
* **The `.docx` report had never run**: `citations=("blast", "aragorn")` with neither key in
  `data/citations.yaml` → `CitationError`. Both added, and the engine had **no report template** at
  all — `report/templates/plastid_annotate.en.md` written.
* **A directory parameter is content-addressed** (P3B-05's promise, one level up): `--set
  reference=<dir>` was hashed as a string, so editing a reference record served stale tables.

**Validation.** `make check` green: **797 fast · 76 golden passed / 4 skipped** (the four are the
missing external tools: PGA.pl, chloe · julia, table2asn). New tests: `tests/unit/test_organelle.py`
(15 — the declaration, both directions, a mixed set, the registry's groups) ·
`tests/unit/test_cli_annotate.py` (6 — the mito sample refused by name, a mito FASTA refused through
its reference set, a mixed batch stopped before anything is annotated, `--audit-only` needs no
sample, a failed sample marked without killing the batch, one artifact per sample) · three golden
tests (the gate stops a real run on the real blocked record · a drop is named · an unmatched
exclusion stops). The three existing golden checks now run with `--set exclude=NC_024936.1` (the gate
forces the drop) and **still hold** — `ycf1` 122.984–128.404 unchanged.

**Residual risks.** (1) **H-20** — dropping the blocked record changes exactly one call of the
golden: `rps19` 279 bp → **225 bp**, because three pinned records carry a truncated `rps19`
(225 bp/74 aa) and one of them wins the tie. Structure-first *between* candidates is the piece H-19
did not implement; it is P3B-10's to classify and a fix of its own. (2) **H-21** — the audit's
partial-coordinate BLOCK refuses `NC_024936.1`, whose `rps19` protein is **complete** (279 bp/92 aa);
that very record produced the published call. The rule's stated reason does not hold for this shape.
(3) The `--exclude` mechanism is a way to *drop*, not to *waive* — but it is the user's act, so a
run's fidelity depends on it being written down; it is, in the manifest and the notes.

**Blocked.** Nothing. The ★ pick is still the end user's (§8.44's provisional override stands).

**Next.** **P3B-10** — the phase's Done-when: the golden to exactly **132 = 86 + 38 + 8** with every
difference listed, and all **34 records leave-one-out**, diffed against the published record *and*
against PGA's reading, every difference classified, with fewer conclusion-level differences than
PGA. It starts from two named differences (the tRNA count — 30 called vs 38 published, the export
mirrors CDS and rRNA but not tRNA — and **H-20**), and its corpus runs carry `--exclude
NC_024936.1` because the gate now stops on it.

### 2026-10-09 · Claude · dev WSL · session 29 (H-19: a locus a record annotated twice is searched once)

**Did.** **H-19**, the conclusion-level difference P3B-08 named: arbitration was identity-first
where D5-08 §3 says structure first, and the golden's `ycf1` came out 900 bp from a reference whose
own `ycf1` is a truncated copy of the gene.

* **The rule is structural, and the fix is not in `_choose`.** A candidate that is never indexed has
  no score to win on, so it lives in the reference index: a locus **one record** annotates twice,
  where the short feature is under 50 % of the long one, is a truncated duplicate — excluded, and
  **named** (`truncated_annotations`, one run note per row; 26 on the golden).
* **`0.5` is the middle of a measured gap, and the comparison stops at one record because that is
  where the measurement works.** CDS only, the pinned 34 records: 243 loci annotated twice — 216 at
  1.000, one at 0.9927 (`rpl2` 271 against 273, a real difference and **kept**), 26 below 0.30 — and
  **nothing at all between 0.2992 and 0.9927**. Across records the same test separates nothing: only
  3 of 80 CDS loci fall below 0.5 and all three are the within-record cases again, while `rpoB`
  runs 614–1,070 residues and `ndhD` 383–523. Within one record the two features are the same gene;
  across records a short ortholog and a truncated annotation are one and the same number.
* **Measured on the golden: exactly one line of the called set changes** — `ycf1` from
  108,448–109,347 (`+`, 900 bp, 98.83 %, `NC_088562.1`) to **122,984–128,404 (`−`, 5,421 bp,
  97.68 %, `NC_081929.1`)**, which is the published `ycf1 [122983:128404](-)` base for base. The
  artifact still reads **86 CDS + 8 rRNA + 30 tRNA** with `/organelle="plastid"`, and the 33 `ycf1`
  candidates left are all complete (1,783–1,862 residues).
* **The Done-when's second clause was not implemented, and a measurement is why.** "an alignment
  that does not reach the reference's ends" would not have fixed this case — the truncated
  reference's alignment *does* reach both of its ends (coverage 1.000, the very number that misled
  the arbiter) — and it is the wrong mechanism: a partial query gene is already expressed by the
  call's own `partial_start`/`partial_end` flags, so preferring a reference in order to avoid one
  would hide a finding. Recorded in §8.48 so P3B-10 reads the rule's real scope.
* **A gap found while measuring, named into P3B-09.** `plastid_annotate` never calls the reference
  audit: `plastid reference audit` **BLOCKs** `NC_024936.1` for a partial `rps19`, and that record
  is still a candidate in a run on the golden. The phase's goal is "an audit gate in front of every
  run"; the wiring is now a line in P3B-09's checklist.

**Validation.** `make check` green: **775 fast · 73 golden passed / 4 skipped** (the four skipped are
the missing external tools: PGA.pl, chloe · julia, table2asn). New: four unit tests on the rule (it fires · a
2-residue difference is not a truncation · the IR's equal copies stay two entries · **a lone short
feature is still searched**) and one golden negative test asserting both the coordinates and that
the run names the feature it did not search. **The 65 mito tests pass** and no mito file is touched.
Recorded as §8.48 · §9.L-4.

**Residual risks.** (1) A record that annotates a locus **once**, briefly, is not caught by this
rule — measured: no such case in the pinned set, and the audit is the gate that catches the class
by coordinates (`<`/`>`) once P3B-09 wires it in. (2) The run's notes grow 30 → 56 lines, 26 of
them the new rows; the `rps12` envelope lines are the pre-existing bulk (P3B-08's envelope-vs-exon
issue, in the *notes* rather than in the table).

**Blocked.** Nothing new.

**Next.** **P3B-09** — CLI & batch, with the audit wiring.

### 2026-10-08 · Claude · dev WSL · session 28 (P3B-07 and P3B-08: the calling core, and the annotator runs)

**Did.** Two tasks, in order.

* **P3B-07 — the shared calling core (§8.46).** The second real copy arrived, so the machinery was
  **extracted, not copied** (the project's own DRY rule): `core/caller.py` holds the rules,
  `core/mito_caller.py` shrinks to what is measured about the mitochondrion (floors 50 %/0.5, table
  1), `core/plastid_caller.py` is the plastid half (reference index from GenBank records, floors
  50 %/**0.7**, table **11**, the self-blast IR reading). Floors and table are **required**
  arguments in the shared core — a path that inherited another path's measured floor would quote a
  number it never measured. **43 mito tests pass with not one line changed.**
* **P3B-08 — arbitration and export (§8.47).** `modules/plastid_annotate.py`: the pinned set with
  its leave-one-out rule, `blastx` over the references' own `/translation`, `blastn` for rRNA, a
  self-search for the IR behind the gate, ARAGORN for tRNA. Measured on the golden: **86 CDS +
  8 rRNA + 30 tRNA** — and 86 is a **rule**: one call per locus (79) plus the seven loci inside the
  IR, which the **export** mirrors (`_ir_copies`), the same way the published record writes them.
* **Three silent failures fixed, all found by measurement rather than by a test** (§9.L): a
  `-max_target_seqs 500` constant that truncated the plastid search from 4,382 HSPs / 2,603
  subjects to 929 / 467 **without saying so**; a gene name with a space (`23S ribosomal RNA`)
  cutting the FASTA header below the three-field contract; rRNA product-sentences reaching the
  vocabulary as gene names. And one lesson paid inside P3B-08: **overlap must be measured by exon,
  not by envelope** — `rps12`'s envelope spans 68,976–98,016 and made it "overlap" 27 genes,
  burying the real pairs.

**Validation.** `make check` green twice: **771 fast tests · 72 golden passed / 4 skipped**. New:
`tests/unit/test_plastid_caller.py` (12), `tests/golden/test_plastid_annotate.py` (3, one real run),
`test_blast.py`'s cap test, `test_gene_names.py`'s rRNA-sentence test. Measured on the golden:
`psbC`×`psbD` **53 bp** · `ndhC`×`ndhK` **10** · `atpB`×`atpE` **4** — three of §8.34's four exact;
`competing_references` marks exactly one chosen per locus; the artifact re-reads and carries
`/organelle="plastid"`.

**Residual risks.** (1) **H-19** — arbitration is still identity-first (D5-08 §3 says structure
first): the golden's `ycf1` came out **900 bp** from a reference whose own ycf1 is truncated, where
the published gene is 5,421 bp. This is conclusion-level and P3B-10's bar counts those. (2) The
export mirrors CDS and rRNA but **not tRNA**, so the count is 30 against the published 38. (3) The
coverage floor 0.7 is §8.41's **corpus** number, named: the golden cannot choose it (0.5–0.8 give
the identical called set) — P3B-10 re-checks.

**Blocked.** Nothing new. The ★ pick still waits on the end user.

**Next.** **P3B-09** — CLI & batch, then **P3B-10** — golden + corpus validation.

### 2026-10-08 · Claude · dev WSL · session 27 (P3B-06: Ts/Tv re-measured — the number holds, and now we know why)

**Did.** Re-measured §6's Ts/Tv under the naming direction P3B-04 flipped — the task's own words
were "a prediction until then", and P3B-04's log said "chưa đo" in as many words.

* **Measured in the pipeline order, not by re-running the older test.** The existing CDS-scope
  golden compares the two *published* records, which is §8.27's measurement of the *old*
  direction; under D5-09 a sample inherits vocabulary A **through** the naming step. So: ran
  `gene_normalize` (`mode=normalize`, the golden as reference) to get `NC_024936.1.normalized.gb`
  — the six renames of §8.44 — then `substitutions --set scope=cds` on **that artifact**. Result:
  **79 shared genes · 68,448 columns · 68,382 comparable · 443 substitutions · 323 Ts / 120 Tv =
  2.6917** — identical to §8.27 to the unit (§8.45).
* **And why it could not have moved.** The alias table maps every member of a locus to that locus
  whichever member is canonical, so the flip **relabels** the shared set instead of re-selecting
  it; it could only re-select if a record carried two *different* genes whose names fall in one
  locus. Measured on both Iris records: 86 CDS → 79 canonical names each, **0** loci reached by two
  spellings. The invariance is a measured property, not a hope.
* **Labels move, numbers do not.** `substitutions_by_gene` now prints `ycf3` · `ycf4` · `psbN` ·
  `clpP` · `ndhL` · `atpI` · `psaI` · `psbB` where §8.27 printed the old spellings for the same
  loci — vocabulary, not conclusion.
* **New guard.** `tests/golden/test_tier_b.py::test_the_naming_step_does_not_move_the_golden_ts_tv`
  runs that order and asserts the artifact **is** in A (`atpI` · `psaI` · `psbB` present), so a
  naming step that silently did nothing cannot pass.
* **A finding, named not fixed.** Reading that artifact: **6 `gene` features still spell
  vocabulary B** (`atpL` · `psaL` · `bspB` · `trnk-UUU` · `trnl-GUA` ×2) over the same spans as the
  renamed CDS/tRNA — `GENE_KINDS` excludes `gene` from every renaming decision. Both published
  records agree 132/132 at the same spans, and CGAS (D5-08's model) renames them too. **H-17**,
  §9.K-1 (a new §9 subsection, so the item numbering rule now reads D–I and K).

**Validation.** `make check` green: **756 fast tests · 69 golden passed / 4 skipped** (PGA ·
chloe+julia · table2asn missing; MAFFT, MISA, BLAST+ on `PATH`). The four named goldens —
`substitutions` (both scopes), `kaks`, `popgen` — pass. `tests/unit/test_project_docs.py` updated
in the same commit: entry count 44 → 45, §9.K added to the restart-at-1 guard.

**Residual risks.** (1) The ★ is still the recorded **override** (`provisional: true`, the golden),
so the printed vocabulary is the dev's provisional choice; a different ★ re-runs this measurement
with the same result only if that record's partition is clean too — the property is measurable per
record and H-17's fix is the other place it matters. (2) H-17 is *display* today: `GENE_KINDS`
keeps `gene` features out of every count and every comparison module, so no number is wrong yet.

**Blocked.** Nothing new. The ★ pick still waits on the end user.

**Next.** **P3B-07** — the calling engine (extract the core from `core/mito_caller.py`).

### 2026-10-08 · Claude · dev WSL · session 26 (P3B-04: one vocabulary, and the alias direction flipped)

**Did.** `gene_normalize` — the naming step D5-09 puts in front of the annotator.

* **The module (§8.44).** Arity 2, `requires_ir=False`, modes `normalize` / `verify`. `genomes[0]`
  is the naming reference (A), `genomes[1:]` the records (B); B's spellings appear in exactly one
  table, in their own column. Four tests in `core/gene_names.py` with §8.41's measured thresholds
  (identity ≥ 50 % · coverage ≥ 0.7 · margin ≥ 10 bits), three statuses (`measured` ·
  `not measured` · `outside vocabulary`) and an `outcome` saying what the step did
  (`confirmed` · `renamed` · `conflict` · `reported`).
* **The reference is a declared parameter.** `ModuleSpec.reference_param`; the CLI resolves
  `naming_reference` (or `data/naming_reference.json`, the recorded **provisional override** —
  the golden), **audits it before use**, makes it `genomes[0]` and records
  `{path, accession, sha256}` in `modules[].reference`. `MANIFEST_SCHEMA_VERSION` → **2**; a
  conflict stops with exit **4**.
* **The alias direction flipped** (D5-09 §8): 5/9 loci now print the reference's spelling, and
  `oriented_aliases` ties the direction to the record a run hands it. Eleven guard tests updated
  to the new direction (same claims, other names).
* **Measured on *I. gatesii* against the golden**: 132 genes · 6 renames · 121 confirmations (82
  margin · 38 neighbours · 1 sequence) · 0 conflicts · 5 `outside vocabulary` · 2 `not measured`.
  The five are tRNA names the reference does not use — one of them `trnL-CAU`, **100 % identical**
  to the reference's `trnI-CAU` — and §8.41's rule says a ranking may not settle a tRNA, so they
  are reported and not renamed.

**Validation.** `make check` green: **755 fast tests · 68 golden passed / 4 skipped** (PGA ·
chloe+julia · table2asn; MAFFT and BLAST+ on PATH). New: 29 tests in
`tests/unit/test_gene_normalize.py` (thresholds · orientation · four tests with one negative case
each · the three statuses · IR copies are not duplicates · hits collapsed by locus · verify writes
nothing · the corpus pair end to end · conflict → exit 4 with coordinates · `on_conflict=continue`
· the CLI's resolve → audit → `genomes[0]` → manifest path · one tier-B test running real blastp).

**Residual risks.** (1) **Ts/Tv is still a prediction** — the flip changed five canonical keys and
P3B-06 has to re-measure it; the goldens pass, which is a test result and not a measurement.
(2) The recorded standard is an override, so the printed vocabulary is the dev's provisional
choice until the end user picks from the audit report. (3) Test 4 (neighbours) is only as good as
the IR frame: without one, a tRNA whose ranking disagrees is reported rather than settled.
(4) `duplicate_names` is structural (IRa/IRb); a record whose IR is undetectable falls back to
exact sequence equality, which is stricter than it should be.

**Blocked.** Nothing new. The ★ pick still waits on the end user.

**Next.** **P3B-06** — re-measure Ts/Tv under the new direction — then **P3B-07**, the calling
engine.

### 2026-10-08 · Claude · dev WSL · session 25 (P3B-03: the reference-set audit, and the two readings that have to agree)

**Did.** The gate that stands in front of every annotation run, as a run of its own.

* **The audit (§8.43).** `plastid reference audit --panel DIR` — six BLOCK classes (unreadable ·
  no gene features · a CDS without `/translation` · a partial coordinate · the IR gate · a CDS that
  does not translate) and two REPORT classes (an inverted SSC · a gene set that differs from the
  majority). Every finding prints as `file — reason` and nothing else; a block exits 2 and nothing
  is annotated; the reason lands in `manifest.json`. The report measures both reports over **every
  record that could be read** — §8.34's 32/34 is about the corpus, not about the survivors.
* **The second reading.** MAFFT `--adjustdirection` over the SSC: 34 sequences · 614,658 bp ·
  17.6 s, and it flags **exactly** the two records the gene-strand rule calls inverted. It is
  *injected* into the core (core may not import `external/`): the CLI passes
  `external.mafft.adjust_direction`, or passes nothing and the report says **not measured**.
  A disagreement would be its own finding; there was none.
* **The picture.** `dotplot` (figure #10) draws every 15-mer run of a deviating record against a
  record the set agrees on, both orientations. Measured: 349 anti-diagonal runs for
  `NC_042827.1` against 20 forward, 332 against 13 for `NC_062351.1`. It reads the
  `region_matches` table, so `plastid figures <run-dir>` re-renders it without re-running the audit
  — and it is drawn **without** MAFFT too: a machine with no aligner must not lose the one view
  that lets a person confirm a flip-flop.
* **Exceptions.** The pin is written by the tool from the file's **current** digest
  (`--confirm-exception FILE:CHECK --confirmed-by … --reason …`), never typed; it can waive a
  REPORT and never a BLOCK (`Finding` refuses to build one; `read_exceptions` refuses to read one).
  Exercised end to end on the real set: pin applied → 21 findings instead of 22 and the row marked
  confirmed; a changed byte → **void**, and the finding returns.
* **A data-model change the gate needed.** `Gene` now carries the record's own `/translation`
  (read, never written — D5-08 §8.21), and it is in the cache's input hash: one sequence with two
  annotations *and* two proteins is two inputs.

**Validation.** `make check` green: **727 fast tests · 68 golden passed / 4 skipped** (PGA ·
chloe+julia · table2asn — MAFFT and BLAST+ were on PATH for this run). New: 31 unit tests
(`tests/unit/test_reference_audit.py` — every block class, the N ≥ 3 and N = 2 majority rules, pin
alive/void/stale/unused, the dot plot's runs, the CLI's exit 2 with the reason in the manifest, the
`--confirm-exception` door) and 4 corpus tests (`tests/golden/test_reference_audit_corpus.py`,
`.panel/`-gated; the MAFFT one is tier B). Measured on the pinned set: **1 record blocks**
(`NC_024936.1` — partial `rps19`, two features = the record §8.34 predicted), **32/34 same SSC
orientation · 2/34 inverted and reported**, 20/34 name sets differ from the majority. The real
negative case is vendored and runs everywhere: `tests/data/corpus/Iris_hermona…pga.gb` — **84 CDS,
83 without `/translation`**, IR gate failed (26.36 %), while the published record of the same
sequence passes every check.

**Residual risks.** (1) The SSC rule cannot see a record inverted *as a whole* — LSC and SSC still
agree inside it; it does not matter for `blastx` (both strands are searched) and the cross-record
dot plot is where it would show, but no test measures that case. (2) The audit re-reads all 34
records on every run (~12 s, IR detection dominating) — no cache behind the gate yet. (3) A pin
naming a check the audit did not raise reads `unused` and does nothing; the status is reported, but
nobody is forced to look at it. (4) The panel currently pins **no** exception: the two flip-flops
are reported on every run until a person confirms them.

**Blocked.** Nothing new. The ★ pick still waits on the end user — the report it is read from now
exists.

**Next.** **P3B-04** — `gene_normalize`, with the measured thresholds of §8.41 (margin ≥ 10 bits ·
identity ≥ 50 % + coverage ≥ 0.7) and a recorded override standing in for the ★ pick.

### 2026-10-08 · Claude · dev WSL · session 24 (Phase 3B opens with its preface: the reference set, the three numbers, and the cache that was lying)

**Did.** Three of the board's ten tasks, each committed and pushed; a pre-existing red test fixed
first, because a commit has to stand on a green `make check`.

* **The red test was a real finding (§9.C-37).** `make check` failed in
  `test_second_identical_run_is_cached_and_identical`: this machine's clock stepped **~1 s
  backwards** between two runs, so a later run was stamped earlier, reused a second, and took the
  `-2` suffix — and `sorted()` over run ids then put the third run *before* the second. The test
  read back the directory it had just emptied. The product is fine (`RunLayout.create` still
  guarantees one run = one directory; a new test pins that for a standing **and** a backwards
  clock); what was wrong is the tests' assumption that a name order is a time order. Every run
  directory is now identified by the set difference before/after a run. **§9's items shift**: the
  new item is `§9.C-37` and §9.J's four move to 38–41 (ranges updated in `CLAUDE.md` and the
  guard's docstring).
* **P3B-01 · the reference set, pinned (§8.40).** `.panel/plastid-reference/` — the 34 Iridaceae
  corpus records, re-pinned, 34/34 past the IR gate, 5,204,491 bp, **every sha256 identical to the
  corpus pin**. The taxonomy rule (D-b) is the family rule, and the family is Iridaceae because
  the contract genome and the answer keys both are. The **leave-one-out design is data**, not
  prose: `core/plastid_reference.py`'s `leave_one_out()` and the committed
  `data/plastid_reference.json` (34 targets × the other 33), because the corpus carries the
  golden's own sequence and PGA's reading of it is 85 CDS where the contract says 86. New doors:
  `plastid reference verify` (re-hash, `file — reason`, exit 1) and `plastid reference loo`.
* **P3B-02 · the three numbers, before any `gene_normalize` code (§8.41).** Measured on the
  corpus with BLAST+ 2.17.0 (`blastp`, `-seg no`, E-value 10) plus a fixed-parameter alignment for
  tRNA. **Two of D5-09's own reasons do not survive**: `rpl2`/`rpl23` are 22.0–44.4 % identical and
  the two `trnL` 61.2 % — neither needs test 4. What does: the **serine tRNAs** (26 of 1,182
  cross-record tRNA comparisons rank another locus first, up to 15.05 points), **two genes under
  one name** (the golden annotates two trnG features as `trnG-UCC`, 73.24 % identical), and names
  the vocabulary cannot resolve (72 tRNA name occurrences + `lhbA`). The margin is a **score**
  margin in bits (min 19.2, median 298.0) — in identity points it can go *negative* (`psaB` vs
  `psaA`). The floor: with **coverage ≥ 0.7** the wrong-locus ceiling is 44.93 % and the true
  ortholog floor 67.75 %, so `identity ≥ 50 %` fits — but only with 0.7, because `psbJ` is 36 aa
  and a 21-aa alignment of it reaches 52.38 % against `matK`. D5-09 §4 and §9.3 are corrected in
  place; the numbers are pinned in `tests/golden/test_naming_numbers.py`. **H-04 measured for the
  first time**: `trnG-GCC` and `trnG-UCC` are two genes and the records use the two names in both
  directions — a curator's call, not a spelling to normalise.
* **P3B-05 · the cache stops lying (§8.42).** The key now carries the **annotation** (not just the
  sequence), the input **order** (not sorted), the module's **source**, and what the step was asked
  to produce (`--docx`) — and any parameter naming a file is **content-addressed**, which is how
  D5-09 §7.2's ★ sha256 reaches a key at all. Each old behaviour was shown broken by loading
  HEAD's `cache.py` beside the new one; four regression tests pin them. Existing cache entries are
  unreachable now — the next run recomputes.

**Validation.** `make check` green after every commit: lint · **691** fast tests · golden **64
passed / 4 skipped** with named reasons (PGA · chloe+julia · table2asn absent). New: 13 tests for
the reference set, 4 tier-B naming-number tests, 5 in `test_run.py`. The old cache behaviours were
demonstrated, not assumed.

**Residual risks.** (1) The three naming numbers are measured on Iridaceae records only — a
cross-family sample (the team's Oldenlandia/Coleus) is **outside what was measured**, and the
identity floor may need re-measuring there; the numbers are recorded with that boundary. (2) The
vocabulary gap §8.41 measured (an alias table that cannot resolve `trnL-CAU`, `trnL-GUA`,
`trmA-UGC`, a bare `trnL`, or `lhbA`) is input to P3B-04 and nothing consumes it yet. (3) `gene_aliases.yaml` still
runs D3-02's direction — flipping it belongs to P3B-04. (4) The cache-key change costs one
recomputation of every cached step, by design.

**Blocked.** Nothing. The ★ pick (end user, from the audit report) blocks only the real-DB
normalize run of P3B-04 — a recorded override is allowed meanwhile.

**Next.** **P3B-03** — the reference-set audit — then **P3B-04** with §8.41's thresholds.

### 2026-10-07 · Claude · dev WSL · session 23 (the coverage check re-run at tool level — and Phase 3B opens)

**Did.** Re-ran ROADMAP's coverage check ("How to use" 5) over "The pipeline" — this time **per tool
named in a box**, which is what session 22's box-level pass had not done. Four gaps found, each
closed with a home; the canonical diagram (§13.0) corrected first and the copies re-synced.

* **Long-read assembly had no milestone.** The diagram names `OATK` · `raven` (and `Flye`) in
  box 1 and reads as "Illumina · **ONT/HiFi**", but Phase 7's milestones covered only the
  short-read path (M7.1) plus Bandage/remote/RAM. → new **M7.6 Long-read path** (BLOCKED: no
  ONT/HiFi reads; both tools declared and run in the survey, §7.1) with rows in ROADMAP's and
  TASKS.md's waiting tables.
* **`Flye` left the canonical diagram** — §7.1 measured it dying on RAM here and replaced it with
  `raven-assembler`/`OATK`; it was also the one tool in the diagram declared nowhere in the repo.
* **`Sibelia` left the canonical diagram** — run for real in the 26/9 survey, never adopted: the
  project's synteny route is `nucmer` (D4-02, M4.2) and no milestone needs LCB blocks.
* **`BWA`/`minimap2`** (the map-back in the read-QC line) → named in **M7.1**'s title.
* **The `.sbt`** → a TASKS.md "Waiting on data / people" row: it is the team's artifact (M3.4 and
  `docs/tools/table2asn.md` both name it as outstanding) and had no row anywhere.
* Rule 5 now reads at **tool level** (a tool with no owner is a gap, or is marked *surveyed — not
  adopted* with the measurement that dropped it); the coverage table rewritten tool by tool;
  PROJECT.md §13.0 gained a dated correction note (canonical record) and the three copies were
  re-diffed byte-identical.
* New guard **`tests/unit/test_pipeline_diagram.py`**: the three diagram copies must be
  byte-identical (the copy rule was prose only before), and every milestone id cited in a ROADMAP
  table must exist as a milestone.
* **Phase 3B opened** (the dev's call): the board is `TASKS.md` **P3B-01 … P3B-10**, ordered by
  execution — the reference set → the three measured numbers → the audit → `gene_normalize` → the
  cache → Ts/Tv → the engine → arbitration → CLI → validation — with a "Decisions to make" table
  (D-a · D-b · D-e, the conflict exit code, the parameter name, the organelle group, the schema
  bump, how Chloe retires), a pitfalls list from §8.15 · §8.21 · §8.22 · §8.31 · §8.34, and
  **H-15 · H-08 · H-14 folded into P3B-05**. Branch `feat/phase-3b`.
* **History rewritten on the dev's call — `main` is now one commit per phase** (9 commits: phases
  1–5 · 4a · the session 22–23 docs · this rule change). The 131-commit linear history is archived
  under the `archive/*` tags, and every archived tree was proven **byte-identical** to the commit
  it was taken from (`git diff` old-main vs new-main: empty) before the force-push; the board branch
  was re-cut on the new main. `CLAUDE.md` → Commits carries the rule now (the branch holds the
  detail and is never deleted; `main` holds the phases and takes `git merge --squash` at close;
  task-level `git bisect` runs on the branch), and this file's header says where the live state is.

**Validation.** `make check` green — lint · **670** fast tests (2 skipped: no mafft) · golden 43
passed / 21 skipped with named reasons. Both new guards negative-checked in memory: a planted
`M1.99` in a table and a one-character drift in a copy are each flagged. D5-09's own pipeline table
(a hand translation, not a copy) synced for the two rows it repeats. The history rewrite verified
lossless: 8/8 archived trees equal to their source commits and `git diff` old-main vs new-main empty
before the force-push.

**Residual risks.** (1) D5-09 Appendix B still lists `Sibelia` among the kind-5 (no-reference)
tools — that appendix is a reference-mechanism list, not the pipeline, so it was left alone.
(2) Six declared tools still carry "# §7.8 records no version" (`minimap2` · `samtools` · `wgsim` ·
`mreps` · `trf` · `novoplasty`) — against CLAUDE.md's pinned-tool invariant; none is on a path that
has produced a number yet, and no task owns it — named here rather than silently fixed.

**Blocked.** Nothing. M7.6 waits on ONT/HiFi reads, now recorded in both waiting tables.

**Next.** Phase 3B is open on `feat/phase-3b` — start at **P3B-01** (the reference set, built and
pinned), then P3B-02 (the three naming numbers, before any `gene_normalize` code).

### 2026-10-07 · Claude · dev WSL · session 22 (doc alignment with D5-09 — and the cache defect the review exposed)

**Did.** Read `PROJECT.md` and `CLAUDE.md` in full against the new D5-09 design record (the
reference-driven naming step) and fixed every contradiction the review found.

* **CLAUDE.md** — the pipeline sentence read `assembly → annotate → gene-name unification →
  analyses`; D5-09 reverses it (names are unified **on the reference DB, before annotation**; the
  step after annotation only **verifies**). Rewritten, plus a D5-09 paragraph in the briefing (the
  naming reference is a sixth kind, chosen by the end user from the audit report; D3-02 is
  reversed; `gene_normalize` is a new module) and a "who decides" row for the ★ choice. Also fixed
  **"Current state: Phase 5 in progress"** — stale since 05/10; the file now says phases 1–5 + 4a
  complete, next board 3B, with its slot (and D5-09's) open.
* **PROJECT.md §3** — D3-02 marked `[ĐÃ THAY 07/10/2026]`; the D4a-05 alias row points at the
  reversal; the **D5-09 decision row** added (it did not exist). §13.3's naming bullet marked as
  moved. §13.0's `[5] RNA editing` said **PREPACT** while the table beneath it said "chưa chọn" and
  §7.4/§7.5 measured PREPACT unreachable (`http=000`) — the diagram now says "chưa chọn"; the D5-09
  copy synced and re-diffed **byte-identical**.
* **D5-09** — §5's status column used "the reference" for two different things; `measured` /
  `not measured` now say **in A**, `outside vocabulary` says **in B**. §9 gained the three gaps the
  review found: the **exit code** for a conflict stop (2 is already the IR gate's — and D5-08's
  audit's), the **organelle group** `gene_normalize` must declare in §13's table, and the fact that
  **§7.2's cache promise cannot hold today**.
* **ROADMAP** — a second open-decision row: the naming step's slot, which interacts with 3B's.
* **The cache defect, measured and named H-15.** `Cache.key` sorts input hashes and `_input_hash`
  hashes the sequence only, so the corpus's own "one sequence, two annotations" pair shares a key in
  either order — while six modules take `genomes[0]` as their reference. Measured on
  `tests/data/NC_050833.1.gb` (86 CDS) vs the corpus PGA copy (85 CDS): equal `_input_hash`, equal
  key with the inputs swapped. The run directory's fingerprint **is** order-sensitive (`cli.py` —
  "ordered: order is meaning"); only the cache is not. → TASKS.md **H-15**.
* **Docs corrected as D5-09 §8 asked**: `gene_name_variants.md` §6 (direction note; the two CGAS
  citations fixed to the measured lines — `return candidates[0]` at **211 and 214**, not 208; the
  prefix match at **149**, not 150–153) · `D5-08` line 54 (pointer: the reference set must be
  name-normalized before the annotator uses it) · `gene_aliases.yaml` header (D3-02's direction
  marked as reversing when the naming step lands — values unchanged until then, so code and tests
  stay consistent).

* **Two more false statements the stale sweep found, corrected in place**: `README.md` still said
  "Phase 5 in progress" and "Nineteen analysis modules" (it is 20 — `stats` arrived in Phase 5);
  `core/gene_names.py`'s docstring still stated D3-02's direction as current.
* **Three more places still stated the old pipeline order** (CLAUDE.md had been fixed first):
  ROADMAP's Phase-3B goal, D5-08's pipeline diagram and its "next steps" bullet, and the D5-08
  decision row in PROJECT.md §3.
* **The four dead items removed**, one commit: `docs/PHASE3_BLOCKERS.md` (its "still blocked — a
  decision on the mito annotator" item was false since 03/10; its one live ask — who `Vu,H.-T.` is
  / whether `PX923121.1` is the manuscript's specimen — moved to `TASKS.md` → "Waiting on data /
  people", and `environment.yml`'s pointer to the file was re-pointed) · `TASKS.md` **H-05** (PMGA
  web re-test — nothing waits on it; the mito path is settled by D3A-02 and frozen by D5-07) ·
  PROJECT.md **§13.10's** 26/9 snapshot table (heading kept — a log entry cites it; every row of it
  now lives in §8 or §7.1) · PROJECT.md **§12's** unpinned install block (the §9.B-14 hazard;
  `environment.yml` is the pinned source, and §7.8's pointer was re-pointed). `docs/tools/pmga.md`
  **stays** — it is the marked record of a measured dead end and keeps a later agent from
  re-investigating PMGA.
* **ROADMAP now covers the whole pipeline — the map, the rule, and the two owners that were
  missing.** Added **"The pipeline"**: a copy of the §13.0 diagram (three copies exist now —
  PROJECT.md is canonical; all three verified byte-identical) plus a **box → milestone** table, and
  "How to use" rule **5**: every box names its owner, and a box with no owner is a gap named in Open
  decisions. The check found **two** ownerless boxes (07/10): the naming step → **M3B.0**, and
  **repeats** (`repeat (vmatch)` — in scope since §1/§7.1, tools declared in `environment.yml`, the
  diagram named it, nobody owned it) → **M2.12** + TASKS **H-16**, with the survey's "9 matches
  ≥30 bp" flagged as measured on other data (the §9.H trap). The same pass **settled the slots**:
  Phase 3B opens next, naming inside it as M3B.0 — recorded in ROADMAP (heading + Open decisions),
  `docs/D5-08` (D-c) and `docs/D5-09` (§9.5). Smaller: M3B.1 now records that the audit's report is
  what the ★ choice reads; Phase 7 notes `fastp`/`BWA` are not declared yet; PROJECT.md §13.2 gained
  the repeats bullet; TASKS gained the PI sign-off row (D5-09 §9.4).

**Validation.** `make check` green (lint · 668 fast tests · golden with named skips). The three
copies of the diagram (§13.0 · D5-09 §1 · ROADMAP) are byte-identical, checked by script. One commit
per step, pushed.

**Residual risks.** (1) The D3-02 reversal is **documented, not built** — the alias table and the
three modules that read the canonical key still run the old direction until `gene_normalize` lands.
(2) `gene_normalize`'s exit code and its organelle group are named open in D5-09 §9, not decided.
(3) H-15 is measured but not fixed.

**Blocked.** Nothing.

**Next.** Open the Phase 3B board (`make-tasks`) from ROADMAP's M3B.0–M3B.5 — the slot is decided,
M3B.0 (naming) first.

### 2026-10-05 · Claude · dev WSL · session 21 (Phase 5 — figures, statistics, publication exports; and the queryable-number check)

**Did.** The whole Phase 5 board, task by task, each committed and pushed on `feat/phase-5`, then
merged to `main`. Order followed the board's own advice: the phase's check first (P5-05), then the
figure system, the rest of the figures, the statistics, the exports.

* **P5-05 · the queryable-number audit.** `report/audit.py` requires every number in a generated
  report to be a bounded number inside the rendered value of a cell of that module's own tables.
  The strict rule was run against the Phase-4 deliverable first and failed with **10 findings over
  its 12 reports, in 6 classes** — every one a quantity colliding with a number printed inside a
  *quoted text* column (a block id `…#10`, `mean + 1 sd …`, `Coleus: 6, Oldenlandia: 5`, a note).
  A `str` value containing a letter is prose: set aside by name, never counted as a cell carrying
  a number. After that rule the deliverable run is clean, and planting an orphan in a real report
  is caught by name. `read_result` came with it — the JSON round-trips into a `Result`.
* **P5-01 · the figure system.** A figure is a pure function of a `Result`: `series()` reads cells
  into plain data, `build()` draws through one theme. The plan said `core/figure.py`; the layer map
  lets `core/` import only `plastid.core`, and two figures read tables from **two different
  modules**, so the system lives in `result/` and renders at **run level from the stored
  `result.json`** — which is why a fully cached run still gets its figures and `plastid run
  --figures` and `plastid figures <run-dir>` are the same function. Determinism was **measured, not
  assumed**: SVG carries a per-process salt *and* a `dc:date` (both pinned), PDF has a metadata
  hook, the PostScript backend has none so its `%%CreationDate` is rewritten.
* **P5-02 · six more figures and the N rule.** `identity_synteny` was rendered at N = 2 · 8 · 34:
  at 8 the series read, at 34 the legend covered the panel and 33 lines merged into one mass. The
  rule is a **cap of 8**, above which the panel draws the span the genomes cover — both edges still
  cell values — and the canvas grows with N.
* **P5-03 · statistics.** `core/stats.py` (scipy, `[stats]`) with every function pinned by a vector
  computed **by hand**; writing those vectors caught two of my own errors. The `stats` module runs
  a real comparison over the stand-in panel: 80 loci × 4 measures = 320 tests, 17 not computable,
  **0 significant after Benjamini–Hochberg** — the honest answer at 4 vs 4 observations per locus.
* **P5-04 · the exports.** Byte-identical multi-sheet Excel, `plastid package` (the Supplementary
  package, regenerated from a run directory), the technical QC report, and a PDF-only one-page
  sheet per genome: **34/34** corpus records, none failing, each one page, measured at 114 s for
  the whole corpus.

**Validation.** `make check` green after every task. Final: **668 fast tests**, golden **60 passed,
4 skipped** (PGA · chloe · julia · table2asn absent on this machine) with the tool environment on
PATH. Every figure was rendered from a real run and **looked at**.

**Residual risks.** (1) The figure and export extras are not installed by a bare `make env` —
`uv sync --all-extras` is what this tree uses. (2) The one-page PDF is A4 portrait with a fixed
layout; a genome with very many regions would crowd it, and no such record is in the corpus.
(3) The N cap of 8 is measured on *this* figure at N = 34; another figure may support more series,
and the rule is stated per figure rather than globally.

**Blocked.** Nothing.

**Next.** Phase 3B (the self-built plastid annotator, D5-08) once its slot is confirmed — or the
open housekeeping items, of which **H-14** was found and measured this session.

### 2026-10-05 · Claude · dev WSL · session 20 (doc ↔ code sweep: false statements fixed in place, §9 renumbered, the log split)

**Did.** A consistency review of `PROJECT.md` · `CLAUDE.md` · `ROADMAP.md` · `TASKS.md` · `docs/`
against each other and against the code, then the fixes it found. One rule clarified in CLAUDE.md:
**a statement that has become false is corrected in place** — history stays in §8/§9 (numbers
untouched), but a stale sentence in a spec is a wrong instruction to the next agent.

* **Superseded numbers leave the live text.** §13's Ts/Tv bullet and §13.4's in-phase check now
  carry **2.6917 · 79 genes · 7 blocks · 97.37** (was 2,686 · 76 · 21 blocks); §10.1 and §12
  follow; §7.8's MAFFT row states what was measured without quoting the old golden number; §13.10's
  two stale rows carry inline provenance like §7.1/§7.4 already did.
* **"Applies to both organelles" is true as written again** — §13's chapeau, §13.2/§13.4 and
  ROADMAP's cross-phase rule + Phase 2/4 headings now say the seven transitional analyses are
  plastid-only today (H-10, on hold). §13.3's "(132) recount" reads like ROADMAP's Phase 3 check:
  the *rule* is fixed, the *number* is the golden genome's.
* **§9 renumbered to unique, continuous numbers 1–40** (A 1–8 · B 9–26 · C 27–36 · J 37–40). The
  duplicates were real and had already bitten: `20` appeared three times, and `§9.C-20` was written
  in one place for the NCBI-filter item and in another for the style-pruner item. Every citation
  repo-wide updated (~25 sites incl. code docstrings); `test_project_docs.py` gained a guard that
  fails on a repeat or a gap.
* **Decisions that had no row:** **D2-03** (compiled table 1 for codon/AA; §3) and **D5-05** (§10.11
  as an open PI question — TASKS claimed it was "recorded in §3" and it was not).
* **New §9.J — four measurements** made while auditing the two-organelle rule in code:
  `gene_names` never loads the mito vocabulary (0 renamed / 0 aliased vs `mito_genes` naming 31 of
  51); the translation-table column implied the record's own table (5 of the golden's 86 CDS start
  GTG/ACG/ATC); our mito artifact has **0** source features and **0** `/organelle`; and
  `submission` on that artifact reports **24** first-pass problems (23 = no `/product`) while still
  writing the files. **Only the table label was fixed** (code, no number moved); the rest are named
  — **H-12** (new), plus H-09 and H-11 — and **on hold with H-10/H-06 by D5-07**.
* **Housekeeping:** `docs/handoff-archive.md` holds sessions 1–13 (a pure move, re-concatenation
  verified byte-identical); `pmga.md` gets the 05/10 update note its siblings already had; README's
  "Phase 1" status; `data/…` → `src/plastid/data/…`; the "Nhật ký & tài liệu" pointer → "The run
  log"; ROADMAP's duplicate M1.6/M1.8 lines, three decided rows and M4a.3 (ticked with its named
  exceptions); PHASE3_BLOCKERS' stale "can start immediately" list.

**Validation.** `ruff check` + `ruff format --check` clean · **567 fast** tests · golden with the
tool env on PATH: **57 passed, 4 skipped** (chloe · julia · table2asn absent on this machine).

**Residual risks.** (1) The §9 renumbering touched citations in code, tests and `docs/` — all
updated, and the new guard fails on a duplicate or a gap. (2) §9.J.39 *states* the initiator
convention (V/T/I where a published protein begins with M) rather than changing it; changing it
moves 5 residues of the golden composition — a decision, not a fix. (3) The HANDOFF archive is a
move; a citation by **content** rather than session number should name the archive file.

**Blocked.** Nothing.

**Next.** Phase 5 — **P5-05** (queryable-number check over generated reports) or **P5-01** (figure
system + first figure).

### 2026-10-04 · Claude · dev WSL · session 19 (Phase 4a: the corpus pinned, the fixtures vendored, and two runs that died on real data)

**Did.** P4a-01 and P4a-02 DONE and pushed; P4a-03 (scale) started, and it has already paid for
itself.

* **P4a-01 · the corpus in the tree and pinned.** `.data/thanhdiem_data/` (gitignored, 136 files,
  28 MB, every sha256 matching `D:\`), both annotation sets pinned as panels, the FASTAs as
  `fasta.json`. Recorded in §8.31 with the measurements that came out of pinning.
* **P4a-02 · fixtures + Phase-1 checks.** Five files · 1.3 MB in `tests/data/corpus/` and twelve
  tier-A tests in `tests/golden/test_corpus.py`: the golden re-annotated is one CDS short, *Iris
  hermona* passes the gate as published and fails it at 26.36 % under PGA **on the same sequence**
  (exit 0 / exit 2 through the CLI), the group's own record 87 → 85, the `ndhI`/`ndhL` alias
  resolving both spellings to one key. Plus a new guard: **every file under `tests/data` must be in
  `SHA256SUMS` and match it** — the corpus files were not, and the test found that before commit.
* **Two corrections to my own board, both worth the ink.** (1) D4a-05 said "add the `ndhI`↔`ndhL`
  alias" — it has been in `data/gene_aliases.yaml` since 03/10 (§8.17); today only re-measured it at
  8/34 pairs. (2) The corpus's tool is identified by the **files themselves**
  (`/note="Annotation Method :: PGA-Plastid Genome Annotator"` in 34/34 of set 02, 0/34 of set 01),
  not by the naming style I had inferred it from.

**P4a-03 is where the corpus earned its place — two records killed two whole runs, and neither
failure was a wrong number.**

1. **`structure`, *Chasmanthe floribunda* NC_072562.1.** LSC 81,003 bp ⇒ the 1,000 bp windowing
   ends in a **3 bp** window (`TTA`). No G, no C, `gc_skew` undefined, module raised ⇒ **run dead**,
   1 of 34. The golden genome survives only because its LSC ends in a 127 bp window that happens to
   hold G and C. Fixed: an undefined window is blank (cumulative blank from there on), counted in
   `n_undefined_windows` and named in a note; a *region* without G and C still raises.
2. **`orf`, *Libertia pulchella* NC_072575.1.** The record annotates `rps16` **twice**
   (`5741-5780,4684-4895` and `5737-5780,4684-4891`); both collapse to span `4684-5780`, so both
   took the id `rps16:4684-5780`, and `Table` refuses a repeated key ⇒ **run dead**, 1 of 34. Fixed:
   ids unique per feature (`#2`), both rows checked, `n_duplicated_cds` counted and named — a
   duplicated annotation is a finding about the record, not rubbish to drop.

Both modules → **0.2.0** (H-08: the cache keys on `module_version`). Regression tests for both are
offline and synthetic; PROJECT.md **§9.I** carries the lesson — *neither fixture could have come
from data this project picked for itself, because it picks records that work.*

**P4a-04 · Chloe over the 34 bare sequences — DONE (§8.33).** Running the annotator on a bare FASTA
needed a **FASTA door in `plastid run`** first (`-i` read GenBank only), and opening it exposed a
third bug of the same family: `orf` could not analyse a record with no CDS, though its own comment
says that is a legitimate input — the per-CDS table was built *before* the branch that guards it.

34/34 annotated, none failed. Totals: published 2,924 CDS · PGA 2,887 · **Chloe 2,922**, with
**0 CDS failing the ORF checks** and 201 alternative starts. **No record has all three agreeing on
the gene set (0/34)**, and the differences are not scattered — they cluster on the two hardest loci
of a plastome: **`ycf1`** (a CDS in 24/34 published records that neither other annotator calls — the
IR/SSC boundary gene §8.11 already flagged) and the **trans-spliced `rps12`** (Chloe calls 2 copies
where 22/34 published records carry 1 — the number of IR copies is an annotator's decision, not the
sequence's). The rest is vocabulary (4/34 differ in tRNA spelling only). And **132 = 86 + 38 + 8 is
not a constant**: published and Chloe both range 84–87 CDS here — the P2-10 trap.

**P4a-03 · scale (M8.3) — DONE, `docs/M8.3_scale.md` + §8.32.** Machine stated with every number
(i7-1165G7, 6 vCPU, 8,135,488 kB — the ordinary machine the toolkit targets).

* Nine single-genome modules: **13–26 s for all 34 genomes each**, peak **87 MB**.
* The whole-genome MAFFT (`substitutions`' default scope): 4 genomes 72 s · 8 → 2 min 44 s ·
  **16 → over 7 min** · **34 → over 11 min**. Peak stays under 100 MB throughout.
* The **CDS scope** — the scope the manuscripts quote — does all 34 in **5 min 24 s** at 179 MB, and
  finishes. On a machine with 8 GB, memory is the cheaper currency.
* Per-pair modules: **`synteny` 19.6 s** (99 tool calls) · **`popgen` 6 min 13 s** (80) · **`kaks`
  8 min 26 s** — 5,202 MAFFT + 4,766 `yn00` = **9,968 tool calls**.
* **The finding: call count is not what costs.** 9,968 two-sequence alignments finish in 8.5 min; ONE
  34-plastome alignment does not finish in 11. A module's price cannot be read off its loop — the
  estimate written into the file first said "about 2,600 of each", and the 33 `nucmer` calls it
  expected to be the slow half turned out to be the fastest thing measured.
* `diversity` and `barcodes` at 34 are **named as not measured**, never estimated.

**Late in the session: the run log moved out of `PROJECT.md`** — §8 was **1,531 of 2,379 lines
(118 of 216 KB)**, so reading one entry meant opening a document that no longer even fits a default
Read. It is now five files under `docs/log/` (`phase-1` · `phase-2-3` · `phase-3a` · `phase-4` ·
`phase-4a`), with a 33-line index left in §8. **Entry numbers did not move**, because **271
references** point at them by number (27 in code docstrings, 36 in tests, 23 in `docs/`, 115 across
the four session files, 70 inside PROJECT.md) — renumbering would have broken citations a docstring
cannot re-derive. The split was proven a **pure move**: the five files concatenated are
byte-identical to the §8 that was in git, and §1–§7 / §9–§13 are untouched. New guard:
`tests/unit/test_project_docs.py` fails the build if a cited `§8.N` has no heading, if the index
misses an entry, if one number appears in two files, or if a cited `§9.x` item is gone — with a
negative test, because a guard nobody has watched fail is not a guard. Rule written where it is
read: `CLAUDE.md` → "The run log".

**Validation.** `make check` green after each commit — **564 fast + 42 golden** at the last run
(tier-B skips name their missing tool; this shell has no tool env unless `PATH` carries
`/root/toolbox/env/pts/bin`).

**Residual risks.** (1) The corpus is gitignored and re-fetchable only from anh Điềm's disk; the
durable record is §8.31 plus `PINNED.json`. (2) The two panels are keyed differently by
construction — anything joining them must join on the **file stem**, never the accession. (3) The
scale numbers are this machine's (WSL, 8 GB cap): they bound a *shape*, not a promise for another
machine.

**Blocked.** Nothing.

**Next.** Finish P4a-03 (§8.32 + board), then **P4a-04** — Chloe over the 34 bare sequences, which
needs a FASTA door in `plastid run` (`-i` reads GenBank only, measured) or the de-annotated-GenBank
route P3A-07 used.

### 2026-10-04 · Claude · dev WSL · session 18 (Phase 4a opened — the Điềm corpus becomes the test dataset)

**Did.** Opened **Phase 4a** (board + ROADMAP milestones M4a.0–M4a.4) and closed **P4a-01**.

* **The corpus, in the tree.** 34 Iridaceae plastomes × 2 annotations of the same sequences
  (published NCBI · PGA) + 34 bare FASTAs = 68 `.gb` + 68 `.fasta`, 136 files, 28 MB. Copied from
  `D:\business\pts\thanhdiem_data` into **`.data/thanhdiem_data/`** (working tree, **`.data/`
  gitignored**), 136/136 sha256 verified. Measured why it belongs in WSL: 300 small files
  7 ms vs 700 ms on `/mnt/d` (**~100×**), while one cold 28 MB bulk read is only 0.38 s vs 0.55 s
  (**1.4×**) — a run is thousands of small writes, so that is the number that matters.
* **Both sets pinned, verified.** `.panel/iridaceae-34/` and `.panel/iridaceae-34-pga/`, 34 records
  each, every sha256 · md5 · byte count and measured value re-read (`PINNED.json` is the table).
  Set 01: **0** gate failures, 2,924 CDS. Set 02: **3** — *Chasmanthe* 29.24 % · *Iris hermona*
  26.36 % · *Trimezia* 27.78 %, each named — 2,887 CDS.
* **The FASTA half needed its own route.** Measured: `panel pin --file` on a `.fasta` pins it under
  the file stem with `parse_error: no GenBank records found`, prints UNVERIFIED and **exits 1** — the
  pin path is GenBank-only. Pinned instead as `fasta.json` (68 entries: sha256 · bytes · length ·
  the paired record's sha256), **and every one of the 68 sequences was proven identical to its
  pinned record**, so the manifest records bytes while the sequence's provenance is the panel.
* **Pinning the two sets into one directory would have been wrong** — `panel pin --file` keys a
  record by its own id, and PGA keys one `Aristea_ecklonii_voucher_KK2022_1_NC_072561` where NCBI
  keys it `NC_072561.1`: one genome, two names, nothing pairing them. A directory per set; the
  pairing is the **shared file stem**.

**What the corpus measured (all new, §8.31).**

* **The golden genome is in it**: `01_…/Iris_domestica_NC_050833.1.gb` is byte-identical to the
  vendored fixture (`a946f23e…`) — the corpus contains the contract it is testing.
* **Set 02 = PGA on the same sequence**: the golden reads **85 CDS where the contract says 86**
  (the second `rps12` copy is gone), and the group's own *E. bulbosa* goes 133 → 131 the same way.
* **3/34 PGA records fail the IR gate** at 26–29 % — and de-novo detection on the **very same
  sequence** finds **100.00 %** every time. On *Chasmanthe*, de-novo returns PGA's own 26,440 bp
  span at 100 % ⇒ the span is right and the pairing is not. 8.8 % of the corpus would be
  silently wrong if anyone trusted the annotation's IR.
* **34/34 records spell at least one locus differently** — and the split matters: **vocabulary**
  (`clpP1`/`clpP` 23 · `pbf1`/`psbN` 21 · `ndhL`/`ndhI` 8 · the §8.27 trio 7 · `lhbA`/`psbZ` 1 ·
  bare tRNA 8 · the rRNA dialect 1) versus **a different call** (`ycf1` has no CDS in PGA's
  annotation in **24/34**; `rps12` loses a copy in 12/34). §8.27's lesson at 34× scale.
* **`ndhL` (published) and `ndhI` (PGA) are the same locus** — same coordinates
  (`118156-118698(-)`), identical 543 bp sequence, on the golden genome itself. D4a-05; the golden
  pair is unaffected (both NCBI records say `ndhL`), so **no golden number moves**.
* **A gap that belongs to P4a-04**: `plastid run -i` reads **GenBank only** (`read_fasta` exists
  solely on the `mito annotate` path) — running Chloe over 34 bare sequences needs a FASTA door, or
  the de-annotated-GenBank route P3A-07 already used (the two carry the same sequence, 68/68
  verified).

**Validation.** `make check` green (**533 fast + 30 golden**; tier-B skips name their missing tool —
this shell has no tool env). Both panels re-read clean; `fasta.json` re-reads clean; the corpus
measurements above are all scripts run in this session, not readings.

**Residual risks.** (1) The corpus is **gitignored and re-fetchable only from anh Điềm's disk** —
the durable record is §8.31 plus `PINNED.json`; a re-fetch must re-produce the same sha256s.
(2) The two panels are keyed differently by construction, so anything joining them must join on the
**file stem**, never on the accession. (3) PGA is present only as pinned output, never as a runnable
tool here — unchanged from Phase 3.

**Blocked.** Nothing in Phase 4a.

**Next.** P4a-02 — vendor the three pairs into `tests/data/corpus/` with `SHA256SUMS` + PROVENANCE,
the tier-A gate checks, and the measured `ndhI`↔`ndhL` alias.

### 2026-10-04 · Claude · dev WSL · session 17 (Phase 4 complete — five tasks, and the numbers they moved)

**Did.** The whole Phase 4 board, task by task, each committed and pushed on `feat/phase-4`:

* **P4-01 · synteny.** `external/mummer.py` (nucmer --maxmatch -l 20 -c 65 -t 1 -> delta-filter -1
  -> show-coords -rcl, MUMmer 4.0.1 pinned), `core/delta.py` (the delta format read as measured:
  one header per pair, `|d|-1` paired columns then one side advances), `modules/synteny.py` (block
  table, per-genome summary, 600/200 % identity series, raw coords/delta as artifacts).
  **The board's "21 blocks, %IDY 83.7–98.2 (mean 91.5)" does not reproduce** from the probe's own
  inputs under any parameter set (§8.26 sweep); the check now holds the measured **7 blocks,
  89.66–99.11, mean 97.37, coverage 99.09 %** and the finding is §9.H.
* **P4-02 · Ka/Ks.** PAML `yn00` through `run_tool()` on a per-gene MAFFT codon alignment, with our
  own frame check (yn00 reads codons and never checks, measured) — and the "list the genes one
  genome does not carry" clause found **three loci spelled differently** (atpL/psaL/bspB). Aliases
  added; the shared set of the golden pair went 76 -> 79, so **§6's golden Ts/Tv was re-measured:
  2.6917 = 323 Ts / 120 Tv over 79 genes, 68,448 columns** (CLAUDE.md updated too).
* **P4-03 · population genetics**, pure Python (D4-03): Hd · θ(S) · π · Tajima's D · Fu & Li D*/F*
  · Nei's Gst, each naming its formula. Tajima's D and π match **scikit-allel 1.3.13 to six
  decimals**; Fu & Li are pinned by hand-computed vectors (EggLib's sdist build returns scrambled
  site views here — recorded); Hudson's Fst was implemented first and **rejected with a measured
  reason** (it returns −1 for identical composition at n=2).
* **P4-04 · barcodes**: π reused from the same core function the diversity module calls (a tier-B
  check compares **773/773 windows**), regions by an exported rule, candidate loci from
  `data/barcode_loci.yaml`. `trnH-psbA` is found in the panel's Oldenlandia records and **not** in
  the golden pair — rps19 sits between the two genes there, and the row says so. Two silent naming
  failures fell out (`trnk-UUU`, `trnl-GUA`) — fixed by a syntax rule and one alias.
* **P4-05 · the deliverable**: three runs over the pinned panel (`standin-oldenlandia`,
  `standin-coleus`, `standin-mito`) with `--label` putting the stand-in provenance into every
  table's scope line and every report header; `docs/P4-05_comparison_tables.md` records how to
  reproduce them and what they hold.

**Three more real defects, all found by running the deliverable rather than by a test:**

1. **The CDS scope emptied the whole table when a record had no annotation** (3 of the panel's 11
   records) — it now leaves them out and **names** them; fewer than one annotated genome is a
   raised error, not an empty CSV.
2. **No shared gene between two mitogenomes** — `OX459128.1` carries no `/gene`, so every CDS
   "name" is a `/product` sentence, and the comparison modules were passing the *plastid* alias
   table explicitly, so the mito vocabulary (`data/mito_gene_names.yaml`, present since P3-04) was
   never asked. `core.gene_names.aliases_for` now picks the table by the record's own
   `/organelle`: 0 shared genes -> **19**.
3. **Every tool call overwrote the previous call's logs.** The log index came from
   `glob("*__{name}__*.stdout.log")`, which does not match the names written
   (`001__mafft__stdout.log`) — 616 MAFFT calls kept **one** log, so "the raw output travels with
   the run" was false for every tool but the last, since Phase 1. Fixed, regression test added
   (three calls -> three pairs, first content intact); the mito run went from 1 to **74 MAFFT ·
   24 yn00 logs**.

Also recorded: the run cache keys on `module_version` but not on the module's source, so a
behaviour change without a version bump served a **stale payload** (the CDS-scope fix came back
`cached`); four modules were bumped to 0.2.0 and the deeper fix is **H-08**.

**Validation.** `make check` green — **488 fast + 30 golden**; the full tier-B set with the tool
environment (MAFFT 7.526 · MUMmer 4.0.1 · PAML 4.10.10 · ARAGORN 1.2.41 · BLAST+ 2.17.0 ·
table2asn 1.29.324 · Chloe 0.1.18): **48 passed, 1 skipped** (PGA, not installed). The three
deliverable runs are in `.runs/phase4/` with their table sets, tool logs and manifests.

**Residual risks.** (1) The phase's checks are new golden numbers measured today (§8.26–§8.30) —
as reproducible as the pinned tools and the pinned panel (gitignored; re-fetch per §8.6).
(2) Ka/Ks on **cross-genus** pairs hits frameshifts in ~20/80 genes (MAFFT is not a codon
aligner); every such gene keeps a row with that status. (3) The mito half is a cross-genus sketch
(D4-05): one sequence per group, so every popgen statistic is blank with its reason.

**Blocked.** M4.5 (RNA editing) — needs one RNA-seq sample, unchanged.

**Next.** Merge `feat/phase-4` to `main`, then Phase 5 (figures & reports — the plots the two
manuscripts need) or the team's own data, whichever comes first. Open for the PI (wording only):
`trnG-GCC` vs `trnG-UCC` (H-04) and which gene name the manuscript prints (D3-06).

### 2026-10-03 · Claude · dev WSL · session 16 (P2-07/08 run for real — and the bug that was hiding there)

**Did.** Merged the three-deep branch stack into `main` (fast-forward, pushed). Then closed the
last Phase-2 leftovers, which had been waiting on MAFFT since Phase 2: **P2-08 (Ts/Tv)** and
**P2-07 (π)** now run on real data.

**The run found a real bug, and it is this project's favourite kind.** MAFFT 7.526 prints its
alignment in **lowercase**; the counters classify a column by `_DNA = {"A","C","G","T"}`, so every
column was silently "skipped": the first real Ts/Tv run reported **zero comparable columns and a
blank ratio** — not a wrong number, no number — while every table still rendered. It survived
because the unit tests build `Alignment` objects by hand (uppercase) and the **MAFFT adapter had no
test at all**. Fixed at both ends: the adapter upper-cases on read (with the measurement in a
comment), and `Alignment` now **refuses** a lowercase row — fail loud instead of counting zeros.

**Validation.** `make check` green: **443 fast + 41 golden** (with the tool env; the 4 remaining
skips are tools this shell lacks). New tier-B checks: the CDS-scope Ts/Tv and the whole-alignment
one, plus π.

**The golden number reproduces exactly, through the toolkit's own path** — §6 is 2.686 = 317 Ts /
118 Tv over 76 shared CDS and 66,066 columns, and that is what the run printed, to the base. The
whole-alignment default is a different number (**1.8416**) over 152,338 columns, and every table
labels its scope. π: **0.007181** whole; **LSC 0.009185 · IRa 0.001682 · IRb 0.001681 · SSC
0.014066**; per gene and a 600/200 sliding window.

**Residual risks.** (1) **π per intron/IGS is still missing** (M2.9's last piece) — named as
**H-07**: intron spans are derived inside `modules/introns.py` and IGS spans nowhere, so the
shared layer needs that derivation extracted first (DRY's second real copy). (2) The
whole-alignment Ts/Tv aligns two circular genomes without normalizing their origins first; the
number is honest for the alignment it describes and its label says so, but a per-region view of a
wider scope would be the more useful shape. (3) M2.10 (manuscript draft) still waits on the team's
real samples — unchanged.

**Blocked.** Nothing.

**Next.** **Phase 4 is open** — branch `feat/phase-4`, board built (`TASKS.md`, five tasks with
measured Done-whens). Start with **P4-01 (identity + synteny)**: install and declare `mummer4`, then
reproduce the in-phase check — **21 blocks, %IDY 83.7–98.2 (mean 91.5)** — whose pipeline §7 records
(`nucmer --maxmatch → delta-filter → show-coords`) and whose inputs are still on disk at
`/mnt/d/business/pts/probe/run_local/{target,compare}.fa`. M4.5 (RNA editing) stays blocked on RNA
data; MAFFT 7.526 is installed and pinned, so the alignment path under every one of these tasks is
already proven.

### 2026-10-03 · Claude · dev WSL · session 15 (Phase 3A closed — the CLI, the golden check, the first real sample)

**Did.** P3A-05 (the `plastid mito annotate` CLI), P3A-06 (the itemized golden check) and P3A-07
(the unannotated sample, pinned). `core/genbank.py` gained `read_fasta` — an assembled sequence is
the input door this phase needed; `cli.py` grew a `mito annotate` command that shares
`_execute_run` with `plastid run`, so the mito path inherits the manifest, the cache and the
provenance instead of re-implementing them (seam 3). The boards and ROADMAP now show Phase 3A
complete.

**Validation.** `make check` green (441 fast + 30 golden, tier-B skips named); with `PTS_BLAST_DIR`
all **8** golden checks pass. Measured: the de-annotated `OX459128.1` (102 features stripped) →
`mito_annotation.gb` → re-read **51 = 31 CDS + 17 tRNA + 3 rRNA** → `orf` + `mito_genes` both run
unchanged. The consolidated count check prints `CDS 31/31 · tRNA 17/17 · rRNA 3/3` beside the
published numbers and lists every difference.

**P3A-07 is the honest one.** `AP045469.1` (*Coleus hadiensis*, 313,585 bp, 0 features) →
**45 genes = 23 CDS + 19 tRNA + 3 rRNA**, pinned `3e7c7114ae5d…`, read back offline from the pin.
**8 vocabulary loci are not called**, and a permissive re-search (E-value 10 against 1e-5) shows
four of them — `mttB` · `rps1` · `rps13` · `rps4` — have **no alignment at all**, so the floors are
not what excluded them. Whether those genes are missing from the genome or from this assembly is
not a question a reference-transfer annotator can answer; that is recorded rather than papered
over by loosening a threshold.

**Residual risks.** (1) The Coleus run is a prototype annotation with 13 of 23 CDS failing the ORF
checks (9 `no_start` — alignments that stop short of the reference's N-terminus); it is **not** a
submission and every table says so. (2) The three tRNA differences on the answer key (one missing,
one extra, one anticodon disagreement) are still awaiting a curator's ruling. (3) The floors are
calibrated on one record; the Coleus run is the first evidence of how they behave on a genome with
no answer key, and it says: they behave fine, the genes really are not findable.

**Blocked.** Nothing.

**Next.** Phase 4 — M4.1 Ka/Ks · M4.2 %identity + synteny · M4.3 barcodes · M4.4 population
genetics · M4.6 Ts/Tv (**scope fixed: CDS only · 76 shared genes · MAFFT · IR excluded**). MAFFT
7.526 is installed and pinned, so the Phase-2 numbers that were waiting on it (P2-07/P2-08) can
also be measured now.

### 2026-10-03 · Claude · dev WSL · session 14 (P3A-04 — the tRNA caller, and 17 is not one set)

**Did.** Installed and declared **ARAGORN 1.2.41** (`environment.yml` + `docs/tools/aragorn.md`,
same session as the install), wrote `external/aragorn.py`, wired it into `modules/mito_annotate.py`
(a `trna_calls` table, tRNA columns in the summary, tRNA genes in the artifact), and added the
golden comparison. 12 new unit tests; the parameter that points at the binaries is now one name
(`--set tools=<dir>`) for both tools.

**Validation.** `make check` green (438 fast + 30 golden, 9 tier-B skips with named reasons); with
`PTS_BLAST_DIR` the five golden tests pass, including the new one. Measured on `OX459128.1`:
**ARAGORN calls 17 tRNA, the record carries 17, and they are not the same 17** — `trnA-UGC
13244-13348` is in the record and not called (105 bp, longer than a tRNA), `trnL-AAG 36355-36427`
is called and not in the record, and the tRNA-His at 11741-11819 is **`trnH-AUG` to ARAGORN and
`trnH-GUG` to the record** — an anticodon disagreement, printed rather than resolved. 9 spans
differ by 1–6 bp; the largest is the one tRNA that crosses the origin (`[258274,73]`, written as
the two parts a `join()` needs).

**Two traps that changed the code.** (1) `aragorn -V` is not a version flag: it exits **0** and
prints "No sequence file specified", which `run_tool` was recording as the version of the tool
that calls tRNA genes — the probe now requires a digit in the answer, and the adapter verifies the
release from `-h` itself. Same family as the BLAST+ `-version` bug from the previous session: the
`--version` convention does not exist. (2) A name is not an identity: the record's 17 tRNAs carry
only 14 distinct names, so any comparison that deduplicates by name reports 14.

**Residual risks.** (1) The three tRNA differences above are **findings, not errors** — nobody has
ruled on whether ARAGORN's `trnH-AUG` or the record's `trnH-GUG` is right for that locus; it is
recorded in PROJECT.md §8.23/§9.E and flagged for review. (2) `-i` (intron search) is off, so an
intron-containing tRNA would be reported missing rather than modelled. (3) ARAGORN's
`trnL-AAG` extra call may be a genuine plastid-derived tRNA or a false positive — not yet
adjudicated.

**Blocked.** Nothing.

**Next.** P3A-05 — `plastid mito annotate --input <fasta> --reference <dir> --out <dir>`, with the
output going through `orf` + `mito_genes` unchanged; then P3A-06's itemized golden check
(31 · 17 · 3, differences listed) and P3A-07 (`AP045469.1`, then pin it).


### Older sessions — moved 05/10/2026

Sessions 1–13 (02/10 – 03/10/2026) are in [`docs/handoff-archive.md`](docs/handoff-archive.md) —
a pure move of the text; no entry was retyped or renumbered. This file keeps the resume point and
the sessions from Phase 3A onward.
