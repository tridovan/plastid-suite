# The Crataegus answer-key corpus — what it is, where it came from, how to prove a copy is the same

**P8B-09.** The corpus behind Phase 8B is *third-party published data*: it is never reported as
the user's own data, and every output made from it carries that label (D-k). This note is the record that lets
a fresh clone (or a reviewer) re-fetch the same bytes and check them — the `docs/tools/*.md`
pattern applied to a dataset instead of a tool (D-m).

## The paper

| | |
|---|---|
| Title | *Comparative genomic and phylogenetic analyses of Crataegus chloroplast genomes: insights for evolution and identification* |
| Authors | Xinyu Sun, Mingqi Cui, Baipeng Zhao, Yu Wang, Xiao Zhang, Yuexue Liu (Shenyang Agricultural University) |
| Journal | *Frontiers in Plant Science* **17** (2026), article 1767012 |
| DOI | `10.3389/fpls.2026.1767012` — Open Access, **CC-BY** |
| Published | 11 Feb 2026 |
| Delivered to us | 10/10/2026, by anh Điềm, as `D:\business\pts\thanhdiem_data\db_test_new` |

The 22 MB PDF and the publisher's supplements are in the drop but **not in git** — the same call as
for tools: the note is durable, the 22 MB of third-party PDF is not ours to relicense.

## What is in the drop

117 files, 39 MB, in three parts:

| Part | What |
|---|---|
| `paper/` | the article PDF, 19 pages, with a text layer |
| `supplement/` | `Table_S1…S7.xlsx` (sampling, barcodes, published plastomes, biogeography models, RSCU, SSR, Ka/Ks) + `Data_Sheet_1…5.pdf` (Figures S1–S5) |
| `gb_data/` | `gb_files/*.gb` — **100 GenBank records** — plus `manifest.tsv`, `groups.json` and `fetch_gb.py` |

The 100 records in four groups (`groups.json`):

| Group | n | Accessions | Role in the paper |
|---|---|---|---|
| sequenced for this study | **18** | `PX413282–PX413299` | Table 1, Table S1, and the columns of S5/S7 |
| already-published *Crataegus* | **26** | see `groups.json` (`ON032469` … `OP964008`) | Table S3 |
| barcode amplicons (`ndhC~trnV-UAC`) | **50** | `PX676094–PX676143` | Table S2, Figure S4 |
| outgroups named in the text | **6** | `MK920292` · `MK920297` · `MN068255` · `MN068257` · `MN068262` (`Amelanchier`) · `NC065485` | Figure S5 |

**The 51st barcode sample does not exist online.** Table S2's `seq45` (*Crataegus scabrifolia*,
isolate YNSZ3) carries no GenBank ID, and the data-availability statement covers
`PX676094–PX676143` — 50 records. That row is why the barcode half of P8B-06 compares 50 of 51, and
50 is the whole of what is available.

## How a copy is proven the same

The working copy is `.data/thanhdiem_data/db_test_new/` (gitignored — `.data/` is the working-tree
convention for team drops). Verified 10/10/2026:

| | |
|---|---|
| Files | **117**, each `sha256` compared against the `D:` original — **all identical** |
| Manifest of the manifests | `sha256(gb_data/manifest.tsv)` = `0a34e67936b60f86503ad5faeadcd4ba16c9e6335e48001d3e786e4a059e9ca4` |
| The whole tree in one line | `sha256` of the sorted `"<sha256>  <relative path>"` lines for all 117 files = `4ee47625c1529c5c83922b3566e802c532462a653ba208de44141c2ba519deb` |
| Re-fetch | `python3 gb_data/fetch_gb.py` (the drop's own script, hits NCBI). **No test runs it** — this project's tests never touch the network; it is provenance, not a step |

To check a copy: recompute the third row. A different result means a different corpus, whatever the
file names say; a difference in one file is visible in the second-to-last row of the list.

## What this corpus is not

- **Not anyone's own data.** It is a published Rosaceae data set used as an answer key. Every
  table made from it carries a "third-party published data" label in its scope.
- **Not a benchmark to be tuned against.** D-o: run ours first with every parameter printed, use
  their tables as an answer key to *infer* an unstated design, and never tune a parameter until a
  table matches — a setting chosen to make a number agree is fitting, not reproduction.
- **Not internally consistent, and that is a finding, not a defect of the copy.** Measured on it:
  their Methods' SSR thresholds disagree with their own Table S6 · the text's IR minimum
  (26,311) appears in no genome they deposited, nor in their own Table 1 · one GC cell does
  not reproduce · Table S7's Ka/Ks is computed against *Amelanchier cusickii* `MN068257`,
  not the *C. kansuensis* its text names · **12 of the 18 deposited records lack the whole
  IR protein-coding set**, while the paper's tables carry values for those genes.

## What the corpus found in *our* code

The point of an answer key is that it was not chosen by us. Running the corpus found, so far:

| Defect | Where | Status |
|---|---|---|
| MISA marks a compound cluster two ways; `_parse_gff` read one — 0 of 16 components flagged on `PX413282` | `modules/ssr.py` | **fixed test-first**, module `ssr` 0.3.0 |
| `ssr`'s `Context` column read a trans-spliced gene's outer span (11 of the golden's 37 rows) | `modules/ssr.py` | **fixed test-first** in Phase 8 |
| `trnV-UAC` is never called — ARAGORN reads that locus as `tRNA-Stop(tta)` here and as `tRNA-Cys(aca)` on the golden, so an intron-containing tRNA is missing from every run | `modules/plastid_annotate.py` + ARAGORN | **open, named H-24**; the run reports it as "Vocabulary loci not called" |

## Re-reading the labels

Two naming traps this corpus sets, both measured:

- **`PX413286` is *C. chlorosarca* in GenBank and *C. jozana* in Table S1** — the paper's own
  Tables S1 and S7 use the second. Column-to-accession maps must come from Table S1, never from
  the `ORGANISM` line.
- **Five records are all named *Crataegus sp.*** and several pairs share a length; a comparison
  must pair by coordinates or by sequence digest, never by name or by length.
