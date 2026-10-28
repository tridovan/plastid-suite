# CGAS 1.0.1 — a comparison, never a standard (P8-03, M8.1)

**Why this file exists.** CLAUDE.md's invariant: a tool that is used is declared in the repo. CGAS
is **not used by any module** — it was read and run once so the maintainer can say what this
toolkit does better than an existing chloroplast suite, and where CGAS is ahead. Minimum install
stays pip-only; the release notes name CGAS because the comparison is in the evidence, not because
anything calls it.

**What it is not.** CGAS is **not a benchmark PTS is graded against, and not a definition of
correctness.** Its defaults are its own choices, several of them different from MISA's, PAML's or
this toolkit's; and where it cannot detect something it prints a zero. A CGAS number therefore
proves nothing about a PTS number. What the comparison is good for is the other question — *what
did PTS improve, and what does PTS still lack* — and it is read in both directions below.

**What it is.** Chloroplast Genome Analysis Suite (Abdullah, MIT) — 14 modules, 19 `.py` files,
~21,400 lines. Source: `.data/cgas-review/` (`github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS`,
read 02/10/2026; defaults re-read from source 28/10/2026).

**Install line** (its own environment; nothing here may move the pinned builds in the shared
toolbox env):

```sh
uv venv /root/toolbox/env/cgas --python 3.12
uv pip install --python /root/toolbox/env/cgas/bin/python -e .data/cgas-review
```

**Measured 10/10/2026 on the dev WSL machine:**

| | |
|---|---|
| Version | `cgas --version` → **CGAS v1.0.1** |
| Python-only modules | 3, 4, 5, 6, 7, 10, 11 — these ran |
| Modules 8, 9, 12, 13 with R absent | they **compute every number in Python and skip only the plots** (`⚠ R is not installed or not in PATH → Skipping visualization`) — measured, not assumed |
| R on this machine | **absent** — so no CGAS figure was produced |
| Not run | 1 (needs FASTQ reads), 2 (annotation through PGA, retired from this project), 14 (needs IQ-TREE and MACSE) |
| Licence | MIT |

## What CGAS has that this toolkit does not

| CGAS | PTS today | Status |
|---|---|---|
| **1 — assembly from reads**: fastp → GetOrganelle (`embplant_pt`, kmers 21–105) → bwa/samtools map-back, 15 rounds, mean quality 20 | no assembly module. GetOrganelle and samtools are declared `future` in [tools](README.md); fastp and bwa are not declared at all, because no PTS module can call them | **CGAS ahead.** Assembly is Phase 7 work here |
| **2 — annotation through PGA** (identity 40 %, coverage `0.5,2`, IR 1,000 bp) | `annotate` transfers annotation with its own BLAST/ARAGORN path; PGA is **retired** here | different engine, not a gap — but CGAS's wrapper expects PGA, which this project does not use |
| **5 / 6 / 8 / 9 / 12 — one multi-genome table and figure across every input GenBank** (merged RSCU/amino-acid/SSR wide table, heatmaps, Word gene-content document) | `batch` writes a long-format `aggregate.csv` plus per-sample runs; the multi-sample R figures (F3/F4/F5) pivot it | **mostly closed.** Still missing natively: the publication **gene-content matrix across samples** (module 6's shape). CGAS's own assembly of the manuscript's 34-species panels is the reason that shape exists |
| **11 — tRNA introns** (`trnK-UUU` 2,617 bp, 8 on the golden) | `introns` publishes CDS and tRNA introns separately with a `kind` column since M11.15 | **closed** — this comparison is what found the gap (H-23) |
| **14 — partitioned supermatrix, MACSE codon-aware alignment, separate gene/intron/IGS trees** | `phylogeny` is CDS-only, one engine (MAFFT per gene → IQ-TREE), with opt-in `partition=gene` or `partition=codon` files; no MACSE, no separate intron/IGS matrices | **CGAS ahead on breadth.** MACSE and the separate intron/IGS matrices remain open on M11.21 |
| **4 — GenBank → FASTA/TBL converter plus a validator** | `submission` emits `.fsa`/`.tbl` and first-pass checks, optionally `.sqn` via table2asn | roughly equal; CGAS has no `.sqn` path |

Module 3 (gene comparison against a reference: normalizes names, lists missing/extra genes, flags
CDS-length differences) overlaps `gene_normalize` and `gene_names`. The difference is the mode:
CGAS **writes the renamed record** as its product, while `gene_normalize mode=verify` is expected
to change nothing on an annotated sample and reports a rename as a finding.

`-og/--outgroup` in module 14 **only reorders alignment rows** — it does not root the tree. Neither
tool roots a tree; the `.treefile` from either is unrooted.

## What this toolkit has that CGAS does not

| PTS | CGAS |
|---|---|
| **IR gate**: the two repeat copies must measure ≥ 99.9 % identical, and position-based modules refuse until they do; a de-novo path finds the pair from the sequence | module 7 works **only from annotation** (`misc_feature`/`repeat_region` notes); it never compares the two copies' sequence and has no identity number. A record whose annotation lacks the IR gets `LSC 0 · SSC 0 · IR 0` |
| **Fail loud**: a refusal is an exception with exit 3 and the reason named | module 7 prints `LSC 0 · SSC 0 · IR 0`, **no warning**, and still counts the record under "Successfully analyzed: 2/2" (measured on `NC_024936.1`; PTS finds 26,184 bp @ 100.00 % there, de-novo included) |
| **Provenance**: every tool call records name · version · argv · exit code (`run_tool`) | 51 `subprocess` call sites across nine modules, each assembling its own argv; CGAS prints versions but stores no run record |
| **A number carries its scope**, and an unmeasured value is written "not measured" | no scope label; a zero can mean "not detected", "not applicable" or "failed" |
| **Ka/Ks** (`kaks`, PAML yn00 per shared CDS, Ks = 0 labelled undefined), **population genetics** (`popgen`), **group statistics with FDR** (`stats`), **synteny/identity** (`synteny`), **barcode candidates** (`barcodes`), **ORF checks** (`orf`) | none of these modules exists |
| **Manifest, cache, digest and run directories** — `<out>/<dataset>/<run_id>/` | writes fixed output folders (`Module7_Comparative_Analysis/`, …) next to the input |
| **Exact repeats** (`repeats`, MUMmer `repeat-match` with Vmatch as the second engine) | no exact-repeat module; SSR only |
| **Minimum install = biopython + PyYAML** | biopython, pandas, numpy, openpyxl, python-docx, matplotlib, seaborn up front |

## Where the numbers differ because the definition differs

These are **not** disagreements to settle and not a score for either side: a cell-for-cell claim
needs the flags below stated next to it.

| Quantity | CGAS | PTS | To compare |
|---|---|---|---|
| **SSR (m12)** | `{1:10, 2:5, 3:4, 4:3, 5:3, 6:3}` — golden **54** | MISA's shipped `{1:10, 2:6, 3:5, 4:5, 5:5, 6:5}` — golden **37**; thresholds are refused as a `--set` key, not silently honoured | run CGAS with `-t "10,6,5,5,5,5"`; then both give **37 = 27·10·0 · 31/4/2**, IR zero included |
| **RSCU (m8)** | every `CDS` via `feature.extract()`, concatenated, then split in triplets — no `/codon_start`, no `/transl_table`, no stop trimming; stops printed as codons | `coding_sequence` reads each part on its own strand; a frame or table outside {1, 11} is **refused by name**; stops and partial codons excluded; the initiator is read as its own amino acid (5 of the golden's 86 CDS start GTG/ACG/ATC) | differs on a non-1 `/codon_start`, a foreign table, a terminal stop, or a non-ATG/GTG start; on the golden: 26,263 sense codons, 0 codons different |
| **Gene count (m6)** | counts `gene` features; a gene with no CDS inside it is marked Ψ and dropped | CDS + tRNA + rRNA; `gene` features are a diagnostic only (CLAUDE.md) | differs by construction on records whose `gene` features lack a CDS |
| **IR floor (m7)** | an IR-like region < 1,000 bp is ignored as "tiny", with a warning | no floor; the gate measures identity | CGAS's floor is a setting; its no-IR zero is the defect above, not a floor |
| **Intron ceiling (m11)** | a gap > 15,000 bp (or a `/trans_splicing` qualifier) becomes a trans-splicing junction | CDS and tRNA introns each published with `kind`; a junction reported, never counted; a tRNA gap beyond the caller's documented 3,000 bp is "not assessed" | the same reclassification, with different ceilings (15,000 vs the tRNA caller's 3,000) |
| **Ts/Tv (m10)** | exactly two sequences; gaps skipped, **`N` is not skipped** — it falls outside the twelve substitution types, so it is neither a transition nor a transversion, while two different ambiguous bases are counted as a substitution and then dropped from both sums; prints `Inf` when Tv = 0 | N genomes; a column with any gap **or any base outside ACGT** is dropped before scoring and counted; the ratio is left **blank** when undefined | **measured 28/10/2026 on the same CDS alignment: 323/120 = 2.692 vs 2.6917, every type equal.** An `A`-vs-`R` case also agreed (both 15 substitutions, Ts/Tv 2.0); the ambiguity difference above was not reproduced. The remaining non-comparable case is `Tv = 0`: CGAS prints `Inf`, PTS leaves the ratio blank |
| **π (m13)** | mean pairwise differences per site over columns with no gap and no `N`; no window/step | the same pairwise-difference idea, published with windows, regions, genes, introns, spacers and the usable-column count beside every value | **measured 28/10/2026: per-gene π identical at 6 dp** (`psbA` 0.005650, `matK` 0.010835), same usable sites. Not comparable: CGAS's **plotted** π is rounded to 1 decimal, its regions are one 218-row list, and it prints no whole-alignment π |
| **Phylogeny (m14)** | MAFFT (or MACSE) per gene, partition files, genes/introns/IGS matrices, `-m MFP -bb 1000 -alrt 1000 -bnni`, `-nt AUTO`, **no seed** | MAFFT per gene → one concatenated matrix → IQ-TREE with the same model/bootstrap flags, opt-in gene/codon partition files (`-p`), `-seed 1`, per-gene inclusion rules and a `matrix_genes` exclusion table | same engine family, different matrix and different reproducibility: a CGAS tree need not repeat byte for byte |

## What the measurement actually found

- **It agrees where it detects.** All four faces of the comparison matched the golden: the four
  regions and their GC values, the 61 sense-codon counts, the 17 CDS introns line by line, and —
  after passing `-t "10,6,5,5,5,5"` — the 37 SSRs. That is the useful half of the comparison: it
  says PTS's arithmetic is not inventing differences.
- **The same-alignment comparison was extended on 28/10/2026.** Module 10, given the exact shared-CDS
  alignment PTS builds (79 CDS, 68,448 sites), reports **323 Ts / 120 Tv = 2.692** where PTS reports
  **2.6917**, with every one of the twelve substitution types equal (A→G 92, A→C 28, …). Module 13,
  on the golden pair, gives the **same per-gene π at 6 decimal places** for `psbA`, `matK`, `rps16`
  and `trnK-UUU`, and the same usable-site count for `psbA` (1,062). An `A`-vs-`R` synthetic case
  also agreed: both sides drop the ambiguous columns and report 15 substitutions, Ts/Tv 2.0 — the
  earlier note that CGAS counts ambiguous bases as substitutions was **not reproduced**.
- **Where a cell cannot be compared, it is named.** CGAS's **plotted** π is rounded to 1 decimal
  (`0.0` for a per-gene 0.0056), so only its text report is comparable; it merges genes, introns
  and intergenic spacers into one 218-row list where PTS publishes three tables; and it prints no
  whole-alignment π to put beside PTS's. Module 3 writes a renamed record — a different product, not
  a cell. Module 14 has **not** been run with the now-installed MACSE/IQ-TREE, so no m14 cell is
  claimed.
- **It fails silently where it cannot detect**, which is the improvement this project exists to
  make: on `NC_024936.1`, module 7 finds no region, prints no warning, writes `LSC 0 · SSC 0 ·
  IR 0`, and counts the record as analysed. The record is not wrong; the zero is.
- **It found one gap in PTS**, since closed: the 8 tRNA introns (H-23, M11.15).

**Not vendored, deliberately:** a Python package with its own environment — the same reasoning as
MUMmer, BLAST+ and ARAGORN.
