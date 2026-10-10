# Platforms, modules and memory (M9.4)

**Why this file exists.** A researcher deciding whether to install this toolkit asks three
questions: *does it run on my machine · what does it need installed · will it fit in my RAM*. The
answers are here, per module, with the measurement each number comes from and the honest "not
measured" where there is none. The cost figures are `docs/M8.3_scale.md`'s, measured on **one
stated machine** — an **Intel i7-1165G7 @ 2.80 GHz, 6 vCPU under WSL2, 8,135,488 kB RAM (~7.8 GB)**,
the ordinary-laptop class this project targets. A RAM number without its machine is not a number.

## The three support levels

| Level | What it means | What it is verified by |
|---|---|---|
| **Linux — official** | The supported path. Native install, `pip install plastid-suite` (or `uv sync` in a checkout). | The dev machine (WSL2 Ubuntu) **and** CI (`ubuntu-latest`): ruff + the whole fast suite (860+ tests) + tier-A golden + tier-B golden with the pinned tools — green on every push (`.github/workflows/ci.yml`) |
| **Windows — official, through Docker** | The toolkit does **not** run natively on Windows (`docs/things_we_do_not_do.md`). Windows users build `docker/Dockerfile` and run everything inside it. | The image builds and the golden suite runs **inside** it, with `PTS_REQUIRE_TOOLS=1` so a missing tool fails rather than skipping — the same numbers as WSL and CI |
| **macOS — experimental, unverified** | Nothing here has ever run on a Mac. The pure-Python path should work (no platform-specific code), and every pinned tool is available from conda-forge/bioconda for osx-64/arm64 — but *should* is not *measured*. | **Nothing.** No machine, so the column stays empty rather than guessed (ROADMAP → Waiting on data/people) |

## The modules — all 22, and what each costs

`arity` = how many genomes one run takes. **IR** = the module needs the IR gate (the two IR copies
≥ 99.9 % identical), which the framework enforces for every position-based analysis. The cost
column is one row of `docs/M8.3_scale.md` unless it says otherwise; the scope of each measurement
is named, because a 34-genome number and a 1-genome number are different numbers.

| Module | What it does | Group | arity | IR | External tool | Measured cost |
|---|---|---|---|---|---|---|
| `genome_stats` | Length, GC, gene counts, regions | plastid | 1 | yes | — | 34 genomes: **20.5 s · 87 MB** |
| `structure` | LSC/IRa/SSC/IRb boundaries, per-region GC, GC-skew | plastid | 1 | yes | — | 34 genomes: **21.6 s · 87 MB** |
| `gene_table` | Publication gene table (+ `.docx`) | plastid | 1 | no | — | 34 genomes: **15.3 s · 51 MB** |
| `codon_usage` | RSCU · ENC · CAI · CBI · Fop | shared | 1 | yes | — | 34 genomes: **22.8 s · 87 MB** |
| `amino_acid` | Composition, GRAVY, aromaticity | shared | 1 | yes | — | 34 genomes: **23.9 s · 87 MB** |
| `introns` | Count, length, position (CDS introns — H-23) | shared | 1 | yes | — | 34 genomes: **21.9 s · 87 MB** |
| `orf` | ORF checks | shared | 1 | no | — | 34 genomes: **12.9 s · 51 MB** |
| `ssr` | Microsatellites by region and motif | shared | 1 | yes | **MISA 2.1** + perl (obtained, not bundled: `docs/tools/misa.md`) | 34 genomes: **25.5 s · 87 MB** |
| `submission` | NCBI `.fsa` / `.tbl` (+ `.sqn` via table2asn) | shared | 1 | no | table2asn (optional step) | 34 genomes: **12.7 s · 52 MB** |
| `repeats` | Exact repeats ≥ 30 bp, IR pair included | shared | 1 | yes | MUMmer's `repeat-match` | golden genome, 1 run: **0.8 s · 87 MB** (measured 11/10/2026; M8.3 has no row for it) |
| `substitutions` | Ts/Tv, by position; `scope=cds` or `whole` | shared | 2 | yes | MAFFT | `cds`, 34 genomes: **5 min 24 s · 179 MB** · `whole`, 8 genomes: 2 min 44 s · 97 MB — **16 did not finish in 7 min** |
| `diversity` | π, whole / by region / per intron and IGS | shared | 2 | yes | MAFFT | 34 genomes: **13 min 14 s · 138 MB** |
| `barcodes` | Barcode loci and hypervariable regions | plastid | 2 | yes | MAFFT | 34 genomes: **12 min 27 s · 138 MB** |
| `popgen` | Fst, gene flow, diversity per gene | shared | 2 | no | MAFFT | 34 genomes: **6 min 13 s · 179 MB** (80 alignments) |
| `kaks` | Ka/Ks per shared gene (PAML's `yn00`) | shared | 2 | no | MAFFT + PAML | 34 genomes: **8 min 26 s · 98 MB** (9,968 tool calls — call count is not what costs) |
| `synteny` | Whole-genome alignment blocks and identity | shared | 2 | yes | MUMmer (`nucmer` → `delta-filter` → `show-coords`) | 34 genomes: **19.6 s · 145 MB** (99 calls) |
| `plastid_annotate` | Plastid annotation by reference transfer | plastid | 1 | yes | BLAST+ and ARAGORN | **not measured** — M8.3 names the annotation path as having no scale row, and this table does not invent one |
| `gene_names` | Gene-name check and normalization | shared | 2 | no | — | **not measured** |
| `gene_normalize` | Reference-driven naming (A/B pairing) | plastid | 2 | no | BLAST | **not measured** |
| `stats` | Group-comparison statistics (`--groups`) | shared | 2 | no | — | **not measured** |
| `mito_annotate` | Mitochondrial annotation by reference transfer (prototype) | mito | 1 | no | BLAST+ and ARAGORN | **not measured** |
| `mito_genes` | Mitochondrial gene set | mito | 1 | no | — | **not measured** |

**What the table says as a whole** — the boring number worth remembering, from M8.3: single-genome
work over a 34-genome panel is **under four minutes and under 100 MB**; the whole panel costs about
**a quarter of an hour and under 150 MB** for any whole-genome module. What costs is the size of a
*single* alignment, not the number of tool calls: 9,968 calls of two short sequences finish in
8 min 26 s, while one alignment of 34 whole plastomes had not finished in 11 minutes.

**Not measured, named rather than implied** (M8.3's own note, kept here because this is the table
users read): the annotation and naming path and the mito pair have no scale row. That is the honest
alternative to letting "everything is measured" be inferred from a table that lists every module —
and it is the first thing to measure when the real datasets arrive.

**One thing this table could not take from the code** (found 11/10/2026): the registry's
`citations` field is **stored but never read**, and it has drifted — `plastid_annotate`,
`mito_annotate` and `gene_normalize` declare no tool there while their runs cite BLAST+ and
ARAGORN through the `Result`. The "External tool" column above is therefore compiled from the
adapters each module actually calls (`src/plastid/external/`), not from that field. The field is
either fed from the module's own citations or dropped; it must not become a second, quieter source
of truth (recorded as a finding, not silently fixed here).

*The module list is guarded: `tests/unit/test_project_docs.py` fails if a registered module is
missing from this table, so a module added later cannot quietly go undocumented.*
