# Platforms, modules and memory (M9.4)

**Why this file exists.** A researcher deciding whether to install this toolkit asks three
questions: *does it run on my machine · what does it need installed · will it fit in my RAM*. The
answers are here, per module, with the measurement each number comes from and the honest "not
measured" where there is none. The cost figures were measured on **one
stated machine** — an **Intel i7-1165G7 @ 2.80 GHz, 6 vCPU under WSL2, 8,135,488 kB RAM (~7.8 GB)**,
the ordinary-laptop class this project targets. A RAM number without its machine is not a number.

## The three support levels

| Level | What it means | What it is verified by |
|---|---|---|
| **Linux — supported development path** | Native install from a checkout: `uv sync` or `pip install .`. | Tested on WSL2 Ubuntu; CI runs fast and golden tests with its declared tools. This is not complete product acceptance; see [known correctness issues](known_issues.md) |
| **Windows — official, through Docker** | The toolkit does **not** run natively on Windows (`docs/things_we_do_not_do.md`). Windows users build `docker/Dockerfile` and run everything inside it. | The image builds and the golden suite runs **inside** it, with `PTS_REQUIRE_TOOLS=1` so a missing tool fails rather than skipping — the same numbers as WSL and CI. Measured 28/10/2026: **64 passed / 61 skipped** before the image carried R and the optional extras; **70 passed / 55 skipped** once it did — and all 55 remaining skips are "`.panel/` is not in this tree", because the panel is third-party data and is mounted at run time, not built in. With the repository mounted (`-v "$PWD:/work"`) the image runs the whole golden |
| **macOS — experimental, unverified** | Nothing here has ever run on a Mac. The pure-Python path should work (no platform-specific code), and every pinned tool is available from conda-forge/bioconda for osx-64/arm64 — but *should* is not *measured*. | **Nothing.** No machine, so the column stays empty rather than guessed (ROADMAP → Waiting on data/people) |

## The modules — all 23, and what each costs

`arity` = how many genomes one run takes. **IR** = the module needs the IR gate (the two IR copies
≥ 99.9 % identical), which the framework enforces for every position-based analysis. The cost
column is one measured row unless it says otherwise; the scope of each measurement
is named, because a 34-genome number and a 1-genome number are different numbers.

| Module | What it does | Group | arity | IR | External tool | Measured cost |
|---|---|---|---|---|---|---|
| `genome_stats` | Length, GC, gene counts, regions | plastid | 1 | yes | — | 34 genomes: **20.5 s · 87 MB** |
| `structure` | LSC/IRa/SSC/IRb boundaries, per-region GC, GC-skew | plastid | 1 | yes | — | 34 genomes: **21.6 s · 87 MB** |
| `gene_table` | Publication gene table (+ `.docx`) | plastid | 1 | no | — | 34 genomes: **15.3 s · 51 MB** |
| `codon_usage` | RSCU · ENC · CAI · CBI · Fop | shared | 1 | no | — | 34 genomes: **22.8 s · 87 MB** |
| `amino_acid` | Composition, GRAVY, aromaticity | shared | 1 | no | — | 34 genomes: **23.9 s · 87 MB** |
| `introns` | Count, length, position and kind — CDS and tRNA introns, trans-splicing junctions (H-23) | shared | 1 | yes | — | 34 genomes: **21.9 s · 87 MB** |
| `orf` | ORF checks | shared | 1 | no | — | 34 genomes: **12.9 s · 51 MB** |
| `ssr` | Microsatellites by region and motif | shared | 1 | yes | **MISA 2.1** + perl (obtained, not bundled: `docs/tools/misa.md`) | 34 genomes: **25.5 s · 87 MB** |
| `submission` | NCBI `.fsa` / `.tbl` (+ `.sqn` via table2asn) | shared | 1 | no | table2asn (optional step) | 34 genomes: **12.7 s · 52 MB** |
| `repeats` | Exact repeats ≥ 30 bp, IR pair included | shared | 1 | yes | MUMmer's `repeat-match` | golden genome, 1 run: **0.8 s · 87 MB** (measured 11/10/2026; M8.3 has no row for it) |
| `substitutions` | Ts/Tv, by position; `scope=cds` or `whole` | shared | 2 | yes | MAFFT | `cds`, 34 genomes: **5 min 24 s · 179 MB** · `whole`, 8 genomes: 2 min 44 s · 97 MB — **16 did not finish in 7 min** |
| `diversity` | π, whole / by region / per intron and IGS | shared | 2 | yes | MAFFT | 34 genomes: **13 min 14 s · 138 MB** |
| `barcodes` | Barcode loci and hypervariable regions | plastid | 2 | yes | MAFFT | 34 genomes: **12 min 27 s · 138 MB** |
| `popgen` | Fst, gene flow, diversity per gene | shared | 2 | no | MAFFT | 34 genomes: **6 min 13 s · 179 MB** (80 alignments) |
| `kaks` | Ka/Ks per shared gene (PAML's `yn00`) | shared | 2 | no | MAFFT + PAML | 34 genomes: **8 min 26 s · 98 MB** (9,968 tool calls — call count is not what costs) |
| `phylogeny` | ML tree from a per-locus supermatrix — CDS (default), intron or IGS blocks (`matrix=`) | shared | 2 (at least 3 genomes) | no | MAFFT + IQ-TREE | 4 genomes: **~50 s · 179 MB** (78 blocks, 70,079 sites) · 34 genomes: **10 min · 179 MB** (78 blocks, 70,091 sites, IQ-TREE `-bb 1000 -alrt 1000`) |
| `synteny` | Whole-genome alignment blocks and identity | shared | 2 | yes | MUMmer (`nucmer` → `delta-filter` → `show-coords`) | 34 genomes: **19.6 s · 145 MB** (99 calls) |
| `plastid_annotate` | Plastid annotation by reference transfer | plastid | 1 | yes | BLAST+ and ARAGORN | **not measured** — M8.3 names the annotation path as having no scale row, and this table does not invent one |
| `gene_names` | Gene-name check and normalization | shared | 2 | no | — | **not measured** |
| `gene_normalize` | Reference-driven naming (A/B pairing) | plastid | 2 | no | BLAST | **not measured** |
| `stats` | Group-comparison statistics (`--groups`) | shared | 2 | no | — | **not measured** |
| `mito_annotate` | Mitochondrial annotation by reference transfer (prototype) | mito | 1 | no | BLAST+ and ARAGORN | **not measured** |
| `mito_genes` | Mitochondrial gene set | mito | 1 | no | — | **not measured** |

**The IR column is a dependency, not a quality grade.** Since M11.5, `codon_usage` and
`amino_acid` declare `requires_ir=False`: every number they print comes from CDS sequences, so a
record with readable CDS runs even when its repeat fails the gate, is linear or is N-rich. They
refuse a CDS whose `/codon_start` is not 1 or whose `/transl_table` is neither 1 nor 11, by name —
reading the wrong frame or another code silently would produce a different measurement. Position-
based modules keep the gate. The panel check is
`tests/golden/test_module_ir_requirements_panel.py` (5 gate-fail, 31 ambiguous, 3 linear, 2
reduced records).

**What the table says as a whole** — the boring number worth remembering, from M8.3: single-genome
work over a 34-genome panel is **under four minutes and under 100 MB**; the whole panel costs about
**a quarter of an hour and under 150 MB** for any whole-genome module. What costs is the size of a
*single* alignment, not the number of tool calls: 9,968 calls of two short sequences finish in
8 min 26 s, while one alignment of 34 whole plastomes had not finished in 11 minutes.

**What a whole pipeline run costs** (M11.8, measured 28/10/2026, same machine). One command over
ten modules on the golden record — `genome_stats structure gene_table codon_usage amino_acid
introns orf ssr repeats submission`, with `--figures --docx` — took **26.3 s at 1.06 GB peak**. The
figure pass is what costs the memory, not the analyses: `genome_stats` alone is 87 MB. Annotating a
bare FASTA against the 34-record reference set took **44.5 s at 125 MB**, and four pairwise modules
over two genomes (`substitutions diversity kaks popgen`, CDS scope) **56.5 s at 111 MB**. All three
fit an 8 GB laptop; quote these when a user asks whether the whole flow fits, rather than the
per-module rows above, which are analyses without a figure pass.

**What a fresh install actually does** (M11.19, measured 28/10/2026 on this machine). A wheel built
from the tree and installed into an empty venv, `pip install
plastid_suite-1.0.0-py3-none-any.whl`, resolves **four packages** — `plastid-suite`, `biopython`,
`numpy` (Biopython's own dependency) and `PyYAML` — and with only those, one command runs seven
modules on the golden record:

```bash
plastid run -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid \
            -m introns -m orf -i NC_050833.1.gb -o runs --dataset minimum
```

Every module reports `ok` and 35 files are written. What needs more says so and nothing else
happens:

| Asked for | Result in the minimum install | Exit |
|---|---|---|
| `-m ssr` (needs `misa.pl`) | names the download line and the `--set misa=` option | **3** |
| `--docx` | `install it with: pip install 'plastid-suite[docx]'` | **3** |
| `--figures` | `install it with: pip install 'plastid-suite[figures]'` | **3** |
| `figures <run>` with R present, matplotlib absent | draws the R half, records `matplotlib figures unavailable: …` | **0** |

The last row is the deliberate difference: `--figures` is the user asking for figures, so a missing
plotting library stops the command; `figures` is a collection step and carries what it can, with
the absence written into the manifest and `qc` rather than dropped silently. A second fresh venv
with `[figures,docx,excel,stats]` produced the PDF/SVG/EPS figures and the Word report from the
same record, so the extras are independent and none is needed for the others.

**The tool files are split by role** (M11.19): `environment.yml` is the runtime analyses;
`environment-check.yml` holds the cross-check engines no module calls; `environment-future.yml`
holds the deferred assembly and survey tools. The role of every tool is written in
[docs/tools/README.md](tools/README.md), and `tests/unit/test_tool_inventory.py` fails if any of
the three files and that table disagree.

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

The previous module-list check was removed with `tests/unit/test_project_docs.py` during the
documentation cleanup. It is restored in `tests/unit/test_product_docs.py` (M11.14), independently
of private session documents: that test asserts every registered module has a row here and that the
README status line names the version in `pyproject.toml`, and it includes negative controls that
prove both checks can fail.
