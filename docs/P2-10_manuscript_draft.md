# Manuscript draft — how it is produced, and what it is not

P2-10 asks for *tables + figures + draft Results for the two species*, with every number traced
to a cell in a run's tables. This file is the recipe. The status — what is present, what is
missing and why — is in `ROADMAP.md` (M2.10) and PROJECT.md §8.8.

## The routine

One dataset, one command, one document. `-m` is repeatable and the run is the deliverable. Since
P2B-03 it is **two runs per species**, because the toolkit's arity rule forbids mixing them: an
arity-1 module takes *exactly one* genome, so the single-genome tables and the comparison tables
cannot travel in one command. Measured 10/10/2026: the mixture exits 2 with *"module 'genome_stats'
takes exactly one genome, got 2"* — a refusal, not a silently partial run.

```bash
# 1. the species' own genome — every single-genome table and figure
plastid run \
  -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid -m introns -m ssr \
  -m repeats \
  -i <genome>.gb --docx --figures --out .runs/phase2b --dataset <label> --label "<label text>"

# 2. the same species against its panel — π and Ts/Tv
plastid run \
  -m substitutions -m diversity \
  -i <species>.gb -i <congener>.gb [-i …] \
  --set scope=cds --docx --figures --out .runs/phase2b --dataset <label> --label "<label text>"
```

`--set scope=cds` matters for `substitutions`: its default is the **whole-alignment** scope and the
manuscripts quote the **CDS scope** (§6/§8.25/§8.27). `diversity` needs no scope switch — every π
carries its own label, and its summary prints the reference's own `n_introns` and `n_spacers`, so
the π tables join the `introns` table row for row (P2B-01).

The comparison run's reference is the species record (row 0 of the alignment is the reference for
every per-region, per-gene, per-intron and per-spacer number), and its panel is **the species'
own genus** — not the whole 11-record panel. Measured 10/10/2026 (P2B-03, §8.59): over the mixed
two-genus panel the whole-plastome alignment spends most of its columns on indels (42,002 of
273,866 usable, 15 %) and the two IR copies come out asymmetric (IRa 3,008 usable columns at
π 0.0866 against IRb 25,215 at 0.0252) — the two copies of one genome must mirror. Over the
genus set: 148,996 of 157,463 columns usable (95 %) and IRa 25,245 / IRb 26,124 at π 0.0058 /
0.0065.

The exact commands, as run (Oldenlandia's panel is every annotated *Scleromitrion*/*Oldenlandia*
record in `.panel/oldenlandia-coleus/`; Coleus's is every annotated *Coleus*):

```sh
P=.panel/oldenlandia-coleus
OL="stand-in (public records, not the team's samples)"
CO="stand-in panel (public records) — the reference PX923121.1 is the group's own GenBank record of C. amboinicus"

plastid run -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid \
  -m introns -m ssr -m repeats -i $P/PP078674.1.gb -o .runs/phase2b \
  --dataset standin-oldenlandia-diffusa --label "$OL" --docx --figures

plastid run -m substitutions -m diversity \
  -i $P/PP078674.1.gb -i $P/NC_057985.1.gb -i $P/NC_057983.1.gb -i $P/PX977793.1.gb \
  -o .runs/phase2b --dataset standin-oldenlandia-diffusa --label "$OL" \
  --set scope=cds --docx --figures

plastid run -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid \
  -m introns -m ssr -m repeats -i $P/PX923121.1.gb -o .runs/phase2b \
  --dataset own-record-coleus-amboinicus --label "$CO" --docx --figures

plastid run -m substitutions -m diversity \
  -i $P/PX923121.1.gb -i $P/NC_057679.1.gb -i $P/NC_058328.1.gb -i $P/NC_069583.1.gb \
  -o .runs/phase2b --dataset own-record-coleus-amboinicus --label "$CO" \
  --set scope=cds --docx --figures
```

What comes out, under `.runs/<label>/<run-id>/`:

| Path | What it is |
|---|---|
| `NN_<module>/tables/*.csv` | one file per table, with its scope and coordinate frame in the header |
| `NN_<module>/result.json` | the same tables with raw values, for recomputation |
| `NN_<module>/report.md` | the draft Results prose, every number a table-cell reference |
| `NN_<module>/report.docx` | that module's own Word document |
| `report/report.docx` | **the draft manuscript** — every module's prose and tables in one file |
| `manifest.json` | what ran, with which parameters, on which inputs (sha256), producing which files |
| `00_run_info/run_info.md` | the tool + parameter report |

`report/report.docx` is assembled by joining the modules' own documents rather than re-rendering
from their results, because a cache hit restores *files*, not results — re-rendering would
produce an empty document on the second run of an unchanged input. Each module's page setup,
tables and section breaks survive the join, and the whole thing is byte-reproducible.

## Labels — the rule that does not bend

A draft built from public records is **not** the team's data. `--dataset` is what carries that
into every output: the document header prints it, the footer repeats it, and the manifest records
it. There are three cases, and they are not the same:

| Where the record came from | Label | Example |
|---|---|---|
| Someone else's submission, used as a stand-in | `standin-<species>` | `standin-oldenlandia-diffusa` (`PP078674.1`, Shen C.) |
| **The group's own published record** | `own-record-<species>` | `own-record-coleus-amboinicus` (`PX923121.1`) |
| The samples this manuscript reports | the dataset's own name | to be set when anh Điềm delivers them |

`own-record-` means the plastome is the group's, so the numbers are the group's — but *published*
is still not the same as *the specimen this manuscript reports*. Ask before quoting one as the
other. When the manuscript's own samples arrive, re-run the same command with the real files and a
real label. Do not patch numbers into a document.

## What the draft contains

Everything a module measured, with the scope label attached:

- genome statistics · regions · the full gene list (every copy, both IRs)
- the publication gene table, grouped by function, with its own reconciliation
- codon usage (RSCU · ENC · CAI · CBI · Fop), amino-acid composition (GRAVY · aromaticity)
- introns (count · length · phase · trans-splicing junctions) and SSR by region and motif
- **repeats** — maximal exact matches ≥30 bp, dispersed and palindromic, each with both copies'
  coordinates and plastid region (M2.12, P2B-02)
- **π and Ts/Tv** — from the comparison run: whole-alignment π, π per region · gene · intron ·
  intergenic spacer, and the CDS-scope Ts/Tv (pooled and per pair)

### What the two runs measured (10/10/2026)

Recorded here as the draft's own numbers — every one of them is a cell in the runs' tables under
`.runs/phase2b/`, and the run directories are what the manuscript tables are copied from.

| | *Scleromitrion diffusum* `PP078674.1` | *Coleus amboinicus* `PX923121.1` |
|---|---|---|
| Length · GC | 153,676 bp · 37.61 % | 152,544 bp · 37.75 % |
| Genes | 131 = 86 CDS + 37 tRNA + 8 rRNA | 133 = 88 CDS + 37 tRNA + 8 rRNA |
| Introns · junctions | 16 · 11 genes · 12,287 bp · 2 (`agrees`) | 17 · 12 genes · 12,134 bp · 2 (`agrees`) |
| SSR | 32 (`agrees`) | 32 (`agrees`) |
| Repeats ≥30 bp | 4 (0 dispersed · 4 palindromic), longest 26,323 bp = the IR pair | 4 (2 · 2), longest 25,711 bp = the IR pair |
| π (whole alignment) | 0.023522 (148,996 / 157,463 usable) | 0.008462 (151,050 / 153,590) |
| π by region | LSC 0.030647 · IRa 0.005816 · IRb 0.006450 · SSC 0.042325 | LSC 0.010759 · IRa 0.001486 · IRb 0.001486 · SSC 0.018165 |
| π per intron (outside the IR) | 10 loci · 0.022237–0.053883 | 11 loci · 0.005594–0.023922 |
| π per intergenic spacer | 91 · 0–0.164751 (`rpl22–rps19`) | 93 · 0–0.080330 (`ccsA–ndhD`) |
| Ts/Tv, CDS scope (pooled) | 1.4616 = 2,017 Ts / 1,380 Tv · 77 shared genes | 1.5134 = 619 Ts / 409 Tv · 79 shared genes |
| Ts/Tv per pair | `NC_057985.1` 0.9583 · `NC_057983.1` 1.5290 · `PX977793.1` 1.5337 | `NC_057679.1` 1.6240 · `NC_058328.1` 1.5000 · `NC_069583.1` 1.3387 |

Two of those numbers are worth a sentence in a Results section, and both are cells rather than
claims: the *S. diffusum* pair `PP078674.1` × `NC_057985.1` — two records of the same species —
gives **0.9583**, barely a thousandth of the divergence the congeners show, and the inverted
repeat pair is the longest exact match in every one of these genomes, which is what the
repeat table says without being told.

### Figures

`--figures` renders, from the tables the run stored, every figure whose source tables exist —
in `svg` (editable), `pdf` (print) and `eps` where the renderer supports it:

| Figure | From | For |
|---|---|---|
| `circular_map` | genome_stats + structure | the genome map every plastome paper opens with |
| `genome_page` | the same, plus codon/gene tables | the one-page genome summary panel |
| `ir_junction` | structure | the LSC/IR · SSC/IR boundaries |
| `rscu_heatmap` | codon_usage | codon usage (RSCU) |
| `ssr_distribution` | ssr | SSR by region and motif |
| `pi_window` | diversity | the mVISTA-style π curve along the genome |

The phase-4 runs (`.runs/phase4/`, `docs/P4-05_comparison_tables.md`) carry the comparison
figures — `dotplot`, `identity_synteny`, `barcode_assessment` — for the same two panels.

### π and Ts/Tv — what the sections say, and what they do not

**π** comes from `diversity` (MAFFT 7.526, §7.8): gapped columns are dropped, never scored, and
every number is printed next to the count of columns it was computed from. The per-intron and
per-spacer tables are the spans of **P2B-01**, derived in `core/spans.py` — the same derivation
the `introns` module publishes, so the two tables join row for row.

**Ts/Tv** comes from `substitutions` in the **CDS scope** — the shared genes, one copy each, IR
excluded — which is the only scope in which the number is comparable to §6's golden ratio. The
whole-alignment scope is a *different* number and every table says which one it is.

Neither number is the golden genome's: §6's **2.6917 over 79 shared genes** belongs to
*Iris domestica* × *I. gatesii*, and a species table printing 2.6917 for another pair would be
quoting a different genome's measurement. What the golden value is for is the **toolkit**: the
golden suite reproduces it on every run (`tests/golden/test_tier_b.py`), which is what makes any
other pair's number trustworthy.

## What it does not contain yet, and why

| Missing | Why | Unblocked by |
|---|---|---|
| The team's own sample numbers | the real Oldenlandia / Coleus records from anh Điềm / cô Trang | the real dataset |
| The manuscript's own comparison set | the tables are built against the **pinned public panel** (`.panel/oldenlandia-coleus/`); the congeners a paper actually cites are a taxonomy choice — the PI's, not the run's | cô Trang's call, then re-run |
| Tandem repeats | `repeat-match` reports maximal exact matches; the **tandem** class is `mreps`/`trf` (M2.12's other half, still open) | a named task, when a manuscript needs it |

The draft says so itself: every report carries the `Draft — needs human review.` banner, and an
absent module is an absent section rather than a silently short one.
