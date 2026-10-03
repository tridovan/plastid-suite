# Manuscript draft — how it is produced, and what it is not

P2-10 asks for *tables + figures + draft Results for the two species*, with every number traced
to a cell in a run's tables. This file is the recipe. The status — what is present, what is
missing and why — is in `TASKS.md` (P2-10) and PROJECT.md §8.8.

## The routine

One dataset, one command, one document. `-m` is repeatable and the run is the deliverable:

```bash
plastid run \
  -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid -m introns -m ssr \
  -i <genome>.gb --docx --out .runs --dataset <label>
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

## What it does not contain yet, and why

| Missing | Why | Unblocked by |
|---|---|---|
| π (nucleotide diversity) and Ts/Tv sections | need MAFFT **7.526**; installing a different version would measure a different aligner (PROJECT.md §7.8, §9.B-13) | a machine with the pinned MAFFT |
| Figures | M5.1 — the figure layer is Phase 5, not Phase 2 | Phase 5 |
| The team's own sample numbers | the real Oldenlandia / Coleus records from anh Điềm / cô Trang | the real dataset |

The draft says so itself: every report carries the `Draft — needs human review.` banner, and an
absent module is an absent section rather than a silently short one.
