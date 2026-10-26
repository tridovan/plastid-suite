# PlastidToolSuite

[![CI](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml/badge.svg)](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml)

Toolkit for **plastid (chloroplast) genome analysis**, developed for academic research and
**MIT licensed**. Twenty-two analysis modules: genome statistics, structure (LSC/SSC/IR), gene
tables, codon usage, amino-acid composition, introns, SSR, repeats, substitutions (Ts/Tv),
diversity (π), barcodes, Ka/Ks, population genetics, synteny, group statistics, ORF checks,
gene-name unification, annotation and NCBI submission files. The existing mitochondrial path
(`mito annotate` · `mito genes`) is an experimental prototype; current development targets plastids.

Every run writes its own directory with tables and a draft report. External tool calls record
their name, version and command line. Output from a new analysis run does not overwrite an older run.

**Status: v1.0.0.** This is the code version, not a claim that acceptance is complete. Open
[correctness issues](docs/known_issues.md) affect IR detection, caching and annotation handling;
the toolkit is under corrective development. What it currently does **not** do is in
[docs/things_we_do_not_do.md](docs/things_we_do_not_do.md); which modules need which tool, and how
much memory they take, is in [docs/platforms.md](docs/platforms.md).

The [canonical pipeline](docs/pipeline.md) shows reference preparation, annotation, analysis and
outputs, including deferred assembly, the open RNA-editing work and the integration still missing.

## Install

**Linux — install from a source checkout.** The minimum install supports the pure-Python analyses;
the open correctness issues still apply. These instructions do not assume a published PyPI release:

```bash
git clone https://github.com/tridovan/plastid-suite && cd plastid-suite
uv sync                 # or: pip install .
```

**Windows** — the toolkit does not run natively on Windows; the supported route is Docker (build
once, then use `docker run` for every command in this file):

```bash
docker build -f docker/Dockerfile -t plastid-suite .
docker run --rm -v "$PWD:/work" plastid-suite plastid --version
```

Some analyses need external tools (MAFFT, MISA, BLAST+, MUMmer, PAML…). The toolkit never bundles
them, and never silently skips one: a missing tool stops the run and names it. What each tool is
and how to install it: [docs/tools/](docs/tools/). The tool-equipped environment for everything at
once is `environment.yml` (`conda env create -f environment.yml`).

## Use

```bash
plastid run --module genome_stats --input my_plastome.gb --out .runs --dataset my_sample
```

Two genomes to compare, and a module that reads them both:

```bash
plastid run --module substitutions --input reference.gb --input sample.gb \
            --out .runs --dataset my_comparison --set scope=cds
```

Every run gets its own directory, named by time and input, so two runs never mix:

```
.runs/my_sample/20261011-090000-1a2b3c4d/
  manifest.json          what ran, with which parameters, on which input, and where every file went
  00_run_info/           run_info.md (Materials & Methods) + the logs of every tool call
  01_genome_stats/       tables/*.csv · result.json · report.md (draft — needs human review)
```

**Exit codes**, so a script can tell the difference: `0` fine · `1` an analysis failed · `2` the run
was refused *before* any position-based work — a wrong command line, the **IR gate** (the two
inverted-repeat copies are not ≥ 99.9 % identical: see
[docs/things_we_do_not_do.md](docs/things_we_do_not_do.md)) or the reference-set audit · `3` a tool
or an optional package is missing, and the message names it · `4` a naming conflict — the data
disagrees with the chosen naming reference. A refused run still leaves its directory, saying why.

## Verify

In a checkout, with the tools you have installed:

```bash
make test      # fast tests — no external tools, no network
make golden    # the golden-number set; names every missing tool
make check     # lint + test + golden
```

The corpus half of `make golden` needs the pinned panels under `.panel/`, which are not in git by
design — those tests skip with a named reason in a fresh clone. Record skips explicitly: skipped
cases are unverified, and passing golden tests alone do not establish correctness on other records.

**Golden sample:** *Iris domestica* NC_050833.1 — **153,724 bp · 37.85 % GC · 132 genes = 86 CDS +
38 tRNA + 8 rRNA** · IR 26,214 bp at 100.00 % identity · LSC 83,127 · SSC 18,169 · 37 SSRs
(31 mono · 4 di · 2 tri). Provenance of every pinned record: `tests/data/PROVENANCE.md`.

## Licence and citation

MIT — see [LICENSE](LICENSE); what that means for the tools it calls is in
[docs/licensing.md](docs/licensing.md). If you use this toolkit in published work, cite it with the
metadata in `CITATION.cff`.
