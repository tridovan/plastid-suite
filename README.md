# PlastidToolSuite

[![CI](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml/badge.svg)](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml)

Toolkit for **plastid (chloroplast) genome analysis**, developed for academic research and
**MIT licensed**. Twenty-four analysis modules: genome statistics, structure (LSC/SSC/IR), gene
tables, codon usage, amino-acid composition, introns, SSR, repeats, substitutions (Ts/Tv),
diversity (π), barcodes, Ka/Ks, population genetics, synteny, group statistics, ORF checks,
gene-name unification, annotation, maximum-likelihood phylogeny, predicted RNA-editing sites and
NCBI submission files, plus assembly from short or HiFi reads. The existing mitochondrial path
(`mito annotate` · `mito genes`) is an experimental prototype; current development targets
plastids.

Every run writes its own directory with tables and a draft report. External tool calls record
their name, version and command line. Output from a new analysis run does not overwrite an older run.

**Status: v1.0.0.** This is the code version — the package metadata says
`Development Status :: 4 - Beta`, and that is the accurate label: no one outside the project has yet
installed it and reproduced its numbers. Every reproduced defect in
[known issues](docs/known_issues.md) is fixed and held by a regression test, and that file also keeps
the limits that remain. What the toolkit does **not** do is in
[docs/things_we_do_not_do.md](docs/things_we_do_not_do.md); which modules need which tool, and how
much memory they take, is in [docs/platforms.md](docs/platforms.md). Every external tool the
toolkit can call, with its role — runtime, cross-check, deferred or survey-only — is listed in
[docs/tools/](docs/tools/). What changed, and what to expect from a number measured with an older
checkout, is in the [release notes](docs/release_notes.md).

The [canonical pipeline](docs/pipeline.md) shows reference preparation, annotation, analysis and
outputs, including assembly from short or HiFi reads, RNA editing (predicted, PREP-Cp-compatible)
and the integration still missing.

## Install

**Linux — install from a source checkout.** The minimum install supports the pure-Python analyses.
These instructions do not assume a published PyPI release:

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
and how to install it: [docs/tools/](docs/tools/) — its front page lists every tool with its role
(runtime · assembly · cross-check · deferred · survey-only). The runtime environment is
`environment.yml` (`conda env create -f environment.yml`). The assembly tools are in
`environment-assembly.yml`, and the cross-check engines and deferred tools in
`environment-check.yml` and `environment-future.yml`, so installing the analyses never pulls in an
assembler.

## The local interface

`plastid ui` opens the toolkit in your browser, on your own computer. Nothing is uploaded: the page
talks to a small server on `127.0.0.1`, and that server starts the same `plastid` command line this
file describes, so a run made there is a run made here.

```bash
plastid ui --workspace my_project
```

The first time, this makes four folders in `my_project`: `records`, `reads`, `references` and `runs`.
The link it prints holds a key for this session, so open that link rather than typing the address.
Closing the tab stops nothing. Ctrl-C closes the interface, and runs already started go on. A busy
port stops the command by name; `--port` chooses another.

### On Windows

The toolkit runs on Windows through Docker. You need no command line:

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and start it.
2. Download this repository as a ZIP and unpack it.
3. Double-click `docker\plastid-ui.bat`. The first time it builds the toolkit, which needs the network
   and takes a long while. After that it starts in seconds, and your browser opens by itself.

Your work goes in `Documents\PlastidToolSuite`. Keep the black window open while you work; closing
it stops the interface. For IR junctions, put your own patched CPJSdraw in that folder and choose
it on the Export page: it is not bundled, because its source states no licence.

### A first run, step by step

1. **Add a record.** Open **Records**, press **Add records** and choose a GenBank file, or copy it into
   the `records` folder. The page reads it and asks the IR gate at once.
2. **Read the gate.** A green mark means the two inverted-repeat copies agree. If the record is refused,
   the page says what it measured and which ways forward keep the check in place. The gate has no off
   switch.
3. **Choose the analyses.** Tick the record and press **Analyse**. Pick the analyses on the left. A
   form appears for any that take parameters, already filled with their defaults. The plan on the
   right says what will run and what is missing, and shows the same run as a command line.
4. **Say where the data came from.** It is printed with every table.
5. **Press Run.** The page follows the work and has a **Cancel run** button. You can close the tab
   and come back.
6. **Read the results.** **Results** lists every run, including runs started at a prompt. A run shows
   its tables with their scope, its figures, and how it was made.
7. **Export.** **Export** draws figures, builds the Supplementary package and the QC report, and draws IR
   junctions across species.

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

From paired short reads to a checked sequence (`conda env create -f environment-assembly.yml`, and
[docs/tools/getorganelle.md](docs/tools/getorganelle.md) for the database it needs once):

```bash
plastid assemble run -1 reads_1.fastq.gz -2 reads_2.fastq.gz --sample my_sample --out .runs \
                     --compare-to a_published_record.gb
plastid annotate -i .runs/my_sample/<run>/01_assembly/assembly.fasta --out .runs \
                 --reference-set a_directory_of_genbank_records
```

The toolkit ships no plastid records: `annotate` transfers from a directory of GenBank records you
supply, same family first, and `plastid panel pin --accessions A,B,C --to that_directory` fetches
and pins one.

PacBio HiFi reads take `--hifi reads.fastq.gz` in place of `-1`/`-2`
([docs/tools/oatk.md](docs/tools/oatk.md)).

The reads are mapped back and every join of the circle is read: `assembly.fasta` is declared for
the next step only when the sequence closed and the reads cross every join, and exit `1` says it
was not. When this machine cannot run the assembly, the command stops with the reason and the
route: `plastid assemble prepare`, the same run on a service or a cluster of your choice, and
`plastid assemble import` to bring the result back. The toolkit never sends reads anywhere by itself.

The command line can tell you what it has, without the documents:

```bash
plastid modules          # every analysis; `plastid modules ssr` shows one, with its parameters
plastid tools            # which external tools this machine has, and their versions
plastid runs             # what has been run, and how each run ended
```

Add `--json` to any of the three for a script. A run is on disk from its first moment: its
`manifest.json` says `running`, moves with each step and holds the command that was typed, so the
same run can be repeated from it. Ctrl-C or a terminate signal stops the run and every tool it
started, and the directory says `cancelled` and keeps what had finished. A run that died without
warning reads as `interrupted` in `plastid runs`.

**Exit codes**, so a script can tell the difference. The status is the word in `manifest.json`:

| Code | Status | Meaning |
|---|---|---|
| `0` | `ok` | finished |
| `1` | `failed` | the work started and did not finish, or a check found a problem |
| `2` | `refused` | nothing was computed: a wrong command line, an input that cannot be read, the **IR gate** (the two inverted-repeat copies are not ≥ 99.9 % identical: see [docs/things_we_do_not_do.md](docs/things_we_do_not_do.md)) or the reference-set audit |
| `3` | `tool_error` | a tool or an optional package is missing or failed, and the message names it |
| `4` | `conflict` | the data disagrees with the chosen naming reference |
| `130` | `cancelled` | stopped on request; what had finished stays in the run directory |

A refused run that reached its directory leaves it, saying why. An assembly that leaves no usable
sequence exits `1` with the status `incomplete`.

## Verify

In a checkout, with the tools you have installed:

```bash
make test        # fast tests — no external tools, no network
make golden      # the golden numbers and pinned-tool checks; names every missing tool
make acceptance  # the tests that read the pinned panels under .panel/
make check       # lint + test + golden
```

The tests that read `.panel/` are marked `panel` and run only under `make acceptance`. The panels
are not in git by design, so a fresh clone cannot run that tier. Record skips explicitly: skipped
cases are unverified, and passing golden tests alone do not establish correctness on other records.
What each kind of golden check proves: `tests/golden/README.md`.

**Golden sample:** *Iris domestica* NC_050833.1 — **153,724 bp · 37.85 % GC · 132 genes = 86 CDS +
38 tRNA + 8 rRNA** · IR 26,214 bp at 100.00 % identity · LSC 83,127 · SSC 18,169 · 37 SSRs
(31 mono · 4 di · 2 tri). Provenance of every pinned record: `tests/data/PROVENANCE.md`.

## Licence and citation

MIT — see [LICENSE](LICENSE); what that means for the tools it calls is in
[docs/licensing.md](docs/licensing.md). If you use this toolkit in published work, cite it with the
metadata in `CITATION.cff`.
