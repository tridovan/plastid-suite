# PlastidToolSuite

[![CI](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml/badge.svg)](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml)

Toolkit for **plastid (chloroplast) genome analysis**, developed for academic research and
**MIT licensed**. Twenty-four analysis modules: genome statistics, structure (LSC/SSC/IR), gene
tables, codon usage, amino-acid composition, introns, SSR, repeats, substitutions (Ts/Tv),
diversity (π), barcodes, Ka/Ks, population genetics, synteny, group statistics, ORF checks,
gene-name unification, annotation, maximum-likelihood phylogeny, predicted RNA-editing sites and
NCBI submission files, plus assembly from short or HiFi reads. The existing mitochondrial path
(`plastid mito annotate` and the `mito_genes` module) is an experimental prototype; current
development targets plastids.

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

Nothing is on PyPI yet, so the toolkit comes from this repository. Windows runs it through Docker,
and Linux through a conda environment. Either way you get the local interface and every analysis.

### On Windows

You need no command line:

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/), start it, and wait until
   it says it is running.
2. Download this repository as a ZIP (**Code**, then **Download ZIP**) and unpack it.
3. Double-click `docker\plastid-ui.bat`. If Windows warns about the file, press **More info**, then
   **Run anyway**. The first time it builds the toolkit, which needs the network, about 15 minutes on
   a laptop and about 10 GB of free disk. After that it starts in seconds, and your browser opens by
   itself.

Your work goes in `Documents\PlastidToolSuite`. Keep the black window open while you work. If it was
closed and the toolkit is still running, double-click the file again: it opens the running one. To
stop the toolkit, stop `plastid-ui` in Docker Desktop, under Containers. A run still going stops with
it, and reads as interrupted.

Every tool is in the image you built except CPJSdraw, which only the IR junction figure needs: its
source states no licence. [docs/tools/cpjsdraw.md](docs/tools/cpjsdraw.md) has the clone line and the
two patches. Copy the patched `CPJSdraw` folder into the `tools` folder of your work folder and press
**Check again** on **Tools**.

To take a newer version, stop `plastid-ui`, delete the image `plastid-suite` in Docker Desktop, under
Images, and double-click the file in the new folder. It builds again, in about 5 minutes if the old
build is still cached.

The same image takes every command in this file. In PowerShell, from a folder you want to work in:

```powershell
docker run --rm -v "${PWD}:/work" plastid-suite plastid --version
```

### On Linux

You need [conda](https://github.com/conda-forge/miniforge), or micromamba with the same arguments. The
environment holds the analysis tools, about 3 GB:

```bash
git clone https://github.com/tridovan/plastid-suite && cd plastid-suite
conda env create -f environment.yml
conda activate plastid-suite
pip install '.[figures,docx,excel,stats]'
```

MISA, which the SSR analysis runs, is a download of its own, because its zip states no licence:

```bash
curl -fsSLO https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip
unzip -j misa_sourcecode_25082020.zip misa.pl -d "$CONDA_PREFIX/bin"
```

`plastid tools` then says what this machine has, and how to get what it lacks. Activate the
environment in each new shell.

The assembly tools take an environment of their own, about 3 GB more. Put it last on `PATH`, so the
analyses keep their pinned BLAST, and fetch GetOrganelle's database once:

```bash
conda env create -f environment-assembly.yml
export PATH="$PATH:$(conda info --base)/envs/plastid-suite-assembly/bin"
conda run -n plastid-suite-assembly get_organelle_config.py -a embplant_pt,embplant_mt --config-dir ~/getorganelle
export GETORG_PATH=~/getorganelle
```

HiFi reads also need Oatk's profiles ([docs/tools/oatk.md](docs/tools/oatk.md)). The `.sqn` step needs
table2asn ([docs/tools/table2asn.md](docs/tools/table2asn.md)), and the IR junction figure needs
CPJSdraw ([docs/tools/cpjsdraw.md](docs/tools/cpjsdraw.md)).

### The pure-Python analyses alone

`pip install .` in a clone needs nothing but Python 3.10 or newer, and runs every analysis that calls
no external tool. A missing tool or extra stops the run and names it; nothing is skipped silently.

What each tool is and how to install it: [docs/tools/](docs/tools/), whose front page lists every tool
with its role (runtime · assembly · cross-check · deferred · survey-only). The cross-check engines and
deferred tools are in `environment-check.yml` and `environment-future.yml`.

## The local interface

`plastid ui` opens the toolkit in your browser, on your own computer. Nothing is uploaded: the page
talks to a small server on `127.0.0.1`, and that server starts the same `plastid` command line this
file describes, so a run made there is a run made here. On Windows, `plastid-ui.bat` does this for
you. On Linux:

```bash
plastid ui --workspace my_project
```

The first time, this makes `my_project` with five folders in it: `records`, `reads`, `references`,
`runs` and `tools`, where a tool you keep yourself goes. The link it prints holds a key for this
session, so open that link rather than typing the address. Closing the tab stops nothing. Ctrl-C
closes the interface, and runs already started go on. A busy port stops the command by name; `--port`
chooses another.

### A first run, step by step

1. **Add a record.** Open **Records**, press **Add records** and choose a GenBank file, or copy it into
   the `records` folder. To try it on a known record, use `tests/data/NC_050833.1.gb` from this
   repository. The page reads it and asks the IR gate at once.
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
   its tables with their scope, its figures, and how it was made. On the known record, **Genome
   statistics** gives 153,724 bp, 37.85 % GC and 132 genes.
7. **Export.** **Export** draws figures again with your choices, such as the root of a tree, builds the
   Supplementary package and the QC report, and draws IR junctions across species.

From reads, start at **Assemble**. From a sequence with no genes, start at **Annotate**, with reference
records you found and audited on **References**, where you also choose the record that names the genes.
**Tools** says what this computer has and how to get what it lacks.

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

From paired short reads to a checked sequence, with the assembly tools of [On Linux](#on-linux)
([docs/tools/getorganelle.md](docs/tools/getorganelle.md) has the method):

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

In a checkout, with the tools you have installed. `make` runs through
[uv](https://docs.astral.sh/uv/), so run `uv sync --all-extras` once first:

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
