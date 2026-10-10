# PlastidToolSuite

[![CI](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml/badge.svg)](https://github.com/tridovan/plastid-suite/actions/workflows/ci.yml)

Toolkit for **plastid (chloroplast) and mitochondrial genome analysis** — academic, non-commercial.

**Status: Phases 1–5, 4a, 3B and 2B complete; Phase 8 (cross-checks, scale, CI) open on
`feat/phase-8`.** Twenty analysis modules run end-to-end —
genome statistics, structure (LSC/SSC/IR), gene tables, codon usage, amino-acid composition,
introns, SSR, substitutions (Ts/Tv), diversity (π), barcodes, Ka/Ks, population genetics, synteny,
statistics, ORF checks, gene-name unification, annotation, NCBI submission, and the mitochondrial
path (`mito annotate` · `mito genes`) — each writing a run directory with its tables, a report and
every tool's provenance. Progress: `ROADMAP.md` · current work: `TASKS.md` · rules: `CLAUDE.md`.

## Install

```bash
uv sync              # development install (Python >= 3.10)
# plain pip: pip install -e .
```

The minimum install runs the pure-Python analyses. External tools (MISA, MAFFT, …) are only needed
for the tier-B checks and later phases — a missing tool is always named, never silently skipped.

## Use

```bash
plastid run --module genome_stats --input tests/data/NC_050833.1.gb --out .runs --dataset iris
plastid run --module ssr          --input tests/data/NC_050833.1.gb --out .runs --dataset iris
```

Every run gets its own directory:

```
.runs/iris/<timestamp>-<hash>/
  manifest.json          what ran, with what parameters, on what input, and where every file went
  00_run_info/           run_info.md (Materials & Methods) + tool logs
  01_genome_stats/       tables/*.csv · result.json · report.md (draft — needs human review)
```

A failed gate stops the run with exit code 2 and still leaves a traceable directory explaining why.

## Verify

```bash
make test      # fast tests — tier A, no external tools
make golden    # the golden-number set; names every missing tool
make check     # lint + test + golden
```

Golden sample: *Iris domestica* NC_050833.1 — 153,724 bp · 132 = 86 CDS + 38 tRNA + 8 rRNA ·
IR 26,214 bp @ 100.00 % · 37 SSR (mono 31 · di 4 · tri 2). Full provenance: `PROJECT.md` §6.
