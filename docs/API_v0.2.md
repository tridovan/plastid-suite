# API v0.2 — the multi-genome contract

Extends `docs/API_v0.1.md`. Nothing in v0.1 changes: an arity-one module is untouched.

## What was decided (D2-01, 2026-10-03)

π, SNP and later Ka/Ks take more than one genome. The v0.1 contract `run(genome, params, ctx)`
could not express that, and widening it to `run(genomes, …)` for *every* module would hand each
single-genome module a sequence it must index at `[0]` — a contract wider than the module needs,
which code seam 5 (small explicit contracts) exists to prevent.

So the two contracts stay separate and a module declares which one it follows:

```python
@module("diversity", version="0.1.0", title="Nucleotide diversity (π)", requires_ir=True, arity=2)
def run(genomes: Sequence[Genome], params: dict[str, Any], ctx: RunContext) -> Result: ...
```

- `arity=1` (the default) — `run(genome, params, ctx) -> Result`, exactly as in v0.1.
- `arity≥2` — `run(genomes, params, ctx) -> Result`, with **at least** `arity` genomes; more is
  allowed because a diversity estimate over N accessions wants all of them.

The runner enforces the count before anything runs: a mismatch raises, and the CLI turns it into
a usage error (exit 2). A module never receives genomes it does not use.

## The command line

`--input` is repeatable, and the order is meaningful (row 0 is the reference):

```
plastid run --module substitutions \
  --input ref.gb --input other.gb --out .runs --dataset iris-pair
```

`arity=1` modules still take exactly one `--input`; passing two is a usage error.

## The IR gate under arity ≥ 2

Every genome must pass the gate, not only the first: one diverged IR corrupts every position-based
comparison that follows. The manifest records the primary genome's IR at `ir` and, when more than
one genome was given, every genome's measured length, identity and source at `ir.genomes`. A
failing gate anywhere stops the run with exit 2 and no module output.

## The alignment

The comparison modules align through MAFFT (`src/plastid/external/mafft.py`, pinned 7.526 per
PROJECT.md §7.8) and never fall back to a hand-rolled aligner — a measured fallback once produced
Ts/Tv 2.667 where MAFFT gave 2.686 (PROJECT.md §9.B). A missing aligner is a named `ToolError`
and exit 3.

All alignment maths lives in `src/plastid/core/align.py` as pure functions, so a number in a
manuscript can be recomputed from the stored alignment without rerunning the aligner.
