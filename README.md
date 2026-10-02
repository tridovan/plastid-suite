# PlastidToolSuite

Toolkit for **plastid (chloroplast) and mitochondrial genome analysis** — academic, non-commercial.

**Status: Phase 1 (framework) — no analysis modules yet.**
Progress: `ROADMAP.md` · current work: `TASKS.md` · working rules: `CLAUDE.md`.

## Develop

```bash
uv sync        # create .venv and install (Python >= 3.10)
make test      # fast tests, no external tools
make golden    # golden-number set; prints a skip reason for every missing tool
make check     # lint + test + golden
```

The minimum install (pip only) runs the pure-Python analyses; external tools (MISA, MAFFT, …)
are only needed for the tier-B checks and will be reported as missing, never silently skipped.
