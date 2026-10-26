# API v0.1 — frozen contracts

Frozen 2026-10-02 (end of Phase 1's framework work). Changing anything here is a breaking change:
bump the version and update this file first.

## 1. Module contract

```python
@module(
    "ssr", version="0.1.0", title="SSR (microsatellites)", requires_ir=True, citations=("misa",)
)
def run(genome: Genome, params: dict[str, Any], ctx: RunContext) -> Result: ...
```

- Pure: no cwd, no globals, no `sys.exit` — the function receives everything it needs.
- `ctx.out(logical, ext)` returns the path to write; modules never name directories.
- `ctx.require_ir()` must be called before any position-based work; it raises `IRGateError`.
- Every external program goes through `ctx.tool(name, argv, …)`; nothing else calls `subprocess`.
- A module never imports another module. Adding one = a new file + one line in `modules/__init__.py`.

## 2. Result schema

`Scope(label, coord_system, include, exclude_regions)` · `Column(name, label, dtype, unit, precision)`
· `Table(id, title, columns, rows, scope, key=None)` · `Result(module, module_version, tables, notes,
citations)`.

- References, literally: `table.column` (single-row table) and `table.rowkey.column`
  (multi-row, `key` set) — e.g. `regions.LSC.length_bp`, `ssr_by_unit.mono.n`.
- An unresolvable reference raises `ResultRefError` — it never prints blank.
- Names and ids may not contain a dot. Multi-row tables without a `key` cannot be addressed.

## 3. Coordinates

Internal coordinates are **0-based half-open**; GenBank and every export are **1-based inclusive**.
Per-genome tables use the GenBank frame; cross-genome comparisons use the normalized frame
(LSC → IRa → SSC → IRb). Every table states its frame in `Scope.coord_system`; the rotation offset
lives in `manifest.json → coordinate_system.rotation_offset`.

## 4. IR gate

`detect_ir(genome)` → annotation first, de-novo (mirror-diagonal) otherwise; both are cross-checked
and any disagreement is recorded in `IRResult.denovo_length`. `enforce_ir_gate(ir)` raises
`IRGateError` below 99.9 % identity. Position-based modules must not run without a passing gate.

## 5. Run layout and manifest

```
<out>/<dataset_slug>/<run_id>/          run_id = YYYYMMDD-HHMMSS-<hash8>, never reused
  00_run_info/run_info.md               the Materials & Methods page
  00_run_info/logs/                     stdout/stderr of every tool call
  NN_<module>/tables/*.csv              one file per table, provenance comments on top
  NN_<module>/result.json               raw values + scopes
  NN_<module>/report.md                 draft report (banner appended by the engine)
  manifest.json
```

Manifest keys (schema_version 1): `run_id` · `command` · `package_version` · `started_at` ·
`finished_at` · `status` · `params` · `inputs[{path,sha256,bytes}]` · `ir{…}` ·
`coordinate_system{name,rotation_offset,note}` · `modules[{id,version,status,key,tools,files}]` ·
`files[{path,module,sha256,bytes}]`. `validate()` is hand-written; no schema library.

## 6. Report references

Templates live in `report/templates/<module>.en.md` and reference cells (`{{…}}`) only — a scan test
rejects literal digits. The engine appends the draft banner and resolves citations from
`src/plastid/data/citations.yaml` (an unknown key raises `CitationError`).

## 7. Exit codes

`0` success · `1` general error · `2` usage error or failed gate · `3` missing or failed tool.
