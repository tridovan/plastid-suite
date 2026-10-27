# Development handoff

The session of 2026-10-27 implemented the Đợt 3 milestone **M11.18** (outputs, QC and submission
as one checkable product) and pushed it on `feat/m11-ir-rework` as **`4e8401d`**; a test-only CI
fix followed in **`ff96417`**. Read `CLAUDE.md`, then `docs/pipeline.md`; in the maintainer
workspace `.dev/TASKS.md` carries the progress notes and
`.dev/records/M11.18_family_crosscheck_2026-10-27.md` the batch evidence.

## What changed

- **One figure catalogue for four commands.** `run --figures`, `figures` and `package` use one
  pass over the matplotlib and R catalogues; `qc` inventories both without drawing. Each figure
  gets a `FigureRecord` (renderer, files, provenance, reason) and a caption file beside the image.
  A missing plotting extra is still exit 3 for `run --figures`; a missing R is always a recorded
  skip; `package` keeps its tables and records either absence in `index.md`.
- **Refresh and staleness.** A re-render deletes the figure files it no longer draws; `figures`
  and `package` drop manifest entries for deleted files and update digests of rewritten ones.
  `package` builds in a staging folder and swaps it in, so a failed rebuild leaves the previous
  package untouched and nothing stale survives a successful one.
- **Package/QC/parity.** `index.md` links every file with its sha256, cites the registry, and
  names every output not produced; the Excel workbook has a `plastid.scope` sheet; the package
  copies the stored draft reports. `tests/unit/test_output_contract.py` checks CSV/JSON/Excel/
  Word/Markdown parity against the same formatted cells, the package manifest and links, the
  failed-rebuild path, and the QC catalogue.
- **Submission.** A record with no gene/CDS/tRNA/rRNA features is refused (a checked GenBank is
  the input, not a bare sequence or an analysis table). `validator_status` distinguishes
  `not requested` / `not run` from a validator that ran; its counts are blank rather than 0; a
  `.sqn` is described as a conversion, not a validation. `.fsa` and `.tbl` are tested to carry the
  same sequence id.
- **Third-party panel provenance.** The batch aggregate repeats `sample_accession`,
  `sample_data_label`, `sample_ir_source` and `sample_ir_gate` on every row (the figure pivot
  skips these metadata columns), and a batch R caption names the data label, accession set and
  IR status counts, pointing at `samples.csv` for the per-sample mapping.
- `docs/known_issues.md` gained a fixed G-03 entry; `docs/pipeline.md` steps 6a/6b state the
  current behaviour and tool policy.

## Checks actually run

- Full fast suite (`make test`, `PTS_REQUIRE_TOOLS=1`, tool env): **1085 passed, 1 failed,
  1 skipped, 66 tier-B deselected** (479 s). The one failure was the K-06 test calling
  `_render_r_figures` with its old positional signature; only that test call was changed, and
  `tests/unit/test_known_issues.py` then **27 passed**. No product code changed after the full
  run, so the full suite was not repeated (user instruction, 27/10/2026).
- Targeted suite after the last product changes: **123 passed, 1 skipped** (the intentional
  machine-without-table2asn case).
- Golden figure path with real R: `tests/golden/test_e2e.py -k "figure or extra"` **2 passed**.
- Private panel cross-check: Rosaceae, 8 records × 5 modules, batch vs single, **31,808 cells,
  0 mismatches** — one family; the other nine and the CGAS/paper definition comparison remain
  open (see `.dev/records/M11.18_family_crosscheck_2026-10-27.md`).
- Not run: the local full golden (`make golden`/`make check`). By the batch plan it belongs to
  Đợt 4, after M11.19/M11.9/M11.8. CI on `ff96417` is green: lint, fast tests, tier-A and tier-B
  golden (run 40529137833).

## Next

Đợt 4: M11.19 → M11.9 → M11.8 (+M9.6), then the single full golden. Still open from Đợt 2: the
batch-vs-single cross-check on the other nine families and the CGAS/paper definition comparison
before any cell-for-cell claim. Do not treat this as a release acceptance.
