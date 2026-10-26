# Development handoff

Updated: 2026-10-26. Branch at review: `feat/phase-11`; reviewed HEAD: `b005845`.
This is the sole current development handoff, not a user guide or a record of scientific decisions.
There is no second HANDOFF in the maintainer workspace; earlier versions remain in Git history.

## Current work

The latest session changed documentation only: direction, acceptance criteria and the next task.
The coverage follow-up added the full canonical pipeline and explicit contracts for stage data,
each analysis, output packaging/submission and the runtime environment. Assembly is deferred.
Plastid RNA editing (M4.5) is open and unimplemented. **The user chose original PREP-Cp**. The
integration uses cutoff 0.8, as in the selected benchmark; PTS only prepares CDS inputs and
imports/maps original results. Do not reimplement the
predictor or silently substitute another engine. See [the selected integration](docs/tools/prep_cp.md).
The Aegilops KJ614412.1/SRX209402–SRX209425 benchmark is selected but not pinned/downloaded yet.
RNA-seq is benchmark evidence, not required prediction input. No official local CLI/API has been
verified and the service could not be accessed by the review browser; do not invent commands or
mark execution accepted. Method selection is complete; prepare/import/engine execution remain open.
Source, tests, configuration, dependencies and data were not changed in this session. The working
tree already contained a large cleanup from an earlier session; review staged and unstaged
changes together before committing. HANDOFF is retained; its previously staged deletion has
been removed from the index. No changes were committed or pushed in this session.

The active scope is plastid analysis from assembled genomes. Mito stays a tested prototype.
Assembly, new mito features and a web UI are deferred.

## Next coding session

Read in this order:

1. [CLAUDE.md](CLAUDE.md) for code rules, then [pipeline](docs/pipeline.md) for the product contract.
2. If this maintainer workspace has `.dev/`, open `.dev/TASKS.md` for the active task, then
   `.dev/ROADMAP.md` for that task's dependencies and Definition of Done. These private files are
   not required in other clones; use the user's requested task and public docs there.
3. Read only the design notes named by the active task. The detailed pipeline has one canonical
   copy in `docs/pipeline.md`; roadmap and historical diagrams do not redefine it.

The coding sequence remains:

1. Reproduce [known issues](docs/known_issues.md). Start with cache correctness and annotation
   preservation during rotation, with failing regression tests before fixes.
2. Restore product documentation checks lost during the earlier cleanup: platform-table coverage
   of registered modules and README/package version consistency.
3. Repair circular IR representation, then detection and identity measurement. The junction
   reader is also defective and must not be used as an independent oracle.
4. Finish reference/annotation integration and the stage data contracts before batch summaries.
   Also cover per-analysis correctness, output packaging/submission and the tool environment;
   the latest coverage review adds those explicitly to the local roadmap.
5. For M4.5 follow the selected PREP-Cp prepare/external-run/import contracts and per-step DoD in
   the roadmap; do not restart algorithm selection. This decision does not change the first cache task.

A task is complete when expected output has an independent basis, its regression passes, applicable
checks actually ran, and remaining skips/limitations are recorded. Do not close a defect while its
regression remains xfail. A changed golden baseline needs an explanation and user review;
preserving a known-wrong value is not an acceptance criterion.

## Verification and remaining work

The review ran 74 related tests successfully, then reproduced further defects using separate
probes. Those passing tests do not establish that defects are fixed. No full tool-equipped
acceptance run was performed for this documentation-only session.

All runtime issues in [known issues](docs/known_issues.md) remain open. Truncated IRs can pass the
gate; cache keys omit annotation fields; rotation drops partial flags; tree options are not
forwarded by the CLI. Keep these limitations visible when interpreting results.
