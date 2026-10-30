# The PlastidToolSuite pipeline

This is the **product direction and the canonical data flow**. The diagram includes both what
exists and what is planned; it does not mean every box is joined into one command. The comparison
table below states the current state as read from the code. Per-module definitions — unit, input,
tables, control, copy/pseudo/partial policy and the missing rule — are in
[module contracts](module_contracts.md). Reproduced defects are in
[known issues](known_issues.md). Session logs and coding order are not sources that define the
pipeline.

## The canonical diagram

```text
   REFERENCE SOURCES
   NCBI RefSeq  ·  the group's own submitted records
                │
                │  plastid panel refresh/propose/pin
                │  accession · source · date · sha256
                ▼
   ┌──────────────────────────────────────┐
   │  RAW REFERENCE DB  (many .gb files)  │
   └──────────────────┬───────────────────┘
                      ▼
   ┌──────────────────────────────────────┐
   │  QUALITY ASSESSMENT                  │
   │  · which file DIFFERS (REPORT)       │
   │  · which file is unusable (BLOCK)    │
   └──────────────────┬───────────────────┘
                      │  report · fix or drop a bad file, recorded
                      ▼
   ┌──────────────────────────────────────┐
   │  USER READS AND CHOOSES              │
   │  1 .gb file as the naming standard   │
   └──────────────────┬───────────────────┘
                      ▼
                ★ NAMING STANDARD    →    vocabulary A
                      │
                      ▼
   ┌──────────────────────────────────────┐
   │  UNIFY NAMES ACROSS REFERENCE DB     │  gene_normalize mode=normalize
   └──────────────────┬───────────────────┘
                      ▼
        NORMALISED REFERENCE DB (A)
        new artifact · new checksum · traceable to the raw DB and ★
                      │
                      │             reads (Illumina · ONT/HiFi)
                      │                         │
                      │                         ▼
                      │       ┌─────────────────────────────────────┐
                      │       │ [1] ASSEMBLY + QC  [DEFERRED]       │
                      │       │ short: GetOrganelle + SPAdes         │
                      │       │ long: OATK / raven, to validate      │
                      │       │ fastp: QC short reads                │
                      │       │ BWA/minimap2: map-back               │
                      │       │ Bandage: graph; check joins/reads    │
                      │       └─────────────────┬───────────────────┘
                      │                         │ FASTA + QC results
                      │                         │
                      │       assembled FASTA ──┤  current entry point
                      │                         ▼
                      └─────────────►┌───────────────────────────────┐
                                     │ [2] ANNOTATION                │
                                     │ in-house engine               │
                                     │ blastx/blastn · ARAGORN · IR  │
                                     └──────────────┬────────────────┘
                                                    │ GenBank — names per A
   externally annotated GenBank ─────────────────────┤
   check format/feature/organelle; keep provenance    ▼
                                     ┌───────────────────────────────┐
                                     │ [3] CHECK NAMES + ANNOTATION  │
                                     │ gene_normalize mode=verify    │
                                     │ ORF · position · qualifier    │
                                     └──────────────┬────────────────┘
                                                    ▼
                         ╔══════════════════════════════════════════╗
                         ║ CHECKED INPUT SET (vocabulary A)         ║
                         ║ sample ID · hash · frame · sample status ║
                         ╚══════════════════╤═══════════════════════╝
                                            │
                   ┌────────────────────────┼─────────────────────────┐
                   ▼                        ▼                         ▼
        [4] ANALYSIS              [5] RNA EDITING              [6b] NCBI SUBMISSION
        structure · gene · codon  PREP-Cp-compatible · C 0.8    checked GenBank
        intron · SSR (MISA)       published PREP method in PTS  → .fsa + .tbl
        exact repeat              CDS + B′ (7 taxa) → MAFFT     + .sbt template
          repeat-match            → sites · score · coverage   → table2asn
          (Vmatch cross-check)    module rna_editing           → .sqn + validator
        synteny (MUMmer4)
        Ka/Ks (PAML yn00)
        Ts/Tv · π · barcode
        popgen · group statistics
        MAFFT when alignment is needed
                   │                        │
                   └──────────────┬─────────┘
                                  ▼
                    [6a] TABLES · FIGURES · DRAFT REPORT
                    CSV/JSON · Excel · Word · PDF for one genome
                    matplotlib · R · user-supplied tree
                    Supplementary · QC · provenance and checksum
```

**How to read it, and the important boundaries:**

- The database-level reference path feeds annotation; it is not a step that replaces assembly. Name
  normalisation on the DB runs **before** annotation. The ★ standard is the user's choice, not a
  fixed genome of the toolkit.
- A sample receives name A through annotation. The verify step is not a step that silently renames
  a second time: a difference is a finding that must be resolved/recorded before the pipeline
  accepts its input as normalised.
- An external GenBank may need an evidenced normalise step if it uses another vocabulary; do not
  equate that difference with an annotator error. Keep the original, the name mapping and the hash
  of the converted copy.
- The IR gate applies only to the part that actually needs position/IR region. A gate-pass does not
  by itself prove the IR boundaries are right; the current defects are blocking acceptance of this
  part.
- RNA editing (M4.5) is the module `rna_editing`. The original PREP-Cp can no longer be run: its
  web service ended in 2022, and its reference alignments were never published. PTS therefore runs
  the **published PREP method itself**, with a rebuilt reference B′ and cutoff 0.8, and labels
  every site **PREP-Cp-compatible**. It reproduces both archived genuine PREP-Cp outputs exactly,
  32 of the original's 34 *Aegilops* sites, and 31 of 33 known sites of the unseen *N.
  sylvestris*. RNA-seq is not required to run a prediction; the public set is used for
  independent validation. The method, the evidence and the acceptance results are in the
  [specification](tools/prep_cp.md).
- The NCBI submission goes from **checked sequence and annotation**, and must not be built from
  statistics tables. The `.sbt` is supplied by the submitter. Producing a `.sqn` does not mean the
  validator passes or that NCBI accepts the submission.
- The circular map already has a matplotlib renderer (`circular_map`); there is no need to restore
  `pycirclize` just because it is in the old diagram. The main repeat uses `repeat-match`; Vmatch is
  a cross-check, not the default engine.
- This is the target flow. Today many steps must be joined by hand; an end-to-end command must not
  be advertised. `plastid run` runs several modules over the same input set, and does not
  automatically pass a module's GenBank artifact into the next module. Listing modules in one
  command does not prove the pipeline is joined.

## Reference roles that must not be mixed

| Role | Data | Choice/origin to record |
|---|---|---|
| Annotation DB | Many audited and normalised GenBanks | File list, original/converted hash, exclusions, audit verdict |
| Naming reference ★ | One GenBank that decides vocabulary A | Chosen by the user; file/hash and name mapping; does not change with taxonomy on every run |
| Comparison reference | Genome used as the first row/coordinate anchor for an analysis | Today many modules take the first input; state the record and input order, do not assume it is ★ |
| Assembly seed/HMM | GetOrganelle seed DB or OATK HMM | DB version/checksum and program used; do not use in place of the annotation DB |
| Mapping reads | Sequence of the sample under test | Hash of the sequence copy and the corresponding frame; do not map to ★ and then treat it as support for the sample |
| RNA editing reference B′ | Mature proteins of the seven PREP-Cp training taxa, rebuilt by PTS because the original alignments were never published | Pinned records, sourced edit lists and build recipe under one digest in the manifest and scope label; never presented as the original PREP-Cp reference ([specification](tools/prep_cp.md)) |
| Submission template | `.sbt` | Submitter metadata, not a genome reference |
| Tree/outgroup | Tree built by the user and declared taxa | Tree hash, root/drop/order performed; do not infer the outgroup from the taxonomy level |

Mode A/B of **panel proposal** are two ways of selecting the candidate set, different from vocabulary
A/B in the naming step. `panel propose` currently suggests by taxonomy, and does not choose ★ itself
or prove the reference is of good quality.

## Current state and gaps at each step

"Has a component" only confirms that a code path/check exists; it is not acceptance of the whole
pipeline.

| Step | Existing component | Gap to accept |
|---|---|---|
| 0 — reference level | `panel refresh/propose/pin`, `reference verify/loo` | Taxonomy/source cache; external data; checksum and provenance across conversion steps |
| 0b — audit/choose ★ | `reference audit`, `annotate --audit-only`, audit exceptions by digest | Default set audits with 0 blocked (H-21 fixed: partial flags reported, annotation-first IR reported with the sequence path named); the **IR gene-mirror rule** is reported per record (a gene inside one IR copy must have its copy in the other on the opposite strand; 0/34 reference records violate it, 23/80 validation-panel records do — measured 29/10/2026); the ★ choice still enters the artifact/manifest only when the user passes it |
| 0c — normalise DB | `gene_normalize mode=normalize`, blastp for protein ranking, GenBank artifact and pairing tables | Artifact retains sequence/locations/translations/qualifiers outside naming, has its own digest in the manifest, and leaves the pin untouched (`tests/golden/test_reference_flow.py`); audit exceptions are bound to the old file digest and do not carry over |
| 1 — assembly/QC | Only tool survey and declaration of some dependencies | No module yet; real reads, seed/HMM, map-back, graph, junction support and resource budget all need integration |
| 2 — annotate | `plastid annotate` calls `plastid_annotate`, batch one run per sample | The normalised DB is not yet joined automatically; feature/locus/IR length discrepancies and validation outside the family need fixing |
| 3 — verify/ORF | `gene_normalize mode=verify`, `orf` | Verify currently reports would-change, does not write a new GenBank itself; needs stop/accept rules in the workflow and annotation preservation |
| 4 — analyses | The modules listed below | Need a data contract per module and validation before aggregation; do not accept only through IR/gene count |
| 5 — RNA editing | `rna_editing`: the PREP method against reference B′ (`plastid editing-reference build/verify`); both archived *Vitis* outputs reproduced; *Aegilops* KJ614412.1 compared with its three evidence sets; *N. sylvestris* scored as Mower 2009 Table 3 did ([specification](tools/prep_cp.md)) | The maize and rice edit lists come from a comparison table, not their primary surveys; *Pisum* and *Saccharum* (Table 3) are not measured; the module does not yet read a GenBank artifact produced earlier in the same pipeline run |
| 6a — export | CSV/JSON, Word, Excel, `figures`, `package`, `qc`, `genome_page` | One figure pass over the matplotlib and R catalogues for `run --figures`, `figures` and `package`; `qc` inventories both catalogues without drawing. A missing plotting extra or R is a named status, never a silent absence. The package is replaced as a whole, its index links every file with a digest, and the workbook names each table's scope; a batch's multi-sample figures use `plastid figures <batch-dir>` |
| 6b — submission | `submission`: `.fsa/.tbl`, optional `.sqn` via table2asn + `.sbt` | Reads the checked GenBank directly (a bare sequence is refused); `.fsa` and `.tbl` share the sequence id; the validator's `not run` state is distinct from zero errors, and a `.sqn` is not ready-to-submit evidence |

Code sources to compare: [CLI](../src/plastid/cli.py),
[panel](../src/plastid/core/panel.py), [audit](../src/plastid/core/reference_audit.py),
[naming](../src/plastid/modules/gene_normalize.py),
[annotator](../src/plastid/modules/plastid_annotate.py),
[Supplementary](../src/plastid/report/package.py), [QC](../src/plastid/report/qc.py),
[submission](../src/plastid/modules/submission.py).

## Coverage of analysis modules and outputs

| Module/group | Output/definition to check | Must not be inferred from the feature name |
|---|---|---|
| `genome_stats`, `structure` | Length/GC, feature count and locus count, four regions, original/normalised frame | A fixed gene total or one gate-pass does not prove every position is right |
| `gene_table`, `gene_names`, `gene_normalize` | Gene groups, name A, synonym/pairing, information lost/changed on serialise | Identical names do not by themselves prove the same locus; a missing gene feature does not mean a missing CDS |
| `orf`, `codon_usage`, `amino_acid` | Transcript strand/order/frame/table, stop/pseudo/partial, RSCU/ENC/CAI/CBI/Fop and reference | The bias indices' self-reference is not a chosen set of highly expressed genes |
| `introns` | CDS/tRNA/trans-splicing, counts and correct genomic spans | CDS and tRNA introns are counted separately with `kind`; a missing feature class is "not measured"; a tRNA gap the caller cannot classify is named "not assessed" |
| `ssr`, `repeats` | MISA + Perl thresholds, motif/context, each exact-repeat pair, coordinates/strand/class | SSR and exact repeats do not mean a general tandem-repeat module exists |
| `substitutions` | Reference/order, CDS/whole scope, shared loci, IR/gap/ambiguity handling, Ts/Tv and substitution types | The current default whole scope differs from the golden's CDS scope |
| `diversity`, `barcodes` | Alignment columns↔reference↔frame, π/window/step/usable sites, gene/intron/IGS, candidate-calling rule | "Barcode candidate" does not by itself prove species-discriminating performance |
| `kaks` | Codon alignment, frame/stop, yn00 estimator, Ks=0, missing and how they are pooled | Do not treat an uncomputable ratio as 0, and do not treat the first input as naming ★ |
| `synteny` | MUMmer4: nucmer → delta-filter → show-coords, block/identity/coverage and parameter scope | Do not equate every genome comparison with another LCB method just because the tool is in the environment |
| `popgen`, `stats` | Sample groups, n, usable sites, formula/assumptions, unit of observation, FDR; correct strand/frame extraction | The Fst column currently uses Nei's Gst; do not promise equivalence with every estimator. Having a haplotype count is not the same as having a network figure |
| `submission`, export, figures | Artifact matching sample/hash/scope, validator, table↔figure↔report↔package | A file existing does not mean QC passed; a CSV does not prove Excel/Word keep the same numbers |

Full sources: [module registry](../src/plastid/modules/__init__.py),
[figure catalogues](../src/plastid/result/figures/__init__.py),
[platform/tool table](platforms.md). `mito_annotate`/`mito_genes` remain a prototype outside the
current plastid development scope; their existing checks are kept.

## Whole-pipeline acceptance conditions

1. There is an input/output copy for each step, sample ID/hash/frame/schema version and the
   provenance of each conversion; the normalised DB is a derived artifact and does not overwrite the
   raw pin. Know exactly which file the next step reads.
2. Run from both assembled FASTA and annotated GenBank on at least two families; include
   single/multiple sample cases and a failing sample. A failed step must not create a "checked"
   state for downstream; the unsupported count must be non-zero.
3. Every metric has an independent check with the same definition, not just a cell matching a paper.
   In particular check sequence extraction, IR copy policy, alignment mapping and missing/unknown.
4. Complete branches for tables/figures/report/QC/Supplementary and submission; R/matplotlib are
   listed within their correct support scope. Re-render does not re-run analysis, and does not keep
   stale output that would mislead.
5. A fresh wheel/checkout has a usable current-scope environment; a tool installed in the
   environment is not automatically an integrated feature. Tools/DBs of deferred branches must be
   clearly classified. **M11.19 has classified them**: [the tool inventory](tools/README.md) gives
   every tool one role, `environment.yml` is the runtime scope only, and the minimum install
   (biopython + PyYAML) was measured running the pure-Python modules from a freshly built wheel
   ([platforms](platforms.md)).
6. RNA editing (M4.5) is built and checked on its specification's evidence. Assembly and the
   already-deferred branches still have their own DoD.

These are target criteria; the documentation update session has not implemented or run acceptance
for the missing parts.
