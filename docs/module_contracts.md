# Module contracts

One block per analysis module: what a row is, what it is computed from, what an independent check
can be, and how the copies, pseudogenes, partial features, naming and missing values are treated.
`docs/pipeline.md` names the modules; this file is the definition a reader can hold a result
against. The registry is the list of modules — `tests/unit/test_module_contracts.py` fails if a
registered module has no block here or a block names no registered module.

Terms used in every block:

- **Copies.** Both inverted-repeat copies are features; a *gene* is a name. A module that counts
  features (gene counts, ORF checks) counts both copies; a module that compares genes counts the
  name once.
- **Pseudo.** A feature the record declares `/pseudo`. Gene-set comparisons exclude pseudogenes;
  a module that checks every annotated CDS keeps them and marks them.
- **Partial.** A feature the record declares partial (`<` / `>` or `/partial`). A partial CDS's
  measurement is published and marked, never silently completed; `phylogeny` excludes it unless
  `include_partial=yes`.
- **Missing.** A value that was not measured is blank with a reason, never `0`. A measured zero is
  written as `0`.
- **Query frame.** Only the reference is normalized. Queries are aligned as deposited, so a query
  written from another origin — or on the other strand from another start — is aligned with end
  gaps; the measurement is still taken over the reference's columns and states its own frame.
- **IR gate.** Modules marked *IR gate: yes* refuse to run until the two copies pass
  `core/ir.enforce_ir_gate` (≥ 99.9 % identity and a minimum known-column count). The gate is never
  bypassed; `--set ir_source=sequence` chooses which copy pair is measured, not whether the gate
  runs.
- **IR scoring policy.** Identity is measured under one fixed policy and its version is recorded
  with every result (match +1, mismatch −1, gap open −2, extend −1). One panel record,
  `NC_035154.1`, differs from a second aligner only because that aligner prices gaps differently:
  the 14 differences are SNPs clustered in ~140 bp, so this toolkit reads 99.9462 % (pass) where a
  cheap-gap alignment inserts 14 gaps and reports 99.8923 % (fail). The case is recorded as a
  **scoring-policy case** — not a pass that was rounded up, and not a defect. The threshold and
  scorer are unchanged.

## Structure and counts

### `structure` — region boundaries, %GC and GC-skew (arity 1, IR gate: yes)

- **Unit.** One row per region (LSC, IRa, SSC, IRb) in the normalized frame; the summary row is
  the genome.
- **Input.** One genome and the IR pair from the gate.
- **Tables.** `structure_summary`, `regions`, `gc_skew_windows`.
- **Control.** The four region lengths tile the chromosome exactly; IRa and IRb have the same
  length because they are the gated pair; GC and skew are computed from the sequence. A rotation
  of the record changes no region value (rotation cases).
- **Copies / pseudo / partial.** Both IR copies are regions; genes do not enter.
- **Naming / comparison reference.** none; the frame is normalized (LSC → IRa → SSC → IRb) and
  `rotation_offset` maps back to GenBank.
- **Missing rule.** A skew window with no bases is absent, not 0.

### `genome_stats` — genome, region and gene statistics (arity 1, IR gate: yes)

- **Unit.** One genome row, one row per region, one row per gene name.
- **Input.** One genome and the gate.
- **Tables.** `genome_stats`, `regions`, `genes`.
- **Control.** Length is `len(sequence)`; GC is counted from the sequence; region lengths come
  from the same gate `structure` uses; gene counts follow `count_genes()`.
- **Copies / pseudo / partial.** Both copies counted; pseudogenes are counted features and marked;
  partial features are counted and marked.
- **Naming / comparison reference.** Region frame normalized; gene names as the record prints
  them.
- **Missing rule.** No value is inferred for an ambiguity the record does not resolve.

### `gene_table` — genes by functional group (arity 1, IR gate: no)

- **Unit.** One row per functional group, plus a `Total` row and a reconciliation table.
- **Input.** One genome's annotation.
- **Tables.** `gene_table_summary`, `gene_groups`, `gene_reconciliation`.
- **Control.** The total is CDS + tRNA + rRNA from `count_genes()`; the group copies add up to the
  total; every kind agrees with the count. A gene no group claims is counted in
  `n_unclassified` and named in the notes.
- **Copies / pseudo / partial.** Both IR copies counted; pseudogenes counted and marked;
  partial features counted.
- **Naming / comparison reference.** The record's own names; canonicalisation only for grouping.
- **Missing rule.** Annotated genes are either grouped or named as unclassified — never dropped.

### `orf` — open-reading-frame check per CDS (arity 1, IR gate: no)

- **Unit.** One annotated CDS feature.
- **Input.** One genome's CDS features.
- **Tables.** `orf_summary`, `orf_status`, `orf_reconciliation`, `annotated_not_coding`,
  `orf_checks`.
- **Control.** Start codon, stop codon, length and internal stops are read from the CDS assembled
  from its own exon intervals; two annotations that assemble the same intervals get the same
  check.
- **Copies / pseudo / partial.** Both copies checked; pseudogenes are checked, marked
  `pseudogene` and excluded from the failure count; a partial CDS is allowed and marked.
- **Naming / comparison reference.** none beyond the record.
- **Missing rule.** Every CDS keeps a row; a status and a detail say why a check did not pass.

## Naming

### `gene_names` — names across a panel (arity 2, IR gate: no)

- **Unit.** One changed or unresolved name; the summary counts the rest.
- **Input.** Two or more records; all provide the sequence evidence for bare tRNA names.
- **Tables.** `panel_summary`, `name_summary`, `name_mappings`, `unresolved_names`.
- **Control.** The canonical table (`data/gene_aliases.yaml`) plus 100 % sequence identity for a
  bare tRNA name against the named copies in the other records; the `rule` column says which test
  decided.
- **Copies / pseudo / partial.** Gene names, one row per name; pseudogenes excluded from the
  comparison; partial features do not change a name.
- **Naming / comparison reference.** The panel itself; names are never guessed — an unresolvable
  name is reported.
- **Missing rule.** A name that needed nothing is counted, not repeated; an unresolved name is a
  row.

### `gene_normalize` — one vocabulary out of many (arity 2, IR gate: no)

- **Unit.** One gene pairing between the naming reference (vocabulary A, `genomes[0]`) and one
  record (vocabulary B).
- **Input.** The chosen naming record first, then the records to rename.
- **Tables.** `normalize_summary`, `pairing`, `differences`, `conflicts`.
- **Control.** Four tests in order — name, sequence ranking, margin, neighbours in the same
  region; a tRNA is never renamed on a ranking. `mode=verify` expects no rename.
- **Copies / pseudo / partial.** One row per gene name; pseudogenes excluded; A's spelling is the
  join key, B's appears only in `pairing`.
- **Naming / comparison reference.** A is the reference; B is compared to it. A and B are never
  mixed.
- **Missing rule.** Every gene has a status: `measured`, `not measured`, `outside vocabulary`.

## Gene-based measurements

### `codon_usage` — codon usage, RSCU, ENC, CAI/CBI/Fop (arity 1, IR gate: no)

- **Unit.** One codon over the counted CDS set; RSCU per codon, ENC per gene set.
- **Input.** One genome's CDS features.
- **Tables.** `codon_usage`, `rscu`, `rscu_families`.
- **Control.** Sense-codon counts from each CDS's own coding sequence; `n_codons` equals their
  sum; internal stops are counted separately.
- **Copies / pseudo / partial.** Both copies contribute; pseudogenes excluded; partial CDS
  included and counted (`n_partial_cds`).
- **Naming / comparison reference.** CAI, CBI and Fop are **self-referenced** — the analysed CDS
  set is its own reference; the table says so.
- **Missing rule.** Stops, partial trailing codons and ambiguous bases do not enter a count and
  are not written as zero-valued observations.

### `amino_acid` — composition, GRAVY, aromaticity (arity 1, IR gate: no)

- **Unit.** One residue count over the counted CDS set.
- **Input.** One genome's CDS features.
- **Tables.** `amino_acid`, `amino_acid_composition`.
- **Control.** The translation of each CDS's own coding sequence; `n_residues` plus
  `ignored_symbols` equals the translation length.
- **Copies / pseudo / partial.** As `codon_usage`.
- **Naming / comparison reference.** The record's declared translation table.
- **Missing rule.** Stops and ambiguous symbols are counted in `ignored_symbols`, never as
  residues.

### `kaks` — Ka/Ks per shared CDS (arity 2, IR gate: no)

- **Unit.** One gene × one compared genome.
- **Input.** Two or more genomes; `genomes[0]` is the reference.
- **Tables.** `kaks_summary`, `kaks_pooled`, `kaks_genes`.
- **Control.** MAFFT aligns the two CDS; the terminal stop is trimmed; yn00 (Yang & Nielsen
  2000) scores the codon alignment; pooled values are weighted by yn00's site counts.
- **Copies / pseudo / partial.** IR second copy counted once; pseudogenes excluded; partial CDS
  scored and marked.
- **Naming / comparison reference.** Canonical CDS names of the reference record; the comparison
  is reference vs each query.
- **Missing rule.** Every gene keeps a row. `Ks = 0` leaves `Ka/Ks` blank (undefined, never 0 or
  1); every other omission carries a status and a detail.

### `introns` — introns, junctions and what was not assessed (arity 1, IR gate: yes)

- **Unit.** One intron, one trans-splicing junction, or one not-assessed gap.
- **Input.** One genome's multi-part features in the normalized frame.
- **Tables.** `intron_summary`, `introns`.
- **Control.** A gap between consecutive parts of one feature inside one region; CDS, tRNA and
  trans-splicing classes are counted separately.
- **Copies / pseudo / partial.** Both IR copies counted; pseudogenes excluded; a gap the caller
  cannot assess is a named `not assessed` row, not a count of 0.
- **Naming / comparison reference.** Gene names as annotated; region frame normalized.
- **Missing rule.** A record with no feature of a kind says `not measured` for that class.

## Position-based measurements

### `substitutions` — Ts/Tv and substitution counts (arity 2, IR gate: yes)

- **Unit.** A substitution count by type, region, genome and gene.
- **Input.** Two or more genomes aligned with MAFFT; `genomes[0]` is the reference.
- **Tables.** `substitutions_summary`, `substitutions`, `substitutions_by_region`,
  `substitutions_by_genome`, `substitutions_by_gene`, `substitutions_by_type`.
- **Control.** Pairwise aligned columns; row 0 is the reference in its normalized frame; the four
  classes (Ts, Tv, gaps, ambiguous) are counted and reported.
- **Copies / pseudo / partial.** `scope=whole` counts every aligned column; `scope=cds` counts
  shared canonical CDS, one copy per gene; pseudogenes excluded from CDS scope; partial CDS
  marked.
- **Naming / comparison reference.** Reference `genomes[0]`, normalized frame; CDS names
  canonical.
- **Missing rule.** Unaligned and ambiguous columns are dropped and counted, not scored; a
  comparison with nothing aligned fails by name.

### `diversity` — π along the genome and per locus (arity 2, IR gate: yes)

- **Unit.** One sliding window, one region, one gene, one intron, one intergenic spacer.
- **Input.** Two or more genomes aligned with MAFFT; `genomes[0]` is the reference.
- **Tables.** `diversity_summary`, `diversity_windows`, `diversity_regions`, `diversity_genes`,
  `diversity_introns`, `diversity_igs`.
- **Control.** π per site over usable alignment columns; reference positions map columns to the
  normalized frame; gene/intron/spacer rows are derived from `core.spans`.
- **Copies / pseudo / partial.** One copy per gene name; pseudogenes excluded; partial CDS marked.
- **Naming / comparison reference.** Reference `genomes[0]`, normalized frame.
- **Missing rule.** A window with no aligned pair has no π (blank), never 0.

### `barcodes` — hypervariable windows and candidate loci (arity 2, IR gate: yes)

- **Unit.** One window, one hypervariable region, one candidate locus.
- **Input.** Two or more genomes aligned with MAFFT; `genomes[0]` is the reference.
- **Tables.** `barcodes_summary`, `barcodes_windows`, `barcodes_regions`, `barcodes_candidates`.
- **Control.** π comes from the same alignment as the windows; a candidate's span is the
  annotation's own feature or spacer length; the threshold rule is printed beside the threshold.
- **Copies / pseudo / partial.** One copy per locus; pseudogenes excluded from candidates; a
  partial locus is not proposed as a candidate.
- **Naming / comparison reference.** Reference `genomes[0]`, normalized frame.
- **Missing rule.** A candidate the annotation does not carry is a row with `found = no` and the
  missing names; it is a candidate, not a proven diagnostic locus.

### `synteny` — blocks, coverage and % identity along the reference (arity 2, IR gate: yes)

- **Unit.** One alignment block, one genome summary, one window.
- **Input.** Two or more genomes; `genomes[0]` is the normalized reference; queries as deposited.
- **Tables.** `synteny_summary`, `synteny_genomes`, `synteny_blocks`, `synteny_identity`.
- **Control.** `nucmer --maxmatch → delta-filter → show-coords`; the delta file and show-coords
  must agree or the adapter refuses; coverage is aligned bp over the reference length; the series
  is read from the same delta file.
- **Copies / pseudo / partial.** Position-based; both IR copies are positions; annotation does not
  select columns.
- **Naming / comparison reference.** Reference normalized frame for target coordinates; query
  coordinates are the query's own.
- **Missing rule.** A window with no aligned base pair is blank; %IDY counts indels as errors,
  the series counts base pairs only — both are labelled.

### `ssr` — microsatellites (arity 1, IR gate: yes)

- **Unit.** One SSR, with its region and gene context.
- **Input.** One genome; MISA 2.1 (`--set misa=` when not on PATH).
- **Tables.** `ssr_summary`, `ssr_by_unit`, `ssrs`, `ssr_by_region`.
- **Control.** MISA's own repeat table; unit thresholds come from `data/misa.ini` and are printed
  in the run.
- **Copies / pseudo / partial.** Both IR copies are sequence positions and both are reported;
  annotation only supplies context.
- **Naming / comparison reference.** Region frame normalized; gene names as annotated.
- **Missing rule.** A run with no SSR reports a measured zero count with its thresholds, not an
  empty table.

### `repeats` — maximal exact matches (arity 1, IR gate: yes)

- **Unit.** One exact match, with both copies' coordinates and regions.
- **Input.** One genome; MUMmer `repeat-match` at or above the length floor.
- **Tables.** `repeat_summary`, `repeats`.
- **Control.** `repeat-match` (Vmatch 2.3.1 is the cross-check engine, not required); the two
  classes are the tool's own: direct and palindromic.
- **Copies / pseudo / partial.** Position-based; annotation does not enter.
- **Naming / comparison reference.** Region frame normalized.
- **Missing rule.** Matches below the floor are not measured; the floor is published. A tandem
  array is not claimed by this module.

## Comparison and aggregation

### `popgen` — π, Tajima's D, Fu & Li, Fst per gene and group (arity 2, IR gate: no)

- **Unit.** One gene × one sample group; one Fst row per gene and group pair.
- **Input.** Two or more genomes and a group label for every one of them.
- **Tables.** `popgen_summary`, `popgen_groups`, `popgen_fst`.
- **Control.** The formulas in `core.popgen`; columns usable under complete deletion; every
  statistic carries its `n` and site count; Fst is **Nei's Gst**, not every estimator called Fst.
- **Copies / pseudo / partial.** IR second copy counted once; pseudogenes excluded; partial CDS
  marked; a gene carried by fewer than two genomes is named in the notes.
- **Naming / comparison reference.** Canonical CDS names; the groups are taxa in this dataset, so
  Fst reads as divergence.
- **Missing rule.** Tajima's D and Fu & Li's statistics are blank at n ≤ 3 and at zero segregating
  sites; π is blank at n < 2 or with no usable column; every blank carries its reason.

### `stats` — group comparisons across genomes (arity 2, IR gate: no)

- **Unit.** One locus × one measure × one group; one test row per locus and measure.
- **Input.** Two or more genomes and their group labels.
- **Tables.** `stats_summary`, `stats_descriptives`, `stats_tests`, `stats_correlations`.
- **Control.** One observation per canonical CDS locus per genome (IR copies collapsed); measures
  are `length_bp`, `gc_pct`, `gc3_pct`, `enc`; Mann–Whitney U for two groups, Kruskal–Wallis for
  more, with Benjamini–Hochberg correction across the run's tests.
- **Copies / pseudo / partial.** IR second copy collapsed; pseudogenes excluded; partial CDS
  counted and flagged.
- **Naming / comparison reference.** Canonical CDS names; groups as declared.
- **Missing rule.** A test that cannot be computed carries its status and reason; descriptives
  carry `n`.

### `phylogeny` — supermatrix and one ML tree (arity 2, IR gate: no for `matrix=cds`; yes for `introns`/`igs`/`all`)

- **Unit.** One locus block; one taxon row; one tree per run.
- **Input.** At least four genomes (IQ-TREE cannot bootstrap fewer).
- **Tables.** `phylogeny_summary`, `matrix_genes`, `matrix_samples`.
- **Control.** Per-gene MAFFT (or MACSE) alignment concatenated in name order; IQ-TREE with a
  fixed seed and `-bb`; support is recorded; the module does not root the tree or choose an
  outgroup.
- **Copies / pseudo / partial.** One copy per gene name; pseudogenes excluded; partial CDS
  excluded unless `include_partial=yes`; a locus carried by fewer than two genomes is excluded and
  named.
- **Naming / comparison reference.** Canonical gene/locus names; taxa are the input accessions.
- **Missing rule.** A genome without a locus gets gaps for that block, and `matrix_genes` reports
  per-locus coverage; the topology is a construction, not a claim of correctness.

## Production and submission

### `plastid_annotate` — annotate a plastome from a reference set (arity 1, IR gate: yes)

- **Unit.** One called locus; one competing reference; one rejected or overlapping call.
- **Input.** A plastome (FASTA or unannotated GenBank) and a reference set.
- **Tables.** `plastid_annotation_summary`, `called_genes`, `not_called`,
  `competing_references`, `rejected_hits`, `overlapping_calls`.
- **Control.** BLAST+/ARAGORN hits scored by the rules in `core/caller.py`; one call per locus;
  overlaps are in the contract and are reported, not resolved by dropping the weaker call.
- **Copies / pseudo / partial.** Genes are named by symbol; no product string is invented; a
  feature is not merged with the input's annotation.
- **Naming / comparison reference.** Vocabulary A from the chosen reference; the naming reference
  is recorded on the artifact.
- **Missing rule.** A locus with no passing hit is `not_called` by name; a search with no hits
  fails the run.

### `submission` — FASTA, five-column table, validator findings (arity 1, IR gate: no)

- **Unit.** One submission problem; one validator finding.
- **Input.** A checked GenBank record.
- **Tables.** `submission_summary`, `submission_problems`, `validator_findings`.
- **Control.** First-pass checks in `core.submission`; `table2asn` validation only when `sqn=true`
  and the tool and `.sbt` template are present.
- **Copies / pseudo / partial.** Both copies written as annotated; a feature past the sequence end
  is a finding, never repaired.
- **Naming / comparison reference.** The checked record's own names.
- **Missing rule.** An unchecked validator is named as unchecked; it is never reported as `ok`.

## Not covered

`mito_annotate` and `mito_genes` are the frozen mitochondrial prototype. They are kept passing
their tests and are deliberately not extended; their contracts are the code as it stands.
