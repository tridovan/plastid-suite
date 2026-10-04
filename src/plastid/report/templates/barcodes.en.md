## Barcode candidates and hypervariable regions — {{barcodes_summary.accession}}

π along *{{barcodes_summary.organism}}* ({{barcodes_summary.accession}}) in
{{barcodes_summary.window_bp}} bp windows with a {{barcodes_summary.step_bp}} bp step, over
{{barcodes_summary.n_genomes}} aligned genome(s) — the same series the diversity module reports,
computed by the same function.

`barcodes_windows` ranks every window by π; `barcodes_regions` calls the hypervariable regions,
which are runs of at least two consecutive windows above **{{barcodes_summary.threshold_rule}}**
({{barcodes_summary.threshold_pi}}) — the rule is exported because a hypervariable region without
its rule is a shape somebody eyeballed. A window with no usable column has no π, no rank and
breaks a run rather than joining it.

`barcodes_candidates` lists the candidate barcode loci of `data/barcode_loci.yaml` — each with the
source of its claim — against this record's own annotation:
{{barcodes_summary.n_found}} of {{barcodes_summary.n_candidates}} were found. A locus the record
does not carry keeps its row and says why: a locus that is silently missing from a table looks
exactly like a locus nobody looked for.
