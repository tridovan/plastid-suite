## Nucleotide diversity — {{diversity_summary.accession}}

π across {{diversity_summary.n_genomes}} aligned genomes, reference
*{{diversity_summary.organism}}* ({{diversity_summary.accession}}), aligned with
{{diversity_summary.tool}}. The alignment is {{diversity_summary.alignment_length}} columns long
and {{diversity_summary.usable_columns}} of those were usable; whole-alignment π is
**{{diversity_summary.pi}}**.

Windows use a {{diversity_summary.window_bp}} bp window and a {{diversity_summary.step_bp}} bp
step over alignment columns — a π curve means nothing without them. Per-window values, their
usable-column counts and each window's plastid region are in `diversity_windows`; per-region and
per-gene values are in `diversity_regions` and `diversity_genes`.

π per intron and per intergenic spacer is in `diversity_introns`
({{diversity_summary.n_introns}} introns, the same spans the `introns` module publishes) and
`diversity_igs` ({{diversity_summary.n_spacers}} spacers — gaps between the exons of neighbouring
called genes inside LSC or SSC). Spacers inside the IR, gaps crossing a region boundary and gene
pairs that leave no gap are named in the notes rather than silently dropped.

Gapped columns are dropped, never scored as matches, and every π is published next to the number
of columns it was computed from.
