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

Gapped columns are dropped, never scored as matches, and every π is published next to the number
of columns it was computed from.
