## Population genetics — {{popgen_summary.accession}}

{{popgen_summary.n_genes}} genes analysed across {{popgen_summary.n_genomes}} genomes in
{{popgen_summary.n_groups}} group(s) ({{popgen_summary.group_sizes}}), each gene aligned with
{{popgen_summary.tool}} and summarised per group in `popgen_groups`; the pairwise Fst between
groups is in `popgen_fst`.

Every statistic travels with its **n** and its usable site count, because they are not
interchangeable: haplotype diversity, π and θ(S) are readable from a handful of sequences, while
**Tajima's D and Fu & Li's D\\* and F\\* need at least four** — with three sequences or fewer
their variance is identically zero, so the cell is blank and the Note column says why. A blank is
never a zero.

`popgen_fst` uses Nei's Gst: two groups with the same haplotype composition give zero, two groups
sharing no haplotype give one, and groups with no diversity left to partition give a blank. In
this dataset the groups are **taxa, not populations of one species**, so an Fst here measures
divergence between them — read it with the n beside it, and do not quote it as a within-species
differentiation.
