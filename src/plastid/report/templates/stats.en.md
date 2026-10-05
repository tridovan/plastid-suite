## Group comparison statistics — {{stats_summary.accession}}

{{stats_summary.n_loci}} loci were measured in {{stats_summary.n_genomes_with_cds}} of the
{{stats_summary.n_genomes}} input genome(s), across
{{stats_summary.n_groups}} group(s) ({{stats_summary.groups}}), once per locus per genome.
{{stats_summary.n_collapsed_cds}} CDS row(s) were collapsed so the IR's second copy is not
counted twice, and {{stats_summary.n_genomes_without_cds}} input(s) contribute no CDS at all
({{stats_summary.genomes_without_cds}}).

Every locus × measure pair was tested with {{stats_summary.test}}:
{{stats_summary.n_tested}} carry a p-value and {{stats_summary.n_not_computable}} are not
computable, each with its own reason in `stats_tests`. With {{stats_summary.correction}} over
the whole family at α = {{stats_summary.alpha}}, {{stats_summary.n_significant}} test(s) fall
below α after correction.

`stats_tests` carries the statistic, the raw p, both corrections and the n behind each test;
`stats_descriptives` carries each measure's distribution per group; `stats_correlations` carries
the association between the measures themselves. A test the data cannot support is reported as
not computable — never as a p-value.
