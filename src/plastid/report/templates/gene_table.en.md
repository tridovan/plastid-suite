## Publication gene table — {{gene_table_summary.accession}}

The plastome of *{{gene_table_summary.organism}}* ({{gene_table_summary.accession}}) carries
**{{gene_table_summary.n_genes_total}} genes — {{gene_reconciliation.CDS.in_table}} protein-coding,
{{gene_reconciliation.tRNA.in_table}} tRNA and {{gene_reconciliation.rRNA.in_table}} rRNA**,
across {{gene_table_summary.n_groups}} functional groups and
{{gene_table_summary.n_unique_names}} distinct gene names. A gene carried in both copies of the
inverted repeat is listed once and marked with its copy number.

The grouping is checked rather than assumed. The table's own totals are reconciled against the
genome-wide count in `gene_reconciliation`: protein-coding
{{gene_reconciliation.CDS.agreement}}, tRNA {{gene_reconciliation.tRNA.agreement}}, rRNA
{{gene_reconciliation.rRNA.agreement}}, grand total {{gene_reconciliation.Total.in_table}}.

{{gene_table_summary.n_unclassified}} genes matched no group in the grouping file. A gene that
matches nothing is reported as unclassified — it is never folded into a neighbouring group, and
never dropped. Per-gene coordinates, both inverted-repeat copies included, come from the
`genome_stats` module's `genes` table.
