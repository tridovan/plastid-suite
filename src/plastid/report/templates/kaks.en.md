## Ka/Ks selection pressure — {{kaks_summary.accession}}

{{kaks_summary.n_genes}} CDS of *{{kaks_summary.organism}}* ({{kaks_summary.accession}}) were
considered; {{kaks_summary.n_shared}} are shared with every genome compared. Each shared gene
was aligned on its own CDS with {{kaks_summary.aligner}} and scored by
{{kaks_summary.estimator}} through {{kaks_summary.tool}} {{kaks_summary.tool_version}}.

`kaks_pooled` carries one row per genome — how many genes produced a value, how many did not,
and the pooled Ka, Ks and Ka/Ks. `kaks_genes` carries one row per gene and genome, with the
reason beside any gene that has no value: a gene the other genome does not carry, a gene whose
annotation is not whole codons, an alignment that breaks the reading frame, a pair the tool
itself refused. Nothing is dropped silently, because a ratio over an unstated subset is not a
result.

**Ka/Ks is undefined when Ks is zero** — the cell is blank, never zero and never one. The pooled
row is not the mean of the per-gene ratios: it weights each gene by the number of synonymous and
nonsynonymous sites the estimator itself counted.
