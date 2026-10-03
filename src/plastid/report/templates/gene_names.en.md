## Gene-name check

Across {{panel_summary.n_records}} records and {{panel_summary.n_genes}} counted genes,
{{panel_summary.n_distinct_names}} distinct names normalize to
{{panel_summary.n_distinct_canonical}} distinct canonical names. Per-record counts are in
`name_summary`. Nothing is renamed in place: the name a record arrived with stays beside the
form it was normalized to in `name_mappings`, and stays in the record itself.

Of those names, {{panel_summary.n_renamed}} were written in a different syntax — a bracket
anticodon, a lowercase prefix, a size suffix on an rRNA — and {{panel_summary.n_aliased}} are
one locus under two nomenclatures. An aliased name is printed with both of its names, so a
reader who knows either name finds the gene, and nobody has to know which dialect a record
happened to use.

{{panel_summary.n_resolved_by_sequence}} bare tRNA name(s) had no anticodon and were resolved
against the named copies elsewhere in this run — by full sequence identity or not at all. A
single differing base is not a match: these tRNAs differ between species, and a near match
would be a different gene or a different sample. Reading the anticodon at a fixed offset
matches only some of the named tRNAs in this dataset, so it is not used to resolve one.

{{panel_summary.n_unresolved}} name(s) could not be resolved. They are listed in
`unresolved_names` as they arrived, with the product the record gives and the reason they
stayed unresolved — reported rather than guessed: a bare tRNA name that matches no reference,
or matches two equally, is not a name this toolkit gets to invent.
