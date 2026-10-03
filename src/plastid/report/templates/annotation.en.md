## Annotation — {{annotation_summary.accession}}

The plastome of *{{annotation_summary.organism}}* ({{annotation_summary.accession}},
{{annotation_summary.length_bp}} bp) was annotated by homology and read back with our own
reader — a file that does not re-read is not a result. The gene set that came out is
**{{annotation_summary.n_genes}}** genes: {{annotation_summary.n_cds}} protein-coding,
{{annotation_summary.n_trna}} tRNA and {{annotation_summary.n_rrna}} rRNA, both inverted-repeat
copies included.

What the annotator produced was **not** taken on trust. Every CDS went through the same ORF
checks as a downloaded record: {{annotation_summary.n_orf_failed}} failed, the rest being
classified in `orf_status`. Every name went through the same normalizer, and
{{annotation_summary.n_aliased_names}} genes carry a second name in the table, printed beside
the canonical one. {{annotation_summary.n_unclassified}} genes matched no functional group;
they are listed in their own table rather than folded into a neighbour.

The raw annotator output is kept beside these tables, so every number above can be traced back
to the file it came from. An annotation is a claim about a sequence; this report says which
program made it and on which reference set.
