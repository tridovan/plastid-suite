## Introns — {{intron_summary.accession}}

*{{intron_summary.organism}}*: {{intron_summary.n_genes_with_introns}} genes carry
{{intron_summary.n_introns}} introns, {{intron_summary.total_length_bp}} bp in total, the longest
being {{intron_summary.longest_bp}} bp. Introns are the exon gaps read from the CDS parts; phase
is the number of bases of the preceding exon(s) modulo three, in transcript order.

{{intron_summary.n_trans_splicing}} junction(s) are **not** introns: their parts sit in different
plastid regions, so they join separately transcribed RNAs rather than being spliced out of one.
The record's own trans-splicing qualifier and this geometry cross-check as:
{{intron_summary.annotation_check}}.

Every intron is listed with its coordinates, length and phase in `introns`.
