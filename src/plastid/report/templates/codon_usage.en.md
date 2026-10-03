## Codon usage — {{codon_usage.accession}}

*{{codon_usage.organism}}*: {{codon_usage.n_cds}} CDS in the set, {{codon_usage.n_codons}} codons
counted, translation table {{codon_usage.translation_table}}. The effective number of codons is
**{{codon_usage.enc}}**, the codon adaptation index **{{codon_usage.cai}}**, the codon bias index
**{{codon_usage.cbi}}** and the frequency of optimal codons **{{codon_usage.fop}}**.

CAI · CBI · Fop are measured against {{codon_usage.reference}} — a bias index means nothing
without the set it was measured against. The per-codon RSCU table is `rscu`; the check that every
amino acid's RSCU sums to its family size is `rscu_families`.

Internal stop codons found inside CDS: {{codon_usage.internal_stops}}.
