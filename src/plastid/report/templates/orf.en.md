## ORF checks — {{orf_summary.accession}}

{{orf_summary.n_cds}} protein-coding sequences of *{{orf_summary.organism}}*
({{orf_summary.accession}}) were checked against the translation table the record itself
declares (table {{orf_summary.translation_tables}}). {{orf_summary.n_failed}} failed:
{{orf_status.not_multiple_of_three.n}} with a length that is not a multiple of three,
{{orf_status.no_start.n}} without an initiation codon, {{orf_status.no_stop.n}} without a stop
codon, and {{orf_status.internal_stop.n}} carrying an in-frame stop codon before the end. A
sequence that disagrees with its annotation is reported as it is — none is rewritten to pass.

Starts are checked against the declared table, not against ATG alone:
{{orf_status.alternative_start.n}} CDS start on an alternative initiation codon. In a plastome
this is expected rather than exceptional — ribosomal-protein and NADH-dehydrogenase genes are
the recurring cases — so each such CDS is flagged with its own codon and start class in
`orf_checks`, and counted as a failure only when the codon cannot initiate at all. The
`start_class` column separates a codon the declared table lists from one it does not: a codon
outside the table is still reported as a start, never silently accepted or excused.

{{orf_summary.n_pseudogene}} CDS are annotated as pseudogenes. They are reported as their own
status and kept out of the failure counts; a broken reading frame in a pseudogene is the
annotation's own statement, not a defect of the record.

{{orf_summary.n_annotated_not_coding}} features are annotated as genes without a CDS of their
own under them. They are listed in `annotated_not_coding` and are not counted as failures: the
annotation says "gene here, not coding" — a truncated or lost copy — and the check's job is to
repeat that, not to repair it.

The check is reconciled against the genome-wide count in `orf_reconciliation`: CDS checked
{{orf_reconciliation.CDS.in_table}}, `count_genes()` {{orf_reconciliation.CDS.counted}},
{{orf_reconciliation.CDS.agreement}}.
