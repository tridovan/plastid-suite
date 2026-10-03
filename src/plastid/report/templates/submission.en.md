## NCBI submission files — {{submission_summary.accession}}

Two files were written for *{{submission_summary.organism}}*
({{submission_summary.accession}}), with the sequence id
`{{submission_summary.sequence_id}}` in both: `submission.fsa`, the sequence, and
`submission.tbl`, its annotation as NCBI's five-column feature table —
{{submission_summary.n_features}} features ({{submission_summary.n_cds}} CDS,
{{submission_summary.n_trna}} tRNA, {{submission_summary.n_rrna}} rRNA). NCBI's web submission
takes exactly this pair, so no external tool stands between this run and a submission.

The table re-reads into the same features it was written from, which the module's own tests
check on the golden genome — a file that does not re-read is not a result. Multi-interval
features are written by stacking their intervals, minus-strand intervals in reverse, and partial
coordinates carry their marker into the table: dropping one would turn "this gene runs off the
contig" into "this gene starts here", a different claim about the same sequence.

Each gene keeps the name it arrived with. A submission describes a record, and renaming a gene
would change what is submitted; the naming decision recorded for the manuscripts concerns what
they print.

The sequence submission file (`.sqn`): {{submission_summary.sqn}}.

{{submission_summary.n_problems}} first-pass problem(s) were found. The checks are the ones
that fail at upload time rather than here: a feature past the end of the sequence, a CDS with no
translation table or no product, two identical features that would generate two identical locus
tags. Findings are listed in `submission_problems` and are reported, never repaired — only a
human knows which side of a disagreement is the wrong one.
