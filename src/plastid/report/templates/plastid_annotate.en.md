## Plastid annotation — {{plastid_annotation_summary.accession}}

The plastome of *{{plastid_annotation_summary.organism}}*
({{plastid_annotation_summary.length_bp}} bp) was annotated by transferring a gene model from a
curated reference set: **{{plastid_annotation_summary.n_calls}}** loci were called
({{plastid_annotation_summary.n_cds}} CDS, {{plastid_annotation_summary.n_rrna}} rRNA and
{{plastid_annotation_summary.n_trna}} tRNA), and
**{{plastid_annotation_summary.n_not_called}}** vocabulary loci were not, each for a reason named
in `not_called`.

These are **calls, one per locus** — the count a reader compares against a gene list, not against
the artifact's feature count. The artifact written alongside this report also carries the inverted
repeat's second copy of every gene inside the repeat, so its CDS and rRNA totals are higher than
the call counts here; a published plastome's gene table counts those copies.

The reference set is {{plastid_annotation_summary.reference}}, searched with
{{plastid_annotation_summary.blast}} (CDS, against each record's own published `/translation`) and
{{plastid_annotation_summary.blast}} (rRNA); tRNA is called de novo by
{{plastid_annotation_summary.aragorn}} and needs no reference. A call had to clear an identity
floor of {{plastid_annotation_summary.min_identity_pct}} % and cover
{{plastid_annotation_summary.min_coverage_pct}} % of its reference sequence at an E-value of at
most {{plastid_annotation_summary.evalue}}, under translation table
{{plastid_annotation_summary.transl_table}}. Those floors are stated rather than implied because
they are the whole difference between a call and a guess:
{{plastid_annotation_summary.n_dropped}} hits did not clear them and are listed with their numbers
and their reason.

**One call per locus, and the losers are listed.** Where several references describe one locus, the
call this run reports is the one that won, and every other candidate is in
`competing_references` with the numbers it lost on — so a reader can disagree with the arbitration
without re-running anything. **Overlap is not an error**: the published records this engine is
measured against carry real overlapping gene pairs, so both calls stand and the pair is listed
({{plastid_annotation_summary.n_overlapping_pairs}} of them). A gene this vocabulary cannot name is
reported in `unnamed_hits` rather than renamed — the count is
{{plastid_annotation_summary.n_unnameable}} — because inventing a symbol is a claim the run has no
evidence for.

The record written alongside this report carries `/organelle` and the inverted repeat's two copies
of every gene inside it, and it is **read back** before the run ends, so the artifact is one the
rest of the pipeline can consume rather than one this report says exists.

Of the called CDS, **{{plastid_annotation_summary.n_cds_failing_orf}}** do not translate the way
their coordinates claim, measured by the same checks a downloaded record gets. A failure here is
not repaired — it is the finding.
