## Mitochondrial annotation — {{mito_annotation_summary.accession}}

The mitogenome of *{{mito_annotation_summary.organism}}* ({{mito_annotation_summary.length_bp}} bp)
was annotated by transferring a gene model from a curated reference set:
**{{mito_annotation_summary.n_calls}}** genes were called ({{mito_annotation_summary.n_cds}} CDS and
{{mito_annotation_summary.n_rrna}} rRNA), and **{{mito_annotation_summary.n_not_called}}**
vocabulary loci were not, each for a reason named in `not_called`.

The reference set is {{mito_annotation_summary.reference}}, searched with
{{mito_annotation_summary.blast}}. A call had to clear an identity floor of
{{mito_annotation_summary.min_identity_pct}} % and cover
{{mito_annotation_summary.min_coverage_pct}} % of its reference protein at an E-value of at most
{{mito_annotation_summary.evalue}}. Those floors are stated rather than implied because they are
the whole difference between a call and a guess: {{mito_annotation_summary.n_dropped}} hits did not
clear them and are listed with their coordinates and their reason.

This is a **prototype**, and the tables say what it does not do. Boundaries come from the alignment
to a reference protein, never from ORF-walking: the stop codon is taken only where the sequence
confirms it, and an end the alignment does not establish is reported as partial — never extended to
the next stop, and never moved to the next start codon. Multi-exon genes appear as the exons the
alignment supports, in the order the reference gives them; no join is claimed. There is no de-novo
gene discovery, no intron model and no RNA-editing prediction, and tRNA is not called here (a
separate step runs ARAGORN). Genes found whose symbol this toolkit does not hold —
{{mito_annotation_summary.n_unnameable}} of them — are reported in `unnamed_hits` rather than
named, because a name is a claim about a locus and the vocabulary that would make it is a
decision of its own.

Of the called CDS, **{{mito_annotation_summary.n_cds_failing_orf}}** do not translate the way
their coordinates claim, measured by the same checks a downloaded record gets. A failure here is
not repaired — it is the finding: the answer-key record this caller is measured against fails
those checks in exactly the same ways, and its failures are RNA-editing-dependent ends rather than
errors.
