## RNA editing — predicted C-to-U sites in {{rna_editing_summary.sample}}

*{{rna_editing_summary.organism}}*: **{{rna_editing_summary.n_sites}}** C-to-U sites are predicted
in {{rna_editing_summary.n_loci_with_sites}} loci at cutoff C = {{rna_editing_summary.cutoff}}.
Every site is **{{rna_editing_summary.method}}**: the published PREP method run by this toolkit,
not output of the original PREP-Cp, whose service no longer exists. A prediction is not a measured
editing frequency, and nothing in the record — sequence, translation or annotation — was changed.

Each supported CDS was translated with the standard code and added to its gene's fixed alignment of
mature reference proteins with {{rna_editing_summary.aligner}}. A codon containing a C is scored
against the reference residues in its alignment column, gaps excluded; the edited state is reported
when it wins and scores at least C. Ties: {{rna_editing_summary.tie_rule}}.

The reference is {{rna_editing_summary.reference}}, digest `{{rna_editing_summary.reference_digest}}`.
Its edit lists are sourced site by site, and the notes below say, taxon by taxon, how complete each
list is against the reference the original PREP-Cp used.

Loci evaluated: {{rna_editing_summary.n_evaluated}}; loci not evaluated:
{{rna_editing_summary.n_not_evaluated}}; supported genes absent from the record:
{{rna_editing_summary.n_absent}}; CDS outside the supported genes:
{{rna_editing_summary.n_not_supported}}. `rna_editing_coverage` gives every CDS its status and
reason — a gene that was not evaluated, not supported or absent is never counted as zero sites.
`rna_editing_sites` places every site on the genome (both strands stated), `rna_editing_query`
records what was extracted from each CDS, and the `alignments/` files let any score be recomputed by
hand.
