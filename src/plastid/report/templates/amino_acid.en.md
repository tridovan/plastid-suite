## Amino-acid composition — {{amino_acid.accession}}

*{{amino_acid.organism}}*: {{amino_acid.n_cds}} CDS, {{amino_acid.n_residues}} residues counted,
translated with {{amino_acid.translation_table}}. GRAVY (Kyte–Doolittle) is
**{{amino_acid.gravy}}** and aromaticity (Phe + Trp + Tyr) is **{{amino_acid.aromaticity}}**.

{{amino_acid.ignored_symbols}} stop or ambiguous symbol(s) were excluded — stops are not residues.
The per-residue composition is in `amino_acid_composition`.

The residue total here equals the sense-codon total of the `codon_usage` run on the same genome:
both modules take their CDS set from `cds_features`, so the two can never disagree.
