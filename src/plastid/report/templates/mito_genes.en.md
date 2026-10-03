## Mitochondrial gene set — {{mito_summary.accession}}

The mitogenome of *{{mito_summary.organism}}* ({{mito_summary.accession}},
{{mito_summary.length_bp}} bp, {{mito_summary.topology}}) carries **{{mito_summary.n_genes}}**
genes: {{mito_summary.n_cds}} protein-coding, {{mito_summary.n_trna}} tRNA and
{{mito_summary.n_rrna}} rRNA.

The inverted-repeat gate does not apply here, deliberately. A plant mitogenome has no canonical
IR — both stand-in mitogenomes fail IR detection, and they are right to — so this module never
asks for one. The gate stays where it belongs, in front of position-based analysis: running a
gated module on the same genome still stops the run.

Naming is the other half of the mitochondrial problem.
{{mito_summary.n_named_by_symbol}} protein-coding genes carry a standard
plant-mitochondrial symbol, mapped from the product text the record itself gives — some
mitogenome records carry no gene qualifier at all, so their "name" is a product sentence. The
mapping lives in its own file, separate from the plastid aliases, because the same sentence
occurs in both organelles and a shared table would rename a plastid gene to a mitochondrial
symbol. Anything the mapping does not know is listed in `unmapped_names` exactly as it arrived:
a name this toolkit cannot map is a name it does not invent. Every row of `gene_set` keeps the
record's own wording beside the symbol, so any name can be traced back to the record it came
from.
