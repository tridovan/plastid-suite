## Substitutions — {{substitutions_summary.accession}}

Every genome compared against *{{substitutions_summary.organism}}*
({{substitutions_summary.accession}}), aligned with {{substitutions_summary.tool}}:
{{substitutions_summary.n_genomes}} genomes over an alignment of
{{substitutions_summary.alignment_length}} columns.

Across all comparisons: {{substitutions_summary.substitutions}} substitutions in
{{substitutions_summary.compared_columns}} comparable columns, Ts/Tv
**{{substitutions_summary.ts_tv}}**.

Ts/Tv is left blank when no transversion was found — it is undefined, not infinite. The
per-genome table is `substitutions`; the per-region breakdown is `substitutions_by_region`.

**Scope:** {{substitutions_summary.scope}}.

The golden Ts/Tv in the project evidence log belongs to one scope — CDS only, shared genes only,
IR excluded — and a number measured in any other scope must never be quoted as if it were that
one.
