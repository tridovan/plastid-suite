## Synteny and identity — {{synteny_summary.accession}}

{{synteny_summary.n_genomes}} genome(s) aligned against *{{synteny_summary.organism}}*
({{synteny_summary.accession}}, {{synteny_summary.reference_length}} bp) with
{{synteny_summary.tool}} — `nucmer --maxmatch`, filtered to a one-to-one alignment, over a
{{synteny_summary.window_bp}} bp window with a {{synteny_summary.step_bp}} bp step.

The per-genome table `synteny_genomes` carries each genome's block count, aligned base pairs and
identity range; `synteny_blocks` lists every alignment; `synteny_identity` is the % identity
series along the reference.

Two numbers in this module are close enough to be confused, so each table names its own:

* the block table's **%IDY** is MUMmer's own — its error count over the alignment length, where
  an indel counts as an error;
* the series compares **aligned base pairs only** — a base the other genome does not carry is a
  gap, not a mismatch — so an indel-rich block reads slightly higher there.

The reference is in its normalized frame (LSC → IRa → SSC → IRb), so window coordinates mean
the same thing in every pairwise plot against it; query coordinates are the query's own. A
window with no aligned base pair has no value: an unaligned stretch is a hole, never a zero.

The raw `show-coords` and delta files travel with the run, so the block table can be re-read
without re-running the tool.
