## Plastome assembly — {{assembly_summary.sample}}

{{assembly_summary.assembler}}, {{assembly_summary.where}}. It ran
{{assembly_summary.assembler_programs}}, with the seed database {{assembly_summary.seed_database}}.
The assembler's result is **{{assembly_summary.status}}** ("{{assembly_summary.assembler_says}}"),
and the read support is **{{assembly_summary.read_support}}**.

The assembled length is {{assembly_summary.length_bp}} bp in {{assembly_summary.n_paths}} path(s),
with GC {{assembly_summary.gc_percent}} %. The sequence's own repeat gives IR
{{assembly_summary.ir_bp}} bp, LSC {{assembly_summary.lsc_bp}} bp and SSC
{{assembly_summary.ssc_bp}} bp. An inverted repeat gives two closed paths because the small
single-copy region recombines; reads shorter than the repeat cannot separate those isomers, so both
are kept and neither is a finding.

The assembler's "complete" is a claim that a path closed. It is tested here with reads mapped back by
an aligner the assembler does not use: {{assembly_summary.mapped_reads}} reads, mean depth
{{assembly_summary.mean_depth}} ± {{assembly_summary.sd_depth}}. The fewest reads crossing any
boundary of the circle is {{assembly_summary.fewest_spanning}}. `junction_support` reads the origin
and the junctions, `weak_regions` lists any boundary too few reads cross, and a join the reads do
not cross is **undetermined**, never verified.

Reads cannot say which arrangement of the pieces is the molecule, because every arrangement crosses
the same joins. That is read from the selected graph, without the assembler's words: `assembly_graph`
counts the structures that close through it and `graph_segments` lists its pieces. The graph is
**{{assembly_summary.graph_call}}**; it leaves open: {{assembly_summary.graph_reasons}}.

Only a closed sequence with every base known, every join supported and a resolved graph is declared
in the manifest as an artifact for the next step.
