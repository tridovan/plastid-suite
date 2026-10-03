## Genome structure — {{structure_summary.accession}}

*{{structure_summary.organism}}* in the {{structure_summary.coord_system}} frame: position one is
the LSC's first base, and the GenBank origin lies {{structure_summary.rotation_offset}} bp into
this frame. Regions follow the order LSC → IRa → SSC → IRb and tile the genome exactly.

| Region | Start | End | Length (bp) | GC (%) | GC-skew (G−C)/(G+C) |
|---|---|---|---|---|---|
| LSC | {{regions.LSC.start_1b}} | {{regions.LSC.end_1b}} | {{regions.LSC.length_bp}} | {{regions.LSC.gc_pct}} | {{regions.LSC.gc_skew}} |
| IRa | {{regions.IRa.start_1b}} | {{regions.IRa.end_1b}} | {{regions.IRa.length_bp}} | {{regions.IRa.gc_pct}} | {{regions.IRa.gc_skew}} |
| SSC | {{regions.SSC.start_1b}} | {{regions.SSC.end_1b}} | {{regions.SSC.length_bp}} | {{regions.SSC.gc_pct}} | {{regions.SSC.gc_skew}} |
| IRb | {{regions.IRb.start_1b}} | {{regions.IRb.end_1b}} | {{regions.IRb.length_bp}} | {{regions.IRb.gc_pct}} | {{regions.IRb.gc_skew}} |

The per-window curve behind the skew figure is tabulated in `gc_skew_windows`:
{{structure_summary.window_bp}} bp window, {{structure_summary.step_bp}} bp step,
{{structure_summary.n_windows}} windows. Positive skew means G outnumbers C on the plus strand.
