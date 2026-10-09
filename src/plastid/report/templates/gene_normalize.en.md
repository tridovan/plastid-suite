## Reference-driven gene naming

The run read {{normalize_summary.genes}} gene(s) from {{normalize_summary.records}} record(s) and
paired each one against **{{normalize_summary.reference}}** — the project's naming reference, the
record whose spelling every table in this project prints. Its names are the left-hand column of
`pairing`; the name a record arrived with sits beside them in its own column, and `pairing` is the
only table in which those names appear at all: one vocabulary per table, always.

Of those genes, **{{normalize_summary.renamed}}** were written a different way and are spelled the
reference's way in the output, and **{{normalize_summary.confirmed}}** needed nothing. The
`rule` column says which of the four tests decided each pairing — the name, the sequence ranking,
the margin between the best match and the runner-up, or the flanking genes inside one region — so
"the name and the sequence agree" can be told apart from "only the gene order could say".

### The three statuses, and why they are not three words for the same thing

{{normalize_summary.measured}} gene(s) are **measured**: a counterpart in the reference was found
and checked. {{normalize_summary.not_measured}} are **not measured** — the check could not run, or
ran and did not settle the pairing (no protein on one side, a ranking too close to trust, flanking
genes that disagree). That is not the same statement as "there is none", which is what the
{{normalize_summary.outside}} gene(s) marked **outside vocabulary** say: the check ran and the
reference has nothing corresponding. Without the distinction, a zero and a blank would look alike
in every later count.

### What differs, and what stopped

{{normalize_summary.conflicts}} conflict(s) were found between a name and its sequence. A conflict
stops the run by default — the names and the sequences disagree, and continuing would write a
record whose genes are named after loci they are not. `conflicts` names every case with its
coordinates; `--set on_conflict=continue` writes the output anyway and keeps that table.

`differences` lists where the record and the reference differ without either being wrong: extra
genes, missing genes, a different `/product`, a length that differs, a gene split into a different
number of parts. Those are reported and never stop the run — a lost gene is a biological result,
and stopping on it would make a real finding look like a broken file.

{{normalize_summary.artifacts}} normalized record(s) were written beside this report, each named
the reference's way and each keeping the record's own `/translation`: the protein a record carries
is what the annotator searches, and a normalization that dropped it would disarm the step standing
behind it.
