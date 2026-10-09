## Reference-set audit

{{audit_summary.records}} record(s) were audited. **{{audit_summary.blocked}}** of them cannot be used
as a reference, and **{{audit_summary.reported}}** further finding(s) are reported without blocking
the run.

A blocked record is one no annotation may be transferred *from*: unreadable, with no gene features,
with a CDS that carries no `/translation`, with a partial coordinate, failing the IR gate, or
holding a CDS whose own sequence does not translate. Every one of those faults would travel into
whatever the record annotates — measured, a reference with one partial coordinate returned zero
genes from a real query while reporting success.

A reported record is one the set disagrees about, and the disagreement is not an error. A plastome's
SSC genuinely inverts (flip-flop) with every one of its own genes translating cleanly, and gene
content genuinely differs between annotations of the same sequence. Blocking those would refuse
published records, so they are named here for a person instead.

Every finding is written as `file — reason`, and the reason is the whole message: what is wrong with
the bytes, and nothing about what to do next. The findings, blocking ones first, are in
`audit_findings` when there are any; every per-record measurement is in `audit_records`.

### The second reading

{{audit_summary.confirmation}}

The gene-strand verdict compares a record's SSC with its own LSC, so it cannot check itself. An
aligner asked to adjust direction reads the region across records instead, and where the two
readings agree the finding is settled; where they disagree, `region_orientation` shows both, and
that is the state in which a person has to look at the dot plot before confirming anything.

### Confirmed exceptions

{{audit_summary.confirmed}} finding(s) are confirmed. A confirmation is pinned beside the reference
set together with the file's digest: it stops applying the moment one byte of that file changes, and
the finding returns. A blocking finding can never be confirmed away — such a file has to be fixed,
replaced, or dropped from the set.

### Choosing the project's naming reference

Every name every record carries is in `name_sets`, with the names the majority of the set carries
marked. Each record's vocabulary is one candidate for the reference the whole project will print, so
this table is what the choice is made from — the record that is picked becomes the spelling every
later table follows, and the analysis that must be re-measured when it changes is re-measured then,
not assumed.
