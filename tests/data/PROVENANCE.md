# Golden test data — provenance

| File | Accession | Organism | Size |
|---|---|---|---|
| `NC_050833.1.gb` | NC_050833.1 | *Iris domestica* | 302,383 B |
| `NC_024936.1.gb` | NC_024936.1 | *Iris gatesii* | 298,692 B |

- **Source:** NCBI Nucleotide, fetched 2026-09-26 via E-utilities
  (`efetch … rettype=gbwithparts&retmode=text`) — PROJECT.md §8.3.
- **Copied into the repo:** 2026-10-02 from the validated probe dataset
  (`/mnt/d/business/pts/probe/data/`, PROJECT.md §11) — the exact files behind every measured
  number in PROJECT.md §8–§9.
- **Pinned:** `sha256sum -c SHA256SUMS` from this directory.
  (`NC_050833.1.gb` = `a946f23e…`, `NC_024936.1.gb` = `fd441727…`)
- **Licence:** NCBI GenBank data is in the public domain.
- **Rule:** tests never touch the network — these files are the only inputs.

## Real five-column tables (P3-05)

| File | Source | Size |
|---|---|---|
| `real_NC_039897.vadr.tbl` | `ncbi/vadr`, `documentation/build-files/NC_039897.vadr.tbl` | 492 B |
| `real_va-dengue-rpn.vadr.tbl` | `ncbi/vadr`, `testfiles/expected-files/va-dengue-rpn.r5/va-dengue-rpn.r5.vadr.pass.tbl` | 1,549 B |
| `real_va-noro.vadr.tbl` | `ncbi/vadr`, `documentation/annotate-files/va-noro.9.vadr.pass.tbl` | 1,907 B |

- **What they are:** feature tables NCBI's own pipeline (VADR) produced — not hand-written
  examples. They are the reference our writer is checked against, field for field, for the
  conventions that fail at *upload* time rather than here: qualifier placement (name in column 4),
  partial markers as prefixes (`<1`, `>515`), and the header/first-line pairing with the FASTA.
- **Fetched:** 2026-10-03 via the GitHub API (`gh api repos/ncbi/vadr/contents/<path>`).
- **Licence:** NCBI VADR carries a PUBLIC DOMAIN NOTICE (United States Government Work, no
  restriction on use or reproduction) — `ncbi/vadr`, `LICENSE`.
- **What they do NOT cover:** none of them has a multi-interval feature, so the *stacking* rule
  for joins is verified against NCBI's documented example (`tests/unit/test_submission.py`)
  rather than against a real submitted table — no such file was found.
- **One of them is a counter-example on purpose:** `real_va-dengue-rpn.vadr.tbl` holds **two**
  records in one table (VADR ran on two sequences). Our reader refuses a multi-record table
  instead of merging the records' features — measured while building the round-trip test, and the
  reason that file is kept.
