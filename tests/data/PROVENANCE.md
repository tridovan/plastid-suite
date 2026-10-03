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

## A real annotator run (P3-03)

| File | Source | Size |
|---|---|---|
| `chloe_NC_050833.1.gff` | Chloe 0.1.18 annotating the golden genome's sequence | 56,975 B |

- **What it is:** what the annotator we now use actually wrote, byte for byte, for the golden
  genome after its own annotation was stripped — the file every P3-03 number comes from
  (PROJECT.md §8.17). It is here so the read-back, the count **132 = 86 + 38 + 8** and the ORF
  profile (80 `ok` + 6 `alternative_start`, none failing) can be checked on any machine, without
  Julia and without the tool.
- **How it was produced:** `julia --project=<chloe> chloe.jl annotate --no-transform <golden
  sequence>.fa`, 2026-10-03, on the checkout at commit `1d7b710` with the reference set
  `Chloe_refs` 1.0.1 (`git-tree-sha1 9163915d…`). The `--no-transform` flag is the point: without
  it Chloe flips the single-copy region and writes coordinates in a frame the input never had.
- **Licence:** Chloe is deposited on Zenodo under the MIT licence (doi:10.5281/zenodo.18194040).
  Its own repository carries no LICENSE file — recorded for the M9.3 review.

## The stand-in submission template (P3-05)

| File | Source | Size |
|---|---|---|
| `standin_template.sbt` | written for this repository, 2026-10-03, from the ASN.1 structure the tool accepts | 1,758 B |

- **What it is:** a submitter template for `table2asn` with **placeholder** names — "Example
  Submitter", `example@example.org`, "Example University". The real one names people and an
  institution, is generated from NCBI's web form behind a login, and is therefore a **team
  artifact**: `--set template=/path/to/your.sbt`.
- **Why it is here:** so the `.sqn` half of P3-05 is testable. Measured 2026-10-03: with it,
  `table2asn` runs the golden genome's `.fsa` + `.tbl` to completion — exit 0, a 361 kB `.sqn`,
  and NCBI's validator report beside it. Without a template it refuses to start.
- **It must never be used to submit anything.** Its own title line says so, and every value in it
  is a placeholder — submitting with it would attribute a sequence to nobody.
- **Licence:** written here, so no third-party claim attaches to it. Its structure follows the
  templates published in public repositories (Cenote-Taker3, enviro-lab/ncbi-submit) and the
  format NCBI's own converter accepts.

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
