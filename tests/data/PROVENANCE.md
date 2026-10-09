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

## The Điềm corpus fixtures (P4a-02)

| File | Accession | Annotation | Size |
|---|---|---|---|
| `corpus/Iris_domestica_NC_050833.1.pga.gb` | `Iris_domestica_NC_050833` | PGA | 230,133 B |
| `corpus/Iris_hermona_NC_081931.1.published.gb` | NC_081931.1 | published NCBI | 298,851 B |
| `corpus/Iris_hermona_NC_081931.1.pga.gb` | `Iris_hermona_voucher_IrisSV1903_NC_081931` | PGA | 229,201 B |
| `corpus/Eleutherine_bulbosa_PX516262.1.published.gb` | PX516262.1 | published NCBI | 272,037 B |
| `corpus/Eleutherine_bulbosa_PX516262.1.pga.gb` | `Eleutherine_bulbosa` | PGA | 228,555 B |

- **What the corpus is:** 34 Iridaceae plastomes that **anh Điềm** provided on 2026-10-04, in two
  annotations of the *same sequences* — the published NCBI record and a **PGA** re-annotation —
  plus the bare sequences (`.fasta`). 136 files · 28 MB. Both annotation sets and the FASTA
  manifest are pinned under `.panel/iridaceae-34{,-pga}/` (gitignored, re-pinnable); the full
  characterisation is PROJECT.md §8.31.
- **Why these five files:** the corpus is far too large to vendor, so what ships is the smallest
  subset in which **every file carries a written-down expected answer** — three pairs, three jobs:
  the golden genome both ways (the contract is already written down, and PGA's annotation of the
  same sequence measures **85 CDS where §6 says 86**), *Iris hermona* both ways (one sequence; the
  IR gate passes on one annotation and fails on the other at **26.36 %**, while de-novo finds
  100.00 % on the same bases), and the group's own *E. bulbosa* (133 = 87 + 38 + 8 published,
  131 = 85 + 38 + 8 by PGA). `tests/golden/test_corpus.py` asserts each one.
- **The golden genome's *published* half is already here** as `NC_050833.1.gb`, so it is not
  duplicated: set 01's own file is **byte-identical** to it (sha256 `a946f23e…`), which is how the
  corpus was checked against this directory in the first place.
- **The tool is named by the files, not inferred:** all 34 records of the second set carry
  `/note="Annotation Method :: PGA-Plastid Genome Annotator"` on their `source` feature; none of
  the first set does. The naming style alone would not have been evidence.
- **Source:** `D:\business\pts\thanhdiem_data\{01_NCBI_tai_truc_tiep,02_da_hieu_chinh}\gb\`
  (anh Điềm, 2026-10-04), copied into `.data/thanhdiem_data/` in the working tree (gitignored) and
  verified file by file — 136/136 sha256 match the source.
- **What these fixtures do NOT cover:** the other 29 genomes, the `.fasta` half (its bytes are in
  `.panel/iridaceae-34/fasta.json`, and every one of its 68 sequences was proven identical to a
  pinned record), and anything about *their* numbers — the corpus is a test dataset, not a
  manuscript's data: no species here is Oldenlandia or Coleus.
- **Licence:** NCBI GenBank data is in the public domain. The PGA annotation is a derived work on
  the same public sequences (PGA is GPL-3.0, PROJECT.md §3 — it is not redistributed here).

## A real delta file (P4-01)

| File | Source | Size |
|---|---|---|
| `iris_pair_sample.delta` | two alignment records of a MUMmer 4.0.1 `nucmer` run on the golden pair | 603 B |

- **What it is:** the first (forward) and the inverted-repeat (reverse) records of the measured
  run, with the file's own header lines. `tests/unit/test_delta.py` parses it and checks the
  spans, the orientation and the walk's totals against what the tool printed — the synthetic
  cases around it were written by the same hand as the parser, so this is the one case that can
  catch a wrong reading of the format.
- **How it was produced:** `nucmer --maxmatch` on `NC_050833.1` (target) vs `NC_024936.1`
  (query), 2026-10-04, MUMmer 4.0.1 from bioconda (`docs/tools/mummer4.md`); records extracted
  from `out.delta` with the block whose query span is `153441..127260` kept as the reverse case.
- **Licence:** generated data, no third-party content — the input sequences are the two public
  GenBank records already vendored here.
