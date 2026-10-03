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
