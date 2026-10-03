# table2asn — declared, not vendored

**Why this file exists.** CLAUDE.md: a tool that exists only on one machine is not part of the
toolkit. The moment a module uses an external tool, that tool must be declared in the repo — in
`environment.yml`, in `vendor/` with a `SOURCE.json`, or here. table2asn fits none of the first
two: it is not on conda and it is a 56 MB NCBI binary, so it is declared here, with its source and
its exact install line, and the module that uses it records name · version · path · argv at every
run.

**What it is.** NCBI's ASN.1 converter: FASTA + five-column feature table in, `.sqn` out. It is
the only step of P3-05 that is not pure Python, and the only one that needs a file the team owns.

**Install line** (verified 2026-10-03 on the dev WSL machine):

```sh
curl -C - -o table2asn.gz \
  https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz
gunzip table2asn.gz && chmod +x table2asn
```

Measured: **20,213,344 B** gzipped → **56,847,504 B** binary, `table2asn -h` → version
**1.29.324** (the tool prints a notice that it is more than a year old; NCBI's FTP serves that
build). `md5sum` of the unpacked binary: `8b5b6cbede9d38b47c7c4cae29a25a0a` — the checksum is recorded so a second
machine can prove it installed the same build.

**Not vendored, deliberately:** 56 MB of third-party binary in a git repository whose golden data
is a 170 kB GenBank file. It is fetched by whoever needs the `.sqn` step, and the step is optional
(D3-04: `.fsa` + `.tbl` are what NCBI's web submission accepts).

**It does not work alone.** table2asn refuses to run without a submitter template (`.sbt`), and a
hand-written one is rejected with `Error loading template file` — measured here, twice, with two
different ASN.1 spellings. The template is generated from NCBI's web form
(<https://submit.ncbi.nlm.nih.gov/genbank/template/submission/>) because it carries the
submitter's name, institution and contact address. That makes it a **team artifact**, not a
toolkit one: `--set template=/path/to/your.sbt`.

So the `.sqn` half of P3-05 has two preconditions, and only one of them is a machine's:
an installed table2asn (this file) **and** the team's `.sbt`.
