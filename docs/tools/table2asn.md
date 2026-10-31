# table2asn — declared, not vendored

**Why this file exists.** CLAUDE.md: a tool that exists only on one machine is not part of the
toolkit. The moment a module uses an external tool, that tool must be declared in the repo — in
`environment.yml`, in `vendor/` with a `SOURCE.json`, or here. table2asn fits neither of the first
two. It is on bioconda, but it cannot join `environment.yml`: 1.29.324 needs `libsqlite` ≥ 3.54.1,
which the solver cannot give beside `r-base=4.3` (measured 31 October 2026). And it is a 56 MB NCBI
binary, too large to vendor. So it is declared here, with its source and its exact install lines, and
the module that uses it records name · version · path · argv at every run.

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

**Or from bioconda, in an environment of its own** (measured 31 October 2026, 62 MB installed):

```sh
micromamba create -n table2asn -c conda-forge -c bioconda table2asn=1.29.324
```

That is build `h190543b_0`, uploaded 28 October 2026. It reports the same version and is not the same
file: 56,851,600 B, md5 `6ae6ffcb80670d42b1ade76af0bf90a4`.
`test_the_sqn_step_runs_and_ncbi_reads_the_table` passes on both.

**Not vendored, deliberately:** 56 MB of third-party binary in a git repository whose golden data
is a 170 kB GenBank file. It is fetched by whoever needs the `.sqn` step, and the step is optional
(D3-04: `.fsa` + `.tbl` are what NCBI's web submission accepts).

**It does not work alone.** table2asn refuses to run without a submitter template (`.sbt`) — a
text ASN.1 ``Submit-block`` naming the people and the institution a sequence is submitted by.
NCBI's web form (<https://submit.ncbi.nlm.nih.gov/genbank/template/submission/>) generates the
real one, and it needs an NCBI login, so a template that carries **Trang's** details is a
**team artifact, not a toolkit one**: `--set template=/path/to/your.sbt`.

**Corrected 2026-10-03, by running it.** An earlier session recorded "a hand-written `.sbt` is
rejected as an unreadable template", after two attempts failed with ``Error loading template
file``. That is not what the tool does: the *structure* has to be right, and with it a template
loads. The repository's stand-in (`tests/data/standin_template.sbt`, placeholder names) runs the
whole `.sqn` step to completion — `exit 0`, a 361 kB ``.sqn``, and NCBI's validator report beside
it. What the two failed attempts had was a wrong ASN.1 spelling, not a wrong idea.

**What the step then reports is the point**: run on the golden genome, NCBI's
validator found **13 errors** in our table at first — a missing source feature, an RNA-edited
start written without its `/exception`, and that same start without a partial marker. Three
fixes later it reports **7**, every one of them the ``GeneXrefStrandProblem`` class that the same
validator produces for a **published, accepted** plastome submission with duplicated IR genes.
Its findings — errors and warnings, by class — are in the run's own tables, and the full report
travels as ``submission.val``.
