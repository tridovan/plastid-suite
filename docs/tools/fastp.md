# fastp 1.4.0 — quality control of short reads

**Why this file exists.** `plastid assemble` runs fastp before the assembler, and a module's tool is
declared in the repository. It lives in `environment-assembly.yml` with the rest of the assembly
tools, because the analyses never call it.

**What it does here.** It removes adapter read-through and low-quality ends from paired reads
before assembly. A read with adapter on it looks like a join to nowhere, in the assembler's graph
and in the map-back. Its report's counts go into the run's `read_qc` table, read from the JSON
report, and the report files themselves stay in the run directory.

Its output goes to the assembler and nowhere else. The map-back reads the raw reads, because a check
that read trimmed ones would lose the evidence trimming removed.

**Install line:**

```sh
conda env create -f environment-assembly.yml
```

**The invocation the adapter uses:**

```
fastp --in1 R1 --in2 R2 --out1 out1 --out2 out2 --detect_adapter_for_pe --thread N --json fastp.json --html fastp.html
```

**Two runs give the same bytes.** Measured 31/10/2026 on `SRR5602589` with `--thread 4`: both output
files of two runs had the same sha256. That does not excuse the import from checking: `plastid
assemble prepare` hashes the files it hands over and the import compares those exact bytes, because
the input of a run elsewhere is whatever was uploaded, not whatever a later fastp would write.

**The version is checked.** A different release is refused by name (`FASTP_VERSION` in
`src/plastid/external/fastp.py`), because trimming decides which reads the assembler sees.

**Skipping it.** `--set trim=no` passes the reads through unchanged. The manifest then records no
fastp step, and the `read_qc` table is absent rather than empty.

Citation: Chen S, Zhou Y, Chen Y, Gu J. *Bioinformatics* 34:i884–i890 (2018).
