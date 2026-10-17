# TASKS.md — the task in flight

One task at a time. When the next one begins this file is rewritten, not ticked: the finished task's
record is its commit. A task is done only when its **Done when** check passes.

**Phase 10 is complete and H-25 is fixed.** Every figure anh Điềm's feedback bundle
(`feedback/goi_tool_hinh_IR_2026-10-16/`) names renders from a run directory, and the IR gate no
longer refuses an exact repeat whose annotated copies differ in length. Read the **whole** bundle
before working anything that cites it — a session that read only its map designed on a premise the
tables do not support, and that has happened once already.

## M7.1 — batch assembly, judged on a published paper

**Goal.** Phase 7 is the only phase with nothing built. It is judged the way every other tool in this
project is judged (CLAUDE.md → How the tools are judged): find a **published paper with public
Illumina reads and its assembled plastome**, assemble those reads with our own chain, and compare
what comes out with the paper's own result. Not a simulated pair — the ROADMAP's 127,679 bp from
20,000 simulated pairs was a feasibility check, not an answer key.

**Steps**
1. **Find the input first.** One paper, one accession set (SRA/ENA), one assembled plastome in
   GenBank. Record accession · date · sha256 with `plastid panel pin`, and label every output made
   from it. If no suitable paper exists, that *is* the finding — say so and stop, do not fall back
   to simulation and call it an evaluation.
2. Declare the tools **in the session they are first used**: `fastp` and `BWA` (or `minimap2`) are in
   neither `environment.yml` nor `docs/tools/` today. Version pinned, install line written, the way
   `docs/tools/r.md` and `docs/tools/misa.md` do it.
3. The chain: `fastp` → `GetOrganelle` → map-back with `BWA`/`minimap2` → the assembled record.
   Every external call through `run_tool()`, so name · version · argv · exit code reach the run's
   provenance.
4. Compare **by input → output**: length, IR boundaries, gene count, and identity to the paper's
   assembled plastome. A number we cannot reproduce is reported as measured, with its cause or
   "unexplained" — never tuned to match.

**Done when:** a run directory holds the assembled record from the paper's own public reads, its
provenance names every tool with its version, a test pins the comparison against the paper's
plastome, and `make check` is green.

**Traps**
- **Peak RAM is part of the result** (M7.5). The user's machine is an ordinary computer; a chain that
  needs 64 GB is a chain most users cannot run, and the Galaxy remote path exists for exactly that.
- `GetOrganelle` and `fastp` are not pure Python: the minimum install must still run everything that
  does not need them, and a missing tool is **named**, never skipped silently.
- The remote path (M7.4) uploads the user's reads to usegalaxy.eu — subsample first, pin the tool
  version, record the history ID, and the data-sharing policy for unpublished reads is **open**
  (PROJECT.md §10.12). Do not take that step without the dev's explicit go.

## After it

**M7.2** (assembly-graph checks with Bandage, after M7.1) and **M7.6** (long reads). The named
defects **H-21 … H-24** remain open in `ROADMAP.md` → Known defects; each is its own task with its
own regression test, never a silent fix inside another.
