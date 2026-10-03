# Phase 3 — blocker clearance

Written 2026-10-03. Purpose: so the next session can run Phase 3 A→Z without stopping to ask
anything. Each blocker names **what exactly is needed** and **from whom** — a blocker that says
"needs the team's call" is not a blocker, it is a delay.

## Cleared here, verified

| Blocker | Was | Now |
|---|---|---|
| **H-01 · MAFFT 7.526** for P2-07/P2-08 (π and Ts/Tv) | absent; `apt` installs a *different* version, and §9.B-13 records that the aligner changes Ts/Tv | **installed and verified against the §7.8 pin**: conda-forge `mafft 7.526` build `h4bc722e_0`, tarball **2,610,460 B**, md5 **`cf3b2acf649bc31e576ed43a2fa4aaa3`** — byte-identical to the recorded pin. `mafft --version` → `v7.526 (2024/Apr/26)` |
| **P3-03 · BLAST+** for PGA | absent | installed, **2.17.0+** (PGA needs ≥ 2.8.1) |
| **P3-03 · Perl** for PGA and MISA | present | present, unchanged |
| **P3-05 · `table2asn`** for `.sqn` | absent, and not on conda | **not finished** — NCBI's FTP is slow and the download timed out twice here (11.6 MB of a partial file sits at `~/toolbox/bin/table2asn.gz`). Source is confirmed live: `ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz` (HTTP 200). One `curl -C - -o table2asn.gz <url> && gunzip table2asn.gz` finishes it. **This blocks only the `.sqn` conversion**, not P3-05's FASTA + TBL half, which is pure Python |

Installed under `~/toolbox/` (env `~/toolbox/env/pts`, micromamba at `~/toolbox/bin/micromamba`).
Nothing was added to the repo: the project pins tool *provenance*, not tool *binaries*.

## Two blockers that turned out to be wrong blockers

**1. MITOS2 cannot annotate a plant mitochondrion.**

M3.1 says "MITOS2/GeSeq (mito)". MITOS2 is a **metazoan** mitogenome annotator: it finds genes by
`blastx` against **metazoan** RefSeq mitochondrial proteins (default reference set `refseq81m`).
On a plant mitogenome that reference set is the wrong kingdom, so its output would be quietly
wrong rather than loudly empty — the exact failure mode this project keeps meeting.

**The tool that fits is PMGA** (Plant Mitochondrial Genome Annotator, Li et al. 2025, *Plant
Communications* 6(3):101191). It is purpose-built for plant mitogenomes, takes one-to-many FASTA,
emits GFF3 + GenBank, and in its own benchmark against *Arabidopsis thaliana* and *Liriodendron
tulipifera* gold standards it produced **1 wrong gene** where GeSeq produced 18 and AGORA 42.
Its RNA-editing module runs **Deepred-mt**, which §7.8 already lists as a known tool.

Cost to note before committing: PMGA wraps MAKER, CPGAVAS2, BLASTn, ARAGORN, tRNAscan-SE, MISA,
TRF and Vmatch. That is a heavy environment — much heavier than PGA's Perl + BLAST+. **Decision
needed at P3-04, not now**, but it should be made with that cost on the table.

*Evidence quality: MITOS2's scope is stated plainly by its own distribution and documentation;
the PMGA benchmark numbers are from a web summary of the paper, not read firsthand. Verify the
paper and PMGA's install requirements before adopting it.*

**2. The naming question blocked nothing.** It looked like a PI decision; it was a dev decision
the moment the four names were measured by locus (PROJECT.md §8.9–8.10). It is now closed.

## Still blocked, with the exact ask

| # | What is needed | From whom | Blocks | The ask, in one line |
|---|---|---|---|---|
| 1 | **Confirmation of who `Vu,H.-T.` is, and whether `PX923121.1` is the specimen this manuscript reports** | cô Trang | the *C. amboinicus* numbers being quoted as the manuscript's | *"Có phải cô nộp PX923121.1 (C. amboinicus, 24/01/2026) không, và mẫu đó có phải mẫu bài báo dùng không?"* |
| 2 | **The team's real Oldenlandia diffusa record** | anh Điềm / cô Trang | the *O. diffusa* half of P2-10 — today it runs on `PP078674.1`, which is **Shen C.'s** submission, a genuine stand-in | *"Nhóm có bản ghi O. diffusa của mình chưa? Nếu chưa thì bản thảo dùng bản của Shen C. hay chờ mẫu?"* |
| 3 | **A decision on the mito annotator** (PMGA vs GeSeq vs MFannot/AGORA), knowing PMGA's dependency weight | dev + anh Điềm (environment) | P3-04 | make it at P3-04 with `docs/D3-01_annotator.md`'s reasoning applied to mito |
| 4 | **The gene name to print in the manuscript and submit to NCBI** | cô Trang | the *wording* of the gene table only — the pipeline is unaffected | recommend `pafI (ycf3)`; see D3-06 |
| 5 | **Item (c) text for the Reference Panel notes** | cô Trang | D-03 (docs) only | the §13.1 paragraph is truncated in the source |
| 6 | **PGA needs a first run by a human** | anh Trí | P3-03 | `git clone https://github.com/quxiaojian/PGA` (GPL-3.0, Perl, 58★, active) into `~/toolbox/src/`, then `perl PGA.pl` once. Executing freshly-cloned third-party code was refused here by the permission classifier — correctly, since it was not a source named in the session. One manual run settles it |
| 7 | **`table2asn` download** | anh Trí, or a run with a longer network window | the `.sqn` half of P3-05 only | `curl -C - -o ~/toolbox/bin/table2asn.gz https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz` |

## Not blocked — the next session can start immediately

- **P3-01 ORF checks** — pure Python, no tool, fully verifiable here
- **P3-02 gene-name normalization** — pure Python; the four semantic aliases are decided (§8.10)
- **P3-05's FASTA + TBL half** — pure Python; only the `.sqn` conversion needs `table2asn`
- **P2-07 / P2-08 real-data runs** — MAFFT is now installed and pinned, so the numbers that were
  "implemented and unit-verified, real run pending" can finally be *measured*
