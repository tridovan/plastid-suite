# Licensing — our licence, and every tool's (M9.3, decision D-p · D-q, 11/10/2026)

**Why this file exists.** A release asks two questions that must be answerable without a lawyer:
*what licence is this toolkit under*, and *what exactly travels inside it*. The second question is
the one that bites — a toolkit is a thin layer over a dozen other people's programs, and the
difference between **calling** a tool and **shipping** it is the whole of licence compliance here.
CLAUDE.md's invariant already says a tool must be declared in the repo; this file says, in one
table, **which of them we redistribute and which we merely run** — and what we may not bundle.

**The rule, stated once.** A work with no licence is **all rights reserved**: public to download is
not free to redistribute. Copying such a file into a wheel, a git repository or a container image
*is* redistribution. Running a program the user installed is not.

## 1. The toolkit's own licence: **MIT** (D-p)

`LICENSE` (MIT, "PlastidToolSuite contributors"), declared in `pyproject.toml` as
`license = "MIT"` + `license-files = ["LICENSE"]`, so every wheel and sdist carries it.

Why MIT: the runtime dependencies are permissive (PyYAML MIT; Biopython's own agreement, a
BSD/MIT-style permissive licence), the ecosystem this toolkit cross-checks against is MIT
(CPStools, CGAS), and the project's stated audience — academic researchers — should be able to
install and modify it without a permissions conversation. "Academic, non-commercial" describes
**the project**, not a restriction placed on the code; a licence that forbade commercial use would
be a different, and here unnecessary, choice.

## 2. What travels inside the package

| What | Licence | Note |
|---|---|---|
| All Python under `src/plastid/` | MIT | ours |
| `src/plastid/data/*.yaml` · `*.json` · `*.md` | MIT | ours — defaults, citations, gene groups, the reference sets' notes |
| `src/plastid/data/misa.ini` | MIT | **ours**: MISA's *default* definitions written down (`1-10 2-6 3-5 4-5 5-5 6-5`, interruptions 100, `GFF: true`). Configuration values describing a tool's documented behaviour; the tool's own file is not copied (D-q) |
| `src/plastid/data/plastid_reference.json` · `mito_reference_*.faa/fa` | derived | built from **public GenBank records** (NCBI — U.S. Government work, no copyright); the derivation is ours |
| `tests/data/**` | derived | GenBank records from NCBI, the corpus's public records, and PGA's *output* on them (a program's output is not covered by the program's licence) |

**Nothing else.** No external tool ships in the wheel, the sdist or the repository: `uv build
--wheel` before this change carried 146 files including `plastid/vendor/misa/misa.pl`; it no longer
carries any third-party program. The guard is `tests/unit/test_licensing.py`.

## 3. Every tool the toolkit touches

Licences below are **measured**, not remembered: conda package metadata read on the dev machine
(`conda-meta/*.json`, 11/10/2026) or the tool's own declaration file. "Bundled" means *redistributed
by us* — which is **none of them**.

| Tool | Version | Called by | Licence | How we know | Bundled |
|---|---|---|---|---|---|
| **MISA** | 2.1 | `external/misa.py` (module `ssr`) | **NOT STATED** — the zip ships no licence file; the IPK page only asks for the citation | the zip itself + the IPK page, 02/10/2026 | **no — left the package 11/10/2026** (D-q) |
| **MAFFT** | 7.526 | `external/mafft.py` (`substitutions`, `diversity`, `barcodes`, `popgen`, `kaks`) | BSD-3-Clause | conda-forge metadata | no |
| **BLAST+** | 2.17.0 | `external/blast.py` (mito annotation, corpus bar) | NCBI-PD (public domain) | bioconda metadata | no |
| **ARAGORN** | 1.2.41 | `external/aragorn.py` (mito tRNA) | GPLv3 | bioconda metadata + `docs/tools/aragorn.md` | no |
| **MUMmer 4** | 4.0.1 | `external/mummer.py`, `external/repeatmatch.py` (`synteny`, `repeats`) | Artistic License 2.0 | bioconda metadata + `docs/tools/mummer4.md` | no |
| **PAML** (`yn00`) | 4.10.10 | `external/paml.py` (module `kaks`) | GPL-3.0-or-later | bioconda metadata (+ academic-free upstream, `docs/tools/paml.md`) | no |
| **Vmatch** | 2.3.1 | `external/vmatch.py` — the *second engine* for the repeats count; never required | **academic use only** (upstream); conda metadata says `Unknown` | `docs/tools/vmatch.md` + bioconda metadata | no |
| **table2asn** | 1.29.324 | module `submission` | NCBI converter — U.S. Government work, as BLAST+ | `docs/tools/table2asn.md`; **not independently measured** (no licence file ships with the binary) | no |
| **perl** | 5.32.1 | runs `misa.pl` | GPL-1.0+ or Artistic-1.0-Perl | conda-forge metadata | no |
| CPStools | 3.0 | **cross-check only** — no module calls it | MIT | `docs/tools/cpstools.md` | no |
| KaKs_Calculator | 2.0.1 | cross-check only (via CPStools) | GPL-3.0-or-later | bioconda metadata | no |
| CGAS | 1.0.1 | cross-check only (dropped from the plan, D5-08) | MIT | `docs/tools/cgas.md` | no |
| PGA | 2.0 | **retired** from the toolkit 10/10/2026 (§8.56); its output survives as a fixture | GPL-3.0 | PROJECT.md §7.8 | no |
| PMGA | — | never ran here | GPL-3.0+ | `docs/tools/pmga.md` | no |

**The obvious question — "PAML and ARAGORN are GPL, how is your package MIT?"** — has one answer,
and it is the reason this table has a "Bundled" column: the toolkit **runs** those programs when
the user has installed them, as separate processes through `run_tool()`; it does not link them,
import their code, or ship them. Invoking a program is not creating a derivative work of it, so
their licences bind their own copies, not ours. The exception proves the rule: MISA *was* shipped,
and its missing licence is why it no longer is.

## 4. What changed on 11/10/2026, and what it costs

* **D-q — MISA left the package.** The trigger was a measurement, not a worry: the built wheel
  carried `misa.pl`, whose licence is not stated. The script is now obtained by whoever needs it,
  with the exact download line and both hashes in **`docs/tools/misa.md`**; the module finds it on
  `PATH` or through `--set misa=/path/to/misa.pl`, and refuses any release but the pinned 2.1 by
  reading the script's own header.
* **Cost, stated plainly.** `ssr` is now a tier-B module in the full sense: without `misa.pl` it
  **fails with a message naming the download line** (and the golden test skips with a named reason,
  as MAFFT's does). Nothing else changed, and §6's golden SSR numbers are unchanged —
  `37 = mono 31 · di 4 · tri 2` reproduced after the move, with the toolkit's own `misa.ini`.
* **The Docker image** (M9.2) obtains MISA at build time from the recorded URL with the md5
  checked — the image a user builds is their own copy. **Publishing a prebuilt image containing
  MISA would be redistribution** and needs IPK's written permission; that ask is on the board.
* **Open, and named:** written permission from IPK Gatersleben (MISA's authors) would let the
  script be bundled again. Until then, the answer above stands.

## 5. For a user who needs to know

If you install this toolkit and run `ssr`, `kaks`, `substitutions`, `synteny`, `repeats` or
`submission`, you are running **other people's programs**, which you installed yourself and whose
licences are yours to honour — the citation for each is in the run's own report
(`01_run_info/run_info.md`) and in `src/plastid/data/citations.yaml`. If you redistribute *this*
toolkit, you redistribute MIT code plus the licence notice; you redistribute none of those tools.
