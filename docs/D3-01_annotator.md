# D3-01 — which annotator, and why the answer is not a preference

**Question:** Phase 3 annotates plastomes. GeSeq produced the group's own published records, but
PROJECT.md §7.8 lists PGA and plann.pl as the offline substitutes. Which one is the pipeline?

**Answer: the offline tool is the pipeline; GeSeq is a one-off cross-check whose output gets pinned.**
The naming problem is handled by D3-02 either way, so it does not decide this.

> **Update, 2026-10-03 — measured, and the offline tool changed.** The constraint below still
> decides *offline vs web*; it never said *which* offline tool, and that part was a hypothesis
> until it was run. PGA, run on the golden genome, returned 147 genes where the record has 132 and
> failed silently three ways (PROJECT.md §8.15). **Chloe replaces it as the pipeline annotator:**
> offline, its own curated references, **132 = 86 + 38 + 8** with every CDS translating (§8.17,
> `docs/tools/chloe.md`). PGA is kept for the genomes Chloe does not cover — it annotates
> flowering plants only. GeSeq's job is unchanged.

> **Update, 2026-10-05 — the annotator changes again, and this time it is ours.** The team's
> decision (D5-08, PROJECT.md §3): the annotate step becomes this project's own tool — CGAS's flow
> with our search-and-decide core — and **Chloe is dropped from the pipeline**. This file stays as
> the record of how the 03/10 decision was measured (PGA 147 genes, §8.15; Chloe 132 = 86 + 38 + 8,
> §8.17). The new engine's design, its reference-set audit and its Done-when live in
> `docs/D5-08_plastid_annotator.md`, its phase is **3B** in `ROADMAP.md`, and `--set
> annotator=chloe|pga` in `modules/annotation.py` remains the current code until that engine passes
> its Done-when.



## Why this is a constraint, not a preference

The toolkit has to run on an ordinary computer, offline, and ship as Linux native + Windows via
Docker (CLAUDE.md, non-negotiable product requirements). Every run must record the tool's **name ·
version · path · argv · exit code** — that is seam 3, and it is why `run_tool()` exists at all.

A web service satisfies none of that. It cannot go in the Docker image, it cannot run without a
network, and it has no version number to record: a run's manifest could only say "GeSeq, sometime".
For a project whose stated failure mode is *numbers that look plausible and are wrong*, an
unversioned step in the middle of the pipeline is the same class of problem.

So the choice is not "which tool is better". It is: **the pipeline must be offline, and then the
only open question is which offline tool.**

## The options, measured

| | GeSeq (CHLOROBOX) | PGA | plann.pl |
|---|---|---|---|
| Runs offline | **no** — web only, no CLI | yes | yes |
| Ships in the Docker image | no | yes | yes |
| Version recordable in the manifest | **no** | yes (git commit / release) | yes |
| Licences | web service | GPL-3.0 | check at P9 |
| Needs, beyond itself | a browser and a network | Perl + BLAST+ | Perl |
| Reference | built-in MPI-MP curated CDS set | **user supplies a GenBank** | user supplies a GenBank |
| Output | GenBank, batch ("GLOBAL") | GenBank + log, batch | GenBank |
| Measured on this machine | not runnable offline | **ABSENT** | **ABSENT** |

**Recommended at the time: PGA primary, plann.pl as the fallback.** PGA is a local command line
tool, batch capable, GPL-3.0, and it already emits GenBank — which our reader handles and which the
phase check requires us to read back. *(Superseded above: PGA was run, and it did not reproduce
the golden record. The table below is what was known before that measurement — Chloe was not in
the survey at all, which is the part worth remembering: the tool that won was one nobody had
listed.)*

### PGA's known traps, and why this repo is already immune to them

PGA's own guide warns about three things. All three are already impossible here:

| PGA's trap | Why it cannot bite us |
|---|---|
| the output folder must be empty, it may delete its contents | every run gets a fresh numbered directory (`RunLayout`), never reused |
| filenames containing spaces cannot be read | inputs go through `slugify()`, which is ASCII-safe by design (§9.30) |
| one sequence per file, a multi-FASTA target can annotate wrongly | one `--input` is one record, enforced by the arity check |

## GeSeq still has a job

Once per manuscript species, someone runs GeSeq by hand and we compare its gene set with ours. That
is cheap (minutes), it is what §13.3's in-phase checks ask for in spirit, and it answers a question a
reviewer can legitimately ask: *does your annotation agree with the tool that produced your own
submitted records?*

Two ways to use the result, and the second is better:

1. Read it as a one-off check and record the date and the agreement in PROJECT.md.
2. **Pin it.** `plastid panel pin` already records accession · date · sha256 and re-reads the file
   with our own reader. A GeSeq GenBank pinned that way becomes an offline, reproducible input —
   the manual web step happens once, and every later run is offline and byte-identical.

Option 2 is the reason this decision is cheap: **the toolkit does not have to reproduce GeSeq, it
has to be able to consume GeSeq's answer.**

## The naming question does not decide this

GeSeq's MPI-MP reference set is what produces `pafI` · `pafII` · `pbf1`, so a GeSeq-annotated genome
arrives already in the group's dialect. PGA inherits whatever names its reference carries. Neither
outcome changes the pipeline: D3-02 (PROJECT.md §8.10) normalises either dialect to the group's
names and keeps the original as `source_name`. **The normaliser is required regardless**, so this
decision can be made on the engineering constraint alone.

## What this costs

- two more tools to pin under §7.8 (Perl — already present for MISA — and BLAST+)
- a **reference GenBank per genome**, which is itself a dataset: pin it with `panel pin`, and prefer
  a reference with at least as many genes as the target (PGA's own guidance)
- one manual GeSeq run per manuscript species, which cannot be automated and cannot go in CI

## What would change this answer

A container image or a documented, versioned CLI for GeSeq. The check was made on 2026-10-03: GeSeq
is a CHLOROBOX web service with no standalone distribution; its "additional command line parameters"
field is a web-form option, not an executable. Re-check before Phase 9 packaging.

## [ĐÃ THAY 10/10/2026] — dòng cập nhật của D3-01

Bảng quyết định ghi *"chloe · pga retires when P3B-10 passes; D3-01 gains its update line"*. **P3B-10
đạt 10/10/2026** (§8.55: bản vàng 132 = 86 + 38 + 8 với 35 khác biệt được liệt kê; corpus 418 so
với PGA 588 — **đọc strict là 1070 so với 699, không đạt**), và **P3B-11 đã thi hành**: hai
annotator ngoài cùng cái bọc `annotation` và `core/gff.py` đã bị xoá khỏi cây làm việc (§8.56).
Đường chú giải lục lạp bây giờ là **`plastid annotate`** — động cơ tự viết theo D5-08.

**Phép đo cũ vẫn nguyên giá trị và vẫn ở đó:** §8.15 (PGA trên bản vàng: 147 gen, ba kiểu hỏng im
lặng) · §8.17 (Chloe qua reader của mình: 132 = 86 + 38 + 8, 0 CDS lỗi ORF — chỉ đúng khi
`--no-transform`). **PGA sống tiếp trong repo** ở `tests/data/corpus/*.pga.gb` + `tools/corpus_bar.py`
— vai "đối chứng đã đo" của D5-08, và là nửa còn lại của thanh chuẩn.

**Nếu một ngày phải chạy lại Chloe** (ví dụ để đo trên một loài mới): tool là
`github.com/…/chloe` (Chloë: Organelle Annotator) **0.1.18**, chạy bằng Julia, cài đặt và bản đã
đo nằm trong lịch sử git của `docs/tools/chloe.md` và `src/plastid/external/chloe.py`; điều duy nhất
phải nhớ là `--no-transform` — không có nó Chloe lật vùng sao chép đơn và ghi toạ độ trong một hệ
quy chiếu mà đầu vào chưa từng có (§8.17). Con số cần tái lập là **132 = 86 + 38 + 8**.
