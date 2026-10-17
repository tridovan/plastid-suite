# Nhật ký đo lường — Giai đoạn 9 — Đóng gói, giấy phép & phát hành (M9.1–M9.6)

**Phần này là `PROJECT.md` §8** (**§8.74–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.74` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.74. Giấy phép: bản kiểm kê, và MISA ra khỏi gói (P9-01, 11/10/2026)

**Việc.** Mở Pha 9 bằng M9.3 — rà soát giấy phép **trước** mọi thứ khác, vì nó quyết định cái gì
được phép nằm trong gói phát hành. Hai câu hỏi phải trả lời được mà không cần luật sư: *bộ công cụ
này theo giấy phép nào*, và *chính xác cái gì đi kèm nó*.

**Đo được — cái gì đang được phát tán.** `uv build --wheel` trên cây trước khi sửa: **146 file**,
trong đó có **`plastid/vendor/misa/misa.pl`** — một tác phẩm của người khác, và `SOURCE.json` của
chính nó ghi `"license": "NOT STATED - the zip ships no licence file and the IPK page only requests
the citation."` **Không có giấy phép nghĩa là bảo lưu mọi quyền**: công khai để tải **không** phải
được phép phát tán lại, và bỏ file vào wheel (hay vào git, hay vào image) *chính là* phát tán lại.
Đây là ranh giới duy nhất của bản rà soát này: **gọi** một công cụ thì không phát tán nó; **đóng
gói** nó thì có.

**Đo được — giấy phép của từng công cụ**, đọc từ metadata gói conda trên máy dev
(`conda-meta/*.json`, 11/10/2026) chứ không từ trí nhớ:

| Công cụ | Giấy phép | Nguồn đo |
|---|---|---|
| MISA 2.1 | **KHÔNG NÊU** | chính file zip + trang IPK (02/10/2026) |
| MAFFT 7.526 | BSD-3-Clause | conda-forge |
| BLAST+ 2.17.0 | NCBI-PD (public domain) | bioconda |
| ARAGORN 1.2.41 | GPLv3 | bioconda (+ `docs/tools/aragorn.md`) |
| MUMmer 4.0.1 | Artistic License 2.0 | bioconda |
| PAML 4.10.10 | GPL-3.0-or-later | bioconda |
| Vmatch 2.3.1 | **`Unknown`** (upstream: chỉ học thuật) | bioconda + `docs/tools/vmatch.md` |
| KaKs_Calculator 2.0.1 | GPL-3.0-or-later | bioconda |
| perl 5.32.1 | GPL-1.0+ hoặc Artistic-1.0-Perl | conda-forge |

**Câu hỏi hiển nhiên — "PAML và ARAGORN là GPL, sao gói của mình MIT được?" — có một câu trả lời,
và nó là lý do bảng trên có cột "đóng gói?"**: không công cụ nào trong bảng **được phát tán**.
Bộ công cụ *chạy* chúng khi người dùng đã cài, như tiến trình riêng qua `run_tool()` — không link,
không import mã của chúng, không ship chúng. Ngoại lệ chứng minh luật: MISA **đã** được ship, và
giấy phép thiếu của nó là lý do nó thôi được ship.

**Đã làm (D-p · D-q, dev quyết, ghi ở §3).**

1. **`LICENSE` = MIT** (+ `license = "MIT"` · `license-files = ["LICENSE"]` trong `pyproject.toml`).
   Đo lại sau khi sửa: wheel **145 file**, `vendor/misa` **không còn**, và `METADATA` mang
   `License-Expression: MIT` + `License-File: LICENSE` (PEP 639 — hatchling chấp nhận cú pháp SPDX).
2. **MISA ra khỏi gói.** Script **không** còn nằm trong cây: nó được **lấy về**, và `docs/tools/misa.md`
   là chỗ khai báo thứ ba của CLAUDE.md — URL · md5 của zip · md5+sha256 của `misa.pl` · dòng tải
   đã kiểm (11/10/2026) · vì sao không vendor. Mô-đun tìm nó trên `PATH`, hoặc qua
   `--set misa=/path/to/misa.pl`, và **đọc release từ chính header của script** — MISA không có cờ
   version — **từ chối mọi bản khác 2.1** (cùng luật `repeat_match_version` áp cho MUMmer).
3. **`misa.ini` là file của mình** (`src/plastid/data/misa.ini`): viết lại **đúng định nghĩa mặc
   định** của MISA (`1-10 2-6 3-5 4-5 5-5 6-5`, interruptions 100, `GFF: true`) — giá trị cấu hình
   mô tả hành vi công khai của công cụ, không phải bản sao file của họ.
4. **`ssr` lên 0.4.0**: nhận `--set misa=`, và cache key theo `module_version` (H-08) là lý do đây
   là một version chứ không chỉ là một sửa đổi.
5. **`docs/licensing.md`** — bản kiểm kê đầy đủ: giấy phép của mình, cái gì nằm trong gói (kể cả
   dữ liệu dẫn xuất từ GenBank công cộng), bảng 14 công cụ, và câu trả lời cho câu hỏi GPL.
6. **Test canh gác** `tests/unit/test_licensing.py`: **đỏ trên cây cũ** — chạy thật, thông điệp
   `bundles ['misa.pl'] but its licence reads 'NOT STATED…'` — và xanh sau khi sửa; cộng một test
   buộc **mọi adapter trong `src/plastid/external/` phải có tên trong bản kiểm kê**.
7. **CI** lấy MISA 2.1 từ IPK ở job tier-B, **kiểm md5** trước khi giải nén (cùng khuôn với
   table2asn) — không có bước ấy thì test SSR vàng sẽ **skip im lặng**, đúng thứ job tier-B sinh ra
   để chặn.

**Không đổi gì về số.** Genome vàng sau khi MISA ra khỏi gói: **37 = mono 31 · di 4 · tri 2**,
compound 3 — đo lại bằng `misa.ini` của mình. `make check` xanh: **856 fast + 85 golden**
(11/10/2026, có đủ công cụ tier-B trên máy).

**Cái giá, nói thẳng.** `ssr` nay là mô-đun tier-B đúng nghĩa: không có `misa.pl` thì nó **dừng với
thông điệp nêu dòng tải** (và test vàng skip **có lý do**, như MAFFT). Đổi lại: gói phát hành chỉ
còn mã MIT của mình. **Việc còn mở, có tên**: xin IPK (tác giả MISA) **sự đồng ý bằng văn bản** —
có thì được đóng gói lại; chưa có thì giữ nguyên như trên. Image Docker (M9.2) lấy MISA lúc build
với md5 đã ghi, nên bản image **người dùng tự build** là bản của chính họ; **phát hành image dựng
sẵn có MISA thì lại là phát tán** và cần đúng văn bản ấy.

---

### 8.75. Gói phát hành: metadata, và **cài từ wheel vào venv sạch** (P9-02, M9.2)

**Việc.** Dựng artifact thật rồi **cài nó ở một chỗ không có gì của dev** — `uv` env của máy này đã
ấm từ Pha 1, nên một lần chạy trong đó không chứng minh được gì về người dùng.

**Đo được (11/10/2026).**

| | |
|---|---|
| `uv build` | `plastid_suite-1.0.0-py3-none-any.whl` **455.720 B** · `plastid_suite-1.0.0.tar.gz` **1.917.958 B** |
| Wheel | **145 file** — không có `vendor/`, có `dist-info/licenses/LICENSE`, METADATA: `License-Expression: MIT` · `License-File: LICENSE` |
| Sdist | **279 file**, mang theo `tests/` và `tests/data` (15 file) |
| Cài wheel vào venv sạch | `python3 -m venv` → **`pip` không có** (Debian thiếu `python3-venv`/`ensurepip`); dùng `uv venv --python 3.12 --seed` ⇒ **cài được, và đúng là pip cài** |
| Phụ thuộc thật của bản cài tối thiểu | `plastid-suite 1.0.0` + **biopython 1.88 · numpy 2.5.3 · PyYAML 6.0.3** — không extra nào |
| `plastid --version` (từ venv sạch) | **`plastid 1.0.0`** |
| **Số vàng từ gói đã cài**, `plastid run -m genome_stats` trên `NC_050833.1.gb` | **153.724 bp · 37,85 % · 132 = 86 CDS + 38 tRNA + 8 rRNA** · IR 26.214 bp @ **100,00 %** · LSC 83.127 · SSC 18.169 — exit 0 |
| Sdist cài được | `plastid 1.0.0` (dựng từ nguồn, cần hatchling — pip tự lấy) |

**Phiên bản: 1.0.0** (D-r). Metadata nay khai đủ: classifiers · keywords · `[project.urls]`. **Cố ý
KHÔNG có `Operating System :: Microsoft :: Windows`** — bộ công cụ không chạy native trên Windows,
đường Windows là Docker; một classifier khai điều ngược lại đúng là kiểu câu nghe hợp lý mà sai mà
dự án này sinh ra để tránh.

**Cái sdist mang, và cái nó không mang** (nói thẳng, vì "cài từ nguồn rồi chạy `make golden`" là
kỳ vọng dễ có): nó mang `tests/` + `tests/data`, nhưng **không** mang `.panel/` — các panel đã pin
nằm ngoài git theo thiết kế (`.gitignore`). Nên **nửa corpus của `make golden` chỉ chạy được trong
một checkout có panel**, còn người cài từ wheel thì chạy *phân tích*, không chạy bộ kiểm tra vàng.
Ghi ở README (P9-04) để không ai đoán.

---

### 8.76. Image phát hành, và **lần chạy golden đầu tiên bên trong container** (P9-03, M9.2 + phép kiểm trong giai đoạn)

**Việc.** Thay `docker/Dockerfile` (689 B skeleton) bằng image thật, rồi chạy **bộ số vàng ngay
lúc build** dưới `PTS_REQUIRE_TOOLS=1`: image không tái hiện được số thì **không build được**.

**Image.** `mambaorg/micromamba:2.0.5` + bộ pin của `environment.yml`/§7.8 — **python 3.12 ·
perl 5.32.1 · mafft 7.526 · mummer4 4.0.1 · paml 4.10.10 · blast 2.17.0 · aragorn 1.2.41 ·
vmatch 2.3.1** — cộng **MISA 2.1 lấy từ IPK** và **table2asn lấy từ NCBI**, cả hai **có kiểm md5**.
Các phiên bản được **in ra ngay trong image**: `mafft v7.526 (2024/Apr/26)` · `nucmer 4.0.1` ·
`blastn 2.17.0+` · `ARAGORN v1.2.41` · `perl 5.32.1` · `plastid 1.0.0` · `misa.pl` 20.148 B
(25/08/2020 — đúng bản đã pin).

**Ba lần build, và mỗi lần dạy một điều — ghi lại vì chúng là lỗi đo được, không phải sự cố vặt:**

1. **Build 1 chết ở `pip install`**: `bash: pip: command not found`. Image micromamba chỉ đưa base
   env vào `PATH` cho các bước `RUN` **khi có `ARG MAMBA_DOCKERFILE_ACTIVATE=1`** — thiếu nó thì
   `pip` không tồn tại dù python đã cài xong.
2. **Build 2 chạy golden và FAIL 3 test**: hai test cần `mkvtree`/`vmatch` và `table2asn` — hai
   công cụ tôi đã *cố ý* bỏ khỏi image ở bản nháp đầu, mà `PTS_REQUIRE_TOOLS=1` thì biến "thiếu
   công cụ" từ skip-có-lý-do thành **fail**. Đúng như thiết kế: một job tier-B im lặng bỏ qua thì
   không chứng minh gì. **Vmatch vào image** (học thuật — image là bản người dùng tự dựng). Và
   **test thứ ba mới là phát hiện thật**:
3. **Build 2 cũng lộ ra lỗi thật đầu tiên mà chỉ container nhìn thấy** — `test_the_naming_step_
   does_not_move_the_golden_ts_tv` fail với `no pyproject.toml above
   /opt/conda/lib/python3.12/site-packages/plastid/core/gene_names.py`. **`gene_normalize` không
   chạy được từ gói đã cài.** Đã sửa, kèm hai test hồi quy — chi tiết ở **§9.N-1**.
4. **Build 3 (đang chạy) — bài học về công cụ tuỳ chọn**: `curl` table2asn (56 MB từ NCBI) **qua
   VM Docker** bò rất chậm; đo từ host cùng lúc thì được ~600 kB/s, tức đường truyền không hỏng,
   chỉ chậm. Một công cụ **tuỳ chọn** (`.sqn`; D3-04 nói `.fsa` + `.tbl` là đủ) **không được phép
   treo build của người dùng**, nên nó thành **best-effort**: tải hỏng vì mạng thì build đi tiếp
   (in ra một dòng), còn **md5 sai thì vẫn chết** — và đúng **một** test của nó được
   `--deselect` **theo tên, có in lý do**, không bao giờ im lặng.

**Bằng chứng độc lập ngoài container** (chạy trên wheel cũ đã cài, trong một thư mục người dùng,
11/10/2026): `plastid run -m gene_normalize …` → `no pyproject.toml above
/tmp/pts-clean/lib/python3.12/site-packages/plastid/core/gene_names.py — cannot locate the
repository root`; sau khi sửa và cài lại wheel, cùng lệnh đó trong cùng thư mục → thông điệp nêu
`--set naming_reference=<file>`, và khi truyền đúng bản ghi thì run **ok**, artifact
`NC_024936.1.normalized.gb` + `NC_050833.1.normalized.gb` được ghi.

**Kết quả cuối (11/10/2026).** `docker build` → **exit 0**, image **1,5 GB**, và bộ golden **chạy
ngay trong lúc build**: **54 passed · 30 skipped · 1 deselected** (110,95 s). Ba mươi skip **đều là
test cần `.panel/`** — panel nằm ngoài git theo thiết kế và bị `.dockerignore` loại — và **mỗi cái
in lý do của nó** (`the pinned reference set lives outside the repo (.panel/)`); **không có skip nào
vì thiếu công cụ**, đúng vì `PTS_REQUIRE_TOOLS=1` biến loại đó thành fail.

**Và cùng một bộ số, đo trong chính image đó** (`docker run --rm plastid-suite:1.0.0 …`):

| | WSL | Docker |
|---|---|---|
| `genome_stats` bản vàng | 153.724 bp · 37,85 % · 132 = 86 + 38 + 8 · IR 26.214 @ 100,00 % · LSC 83.127 · SSC 18.169 | **y hệt** |
| `ssr` (tier B, MISA vừa lấy về) | 37 = mono 31 · di 4 · tri 2, compound 3 | **y hệt** (`MISA 2.1`, `(1/10) (2/6) (3/5) (4/5) (5/5) (6/5)`, `agrees`) |

**Một điều học được khi chạy image, ghi lại vì người dùng sẽ gặp**: `docker run … --out /runs` báo
`PermissionError: '/runs'` — image chạy bằng user `mambauser` (uid 1000), nên thư mục ghi phải là
`/tmp/…` hoặc một thư mục được mount vào. README dùng đúng lối mount (`-v "$PWD:/work"`).

---

### 8.77. Bộ tài liệu phát hành: README chạy thật · bảng mô-đun × nền tảng × RAM · "những thứ không làm" (P9-04 · P9-05 · P9-06)

**Việc.** Ba tài liệu mà người dùng cuối thật sự đọc — và nguyên tắc chung: **tài liệu phải được
chạy**, không phải được viết ra rồi hy vọng.

**P9-04 — README (`README.md`, M9.1).** Viết cho **nhà nghiên cứu, không phải dân CLI**: cài (Linux
· Windows-qua-Docker · từ checkout) → lần chạy đầu → thư mục kết quả → **mã thoát**. Đo được, trong
**venv sạch** (`uv venv --seed`, wheel 1.0.0):

| Lệnh trong README | Kết quả đo |
|---|---|
| `plastid --version` | `plastid 1.0.0` |
| `plastid run -m genome_stats -i my_plastome.gb …` | ok — run dir có `manifest.json` · `00_run_info/` (kèm `run_info.md` + `logs/`) · `01_genome_stats/`, **đúng như README vẽ** |
| `plastid run -m substitutions -i reference.gb -i sample.gb --set scope=cds` | ok (cần MAFFT — tier B) |
| một **amplicon 600 bp** (FASTA) qua `genome_stats` | **exit 2**, `IR gate: amplicon_test: no valid IR pair found` |

Mã thoát trong README nay đúng cả năm: `0` · `1` · **`2` = lỗi cú pháp *hoặc* cổng IR *hoặc* audit
bộ tham chiếu** (`cli.py: EXIT_USAGE`) · `3` thiếu công cụ/gói tuỳ chọn · `4` xung đột tên
(`EXIT_CONFLICT`, D5-09 §5). Dòng **Status** từng ghi "Pha 8 đang mở trên `feat/phase-8`" trong khi
Pha 8 **và** 8B đã đóng — nay nó chỉ nói **phiên bản**, và `tests/unit/test_project_docs.py` fail
nếu nó lệch `pyproject.toml`, nên một lần bump phiên bản phải chạm cả hai trong cùng commit.

**P9-05 — `docs/platforms.md` (M9.4).** **Cả 22 mô-đun** trong một bảng: nhóm · arity · cổng IR ·
công cụ ngoài · **chi phí đo được kèm phạm vi** (34 genome theo giao thức M8.3, hoặc 1 genome với
`repeats` — **0,8 s · 87 MB**, đo 11/10/2026) · ba mức hỗ trợ (Linux chính thức · Windows/Docker
chính thức · **macOS thực nghiệm, chưa kiểm — không có máy**). **Bảy mô-đun không có số đo** được
nêu tên chứ không im lặng (đúng ghi chú của M8.3). Máy đo nêu một lần: i7-1165G7 · ~7,8 GB.

*Phát hiện khi dựng bảng:* trường `citations` của registry **được lưu nhưng không nơi nào đọc**, và
đã trôi — `plastid_annotate`, `mito_annotate`, `gene_normalize` khai **rỗng** ở đó trong khi run
của chúng **có** trích BLAST+/ARAGORN qua `Result`. Nên cột "công cụ" của bảng được lập từ **adapter
mà mô-đun thật sự gọi** (`src/plastid/external/`), và sự trôi này được ghi ra như một phát hiện —
không sửa lặng lẽ, vì nó là *nguồn sự thật thứ hai* đang nói khác.

**P9-06 — `docs/things_we_do_not_do.md` (M9.5).** 12 mục, **mỗi mục kèm lý do đo được** và việc nên
làm thay: không Windows native · không hứa macOS · **không lắp ráp** (Pha 7 chờ reads thật) · không
đoán RNA editing nếu không có RNA-seq · không align cả panel (`scope=whole`: 8 genome 2 ph 44 s,
**16 genome không xong sau 7 phút**) · một bào quan một lần chạy · vòng ty thể hoãn · thời gian phân
kỳ/ancestral state là việc của PI · **không ship công cụ nào** · cổng IR từ chối mảnh amplicon **có
chủ đích** · chưa có web UI · chi phí panel kèm máy đo. README trỏ tới nó ngay từ đầu.

**Canh gác thêm:** `test_every_registered_module_has_a_row_in_the_platform_table` so **tập hợp**
mô-đun trong bảng với registry — thêm mô-đun mà quên tài liệu là fail.
