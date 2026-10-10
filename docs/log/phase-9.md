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
