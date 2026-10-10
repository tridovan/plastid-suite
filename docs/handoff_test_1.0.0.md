# Bàn giao test — PlastidToolSuite 1.0.0 (11/10/2026)

*Tài liệu này viết bằng tiếng Việt vì người đọc là cô Trang và anh Điềm, không phải vì nó là tài
liệu kỹ thuật của repo — phần kỹ thuật của repo vẫn tiếng Anh (`docs/`).*

## Đây là gì

**Bản ứng viên 1.0.0** của bộ công cụ — **chưa phát hành**, chưa có DOI, chưa có trên PyPI. Mục
đích của đợt này: **chạy nó trên máy người khác**, vì tới giờ mọi con số đều đo trên **một máy**
(máy dev). Một bộ công cụ chỉ chạy được trên máy người viết ra nó thì chưa phải là bộ công cụ.

Giấy phép **MIT** (dùng tự do, kể cả sửa); không công cụ ngoài nào được đóng gói kèm — xem §6.

## 1. Ai test phần nào

| | Việc phù hợp |
|---|---|
| **Anh Điềm** | Chạy thật trên WSL/máy anh: cài, chạy phân tích, và nếu được thì chạy **bộ kiểm tra đầy đủ** (mức 3). Anh là người duy nhất ngoài máy dev có đủ công cụ (MAFFT, BLAST, PAML…) |
| **Cô Trang** | Đọc **báo cáo** của một lần chạy (mức 1) và nói xem bảng với câu chữ trong đó **có dùng được cho bài báo** không, thiếu gì, thừa gì. Đây là câu chỉ cô trả lời được — không phải việc cài đặt |

## 2. Mức 1 — 5 phút, không cần công cụ gì (ai cũng làm được)

Trong thư mục này có `plastid_suite-1.0.0-py3-none-any.whl` và `NC_050833.1.gb` (một plastome
công khai, dùng làm mẫu chuẩn của dự án).

```bash
python3 -m venv .venv-test          # nếu báo lỗi ensurepip: uv venv .venv-test --seed
.venv-test/bin/pip install plastid_suite-1.0.0-py3-none-any.whl
.venv-test/bin/plastid --version    # phải in ra: plastid 1.0.0

.venv-test/bin/plastid run -m genome_stats -i NC_050833.1.gb --out .runs --dataset iris
grep -v '^#' .runs/iris/*/01_genome_stats/tables/genome_stats.csv
```

**Số phải ra** (nếu ra khác, đó là lỗi — xin gửi lại nguyên văn):

```
153724 bp · 37.85 % GC · 132 gen = 86 CDS + 38 tRNA + 8 rRNA
IR 26214 bp @ 100.00 % · LSC 83127 · SSC 18169
```

Rồi mở `report.md` trong cùng thư mục run — **cô Trang đọc file này**: nó là bản nháp kết quả
(draft — needs human review), viết theo lối Materials & Methods.

## 3. Mức 2 — có công cụ ngoài (anh Điềm)

Hai phân tích cần công cụ. **MISA không còn được đóng gói kèm** (zip của nó không ghi giấy phép,
nên bên em không được phát tán lại) — lấy bằng đúng hai lệnh này, có kiểm md5:

```bash
curl -sSLO https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip
md5sum misa_sourcecode_25082020.zip     # phải là: cf9d3579f7caacd02f99830b329c7d9c
unzip -j misa_sourcecode_25082020.zip misa.pl -d "$HOME/misa"
md5sum "$HOME/misa/misa.pl"             # phải là: cfebc1d1a36fddfd6e0d897c96babd12
```

```bash
# SSR (cần MISA + perl)
.venv-test/bin/plastid run -m ssr -i NC_050833.1.gb --out .runs --dataset iris \
    --set misa="$HOME/misa/misa.pl"
# phải ra: 37 SSR = 31 mono + 4 di + 2 tri, compound 3

# So sánh hai bộ gen (cần MAFFT) — đây là phân tích bản thảo dùng
.venv-test/bin/plastid run -m substitutions -i NC_050833.1.gb -i NC_024936.1.gb \
    --out .runs --dataset pair --set scope=cds
```

Anh xem luôn `00_run_info/run_info.md` trong thư mục run: nó ghi **tên, phiên bản, đường dẫn và
nguyên dòng lệnh** của mọi công cụ đã chạy — nếu máy anh ra số khác máy em, file đó là chỗ đối
chiếu đầu tiên.

## 4. Mức 3 — bộ kiểm tra đầy đủ (anh Điềm, nếu muốn đi xa hơn)

Gói `plastid_suite-1.0.0.tar.gz` **mang theo cả mã kiểm tra và dữ liệu vàng**:

```bash
tar xzf plastid_suite-1.0.0.tar.gz && cd plastid_suite-1.0.0
python3 -m venv .venv && .venv/bin/pip install -e . pytest
.venv/bin/pytest -q -m "not tier_b"     # 862 test nhanh, không cần công cụ
PTS_REQUIRE_TOOLS=1 .venv/bin/pytest -q tests/golden -rs   # bộ số vàng, cần công cụ
```

(Bộ đầy đủ cần thêm các panel đã pin nằm ngoài git — những test đó tự bỏ qua **có nêu lý do**,
không im lặng.) Nếu anh muốn clone cả repo thì cần được cấp quyền — repo đang private.

## 5. Gửi lại gì

* **Chạy được**: một dòng "mức X xong, số đúng" là đủ.
* **Không chạy được**: **nguyên văn thông báo lỗi**, mã thoát (`echo $?`), phiên bản Python
  (`python3 -V`) và hệ điều hành. Nếu có thư mục run thì gửi cả `manifest.json` +
  `00_run_info/` — nó ghi đủ để bên em dựng lại đúng tình huống.
* **Chậm hay tốn RAM** bất thường: nói giúp máy cấu hình gì, chạy mất bao lâu.
* **Cô Trang**: nhận xét về `report.md` — chỗ nào khó hiểu, thiếu bảng gì, câu chữ có dùng được
  không.

## 6. Những thứ **đã biết là chưa có** (đừng mất thời gian tìm)

| Chưa có | Vì sao |
|---|---|
| **Lắp ráp bộ gen từ reads** | Chưa làm — đang chờ **reads thật** (dữ liệu chưa có). Đầu vào là file GenBank/FASTA đã lắp |
| **Ty thể** ngoài 2 mô-đun thử nghiệm | Vòng ty thể hoãn có chủ đích; vòng này chỉ làm lục lạp |
| **Chạy native trên Windows** | Không hỗ trợ; đường Windows là **Docker** (có `docker/Dockerfile` trong gói .tar.gz) |
| **macOS** | Chưa từng chạy thử — không có máy. Cột này ghi "thực nghiệm, chưa kiểm" |
| **Dự đoán RNA editing** | Cần dữ liệu RNA-seq; gọi từ DNA là đoán |
| **Công cụ ngoài đóng gói kèm** | Không công cụ nào (giấy phép: MISA không nêu, Vmatch chỉ học thuật) — mỗi cái có dòng cài trong `docs/tools/` |
| **Giao diện web** | Hoãn tới khi lõi ổn định |

Danh sách đầy đủ (tiếng Anh): `docs/things_we_do_not_do.md` trong gói .tar.gz.
Số đo hiệu năng theo từng mô-đun (RAM, thời gian, máy đo): `docs/platforms.md`.

## 7. Lưu ý về dữ liệu

Chạy trên **bản ghi công khai** (GenBank) thì kết quả là **dữ liệu thay thế (stand-in)** — dùng để
kiểm tra công cụ, **không phải số của đề tài**, và không đưa vào bản thảo. Khi có mẫu thật thì
chạy lại, không vá số.

## 8. Kiểm tra file nhận được

```bash
sha256sum -c SHA256SUMS
```

| File | sha256 |
|---|---|
| `plastid_suite-1.0.0-py3-none-any.whl` | `a9be916741006ec9a4bb720467e598e5c2d9f94bdbb0910c4d488e783063371c` |
| `plastid_suite-1.0.0.tar.gz` | `28c05e09a673609443cc86e07edbeef22fa6d2e53308e331adfbc4701b1fd563` |
| `NC_050833.1.gb` | `a946f23e6e5f6a886eea93c52d72497209531aaf855ec08853dbef7997638493` |

Bản dựng này lấy từ nhánh `feat/phase-9`, commit `003cfa1`, và đã qua: 862 test nhanh + 85 test
vàng trên máy dev · CI xanh · image Docker dựng được và ra **đúng bộ số trên**.
