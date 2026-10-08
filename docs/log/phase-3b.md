# Nhật ký đo lường — Giai đoạn 3B — annotator lục lạp tự viết (D5-08) & bước chuẩn tên (D5-09)

**Phần này là `PROJECT.md` §8** (**§8.40–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.40` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.40. Bộ tham chiếu lục lạp: 34 bản ghi Iridaceae được pin lại thành một bộ riêng, và luật leave-one-out viết thành dữ liệu (P3B-01, 08/10/2026)

**Việc.** Dựng bộ tham chiếu cho annotator lục lạp (D5-08, D-b) và pin nó: `.panel/plastid-reference/`
— 34 bản ghi, `PINNED.json` (accession · sha256 · md5 · byte · ngày pin) + `REFERENCE_SET.json`
(bộ + luật + thiết kế leave-one-out). Không phải dữ liệu mới: **chính 34 bản ghi của bộ test
anh Điềm** (`.panel/iridaceae-34/`, §8.31), pin lại thành một bộ riêng vì lý do ở dưới.

**Vì sao là họ Iridaceae.** Luật taxonomy của dự án là **cùng họ trước** (bộ `mito-reference` đã
theo đúng luật này). Genome hợp đồng là *Iris domestica* `NC_050833.1` (họ Iridaceae) và **bộ trả
lời** cũng là Iridaceae — 34 plastome có bản chú giải công bố, thứ duy nhất khiến một kết quả
chuyển chú giải kiểm được. Bộ tham chiếu là **dữ liệu theo họ**: mẫu thuộc họ khác sẽ được pin một
bộ theo họ đó, cùng luật, không phải chọn lại mỗi run.

**Số đo.** 34/34 bản ghi qua cổng IR; tổng **5.204.491 bp**; CDS mỗi bản ghi **84–87** — tức
`132 = 86 + 38 + 8` là **số của genome vàng, không phải hằng số**, và một bộ tham chiếu mà các
thành viên bất đồng về sự tồn tại của một gen chính là thứ mà bước kiểm số 3 sinh ra để in ra.
`PX516262.1` (*Eleutherine bulbosa*) là bản nộp GenBank, không phải RefSeq — vẫn nằm trong bộ vì
bộ test có nó. **34/34 sha256 của bộ tham chiếu trùng khít bộ test** (`test_the_reference_set_holds_the_corpus_records_byte_for_byte`):
pin lại không được đổi byte, vì bộ trả lời thuộc về các file đó.

**Luật leave-one-out, viết thành dữ liệu.** Bộ test chứa chính chuỗi của genome vàng, và bản chú
giải PGA của nó đọc **85 CDS** nơi hợp đồng ghi 86 (§8.31) — để bộ trả lời nằm trong bộ tham chiếu
thì annotator chỉ đang tự chấm bài mình. Nên: mỗi bản ghi trong bộ được chuyển chú giải bởi **chính
bộ đó trừ nó ra** (33 bản ghi), và `leave_one_out` là luật ấy ở dạng đọc được bằng máy
(`core/plastid_reference.py`). Bản cam kết `src/plastid/data/plastid_reference.json` để thiết kế đọc
được **mà không cần 28 MB dữ liệu**, kèm `plastid_reference.md` ghi nguồn, luật và cách dựng lại.
Một đích **ngoài** bộ thì được chuyển chú giải bởi cả bộ — không có gì của nó ở trong để phải loại.

**Một điều phải nói rõ vì nó đi vào bài báo nếu bị bỏ quên.** Đây là bộ test của dự án
(D-b: "bộ test là để kiểm, không bao giờ là production"). Nó làm bộ tham chiếu ở đây vì genome hợp
đồng thuộc Iridaceae và không có nguồn Iridaceae nào khác được pin; **không con số nào đo trên bộ
này được vào bảng bản thảo của mẫu thật.**

**Hiện thực.** `core/plastid_reference.py` (luật + kiểm tra pin) · `plastid reference verify`
(băm lại từng bản ghi, in `file — reason`, thoát 1 khi bộ đã đổi) · `plastid reference loo` (sinh
manifest). Một chi tiết của đường dẫn: `pin_files` ghi lại **đúng đường dẫn nó được đưa**, nên một
pin tạo bằng `--out` tương đối mang đường dẫn tương đối — `verify` tìm theo đường dẫn ghi trong pin
trước, rồi tới chính thư mục pin, nên nó chạy đúng từ bất kỳ thư mục nào (có test riêng).

---

### 8.41. Ba con số của D5-09 §9.3 — đo trước khi viết `gene_normalize` (P3B-02, 08/10/2026)

**Việc.** D5-09 §9.3 để mở ba con số mà bước chuẩn tên phải có trước khi thành code: **biên** của
phép thử 3, **phép thử 4 (hàng xóm)** có tách được `rpl2`/`rpl23` và hai `trnL` không, và
**`identity_floor`**. Cả ba đo trên dữ liệu đã có trong repo + bộ test 34 bản ghi.

**Cách đo (ghi lại vì con số chỉ có nghĩa kèm cách đo).** Tham chiếu = genome vàng
`NC_050833.1` (**79 locus CDS** sau khi gộp tên canonical); truy vấn = 33 bản ghi còn lại của bộ
test, **bản chú giải công bố của chính bản ghi** (`/translation`, không tự dịch — §8.21).
Công cụ **BLAST+ 2.17.0**, `blastp`, `-seg no`, E-value 10, `-max_target_seqs 100` — đúng chương
trình và cờ mà động cơ D5-08 dùng. `pident` là độ tương đồng, `coverage` = chiều dài alignment /
chiều dài truy vấn. Riêng tRNA (quá ngắn để BLAST tìm thấy seed — đo được: `blastn` **0 hit** giữa
hai `trnL`) thì so bằng alignment global, match +1 · mismatch 0 · mở gap −2 · nối gap −0,5, và
identity = matches / max(chiều dài hai chuỗi). **Tham số alignment là một phần của con số**: cùng
cặp trnG, bộ tham số rẻ gap (−0,5/−0,1) cho **78,87 %**, bộ trên cho **73,24 %**.

#### 1. Biên của phép thử 3 — đo bằng **bit**, không bằng điểm identity

2.838 lượt truy vấn CDS liên bản ghi. Mọi truy vấn đều có hit. Biên = điểm của hit tốt nhất trừ
điểm của hit tốt nhất **thuộc locus khác**:

| | giá trị |
|---|---|
| nhỏ nhất | **19,2 bit** (`rps19`, `NC_081931`) |
| percentile 5 | 53,6 bit |
| trung vị | **298,0 bit** |
| dưới 5 bit | **0** |
| dưới 20 bit | 2 |

**Biên tính bằng điểm identity là vô nghĩa** — đo được trên chính bộ này: truy vấn `psaB` cho
hit `psaA` với pident **cao hơn** (alignment cụt, bitscore thấp hơn), nên "biên" ra **−418 điểm**.
Điểm xếp hạng phải là **điểm có trọng số chiều dài** (bit), và ngưỡng biên phải đi kèm ngưỡng
coverage. Đề nghị cho P3B-04: **biên ≥ 10 bit** — một nửa giá trị nhỏ nhất đo được, đủ để từ chối
một thế hoà và không từ chối gì khác.

**Ba ca xếp hạng sai, cả ba đều là lỗi TÊN chứ không phải lỗi chuỗi:** `NC_042827` gọi `lhbA` một
gen mà chuỗi của nó **trùng khít `psbZ`** của tham chiếu; `NC_072574`/`NC_072575` gọi `ycf4` một
CDS **132 bp** có `/product` = "photosystem biogenesis factor 1" — tức `psbN` (43 aa). Chuỗi luôn
đúng; tên bản ghi mới là thứ sai. Đây đúng là loại ca mà bước chuẩn tên sinh ra để **in ra**.

**tRNA thì ngược lại — xếp hạng theo chuỗi KHÔNG đủ.** 1.182 lượt so tRNA liên bản ghi:
**26 lượt (2,2 %) xếp một locus khác lên đầu**, tất cả đều là truy vấn họ **`trnS`** (hit tốt
nhất là một tRNA serine khác), lệch tới **15,05 điểm**. Và hai dải **chồng nhau**: locus đúng có
thể chỉ **53,76 %** trong khi locus khác đạt **73,97 %**.

#### 2. Phép thử 4 (hàng xóm) — lý do D5-09 viết ra KHÔNG tái lập được; lý do thật thì có

| Cặp D5-09 §4 nêu | Đo được | Kết luận |
|---|---|---|
| `rpl2` vs `rpl23` | protein **22,03–44,44 %** (273 aa × 93 aa), 34/34 bản ghi | **không** gần giống nhau; xếp hạng tách được, 0 ca xếp sai |
| hai `trnL` (`trnL-UAA` vs `trnL-CAA`) | **61,18 %** (vs `trnL-UAG` 64,71 %) | **không** gần giống nhau |

Ba loại ca **thật sự** cần phép thử 4, đo được:

1. **họ serine** — 26/1.182 lượt như trên;
2. **hai gen mang MỘT tên** — chính genome vàng: nó chú giải **hai** feature tRNA là `trnG-UCC`
   (8889-8911/9586-9633 và 36413-36483), hai chuỗi ghép **73,24 %** giống nhau. Một tên, hai gen:
   bước chuẩn tên không được gộp, và ★ thừa hưởng vấn đề này;
3. **tên không phân giải được** — 33 bản ghi chứa **72 lượt tên tRNA** và **1 tên CDS** (`lhbA`)
   mà từ vựng của tham chiếu không có. **44 trong 72** lượt tRNA có chuỗi **trùng khít** một chuỗi
   tham chiếu dưới tên khác (`trnL-CAU`=`trnI-CAU` 12 · `trnL-GUA`=`trnI-GAU` 10 ·
   `trmA-UGC`=`trnA-UGC` 10 · `trnG-GCC`=`trnG-UCC` 6 · `trnL` trần =`trnL-UAG` 5 ·
   `trnA`=`trnA-UGC` 1); 28 lượt còn lại là chuỗi tham chiếu **không mang**.

**H-04 được đo lần đầu, và nó không phải câu hỏi chính tả.** Hai tên `trnG-GCC`/`trnG-UCC` là
**hai gen**, không phải hai cách viết của một gen — và các bản ghi công bố **dùng hai tên đó ngược
nhau**: chuỗi thứ nhất của genome vàng được *Gladiolus* (`NC_062351`) gọi `trnG-GCC` còn 8 bản ghi
gọi `trnG-UCC`; chuỗi thứ hai được **25 bản ghi** gọi `trnG-GCC` nhưng **genome so sánh
`NC_024936.1`** và 7 bản ghi khác gọi `trnG-UCC`. Chuỗi **không** phân xử được cách viết — đây là
việc phải hỏi người (cô Trang), và là ca `stop` của D5-09 §5, không phải ca tự sửa.

#### 3. `identity_floor` — và vì sao ngưỡng coverage của đường mito không đủ

| Ngưỡng coverage | locus đúng thấp nhất | locus KHÁC cao nhất (bỏ 3 ca lỗi tên) |
|---|---|---|
| ≥ 0,5 | 67,75 % | **52,38 %** (`psbJ` vs `matK`, 21 aa) |
| **≥ 0,7** | **67,75 %** | **44,93 %** (`psaB` vs `psaA`, 750 aa) |
| ≥ 0,8 | 67,75 % | 44,93 % |
| ≥ 0,9 | 69,23 % | 44,93 % |

Đường mito đã đo `identity ≥ 50 %` + `coverage ≥ 0,5`; **với plastid thì ngưỡng coverage 0,5 thủng
một lỗ**: `psbJ` chỉ dài **36 aa**, nên một alignment 21 aa (coverage 0,53) đạt **52,38 %** với
`matK` — 13 trong 16 hit sai trên 50 % là đúng cặp `psbJ`/`matK` này. Gen CDS của plastome ngắn
(`psbJ` 36 · `psbL` 38 · `psbM` 34 aa…), nên **coverage phải cao hơn**:

> **Đề nghị cho P3B-04/P3B-07: `identity ≥ 50 %` VÀ `coverage ≥ 0,7`.** Khoảng trống đo được khi đó
> là **(44,93 % · 67,75 %)** — ngưỡng 50 % nằm trong đó, và ba hit vượt ngưỡng (97,67–100 %) là
> **cùng một locus dưới tên bản ghi viết khác** (`lhbA`=`psbZ`, `ycf4`=`psbN`), tức ca "đề nghị đổi
> tên", không phải ca "dưới ngưỡng".

**tRNA thì KHÔNG có ngưỡng identity nào tách được** — hai dải chồng nhau (đúng: 53,76–100 % · khác:
tới 73,97 %). Quyết định cho tRNA phải đến từ **tên/anticodon/hàng xóm**, không từ một ngưỡng số.

**Kiểm chứng.** Ba phép đo trên được ghim thành test (`tests/golden/test_naming_numbers.py`,
tier-B: cần `blastp` + bộ test đã pin): biên nhỏ nhất 19,2 bit · trung vị 298,0 · locus đúng thấp
nhất 67,75 % · locus khác cao nhất 44,93 % (coverage ≥ 0,7) · 3 ca lỗi tên 97,674/100/100 ·
26 lượt tRNA xếp sai, tất cả họ `trnS` · locus đúng thấp nhất 53,76 % · cặp trnG 73,24 %.

---

### 8.42. Khoá cache mang thêm bốn thứ — và bốn hành vi cũ được chứng minh là sai (P3B-05, 08/10/2026)

**Việc.** Khoá cache của run layer (`run/cache.py` + `run/runner.py`) thiếu đúng những thứ quyết
định kết quả. Sửa bốn chỗ, mỗi chỗ là một việc mở đã đo trước đó: **H-15** (chú giải + thứ tự
input) · **H-08** (mã nguồn module) · **H-14** (`--docx`) · **D5-09 §7.2** (sha256 của bản ghi ★).

**Bốn hành vi cũ, chứng minh bằng cách nạp lại chính `cache.py` của HEAD và chạy cạnh bản mới:**

| Việc | Khoá cũ | Khoá mới |
|---|---|---|
| hai bản chú giải của **một chuỗi** (`Iris_domestica` công bố vs PGA, 86 vs 85 CDS) | `_input_hash` **giống nhau** (chỉ băm chuỗi) | khác nhau (băm cả gen, cả `source`) |
| thứ tự input (`ab` vs `ba`) | khoá **giống nhau** (`sorted`) | khác nhau |
| một run xin `--docx` | khoá **không có** khái niệm này | `payload_kinds` vào khoá |
| module đổi hành vi mà không bump version | khoá **không có** mã nguồn | `module_source` = sha256 file nguồn |
| `--set naming_reference=<file>`, rồi sửa file đó | chuỗi đường dẫn là hằng số ⇒ khoá **không đổi** | tham số trỏ tới file được **băm theo nội dung** |

Hai dòng đầu là **H-15**: sáu module lấy `genomes[0]` làm tham chiếu, nên một run `kaks` với hai
input đảo thứ tự bị phục vụ bảng của run kia. Dòng thứ tư là **H-08**. Dòng thứ ba là **H-14** —
đo 05/10: run thứ hai trả `cached` với **0 file `.docx`**, manifest cũng không ghi cờ đó. Dòng
cuối là **D5-09 §7.2** — điều kiện để bản ghi ★ vào được khoá cache, và giờ nó vào bằng một luật
chung: **tham số nào trỏ tới một file đang tồn tại thì được băm theo nội dung file**, nên sửa file
là đổi khoá; tham số không phải đường dẫn thì băm như chính nó.

**Hệ quả phải nói rõ.** Khoá đổi ⇒ **mọi cache cũ không còn dùng được** (`.runs/.cache`), lần chạy
kế tiếp sẽ tính lại chứ không phục vụ payload cũ — đúng hướng an toàn, nhưng là một lần chậm hơn.

**Kiểm chứng.** Bốn test hồi quy mới trong `tests/unit/test_run.py`, mỗi test đỏ trên hành vi cũ
(chứng minh ở bảng trên bằng cách chạy code HEAD cạnh code mới): hai bản chú giải của một chuỗi cho
hai khoá · khoá giữ thứ tự input · tham số trỏ file được băm theo nội dung · đổi mã nguồn module
thì đổi khoá · run xin `--docx` không được phục vụ payload không có tài liệu.
