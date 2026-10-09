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

---

### 8.43. Cổng kiểm bộ tham chiếu: cái gì chặn, cái gì chỉ báo, và lần đọc thứ hai bằng MAFFT (P3B-03, 08/10/2026)

**Việc.** Dựng cổng audit của D5-08 (`core/reference_audit.py`) và cửa dòng lệnh
`plastid reference audit --panel <bộ> --out <runs>`: một **run** thật (thư mục run · manifest ·
run_info · log công cụ), trong đó mỗi phát hiện in ra đúng dạng `file — reason`, không kèm lời
khuyên nào. Bộ bị chặn ⇒ **exit 2, không chú giải gì**, lý do nằm trong `manifest.json`
(`modules[0].error`). Bộ chỉ bị báo ⇒ exit 0, run đi tiếp.

**Luật, và số đo đã chọn luật ấy.** Sáu lớp **CHẶN**: không đọc được · không có feature gen ·
CDS thiếu `/translation` · toạ độ partial · cổng IR · CDS không dịch được. Hai lớp **BÁO**: chiều
SSC lệch đa số · bộ tên gen lệch đa số. Đo trên chính bộ đã pin (§8.40), 34 bản ghi, 08/10/2026:

| Lớp | Số đo trên 34 bản ghi |
|---|---|
| thiếu `/translation` · ORF · cổng IR | **0 · 0 · 0** |
| toạ độ partial | **1 bản ghi**, 2 feature: `NC_024936.1` — `CDS rps19` (`<153124`, `>83020`) |
| chiều SSC | **32 "same" · 2 "opposite"** — `NC_042827.1` (*I. missouriensis*) · `NC_062351.1` (*Gladiolus*), y hệt §8.34 |
| bộ tên gen | **20/34 lệch** bộ tên phổ biến nhất (bộ ấy do 14 bản ghi dùng) |

Ca `NC_024936.1` là đúng ca §8.34 đã dự báo: **một bản công bố thật** (và là genome so sánh
*Iris gatesii*) bị chặn, vì một CDS mang toạ độ partial sẽ đóng dấu partial lên bất cứ thứ gì nó
chú giải — §8.15 bẫy 1 đo PGA trả **0 gen** từ chính bản ghi này. Bộ tham chiếu vì thế phải **bỏ
hoặc thay** bản ấy; một lệnh ngoại lệ không cứu được (xem dưới).

**Lần đọc thứ hai, và nó độc lập thật.** Luật chiều SSC so SSC với LSC **bên trong một bản ghi**,
nên nó không tự kiểm được mình. MAFFT `--adjustdirection` chạy trên **vùng SSC của cả bộ** (34
chuỗi · 614.658 bp · **17,6 s** trên 7.526) và đánh dấu `_R_` cho chuỗi nào nó phải đảo chiều —
đo 08/10/2026: **đúng hai bản ghi ấy**, không hơn không kém, và cả ba dòng của bảng
`region_orientation` đều `Agree = yes`. Hai phương pháp, một câu trả lời; khi chúng lệch nhau,
audit **báo là lệch** chứ không chọn bên thắng. Không có MAFFT thì dòng ấy ghi **"not measured
(mafft 7.526 not found)"** — không bao giờ tự nhận là đã xác nhận.

**Bức hình.** Dotplot (`result/figures/dotplot.py`, mục thứ mười của danh mục hình) vẽ mọi run
khớp 15-mer của bản ghi lệch chiều với một bản ghi đa số, **cả hai chiều**. Đo trên bộ thật: bản
`NC_042827.1` cho **349 run trên đường chéo nghịch** (đối chiếu: 20 run xuôi), `NC_062351.1` cho
**332** (13 xuôi) — tức hình đọc ra đúng điều luật nói, và đây là thứ người dùng nhìn trước khi
xác nhận một flip-flop. Hình vẽ **từ bảng** `region_matches` (716 dòng), không tính lại gì, nên
`plastid figures <run-dir>` dựng lại được nó từ `result.json` mà không chạy lại audit.

**Ngoại lệ đã xác nhận thì buộc vào byte.** `AUDIT_EXCEPTIONS.json` nằm cạnh bộ tham chiếu; mỗi
mục ghi `file · sha256 · check · reason · confirmed_by · confirmed_on`, và **sha256 do công cụ
viết** (`--confirm-exception FILE:CHECK --confirmed-by … --reason …`), không phải người gõ tay.
Đo end-to-end trên bộ thật: pin cho `NC_042827.1:orientation` ⇒ run sau báo **21** thay vì 22 phát
hiện, dòng ấy thành `confirmed` và **vẫn hiện** trong bảng; sửa một byte của file ⇒ pin **void**,
phát hiện quay lại. Một pin **không bao giờ** miễn được lớp CHẶN — `Finding` từ chối dựng, và
`read_exceptions` từ chối đọc (`translation` không nằm trong `REPORT_CHECKS`). Bài học suýt nữa thì
trả giá: lần chạy đầu, pin được áp **trước khi** phát hiện cấp-bộ ra đời, nên nó thành `unused` và
im lặng hết tác dụng — thứ tự nay là **phát hiện cấp bộ → áp pin → lần đọc thứ hai**.

**Bản ghi PGA của *Iris hermona* là ca âm thật, có sẵn trong repo.** `tests/data/corpus/…pga.gb`:
**84 CDS, 83 thiếu `/translation`**, và cổng IR hỏng (26,36 %, §8.31) ⇒ bị chặn; bản công bố cùng
chuỗi ấy qua sạch. Đây là test âm chạy được trên mọi máy, không cần `.panel/`.

**Một thay đổi ở tầng dữ liệu, vì cổng cần nó.** `Gene` nay mang thêm `/translation` của chính
bản ghi (`core/model.py`), đọc ở `core/genbank.py` (bỏ khoảng trắng, viết hoa), **không bao giờ
ghi ra** — D5-08 §8.21 nói engine tìm trên protein công bố, không tự dịch; và `_input_hash` của
cache nay băm cả trường ấy, vì hai bản ghi cùng chuỗi mà khác protein là hai input khác nhau.

**Chưa đo, và phải nói ra.** Audit **không** phát hiện một bản ghi bị đảo chiều *toàn phần* so với
bộ: bên trong bản ghi ấy LSC và SSC vẫn cùng chiều, nên luật không có gì để so — và với blastx
(tìm cả hai mạch) điều đó không ảnh hưởng kết quả chuyển chú giải. Nó cũng **không** kiểm độ dài
vùng IR của bản ghi chưa chú giải (H-13 vẫn mở). Cả hai là giới hạn đã biết, không phải số ước lượng.

**Kiểm chứng.** `tests/unit/test_reference_audit.py` (31 test, chạy mọi máy: 6 lớp chặn, luật đa
số N ≥ 3 và ca N = 2, ngoại lệ sống/void/stale/unused, dotplot, cửa CLI exit 2 + lý do trong
manifest, cửa `--confirm-exception`) · `tests/golden/test_reference_audit_corpus.py` (4 test, cần
`.panel/`; test MAFFT là tier B): tái lập **32/34 · 2/34** và tên đúng hai bản ghi, MAFFT đánh dấu
đúng hai bản ghi ấy, 0 phát hiện `confirmation`.

---

### 8.44. `gene_normalize` — bốn phép thử, ba trạng thái, ngưỡng đã đo; và chiều alias đảo (P3B-04, 08/10/2026)

**Việc.** Mô-đun mới `gene_normalize` (arity 2, `requires_ir=False`, hai chế độ `normalize` /
`verify`) — bước đứng trước annotator: mọi GenBank đi vào pipeline được **ghi lại theo cách viết
của bản ghi ★** (từ vựng A), kèm bảng nói cái gì đổi và vì sao. `genomes[0]` là **tham chiếu
chuẩn tên** (A), `genomes[1:]` là các bản ghi cần chuẩn hoá (B). **B chỉ xuất hiện ở đúng một
bảng** — `pairing`, ở cột riêng của nó (D5-09 §3).

**Tham chiếu chuẩn là một tham số có ghi lại, không phải một input ngang hàng.** Registry nay có
`ModuleSpec.reference_param`; `gene_normalize` khai `naming_reference`, và `cli.py` làm bốn việc
mà mô-đun không được làm (mô-đun không chạm filesystem): đọc file (từ tham số, hay từ
`data/naming_reference.json` khi run không override) → **chạy audit trước khi dùng** (D5-08: một
tham chiếu hỏng làm mọi cái tên sai theo cùng một kiểu) → đặt nó thành `genomes[0]` → ghi
`{path, accession, sha256}` vào `modules[].reference` của manifest. `MANIFEST_SCHEMA_VERSION` lên
**2**; mã thoát **4** là của riêng ca xung đột (1 lỗi · 2 dùng sai/cổng · 3 thiếu công cụ ·
**4 xung đột tên**).

**Bản ghi chuẩn hiện tại là một override có ghi, không phải quyết định của ★.** Người dùng cuối
chưa chọn (D5-09 §7.1), nên `data/naming_reference.json` ghi `provisional: true`, trỏ vào
`tests/data/NC_050833.1.gb`, kèm sha256 **đo từ chính file đó** — đổi byte là run dừng có tên.
Cùng lúc, **chiều của `data/gene_aliases.yaml` được đảo** (D5-09 §8): 5/9 locus đổi
(`ycf3`·`ycf4`·`psbN`·`clpP`·`ndhL` thay `pafI`·`pafII`·`pbf1`·`clpP1`·`ndhI`; bốn locus kia giữ
nguyên vì bản ghi ★ đồng ý với chúng). Chiều thật sự được **định hướng theo bản ghi mà run nhận**
(`oriented_aliases`), nên một ★ khác đảo lại mà **không phải sửa file** — file chỉ ghi chiều của
bản ghi đã ghi, và một test canh hai thứ đó khớp nhau.

**Số đo trên dữ liệu thật.** *Iris gatesii* `NC_024936.1` (86 CDS, có `/translation`) chuẩn hoá
theo genome vàng: **132 gen · 6 đổi tên · 121 xác nhận · 0 xung đột · 5 `outside vocabulary` ·
2 `not measured`**. Sáu ca đổi tên, và đúng phép thử nào quyết:

| B ghi | A ghi | phép thử | bằng chứng |
|---|---|---|---|
| `atpL` | `atpI` | margin | 100 % · phủ 1,00 (bảng alias biết locus, chuỗi xác nhận) |
| `psaL` | `psaI` | margin | 100 % · phủ 1,00 |
| `bspB` | `psbB` | margin | **99,80 %** · phủ 1,00 |
| `trnk-UUU` | `trnK-UUU` | neighbours | lỗi chữ nhỏ, hàng xóm xác nhận |
| `trnl-GUA` (×2) | `trnI-GAU` | neighbours | bảng alias + hàng xóm |

121 ca xác nhận chia theo phép thử quyết định: **82 margin · 38 neighbours · 1 sequence**. Năm ca
`outside vocabulary` đều là tRNA mà bản ghi ★ không dùng tên đó — và **một trong chúng đáng ghi
lại**: `trnL-CAU` của bản ghi trùng khít 100 % với `trnI-CAU` của bản ghi vàng (anticodon CAU
không phải của Leu), `trmA-UGC` trùng 100 % với `trnA-UGC`, `trnL` trần gần nhất là `trnL-UAG`.
Theo §8.41, với tRNA **không ngưỡng identity nào tách được locus đúng khỏi locus sai**, nên cả
năm ca được **báo**, không được đổi tên — kể cả ca 100 % ấy.

**Hai lỗi thật do test tìm ra, cả hai đều là loại im lặng.**
1. **Margin đọc bằng 0 trên gen IR.** Bản ghi ★ mang **hai bản sao IR** của một gen dưới một tên;
   BLAST trả cả hai với điểm gần bằng nhau, nên "á quân" chính là locus của quán quân ⇒ `ycf2`
   (và `rpl2`, `rpl23`, `rps19`) bị gọi là **xung đột** trong khi tên và chuỗi *đồng ý*. Luật nay:
   hit được **gom theo locus, giữ cái tốt nhất**, trước khi tính margin — "á quân" nghĩa là locus
   *khác* tốt nhất.
2. **`rrn23` của một bản công bố bị coi là hai gen.** Hai bản sao IR của nó khác nhau một hai base
   (cổng IR đòi ≥ 99,9 %, không đòi 100 %), nên phép so *chuỗi* gọi một gen thành hai và **chặn**
   một run không có gì sai. Luật nay dựa vào **cấu trúc**: hai feature cùng tên là một gen khi cả
   hai nằm trong cặp IR (một bên IRa, một bên IRb) — không dùng ngưỡng tương đồng nào. Đo lại:
   `trnG-UCC` (hai gen thật, cùng ở LSC, 73,24 % giống nhau, §8.41) vẫn bị bắt trên cả ba bản ghi
   đã đo; `rrn23` thì không.

**Ba giá trị `status` phải khác nhau, và khác thật.** `measured` = tìm được đối tác trong A **và
đã kiểm**; `not measured` = **không kiểm được** (thiếu `/translation`, ranking quá sát, hàng xóm
không đo được vì thiếu khung IR) — không bao giờ có nghĩa "không có"; `outside vocabulary` = đã
kiểm và A **không có gì** tương ứng. Kèm `outcome` nói bước này làm gì: `confirmed` (không đổi) ·
`renamed` (đầu ra viết theo A) · `conflict` (dừng) · `reported` (không đổi, và người đọc cần biết).
Hai nửa corpus PGA (84 và 85 CDS, `/translation` chỉ có ở **1** feature) chạy được trên mọi máy
mà không cần aligner: mọi CDS báo `not measured` với lý do, không đoán.

**Điều đo được đã đổi sau khi đảo chiều.** `tests/golden/test_naming_numbers.py` (P3B-02) pin ba
giá trị nằm trên sàn ortholog; sau khi đảo chiều, số *dòng* trên sàn đổi 3 → 5 (cùng hai giá trị
97,674 % và 100 %) vì cách gom tên theo locus đổi — nghĩa là **cách đếm cũ phụ thuộc vào chiều
chuẩn hoá**. Test nay pin **tập giá trị** (97,674 % và 100 %) chứ không pin số dòng: điều §8.41
khẳng định là *khoảng trống* giữa 44,93 % và 67,75 %, và điều đó không đổi.

**Chưa đo.** Ts/Tv dưới chiều mới (**P3B-06** — 2,6917 vẫn là dự đoán, không phải kết quả);
`kaks` · `substitutions` · `popgen` chạy lại trên bộ thật (bộ test xanh, nhưng đó là test, không
phải phép đo); và bản ghi **★ thật** — người dùng cuối chọn từ báo cáo audit (§8.43), lúc đó
`test_the_recorded_reference_is_what_the_audit_and_the_flip_agree_about` là chỗ hai thứ gặp nhau.

**Kiểm chứng.** `tests/unit/test_gene_normalize.py` — 29 test: ngưỡng đo được · định hướng alias
theo bản ghi · bốn phép thử, mỗi phép một ca âm · ba trạng thái · IR copy không phải trùng tên ·
locus gom trước khi tính margin · chế độ verify không ghi file · chạy cặp corpus ra `.gb` + bảng
`pairing` · xung đột dừng với toạ độ · `on_conflict=continue` ghi bảng `conflicts` · cửa CLI: giải
tham chiếu → audit → `genomes[0]` → `modules[].reference` (schema 2) · thoát **4** khi xung đột ·
một test tier B chạy blastp thật. 11 test canh chiều cũ trong `test_gene_names.py` ·
`test_annotation.py` · `test_tier_a.py` · `test_corpus.py` được **cập nhật theo chiều mới** (không
nới: cùng khẳng định, tên khác).

### 8.45. Ts/Tv đo lại dưới chiều tên mới — con số không đổi, và vì sao nó *không thể* đổi (P3B-06, 08/10/2026)

**Việc.** P3B-04 lật chiều `gene_aliases.yaml` (5/9 locus) và **ba mô-đun đọc khoá chuẩn**
(`kaks` · `substitutions` · `popgen`) — nên §6 ghi 2,6917 là **dự đoán** cho tới khi đo lại, và
§8.44 ghi lại đúng chữ "chưa đo". Đây là lần đo đó.

**Đo theo đúng thứ tự pipeline, không phải bằng cách chạy lại test cũ.** Một mình
`substitutions` trên hai bản ghi công bố là phép đo của §8.27, không phải phép đo của chiều mới:
trong pipeline thật, mẫu **thừa hưởng** từ vựng A *qua* bước đặt tên. Nên lần đo này chạy
`gene_normalize` (`mode=normalize`, tham chiếu = genome vàng) trước — artifact
`NC_024936.1.normalized.gb`, đúng 6 ca đổi tên của §8.44 — rồi mới đưa **artifact đó** vào
`substitutions --set scope=cds`, tham chiếu là chính `NC_050833.1.gb` công bố (bản ghi vàng không
đổi tên gì, nên nó là A ở cả hai dạng).

**Số đo — trùng khít, không một đơn vị nào lệch.**

| | §6 / §8.27 (04/10) | đo lại 08/10, chiều mới |
|---|---|---|
| gen dùng chung | 79 | **79** |
| cột căn | 68.448 | **68.448** |
| cột so sánh được | 68.382 | **68.382** |
| thay thế | 443 | **443** |
| Ts / Tv | 323 / 120 | **323 / 120** |
| Ts/Tv | 2,6917 | **2,6917** |

**Vì sao nó không thể đổi — đo, không suy luận.** Bảng alias ánh xạ **mọi** thành viên của một
locus về locus đó, bất kể thành viên nào đang là canonical ⇒ lật chiều là phép **đổi nhãn** của
tập gen dùng chung, không phải phép **chọn lại** nó. Nó chỉ có thể chọn lại trong đúng một tình
huống: một bản ghi mang **hai gen khác nhau** mà tên rơi vào cùng một locus — khi đó `setdefault`
giữ một feature và chiều quyết định giữ cái nào. Đo trên cả hai bản ghi Iris: **86 CDS → 79 tên
chuẩn** ở mỗi bên, và **0 locus** có hai cách viết (86 − 79 = 7 gen nằm trong IR, mỗi gen hai
bản). Không có ca nào như vậy ⇒ tập gen, cột căn và mọi con số là **bất biến theo chiều**.

**Nhãn đổi, số không đổi.** `substitutions_by_gene` nay in `ycf3` · `ycf4` · `psbN` · `clpP` ·
`ndhL` · `atpI` · `psaI` · `psbB`; bảng tương ứng của §8.27 in tên cũ cho cùng những locus đó.
Đây là *từ vựng*, không phải *kết luận* — cùng locus, cùng số (cách phân biệt mà §8.31 đặt ra ở
quy mô 34 bản ghi).

**Test canh.** `tests/golden/test_tier_b.py::test_the_naming_step_does_not_move_the_golden_ts_tv`
— chạy đúng thứ tự trên và khẳng định **artifact thật sự ở từ vựng A** (`atpI` · `psaI` · `psbB`
có mặt trong nó), nên một run mà bước đặt tên im lặng không làm gì **không** qua được test này.
Test cũ `test_the_cds_scope_reproduces_the_golden_ts_tv` giữ nguyên (đường công bố, không qua
bước đặt tên) — hai test, hai đường, cùng một con số.

**Phát hiện kèm, không thuộc phép đo này: artifact trộn hai từ vựng trong một file.** Đọc lại
artifact của chính run trên: **6 feature `gene` còn nguyên từ vựng B** (`atpL` · `psaL` · `bspB` ·
`trnk-UUU` · `trnl-GUA` ×2), nằm **cùng span** với CDS/tRNA đã đổi thành `atpI` · `psaI` · `psbB` ·
`trnK-UUU` · `trnI-GAU` — `GENE_KINDS` loại feature `gene` khỏi mọi quyết định đổi tên. Hai bản
ghi công bố thì khớp **132/132** feature với feature `gene` cùng span cùng tên (0 lệch, 0 thiếu),
và CGAS — mô hình D5-08 bám theo — **có** đổi feature `gene` (D5-09 Appendix A: "8 features per
target, 4 tRNA + 4 `gene`"). Ghi tên **H-17**, chi tiết ở §9.K-1; hôm nay nó là lỗi *hiển thị*,
không phải lỗi số.

### 8.46. Lõi gọi gen dùng chung, và ba lỗi im lặng nó lôi ra ánh sáng (P3B-07, 08/10/2026)

**Việc.** Bộ gọi gen plastid (D-e), theo đúng luật DRY của dự án: bản sao thật thứ hai đã tới, nên
máy móc được **tách ra dùng chung** thay vì chép. `core/caller.py` (mới) giữ toàn bộ luật biến HSP
thành `Gene` — chuỗi HSP, exon, tinh chỉnh codon, các finding; `core/mito_caller.py` co lại thành
lớp mỏng chỉ còn **thứ đo được riêng cho ty thể** (ngưỡng 50 %/0,5 và bảng mã 1);
`core/plastid_caller.py` (mới) là nửa plastid: chỉ mục tham chiếu dựng từ bản ghi GenBank, ngưỡng
50 %/**0,7**, bảng mã **11**, và cách đọc inverted repeat từ self-blast.

**Bằng chứng việc tách không đổi hành vi: 43 test mito xanh nguyên vẹn, không sửa một dòng.** Đó
là điều kiện "Done when" của task, và là lý do tách theo cách này chứ không chép.

**Ngưỡng đo trên bản ghi vàng (đáp án).** Chạy thật blastx/blastn: genome vàng × 33 bản ghi còn
lại của bộ tham chiếu đã pin (leave-one-out), BLAST+ 2.17.0, `-seg no`, evalue 1e-5.

| Đo | Kết quả |
|---|---|
| HSP tìm được | **4.382** (CDS) · 528 (rRNA) |
| Call ở ngưỡng 50 % / 0,7 | **80 CDS** (79 locus công bố + `lhbA`) · **4 rRNA** = `rrn16`·`rrn23`·`rrn4.5`·`rrn5` |
| Call đáp xuống một locus công bố | **80/80** — không call nào rơi vào vùng không được chú giải |
| Ngưỡng phủ 0,5 → 0,8 | **cùng một tập call** — bản ghi vàng *không phân biệt được* ngưỡng phủ |

Dòng cuối là điều đáng ghi nhất: bản ghi vàng là họ hàng gần, nên nó **không chọn được** ngưỡng.
Ngưỡng **0,7** là số của §8.41 — đo trên **corpus 34 bản ghi**, nơi khoảng trống thật nằm giữa
44,93 % (locus khác) và 67,75 % (ortholog thật), và nơi 0,5 thủng vì `psbJ` chỉ 36 aa. Nói cách
khác: ngưỡng *không* được chọn hôm nay, nó được **kế thừa kèm tên phép đo**, và P3B-10 sẽ kiểm lại
trên corpus.

**Ba lỗi im lặng, cả ba do chính phép đo này lôi ra — không phải do test:**

1. **`-max_target_seqs 500` cắt cụt bộ tham chiếu trong im lặng.** Hằng số 500 đúng cho 287 chuỗi
   của đường mito; bộ plastid có **2.838** chuỗi. Đo: cùng một lệnh, cap 500 → **929 HSP / 467
   chuỗi**, cap = cỡ bộ → **4.382 HSP / 2.603 chuỗi**. BLAST không nói gì. Sửa: ngưỡng nay **đếm
   từ chính file** bộ tham chiếu (`subject_count`), và không còn là hằng số nào cả.
2. **Một tên có khoảng trắng làm header FASTA bị BLAST cắt.** `NC_042827.1` đặt tên rRNA là
   `23S ribosomal RNA` (không có `/gene`, reader rơi về `/product`) ⇒ header thô về tới nơi là
   `NC_042827.1|23S`, và parser **từ chối có tên** thay vì khớp một khoá không bao giờ khớp. Sửa:
   `ReferenceEntry.key` bỏ khoảng trắng trong từng field, và header FASTA dùng chính `key` ấy.
3. **Tên rRNA dạng câu sản phẩm chưa được chuẩn hoá.** Bốn locus ấy vào bộ từ vựng như *tên gen*
   (`16S ribosomal RNA`…), nên không run nào gọi tên được chúng. Sửa: `core/gene_names` mở rộng
   luật đã có (`rrn16S` → `rrn16`) cho dạng câu — cùng một tri thức, không phải tri thức mới; bộ
   từ vựng rRNA nay đúng **4** ký hiệu.

**Một phát hiện của tầng đặt tên, chưa sửa: `lhbA` là `psbZ`.** Bộ tham chiếu có **80** tên CDS
trong khi số locus thật là **79**: `NC_042827.1` gọi locus 62 aa ấy là `lhbA`, 33 bản ghi còn lại
gọi `psbZ` (bản ghi vàng: `psbZ`, 62 aa, 35.895–36.083). Hôm nay caller gọi **cả hai** — đúng theo
luật "một call cho mỗi locus *của từ vựng*", và sai theo luật "một call cho mỗi locus *thật*".
Gốc là bảng alias (D5-09) thiếu một mục; đặt tên **H-18**, kèm bằng chứng đo được.

### 8.47. Annotator plastid chạy thật trên genome vàng — 86 CDS, và hai khác biệt phải phân loại (P3B-08, 08/10/2026)

**Việc.** `modules/plastid_annotate.py` — lớp wiring + xuất: đọc bộ tham chiếu đã pin (34 bản ghi
Iridaceae, **leave-one-out** lấy từ `plastid_reference.json` nên bản ghi vàng không tự chú giải
chính nó), chạy `blastx` (protein) · `blastn` (rRNA) · `blastn` tự thân (IR) · ARAGORN (tRNA), gọi
gen qua `core/plastid_caller.py`, rồi trọng tài và xuất.

**Số đo trên bản ghi vàng** (BLAST+ 2.17.0, ARAGORN 1.2.41, ngưỡng 50 %/0,7, bảng mã 11):

| | Gọi được | Bản công bố |
|---|---|---|
| CDS | **86** | 86 |
| rRNA | **8** | 8 |
| tRNA | **30** | 38 |

**86 CDS là kết quả của một luật, không phải của may mắn.** Bộ gọi cố ý ra **một call cho mỗi
locus** (79 locus CDS); 86 = 79 + 7 vì **bảy locus nằm trong IR** và bản công bố ghi mỗi locus hai
lần. Việc nhân bản là của **lớp xuất**, không phải phép đo thứ hai: `_ir_copies` ánh xạ một gen nằm
trọn trong một bên IR sang bên kia bằng chính phép đảo (đối xứng qua tâm vùng, đổi mạch), và **báo
tên** một call vắt qua ranh giới IR thay vì đoán. Cùng luật ấy cho rRNA: 4 locus → **8** feature.

**Ba điều Done-when, đo được:**

- **Chồng lấn được giữ và báo, không bị giải quyết bằng cách bỏ.** Đo theo **exon**, không theo
  envelope — bài học đo được ngay trong task này: `rps12` là gen trans-spliced, envelope của nó
  trải 68.976–98.016, nên một phép thử theo envelope làm nó "chồng lấn" **27 gen** và chôn mất bốn
  cặp thật. Theo exon: `psbC`×`psbD` **53 bp** · `ndhC`×`ndhK` **10** · `atpB`×`atpE` **4** — ba
  cặp của §8.34, đúng từng base.
- **Một call cho mỗi locus, và mọi kẻ thua được liệt kê.** `competing_references`: mỗi locus đúng
  **một** dòng `chosen=yes`, kèm mọi bộ tham chiếu khác đã được xét với số của nó.
- **Artifact đọc lại được và khai báo organelle.** `/organelle="plastid"` trên `source`, hai
  `repeat_region` = hai bên IR, và vòng đọc-lại (`read_genbank`) chạy **trước khi** kết quả được
  trả về — file không đọc lại được thì không phải kết quả (H-11).

**Hai khác biệt với bản công bố, đo được, phải phân loại — chưa sửa:**

1. **`ycf1` được gọi thành 900 bp ở 108.448–109.347, trong khi bản công bố là 5.421 bp ở
   122.984–128.404.** Ứng viên thắng là `NC_088562.1` ở **98,83 %** với coverage 1,000 — nghĩa là
   bộ tham chiếu ấy mang một `ycf1` **cụt**, và luật "identity cao nhất thắng" đã chọn nó. Đây
   đúng là ca mà D5-08 §3 gọi tên: **cấu trúc phải quyết trước, điểm chỉ tie-break trong cùng một
   chương trình** — và luật ấy **chưa được thi hành** trong `core/caller.py` (`_choose` vẫn
   identity-first). Hệ quả kèm theo: cặp chồng lấn `ycf1`×`ndhF` 21 bp trong bảng là *của artifact
   này*, không phải của bản công bố. Ghi tên **H-19**.
2. **tRNA 30 so với 38.** ARAGORN tìm 30 feature; bản công bố ghi 38 (29 tên khác nhau, 9 tên có
   hai bản trong IR). Lớp xuất **chỉ nhân bản CDS/rRNA**, chưa nhân bản tRNA — nên phần thiếu
   trước hết là của phép nhân bản, phần còn lại (nếu có) là của ARAGORN. P3B-10 phân loại.

Cả hai đều là **khác biệt mức kết luận**, không phải mức từ vựng — đúng loại mà P3B-10 phải đếm
và so với PGA.

### 8.48. Một locus bị chú giải hai lần: bản cụt không được tìm kiếm, và 0,5 nằm giữa một khoảng trống đo được (H-19, 09/10/2026)

**Việc.** H-19, đặt tên ở §8.47: `ycf1` của bản ghi vàng bị gọi thành 900 bp vì trọng tài
identity-first chọn một bộ tham chiếu mà **chính nó** mang `ycf1` cụt — 98,83 % trên 299 residue,
coverage 1,000 — thay vì các ứng viên đầy đủ (97,2–97,7 %). D5-08 §3 nói **cấu trúc quyết trước**.
Luật thi hành: **một locus mà một bản ghi chú giải hai lần, bản ngắn dưới 50 % bản dài, là bản sao
bị cụt — không được đưa vào tìm kiếm, và được nêu tên** (`truncated_annotations`, in thành note
của run). Chỗ sửa không phải `_choose` mà là **bộ chỉ mục tham chiếu**: một ứng viên không tồn tại
thì không có điểm nào để thắng.

**Vì sao 0,5, và vì sao luật dừng ở trong-một-bản-ghi.** Đo trên bộ đã pin (34 bản ghi, **chỉ
CDS**): **243** locus được chú giải hai lần — **216 ở đúng 1,000** (hai bản của IR), **một ở
0,9927** (`rpl2`, 271 so với 273 residue: khác biệt thật, và nó **được giữ**), **26 nằm dưới 0,30**
(`ycf1` 0,16–0,30 ở 24 bản ghi; `ycf4` 0,23 ở 2 bản ghi). **Không có gì nằm giữa 0,2992 và
0,9927.** 136 locus rRNA trùng nhau đều ở 1,000.

Phép so sánh **chỉ có nghĩa trong một bản ghi**, và đó là chỗ dựa của luật, không phải một giới
hạn kỹ thuật. Trên 80 locus CDS, chỉ **3** locus có dải dài giữa-các-bản-ghi dưới 0,5 — và cả ba
**chính là** ca trong-một-bản-ghi (`ycf1` · `ycf4` · `rps19` 22/22 của `NC_081931.1`) — trong khi
`rpoB` chạy 614–1.070 residue và `ndhD` 383–523. Giữa các bản ghi, một **ortholog ngắn** và một
**chú giải cụt** không phân biệt được bằng độ dài; trong một bản ghi thì hai bản là *cùng một gen*,
nên so được. Luật vì thế không mở rộng ra ngoài, và lý do là hệ quả chứ không chỉ là độ chính xác:
một luật độ dài xuyên-bản-ghi sẽ bỏ những entry như `rpoB` (0,574) hay `ndhD` (0,733) mà **chưa đo
được** là cụt hay thật — và bỏ một reference thì không lấy lại được, mọi dòng `not_called` sau đó
thành lời nói dối. Luật trong-một-bản-ghi chỉ bỏ một feature **đã có bản đầy đủ cùng bản ghi**.

**Đo trên bản ghi vàng, trước/sau: toàn bộ tập call đổi đúng một dòng.**

| | trước | sau |
|---|---|---|
| `ycf1` | 108.448–109.347 (**+**, 900 bp, 98,83 %, `NC_088562.1`) | **122.984–128.404 (−, 5.421 bp, 97,68 %, `NC_081929.1`)** |
| bản công bố | `ycf1 [122983:128404](-)`, 1.806 aa | **khớp chính xác từng base** |

Ứng viên `ycf1` nay còn **33**, tất cả đầy đủ (1.783–1.862 residue); artifact vẫn **86 CDS +
8 rRNA + 30 tRNA** với `/organelle="plastid"`; note của run từ 30 lên **56** dòng, 26 dòng mới nêu
tên từng feature bị bỏ ("`NC_088562.1` annotates `ycf1` twice: the 1.806-residue feature is used
and its 299-residue copy (17 % of it) is not searched"). Lính gác mito: **65 test** mito (`-k mito`, tầng unit) xanh và **không file mito nào bị sửa**
(`git status` chỉ có `plastid_caller.py` · `plastid_annotate.py` · test · tài liệu).

**Vế thứ hai của Done-when đã không làm, và lý do là một phép đo.** Done-when viết "hoặc một ứng
viên mà alignment không chạm hai đầu reference, trong khi một ứng viên dài hơn thì chạm". Đo lại
thì vế ấy **không sửa được ca này**: alignment của bộ tham chiếu cụt **có** chạm cả hai đầu của nó
— coverage 1,000, chính con số ấy là thứ đánh lừa trọng tài. Và nó là **cơ chế sai**: "gen của mẫu
bị cụt" đã có cách biểu diễn đúng, cờ `partial_start`/`partial_end` trên call, nên ưu tiên một
reference để *tránh* một call partial sẽ **giấu một phát hiện thật**. Ghi lại đây để P3B-10 đọc
đúng phạm vi của luật.

**Còn lại, có tên.** (1) Một bản ghi chú giải locus ấy **một lần** thì không có gì để so — luật
không bắn; **đo được: bộ đã pin không có ca nào như thế**. (2) Audit **có** bắt lớp đó bằng toạ độ
(`<`/`>`), đo trên chính bộ này: `NC_024936.1` bị **BLOCK** vì `rps19` partial (`>83020`,
`<153124`) — **nhưng audit chưa được đặt trước run annotate**: `plastid_annotate` không gọi nó, và
bản ghi ấy vẫn là ứng viên trong run vàng. Nối cổng ấy vào lệnh `annotate` là **P3B-09**, đã ghi
vào checklist của task đó.

**Test.** Ba test đơn vị đã có từ trước (luật bắn · khác biệt 2 residue không bị coi là cụt · hai
bản IR bằng nhau vẫn là hai entry) cộng **một mới**: một feature lẻ, ngắn, **vẫn được tìm** — luật
bắn trên một phép *so sánh*, và một feature thì không có phép so nào. Golden:
`test_ycf1_lands_on_the_published_gene_not_on_a_truncated_reference` — khẳng định toạ độ
122.984–128.404 **và** rằng run nói ra feature nó đã không tìm.

### 8.49. `plastid annotate` — cổng audit trước run, batch một mẫu một thư mục, và hai phát hiện đo được (P3B-09, 09/10/2026)

Task **P3B-09** (M3B.4 · H-09). Lệnh `plastid annotate` ra đời: 1..N mẫu × 1..N bản ghi tham
chiếu, **mỗi mẫu một thư mục run riêng** (manifest · report · artifact) — batch là *một loạt run*,
không phải một run chứa N mẫu, vì §2 chốt "output tách theo lần chạy/dataset, không ghi đè tại
chỗ", và một mẫu hỏng **không được** làm hỏng dấu vết của các mẫu còn lại.

**Cổng audit nay nằm trước run — và nằm trong module, không nằm ở CLI.** H-19 gọi tên lỗ này
(§8.48): `plastid_annotate` không hề gọi audit. Chỗ nối là `_reference_set`, tức **đúng tập bản ghi
run sắp tìm**, chạy **trước HSP đầu tiên** — nên "không chú giải gì" là nghĩa đen chứ không phải
hy vọng. Đo trên chính bộ đã pin, `plastid annotate -i tests/data/NC_050833.1.gb --reference
.panel/plastid-reference`: **exit 2** sau **13 giây**, hai dòng `file — reason` và **không artifact
nào**:

```
.panel/plastid-reference/NC_024936.1.gb — CDS rps19: partial coordinates (>83020)
.panel/plastid-reference/NC_024936.1.gb — CDS rps19: partial coordinates (<153124)
```

**Cách đi tiếp là `--exclude`, và nó không phải một lối vòng qua cổng.** D5-08 viết vòng lặp là
*"fix or drop the bad file → run again"*; `--exclude <accession>` là dạng cơ học của "drop", **khai
trên dòng lệnh** (manifest ghi lại nguyên argv) và **run nói ra**: nhãn bộ tham chiếu thành
`32 GenBank record(s) …, audited (0 blocked), 1 left out by leave-one-out (…), 1 dropped by name
(NC_024936.1)`, và `result.json` mang một note riêng. Phán quyết về file đó **không đổi** — pin
không bao giờ waive được một BLOCK (`core/reference_audit`) — chỉ có điều file ấy không được tìm.
Một cái tên không khớp bản ghi nào thì **dừng run** chứ không im lặng: người dùng tin rằng một bản
ghi đã bị bỏ ra trong khi run vẫn tìm nó là đúng kiểu sai-im-lặng của dự án (§8.43, luật của pin).

**Đo được — bỏ bản ghi bị chặn làm đổi đúng một dòng của tập call, và dòng ấy xấu đi.** Vàng
(33 tham chiếu, có `NC_024936.1`) vs vàng (32 tham chiếu, đã bỏ nó): **cùng 84 dòng** (80 CDS + 4 rRNA), khác ở chỗ
bản ghi được chọn cho nhiều locus (cột 9), và **đúng một dòng khác về nội dung** — `rps19`:

| | toạ độ | dài | bản ghi được chọn |
|---|---|---|---|
| có `NC_024936.1` | 83.173–83.451 | **279 bp** | `NC_024936.1` |
| đã bỏ nó | 83.173–83.397 | **225 bp** | (một bản ghi 225 bp) |
| bản công bố | 83.173–83.451 | 279 bp | — |

Nguyên nhân đo trên cả bộ: `rps19` của 34 bản ghi là **27 × 279 bp/92 aa · 2 × 285 bp/94 aa ·
3 × 225 bp/74 aa · 1 × 69 bp/22 aa** — ba bản ghi (`NC_081929.1` · `NC_081930.1` · `NC_081932.1`)
mang một `rps19` **bị chú giải cụt**, và khi bản ghi đầy đủ bị bỏ ra thì một trong ba bản cụt ấy
thắng (cùng 100,00 % identity) rồi **kéo call về độ dài của chính nó**. Đây là đúng vế mà H-19
**chưa làm**: D5-08 §3 "cấu trúc trước, điểm sau" trong *trọng tài giữa các ứng viên* — H-19 mới
sửa được vế trong-một-bản-ghi (bản sao cụt bị loại khỏi chỉ mục), còn vế **giữa các bản ghi** thì
chưa. **Ghi tên: H-20** (§9.L-5), và nó là một trong những khác biệt P3B-10 phải phân loại.

**Phát hiện thứ hai, sửa ngay: báo cáo `.docx` của engine chưa từng chạy được.** `plastid annotate
--docx` dừng với `CitationError: citation key 'blast' is not in the registry` — module khai
`citations=("blast", "aragorn")` mà `data/citations.yaml` không có hai khoá ấy (registry mới có
misa · biopython · mafft · pga · chloe · mummer · paml · cbol · yang_nielsen). Đã thêm hai mục
(BLAST+ Camacho và cs. 2009 · ARAGORN Laslett & Canback 2004) và viết
`report/templates/plastid_annotate.en.md` — engine mới **chưa có template**, nên trước đó run không
sinh `report.md` nào cả. Đo lại: `--docx` chạy hết, `report.md` + `report.docx` ra đời, phần
References in đúng hai bài.

**Bào quan, khai ở đầu vào (H-09).** `core/organelle.py` là chỗ duy nhất trả lời "bản ghi này là
bào quan nào": đọc `/organelle` (đo trên hai bộ đã pin: `plastid:chloroplast` ×29 · `plastid` ×5 ·
`mitochondrion` ×8) và **không đoán** khi bản ghi không khai — D5-06: chuỗi trần tự khai bằng
**đường vào**. Registry nay mang `group` cho cả **22 module** (§13 → "Hai bào quan" là bảng nguồn;
`group` **bắt buộc**, không có mặc định, vì bỏ sót không phải là một lời khai). Hai chiều lệch đều
dừng **có tên**, kèm đường đúng: `mito_genes` trên plastome hôm nay không còn đi qua im lặng (§9.J
đo 05/10). Một **bộ trộn** hai bào quan cũng dừng, và dừng **trước khi** mẫu nào được chú giải.

**Phát hiện thứ ba, ghi lại chứ chưa sửa: luật BLOCK "toạ độ partial" từ chối một bản ghi mà
protein của nó đầy đủ.** `NC_024936.1` có **hai** feature `rps19`, mỗi cái **279 bp / 92 aa** với
`/translation` trọn vẹn; cờ partial (`>83020`, `<153124`) là cách bản công bố mô tả gene nằm ở
ranh giới IR, không phải một protein cụt. Bằng chứng mạnh nhất: **chính bản ghi ấy** là bản ghi đã
cho ra `rps19` khớp bản công bố từng base (83.173–83.451, 100,00 %, không cờ partial nào trên
call). Vậy lý do của luật ("một reference partial đóng dấu partial lên mẫu", §8.15 trap 1) **không
đúng cho dạng này** — nhưng sửa một luật BLOCK đã chốt ở P3B-03 là việc của một task có đo lại
toàn bộ, không phải việc lặng lẽ làm trong P3B-09. **Ghi tên: H-21**, kèm giá đo được của việc bỏ
bản ghi ấy (bảng `rps19` ở trên).

**Cache: một tham số là *thư mục* nay cũng được định địa chỉ theo nội dung.** P3B-05 hứa "một tham
số mang giá trị file thì tính theo nội dung" nhưng chỉ áp cho `is_file()`; `--set reference=<dir>`
được băm như chuỗi đường dẫn, nên sửa một bản ghi trong bộ tham chiếu sẽ trả lại bảng của lần chạy
trước — đúng kiểu H-15, lùi một tầng. Nay thư mục được băm theo **tên file + sha256 từng file**, và
test khẳng định cả ba cách bộ đổi (thêm · bớt · sửa) đều đổi khoá.

**Test.** `tests/unit/test_organelle.py` (15) — khai báo · hai chiều lệch · bộ trộn · registry mang
`group` cho cả 22 module (hai tập tên chính xác là tripwire có chủ ý) · không thể đăng ký module mà
không khai nhóm. `tests/unit/test_cli_annotate.py` (6) — mẫu mito bị từ chối **có tên** · một
FASTA mito bị từ chối qua bộ tham chiếu · batch trộn hai bào quan dừng trước khi chú giải gì ·
`--audit-only` không cần mẫu còn batch thì cần · **mẫu hỏng được đánh dấu và nêu tên mà không giết
batch** · mỗi mẫu một thư mục và một artifact. Golden
`tests/golden/test_plastid_annotate.py`: ba test cũ nay chạy với `--set exclude=NC_024936.1` (cổng
buộc phải bỏ bản ghi ấy) và **vẫn khẳng định đúng** — `ycf1` 122.984–128.404 không đổi — cộng ba
mới: cổng chặn run thật · bỏ theo tên thì không tìm và run nói ra · một cái tên không khớp thì
dừng. Manifest nay mang khoá `dataset` (tuỳ chọn, schema 2 không đổi): một run hỏng trước khi kịp
ghi gì vẫn phải nói nó là mẫu nào — đo được ở chính batch này.

### 8.50. P3B-10 mở màn: tRNA có intron, bản sao IR đặt sai chỗ, và con số corpus đầu tiên (10/10/2026)

**P3B-10 chưa xong.** Mục này ghi lại ba phép đo đã làm và **nói thẳng chỗ thanh chuẩn chưa đạt**.

**(1) ARAGORN bỏ qua tRNA có intron — trong im lặng.** Help của chính nó viết: *"Unless `-i` is
specified, tRNA genes containing introns are not detected"*, và nó không nói gì khi bỏ qua. Adapter
của dự án không truyền `-i`, nên bản ghi vàng ra **30 tRNA** nơi bản công bố có **38**; tám cái
thiếu đúng là tám tRNA có intron (`trnK-UUU` · `trnL-UAA` · `trnV-UAC` · một `trnG` · cả hai bản
`trnI-GAU` và `trnA-UGC`) — đo bằng cách khớp khoảng: 30/30 call khớp một feature công bố, và 8
feature công bố không có call nào. Với `-i`: ARAGORN trả **38**, và artifact của bản ghi vàng nay là
**132 = 86 CDS + 38 tRNA + 8 rRNA** — đúng hợp đồng §6. Span của gen có intron được **tách thành
hai exon** (`i(start,end)` của ARAGORN, đổi theo mạch và theo gốc) — không tách thì một feature
viết trên cả span sẽ nói rằng intron là vùng mã hoá.

**`-i` là tham số của từng đường, không phải mặc định**, vì hai phép đo không đồng ý: trên lục lạp
nó cho đúng 38/38; trên bản ghi mito answer-key nó cho 17 → **21**, bốn cái bản công bố không có.
Đường lục lạp truyền cờ; mở cho mito là một câu hỏi có tên, không phải một mặc định.

**(2) Bản sao IR trong artifact bị đặt sai — lệch đúng bằng chiều dài IR.** `_ir_copies` tính phép
nghịch đảo từ **đầu** bản IR thứ nhất thay vì **cuối** bản ấy, nên mỗi bản sao được ghi ra cách chỗ
đúng đúng một chiều dài IR (26.214 bp trên bản vàng) — tức **ngoài vùng lặp**, rơi vào vùng một
bản. Đo 10/10/2026: `rpl2` bản thứ hai được ghi ở 125.457–126.940, bản công bố ở 151.668–153.154;
`ndhB` ở 114.075–116.306 thay vì 140.289–142.520. **Số lượng gen vẫn đúng** (86) nên mọi test đếm
đều xanh — đúng kiểu "sai mà không báo" mà dự án này tồn tại để chặn. Sửa: `target.start + (source.end
- part.end)`. Đo lại: `rpl23` · `ndhB` · `ycf2` · `rps7` · `rrn16` · `rrn5` khớp bản công bố **từng
base**, `rpl2` lệch 3 bp (ranh giới IR do bộ dò tìm khác bản ghi vài base). Test mới:
`test_every_copy_the_inverted_repeat_is_given_lands_where_the_record_puts_it` — kiểm bản sao **theo
toạ độ của bản công bố**, không theo số học của chính nó.

**(3) Con số corpus đầu tiên — và thanh chuẩn CHƯA đạt.** 34/34 bản ghi chạy leave-one-out xong
(0 run hỏng, `--exclude NC_024936.1` theo cổng audit). Đối chiếu với bản công bố, phân loại theo
luật ghi trong harness (`vocabulary` = cùng locus khác cách viết; còn lại là `conclusion`):

| | kết luận | trong đó: thiếu · thừa · khác độ dài | vocabulary |
|---|---|---|---|
| **bộ gọi của dự án** | **1.579** | 167 · 161 · 1.251 | 376 |
| **PGA** | **637** | 39 · 4 · 594 | 149 |

**PGA đang thắng ở mức kết luận** (637 so với 1.579) — Done-when của P3B-10 *không* đạt, và đây là
con số nói điều đó. Bộ gọi của dự án **thừa/thiếu nhiều hơn hẳn** (161+167 so với 4+39), trong khi
PGA gần như chỉ khác ở độ dài — hợp lý với hai cơ chế: PGA **chép toạ độ** từ bản tham chiếu (nên
độ dài hay đúng) nhưng bỏ sót/thêm locus (đo được: 39 thiếu · 4 thừa), còn bộ gọi ở đây **suy toạ độ
từ alignment** (nên độ dài lệch vài base, và bốn gen có exon đầu 6–40 bp không được gọi). Phân rã
của bộ gọi: CDS 753 kết luận · tRNA 717 · rRNA 109; phần lớn là lệch 1–3 bp.

**Việc còn lại của P3B-10, có tên:** (a) vì sao 161 call thừa / 167 thiếu — chưa đo; (b) exon cụt
của `rps16` · `petB` · `petD` · `rpl16` (exon đầu 6–40 bp, blastx không bắt được — có nên lấy từ
bản tham chiếu không, và cái giá của luật ấy); (c) `rps19` −54 bp (**H-20**) và `lhbA` (**H-18**);
(d) 8 tên tRNA khác bản công bố (**H-04**'s `trnG`, `trnI-GAU` bị ARAGORN đọc thành `trnE-UUC`,
`trnfM-CAU`/`trnI-CAU` thành `trnM-CAU` — ARAGORN đọc anticodon của tRNA có intron); (e) harness
phân loại chưa được commit.

### 8.51. P3B-10: chạy lại corpus sau hai bản sửa, dựng thước đo thành code trong repo — và ba phát hiện mới (10/10/2026)

**(1) Con số của §8.50 là của bộ gọi *trước* hai bản sửa.** Chạy lại 34 bản ghi leave-one-out bằng
code hiện tại (0 run hỏng, `--exclude NC_024936.1`): **thiếu 167 → 47, thừa 161 → 41**. Toàn bộ họ
"sao chép IR" biến mất — đúng như bản sửa §8.50 dự đoán: `rps7` · `rpl23` · `rpl2` · `ndhB` · `ycf2`
· `rps12` · bốn locus rRNA, mỗi thứ 11–12 bản ghi "thiếu" một bản sao, nay khớp. Số khác độ dài
tăng nhẹ (1.251 → 1.278) vì bản sao nay đúng chỗ nên **nhiều cặp được đối chiếu hơn**, không phải
vì có thêm lỗi.

**(2) Thước đo nay nằm trong repo, có test.** P3B-10 hỏi một câu so sánh, nên nó chỉ tranh luận được
nếu cái thước ở trong repo, chạy lại được, và **phát biểu luật của nó ở đúng một chỗ**:
`tools/corpus_bar.py` (chạy `make corpus-bar`) + `tests/unit/test_corpus_bar.py` (19 test cho chính
luật ấy). Khi chuyển luật từ script nháp vào repo, **hai lỗ của luật tự lộ ra** và đều đã sửa, mỗi
lỗ một test: (a) bản đầu chỉ so **envelope** — hai cách đọc `rpoC1` của bản vàng cùng trải
20.179–22.962 mà lệch nhau 27 bp *bên trong*, luật cũ gọi cặp ấy là giống hệt; (b) luật cũ **không
so sợi** — một call trên sợi đối diện với toạ độ y hệt đi qua như "không khác gì" (corpus có đúng
một ca: `trnS-UGA`).

**(3) Luật, và một ngưỡng được khai ra.** Cặp khớp = cùng locus, cùng chỗ (chồng lấn ≥ một nửa
exon ngắn hơn, **theo exon, không theo envelope**). Khác **tên** = `vocabulary`. Khác **exon** =
`conclusion`, **trừ** khi mọi đầu exon lệch ≤ **1 bp**: ấy là `boundary` — cùng gen, cùng chỗ, hai
phương pháp đặt đầu mút lệch nhau một base. Đo được, và là lý do của ngưỡng: **719/1.264** cặp khớp
lệch đúng 1 bp ở một hoặc hai đầu, 634 trong đó bộ gọi của dự án **bao trùm** bản công bố; và chạy
ARAGORN bằng tay trên bản vàng trả về **đúng** toạ độ bộ gọi xuất ra (`tRNA-Asp c[30540,30615]` so
với `complement(30541..30614)` của bản ghi) — tức đây là hai phương pháp cắt đầu tRNA, không phải
hai kết luận về hệ gen. Ngưỡng chọn **1 bp — mức chặt**, không phải mức làm con số đẹp: ở 3 bp
(một codon, quy ước codon dừng của CDS) bộ gọi thắng ngay cả khi chưa sửa gì, nên **thanh chuẩn
được đọc ở 1 bp** và con số 3 bp vẫn được in ra cạnh nó. Ngưỡng áp cho **cả hai** cách đọc như nhau.

**(4) Kết quả, cả hai cách đọc:**

| | kết luận | thiếu · thừa · độ dài | `boundary` | nếu tính mọi đầu mút (không ngưỡng) | vocabulary |
|---|---|---|---|---|---|
| **bộ gọi của dự án** | **811** | 47 · 41 · 723 | 602 | 1.413 | 373 |
| **PGA** | **588** | 39 · 4 · 545 | 111 | 699 | 149 |

**Thanh chuẩn KHÔNG đạt, ở cả hai cách đọc** (811 so với 588; không ngưỡng: 1.413 so với 699), và
**0/34 bản ghi** bộ gọi của dự án tốt hơn PGA (32 kém hơn, 2 hoà). Không có cách đọc nào của luật
đảo được kết luận này — nên §8.51 không phải một lần chỉnh luật để cho qua.

**(5) Ba phát hiện mới, đều đo được, đều chưa sửa.**

- **`lhbA` là một *call trùng*, không chỉ là một cái tên khác (H-18, nay có số):** artifact bản vàng
  gọi **cả `psbZ` lẫn `lhbA` ở đúng 35.895–36.083** — hai feature CDS, một locus, phạm luật "một
  locus một call" của P3B-08. Đây là **33 trong 41 call thừa** của cả corpus. Và nó làm **con số 86
  của bản vàng đúng một cách tình cờ**: bộ gọi có **85 locus + 1 bản trùng `lhbA` = 86**, bản công
  bố có **85 locus + `rps12` viết hai lần = 86**. Sửa H-18 (một mục trong `data/gene_aliases.yaml`)
  mà không sửa `rps12` thì bản vàng **rơi xuống 85 CDS** — hai việc này đi cùng nhau.
- **`rps12` bị gọi thiếu một exon và thiếu cả bản sao IR:** bản công bố viết `rps12` thành **hai
  feature, mỗi cái 3 exon** (exon 1 dùng chung ở LSC, hai exon trong mỗi bản IR); bộ gọi xuất **một
  feature, 2 exon** (68.976–69.089 + 97.786–98.016) — thiếu **exon giữa 26 bp (97.218–97.243)** và
  thiếu hẳn bản sao trong IRa. Đây là **12 call thiếu** trong corpus (PGA cũng thiếu: 11) và 21
  khác độ dài — gen trans-spliced nằm vắt qua ranh giới IR nên luật sao bản sao hiện tại
  (`_ir_side`: chỉ sao gen **nằm trọn** trong một bản IR) không chạm tới nó.
- **Exon đầu 6–40 bp: bốn gen, 136/242 ca lệch lớn nhất.** `petB` (exon đầu **6 bp**), `petD` (8),
  `rpl16` (9), `rps16` (40) — bộ gọi ra **đúng exon thứ hai** và bỏ exon đầu, vì exon 2–13 codon
  nằm dưới ngưỡng seed của blastx. Cộng `rps12` (21) · `rps19` (13) · `psbC` (18) · `rpl2` (8) ·
  `ndhK` (8) · `clpP1` (7), nhóm ">30 bp" là **242 ca**. Đây là việc lớn nhất còn lại: lấy **cấu
  trúc exon từ bản tham chiếu** thay vì suy từ alignment — đúng nguyên tắc "cấu trúc trước" của
  D5-08 và cùng gốc với **H-20**.

**(6) Phân rã độ lệch, để lần sau không phải đo lại:** bộ gọi — `rpl2` 68 · `trnA-UGC` 57 · `rps12`
34 · `rps16` 34 · `rpoC1` 34 · `petB` 34 · `petD` 34 · `rpl16` 34; PGA — `rpl2` 62 · `ycf2` 62 ·
`trnI-GAU` 44 · `rpoC2` 31 · `psbJ` 31 · `rpl16` 31 · `psbK` 29 · `ndhD` 28. Thiếu, bộ gọi: `ycf1`
24 · `rps12` 12; PGA: `ycf1` 23 · `rps12` 11 — **hai bộ gọi cùng thiếu `ycf1` và `rps12`** ở cùng
những bản ghi, tức hai locus ấy khó với cả hai phương pháp, không phải lỗi riêng của bộ gọi này.

### 8.52. ARAGORN's `i(a,b)`: trường thứ hai là **độ dài** intron, không phải toạ độ kết thúc — và cái trần 3.000 bp nói không gì cả (10/10/2026)

**(1) Lỗi, đo được từ chính con số của thước đo.** Bảng "extent theo locus" của §8.51 có ba dòng
đứng đầu là tRNA: `trnA-UGC` 57 · `trnK-UUU` 27 · `trnL-UAA` 34. Đọc thẳng ra: cả ba là tRNA **có
intron**. Bộ đọc `i(a,b)` của ARAGORN được viết hôm 09/10 (khi thêm `-i`) và đọc `b` như **toạ độ
kết thúc** của intron. Đo lại 10/10: đọc như thế, tRNA sau khi nối ra **104–134 nt** — không tRNA
nào dài thế; đọc `b` như **độ dài** intron thì ra **72–95 nt**, và hai ca khớp **chính xác** exon
của bản công bố: `trnL-UAA` 35 + 50 = **85** (bản ghi: `join(47453..47487, 47894..47943)` = 85) và
`trnA-UGC` 38 + 35 = **73** (bản ghi: 73). Thêm một chỗ dựa nữa: help của ARAGORN mô tả intron bằng
**độ dài** (`-i<min>,<max>`), không hề có toạ độ kết thúc.

**(2) Đúng loại "sai mà không báo" của dự án.** Số gen vẫn 38 tRNA, tên vẫn đúng, dòng "N genes
found" vẫn khớp — chỉ có **cấu trúc exon** sai, nên không test đếm nào bắt được; cái bắt được là
**bản ghi công bố** (join của nó) và **sinh học** (một tRNA 112 nt là bất khả). Test cũ
(`test_an_intron_splits_the_hit_into_the_two_exons_a_join_needs`) **khẳng định chính cách đọc sai**
ấy — nó được viết từ cùng một lần đọc nhầm, nên nó được **sửa theo đo lường**, kèm lý do và số, chứ
không phải nới ra; thêm test thứ hai
(`test_a_spliced_trna_comes_out_the_length_a_trna_can_be`) chốt cả ba ca vào **độ dài tRNA hợp lệ**.

**(3) Đo lại toàn bộ corpus sau khi sửa:** **811 → 692** khác biệt mức kết luận (mất đúng 119, bằng
57 + 27 + 34 + 1 của họ tRNA có intron). Bảng xếp hạng extent nay bắt đầu bằng `rpl2` 68 · `rps12`
34 · `rps16` 34 · `rpoC1` 34 · `petB` 34 · `petD` 34 · `rpl16` 34 · `trnP-UGG` 34. PGA (đo đủ 34
bản ghi lần này — lần trước `PX516262.1` chưa khớp được file PGA vì tên file là tên loài, nay khớp
theo **digest chuỗi**): **588**. **Thanh chuẩn vẫn KHÔNG đạt**: 692 so với 588; bộ gọi **tốt hơn ở
9/34** bản ghi (trước: 0), kém hơn 23, hoà 2.

**(4) Một rủi ro mới, có tên và có số — trần intron 3.000 bp của `-i`.** ARAGORN's `-i` chỉ tìm
intron **≤ 3.000 bp** và **không nói gì** khi bỏ qua (cùng hình dạng với `-max_target_seqs` cắt cụt,
§9.L-1). Đo trên 34 bản ghi: **266 lượt tRNA có intron**, intron dài nhất **2.671 bp** (`trnK-UUU`)
— chỉ cách trần **12 %**. Một loài có intron `trnK` dài hơn 3.000 bp sẽ **mất tRNA trong im lặng**.
Chưa đổi trần: đổi mà chưa đo thì lại là một thay đổi không đo được, và câu hỏi đúng là "nâng trần
có thêm call sai nào không" — ghi thành **H-22** kèm cách đo.

### 8.53. Bản vàng: 45 khác biệt so với bản công bố, từng cái một, thành test (10/10/2026)

Done-when của P3B-10 có hai vế; vế thứ nhất là *"mọi khác biệt so với bản công bố được liệt kê
(kind + span, không tóm tắt cái nào)"*. Nay nó là một test: `GOLDEN_DIFFERENCES` trong
`tests/golden/test_plastid_annotate.py` ghi **cả 45 khác biệt** — **18 mức kết luận · 19
`boundary` · 8 `vocabulary`** — và test so **tập hợp**, không so con số, nên một khác biệt đổi
verdict, đổi toạ độ, biến mất hay mới xuất hiện đều **fail kèm hai danh sách** (cái không còn thấy ·
cái chưa có trong danh sách). Đã **xem nó fail** theo cả hai chiều trước khi tin nó.

Bốn khiếm khuyết đã biết **nằm trong danh sách, không bị danh sách che**: exon đầu cụt của `petB`
(6 bp) · `petD` (8) · `rpl16` (9) · `rps16` (40); `rps12` thiếu bản sao thứ hai và exon giữa 26 bp;
`lhbA` (call trùng tại đúng `psbZ` 35.895–36.083 — và là lý do **86 của bản vàng đúng một cách tình
cờ**, xem §8.51); `rps19` −54 bp ở cả hai bản (**H-20**); tám tên tRNA ARAGORN đọc khác (**H-04**'s
họ). Mỗi cái có task riêng, ghi ngay trên danh sách — nên khi một cái được sửa, test sẽ nói **cái
gì đã đổi**, thay vì im lặng đi qua.

Sau bản sửa §8.52, danh sách tRNA có intron đã sạch: `trnL-UAA` và `trnA-UGC` **không còn** trong
danh sách (khớp bản công bố từng exon), `trnK-UUU` chỉ còn **1 bp mỗi đầu** (xếp vào `boundary`).

### 8.54. `lhbA` là call trùng, và `rps12` phải viết hai lần — hai lỗi bù nhau làm con số 86 đúng tình cờ (10/10/2026)

**(1) `lhbA` (H-18) — ba tín hiệu, đo được, khớp tuyệt đối.** Bản vàng `NC_050833.1` gọi locus
62 aa là `psbZ`; `NC_042827.1` gọi **cùng locus** là `lhbA`, và ba tín hiệu của §8.9 trùng khít:
**cùng chuỗi `/product`** ("photosystem II protein Z"), **cùng độ dài** 189 bp / 62 aa, **cùng
hàng xóm** `psbC … rps14`. Đếm trên bộ đã pin: **33 bản ghi viết `psbZ`, 1 bản ghi viết `lhbA`** —
đúng như §8.46 đã ghi. Không có entry alias thì bộ từ vựng tham chiếu mang **80 tên CDS cho 79
locus**, và bộ gọi — vốn "một locus một call" — **viết cả hai tên tại đúng 35.895–36.083** của bản
vàng: một **call trùng**, phạm chính Done-when của P3B-08. Đo sau khi thêm entry: **80 → 79 tên**,
`lhbA` biến mất. Test: `test_the_pinned_set_spells_every_cds_locus_once`.

**(2) `rps12` vắt qua ranh giới IR (luật sao bản sao).** Bản vàng có **đúng một** gen CDS vắt qua
ranh giới IR: `rps12` — exon 1 ở LSC (68.976–69.089), hai exon trong IRa. Luật cũ chỉ sao gen **nằm
trọn** trong một bên, nên gen này không được sao; bản công bố thì viết nó **hai lần**, mỗi feature
mang **cùng exon 1** và exon của bản sao mình. Luật nay: exon nào nằm trong IR thì được soi gương,
exon còn lại giữ nguyên — nên hai feature **dùng chung** exon 1, đúng như bản ghi. Đo lại bản vàng:
`rps12` **có** bản sao thứ hai, khớp bản công bố tới **3 bp** (138.836–139.066 so với
138.839–139.067), và vẫn thiếu **exon thứ ba 26 bp** (139.609–139.634) — cùng gốc với bốn exon đầu
cụt (dưới ngưỡng seed của blastx). Test đơn vị mới `tests/unit/test_ir_copies.py` chốt luật
trên plastome tổng hợp, **kèm ca đối chứng `ycf1`**: `ycf1` chạy từ SSC *vào* IR như **một** call,
và một call không phải bản sao thứ hai — bản vàng công bố có **một** `ycf1` và **hai** `rps12`; một
luật sao mọi gen chạm IR sẽ **phát minh** ra `ycf1` thứ hai.

**(3) Vì sao hai việc phải đi cùng nhau — đo được, không phải suy luận.** Sửa H-18 một mình đưa bản
vàng xuống **85 CDS**: con số 86 của nó đang đúng vì **một cái sai bù một cái sai** — artifact của
ta là **85 locus + `lhbA` gọi trùng = 86**, bản công bố là **85 locus + `rps12` viết hai lần = 86**.
Sau cả hai: **86 CDS = 86**, lần này vì đúng lý do, và test danh sách khác biệt (§8.53) tự báo đúng
ba dòng đổi: mất `lhbA` extra, mất `rps12` missing, thêm `rps12` extent. Danh sách còn **44**.

**(4) Cái giá trên corpus, đo được: 22/34 bản ghi viết `rps12` **một lần**.** Đo phân bố trên bộ
Iridaceae: **22 bản ghi** viết `rps12` một feature (3 exon, gọn trong một bản IR), **12 bản ghi**
viết hai feature — bản vàng thuộc nhóm 12. Luật "viết hai lần" theo hợp đồng của bản vàng
(CLAUDE.md: 132 = 86 + 38 + 8), nên trên 22 bản ghi kia bản sao của ta thành **extra**: corpus
692 → **681** (thiếu 47 → 34, thừa 41 → 30, độ dài 617). Đây là **khác biệt về cách trình bày**, và
nó đáng một câu hỏi riêng cho thước đo — nhưng **không đổi kết luận**: nếu phân loại nó là khác
biệt trình bày (bỏ 22 của ta) và áp y hệt cho PGA (bỏ 12 "missing rps12" của nó), thì **659 so với
576 — vẫn KHÔNG đạt**. Thanh chuẩn không đạt ở mọi cách đọc trung thực của luật, và đó là điều
đáng ghi lại: nó không phải một con số sát nút có thể lật bằng cách chỉnh luật.
