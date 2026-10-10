# Nhật ký đo lường — Giai đoạn 2B — bộ phân tích plastome, phần còn lại (M2.9 · M2.10 · M2.12)

**Phần này là `PROJECT.md` §8** (**§8.57–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.57` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.57. Span intron và IGS vào lớp dùng chung — π theo intron và theo IGS trên cặp vàng, và một cổng kiểm không thể fail (P2B-01, 10/10/2026)

**Việc.** Phép suy span rời `modules/introns.py` sang **`core/spans.py`** — đây là bản sao thật thứ
hai của đúng phép suy ấy (DRY chỉ kích hoạt từ bản thứ hai), vì `diversity` cần chính các ranh giới
đó để đo π theo intron và theo IGS. `core/spans.py` nay là nơi duy nhất định nghĩa: **intron** (khe
giữa hai part liên tiếp của một CDS, hai đầu cùng một vùng) · **junction trans-splicing** (khe ấy
vắt qua ranh giới vùng — báo ra, không bao giờ tính là intron) · **IGS** (khe giữa hai gen *được
gọi* liền kề, nằm gọn trong một vùng). `diversity` thêm hai bảng `diversity_introns` và
`diversity_igs`, và hai cột đếm vào `diversity_summary`.

**Ba lỗi tìm thấy khi đọc code cũ, sửa luôn trong cùng việc:**

1. **Phép tra vùng trộn hai hệ toạ độ.** `introns` và `diversity` so vị trí **hệ GenBank** với
   ranh giới vùng **hệ chuẩn hoá**. Bản vàng có `rotation_offset = 0` nên hai hệ trùng nhau — đúng
   vì may. Đã thêm `NormalizedFrame.region_at()` (chuyển hệ bên trong) và một test **xoay bản vàng
   7.500 bp** rồi kiểm nhãn vùng không đổi.
2. **`diversity` không gọi `ctx.require_ir()`** mà tự `detect_ir(reference)` — tức cổng IR của
   riêng nó không chạy qua cổng chung (seam 4). Nay gọi `ctx.require_ir()`, dùng đúng IR của lượt
   chạy.
3. **`_check_span` — một cổng kiểm không thể fail.** Phép kiểm "exons + gaps lấp đầy span" của
   `introns` (P2-05) **đúng với mọi cặp parts**: khe được tính bằng `min`/`max` nên phần chồng lấn
   tự triệt tiêu trong tổng. Đo 10/10/2026: nó chưa từng chạy, và không thể chạy. Cái nó để lọt
   đúng là ca nó sinh ra để bắt — parts chồng/lồng nhau cho "intron" **âm**, và độ dài âm ấy được
   **xuất ra như một con số**. Thay bằng phép kiểm khe ≥ 1 bp; đo trên **131 bản ghi / 2.850 gen
   nhiều phần** trong cây: đúng **1** ca (OX459128.1, mito thay thế, parts lồng nhau −138 bp), và
   ca đó nay dừng có tên thay vì in ra một intron âm (§9.H-6).

**IGS: định nghĩa D-g, và một phép đo buộc phải sửa định nghĩa.** Bản nháp đầu lấy **span ngoài**
cùng của gen (min→max của các part). Đo trên bản vàng thì nó **nuốt mất spacer thật**: `rps12` bị
ghi thành hai feature (chung exon 1 ở LSC, exon còn lại ở IR) nên span ngoài của nó trùm qua `clpP`,
và `trnK` (tRNA có intron) trùm qua `matK`. Bản chốt: **gen chiếm các exon của nó**, không chiếm
khoảng giữa chúng; khe giữa hai exon của *cùng một* feature không bao giờ là spacer (nó là intron
hoặc junction — việc của `introns_of`). Nhờ vậy ba spacer thật quay lại: `rps12-clpP` 138 bp,
`trnK-UUU-matK` 260 bp, `matK-trnK-UUU` 789 bp.

**Số đo — cặp vàng (`NC_050833.1` × `NC_024936.1`), MAFFT, `plastid run -m diversity`:**

| | Giá trị |
|---|---|
| Intron | **17** = 10 LSC · 3 IRa · 3 IRb · 1 SSC — **khớp từng dòng** với bảng `introns` (id · toạ độ · phase) |
| Junction trans-splicing | **2** (`rps12.1.1` · `rps12.2.1`) |
| IGS | **91** = 79 LSC · 12 SSC · **32.134 bp** · 30.925 cột dùng được |
| π theo IGS | 0,000000 – **0,082353** (`psbT-psbN`) · trung bình 0,014693 · theo cột 0,015618 · **12** spacer π = 0 |
| π theo intron (ngoài IR) | 11 locus · 0,004138 (`rpoC1.1.1`) – 0,012658 (`clpP.1.1`) · trung bình 0,009082 |
| π toàn alignment | 0,007181 — **không đổi** so với §8.25 |

**Ba cặp intron trong IR đo y hệt nhau** — đây là lý do đo được cho quyết định D-g: `rps12.1.2` và
`rps12.2.2` cùng 0,000000 trên 541 cột; `rpl2.1.1` và `rpl2.2.1` cùng 0,001506 trên 664 cột;
`ndhB.1.1` và `ndhB.2.1` cùng 0,000000 trên 699 cột. Hai bản sao IR là **một locus in hai lần**.

**Vì sao IGS loại hẳn IR, còn intron thì không.** Intron vẫn in đủ 17 dòng (cả hai bản sao IR) vì
bảng này phải **khớp từng dòng** với bảng `introns` — cross-check mà board yêu cầu — và cột `region`
nói rõ dòng nào là bản sao. Spacer thì không thể làm vậy: 18 cặp spacer của IR **không phải cặp nào
cũng đối xứng chính xác** (đo trên bản vàng: 3/18 cặp lệch — bản `ycf1` cụt không có CDS nên bên
IRa khoảng ấy là một spacer dài, hai spacer rRNA lệch đúng 1 bp), nên một luật "mỗi locus một lần"
theo phép đối xứng sẽ **im lặng** ở đúng ba chỗ ấy. Chọn cách nói được thành lời: **spacer chỉ in
ngoài IR**, phần bị loại **đếm và nêu tên** trong notes — 38 spacer trong IR · 2 khe vắt ranh giới
vùng (`rpl22–rps19` 108 bp · `rps19–psbA` 116 bp) · 4 cặp gen không có khe để đo (`psbD·psbC` −53 ·
`ndhK·ndhC` −10 · `atpE·atpB` −4 · hai feature `rps12` chung exon 1, −114).

**Kiểm chứng.** `make check` xanh. Golden: test mới chốt **từng dòng** của phép suy so với bảng
`introns` (17 id · toạ độ · phase · 2 junction), chốt số spacer (91 · 79 LSC · 12 SSC · 32.134 bp ·
38 trong IR · 2 vắt ranh giới · 4 không khe), và tier B chốt số π theo cặp vàng (17 intron · 91
spacer · ba cặp IR đo y hệt · `rpoC1.1.1` 0,004138 · `clpP.1.1` 0,012658 · `psbT-psbN` 0,082353).
Test tổng hợp (`tests/unit/test_spans.py`) đặt từng ca một cách có chủ ý: intron trong một vùng ·
junction vắt ranh giới · spacer có tRNA làm mốc · pseudogene **không** chia spacer · cặp chồng lấn ·
exon vắt qua gen khác · và ca âm (parts lồng nhau phải dừng).

---

### 8.58. Module `repeats` — hai engine độc lập cho cùng **9** match, và nguồn thật của `repeat-match` (P2B-02, 10/10/2026)

**Việc.** Mô-đun phân tích cuối cùng của bộ công cụ còn thiếu (**H-16**) nay có: `modules/repeats.py`
— **match chính xác cực đại** ≥ 30 bp, cả hai sợi, hai lớp theo đúng cách công cụ đánh dấu:
**dispersed** (bản sao thứ hai đọc xuôi) và **palindromic** (đọc trên sợi đảo). Engine chạy là
`repeat-match`; **Vmatch 2.3.1** là engine thứ hai độc lập (D-f), và nó **không** là engine của
mô-đun: Vmatch chỉ miễn phí học thuật, nên bản cài tối thiểu vẫn chạy được — `external/vmatch.py`
phục vụ đúng phép kiểm chéo, và bài kiểm chéo **skip có nêu tên** khi máy không có nó.

**Số đo — bản vàng, `repeat-match -n 30` (MUMmer 4.0.1), và Vmatch 2.3.1 (`-l 30 -d -p`):**

| | |
|---|---|
| Match ≥ 30 bp | **9** — **4 dispersed · 5 palindromic** (cả hai engine) |
| Dài nhất | **26.214 bp palindromic, IRa ↔ IRb — chính là cặp IR** (`83128–109341` ↔ `127511–153724`) |
| Còn lại | 42 bp ×4 (cụm trong IR: `88397/88418` ở IRa, `148393/148414` ở IRb, cộng hai match vắt IRa↔IRb) · 52 bp ×2 (kẹp tóc tự đối xứng ở LSC: 5.917 · 29.416) · 32 bp · 31 bp (SSC) |
| Hai engine | **trùng nhau từng dòng**: cùng 9 match, cùng toạ độ (sau khi đổi 1-based ↔ 0-based), cùng hai lớp |

**Nguồn thật của `repeat-match` — dòng §7.8 cũ đã sai, nay sửa tại chỗ.** §7.8 ghi nó "thuộc
MUMmer3 (tarball SourceForge)". Đo 10/10/2026: nhị phân đang chạy thuộc gói **`mummer4-4.0.1`**
(bioconda) — `conda-meta` ghi tên gói, và nó nằm cạnh `nucmer`. Hệ quả kéo theo: **`repeat-match`
không có cờ phiên bản** (`--version` là option không hợp lệ, in usage rồi exit 1) nên bản phát
hành phải đọc từ `mummer` cùng thư mục, và lệch pin thì **dừng** chứ không ghi chú.

**Toạ độ — cái bẫy đã khoá bằng test.** Hậu tố `r` của công cụ nằm ở **đầu của bản sao đọc trên
sợi đảo, trong khung reverse-complement** — không phải toạ độ xuôi. Cặp IR là phép kiểm: dòng
`83128 153724r 26214` phải đổi thành `[83127, 109341)` và `[127510, 153724)` = IRa và IRb. Đọc nó
như một đầu xuôi cho ra hai khoảng **trông vẫn hợp lý** nhưng sai chỗ — đúng kiểu lỗi im lặng mà
dự án này gặp nhiều lần, nên nó thành test riêng (`test_repeatmatch.py`, ví dụ có thật của bản vàng).

**Vì sao hai engine được đổi toạ độ ở **một** chỗ.** Vmatch in **0-based** và lớp bằng chữ `D`/`P`;
`repeat-match` in **1-based** và sợi bằng hậu tố `r`. Cả hai đổ về cùng một
`ExactMatch` 0-based nửa-mở, nên phép so là so trên cùng một tập — nếu mỗi bên tự đổi một kiểu,
"hai engine trùng nhau" sẽ là một câu không kiểm được.

**Cái bị loại, và vì sao.** Đầu ra của `repeat-match` là **match cực đại**, không phải locus: một
mảng tandem hiện ra thành nhiều dòng chồng lấn chứ không một dòng (nên `-t` và `-f` bị bỏ qua có
chủ đích), và lớp **tandem** vẫn là đất của `mreps`/`trf` (nửa còn lại của M2.12, còn mở). Con số
vì thế được phát biểu là "số match", không phải "số locus repeat" — scope label in đúng câu đó.

**Kiểm chứng.** `make check` xanh. Golden chốt: 9 match · 4/5 hai lớp · cặp IR (`83128–109341` ↔
`127511–153724`, IRa–IRb) · mọi dòng ≥ 30 bp · và **hai engine trùng tập** trong cùng một test
(`test_vmatch_agrees_with_repeat_match_on_the_golden`, skip có tên khi thiếu Vmatch). Unit: parser
với ví dụ thật của bản vàng (kể cả ca tự đối xứng `start1 == start2`), ca âm cho dòng không đọc
được và cho match dưới sàn, và module chạy với binary giả để chốt nhãn vùng + thứ tự + sàn.
`repeats` vào **hàng "chuyển tiếp"** của §13 (in vùng plastome + gọi `ctx.require_ir()`), nên câu
"bảy mô-đun" ở §13 và `CLAUDE.md` đã sửa thành **tám**.
