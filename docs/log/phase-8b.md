# Nhật ký đo lường — Giai đoạn 8B — bộ dữ liệu Crataegus có đáp án công bố (M8B.1–M8B.8)

**Phần này là `PROJECT.md` §8** (**§8.65–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.65` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.65. Cờ compound của MISA — một cụm, hai cách viết, và bảng chỉ đọc được một (P8B-01, 10/10/2026)

**Việc.** Bộ **Crataegus answer-key** (anh Điềm đưa 10/10/2026: 100 bản ghi GenBank + bài Sun và
cs. 2026, *Front. Plant Sci.* 17:1767012, CC-BY, cùng **toàn bộ bảng kết quả của bài**) được đưa
vào cây làm việc ngay trong ngày (§8.66). Lần chạy thử đầu tiên — `plastid run -m ssr` trên
`PX413282` (*C. altaica*) — trả về 73 SSR, và cột `In compound formation` gắn **`no` cho cả 16
thành phần**, trong khi chính MISA báo *"Number of SSRs present in compound formation: 8"*.

**Đo được.** MISA 2.1 ghi một cụm compound bằng **hai cách**, và **một file có thể chứa cả hai**:

| Cách viết | Ví dụ (nguyên văn) | Parser cũ |
|---|---|---|
| Một dòng `repeat_region` liệt kê các thành phần | `Note=repeat_region,dinucleotide_repeat_microsatellite_feature,(TA)6,monomeric_repeat,(A)10` — **bản vàng Iris**, 48.659–48.679 | hiểu → `compound=True` |
| **Mỗi thành phần một dòng**, mở đầu bằng `compound_repeat` | `Note=compound_repeat,monomeric_repeat,(T)16` — **PX413282**, 211–226 | **không hiểu** → `compound=False` |

`_parse_gff` suy cờ từ "dòng này chứa >1 thành phần", nên mọi dòng chỉ có một thành phần đều ra
`no`. Hệ quả đo được:

| Bản ghi | Thành phần compound thật | Bảng cũ gắn cờ | MISA tự báo |
|---|---|---|---|
| `PX413282` (*C. altaica*) | **16** (8 cụm × 2) | **0** | 8 |
| `NC_050833.1` (bản vàng) | **6** (3 cụm × 2) | **2** | 3 |

Bản vàng **che mất lỗi**: nó chỉ có **một** cụm viết theo cách thứ nhất, nên bảng vẫn gắn cờ được
2 dòng và không ai nghi ngờ. Và phép đối chiếu GFF-vs-statistics của mô-đun **không bắt được** vì
nó so **số lượng**, không so **cờ** — đúng họ lỗi "một cột không ai kiểm là một cột không ai biết
nó sai" của §9.M-1, lần này ở một cột khác.

**Sửa.** `compound = len(components) > 1 or "compound_repeat" in fields[8]` — nhận cả hai cách
viết. Test hồi quy viết **trước** khi sửa (`tests/unit/test_misa.py`, fixture đã có sẵn dòng
`compound_repeat` từ 03/10 nhưng **chưa từng có assertion nào đọc cờ của nó**); test cũ
`test_a_compound_cluster_yields_one_feature_per_component` được **thu hẹp phạm vi** chứ không nới
lỏng — nó vốn khẳng định "MISA ghi một cụm thành một dòng", câu ấy nay đo được là **chưa đủ**, nên
nó chỉ còn nói về dòng `repeat_region` và trỏ sang test mới. Mô-đun lên **0.3.0** (cache khoá theo
`module_version` — H-08).

**Bản vàng không đổi số nào:** 37 = mono 31 · di 4 · tri 2, MISA báo 3 — chỉ cột cờ đổi từ 2 lên
**6** dòng. Và ghi chú của report nay nói rõ hai con số **đo hai thứ khác nhau**: MISA đếm **một
lần cho mỗi cụm** (3), bảng đếm **mọi thành phần** (6) — không con số nào suy ra từ con số kia.

**Vì sao đáng ghi hơn một lần so khớp trùng nhau:** đây là lỗi **thứ hai** trong hai ngày tìm ra
bởi một phép so với bên ngoài (lỗi thứ nhất: cột `Context` của `ssr` đọc span ngoài của gen
trans-spliced, §9.M-1). Cả hai đều nằm ở **một cột không có test nào đọc**, và cả hai đều đã tồn
tại qua nhiều lần chạy bản vàng mà không ai thấy.

### 8.66. Bộ Crataegus vào cây làm việc — kiểm toàn vẹn, đọc 100/100, và ba phép so đầu tiên (P8B-02, 10/10/2026)

**Việc.** anh Điềm đưa `/mnt/d/business/pts/thanhdiem_data/db_test_new`: bài báo full-text (PDF
22 MB, có text layer) · **7 bảng supplementary** (Table_S1…S7) + 5 Data Sheet · **100 bản ghi
GenBank** chia 4 nhóm (`groups.json`) kèm `manifest.tsv` và `fetch_gb.py` tải lại được. Đây là bộ
**thứ ba** anh ấy gửi và là bộ đầu tiên có **đáp án công bố kèm theo**.

**Đã làm.** Copy vào `.data/thanhdiem_data/db_test_new/` (117 file, 39 MB — chỗ làm việc nhanh,
M4a.0 đã đo `/mnt/d` chậm hơn ~100× ở các file nhỏ), rồi **đối chiếu sha256 từng file** với bản
gốc: **117/117 byte-identical**.

**Đo được, ngay trong lần kiểm đầu.**

| Phép đo | Kết quả |
|---|---|
| Đọc 100/100 bản ghi bằng reader của mình | **0 lỗi** |
| 50 plastome (18 mới + 26 công bố + 6 ngoài nhóm) | 159.638–159.973 bp · **tất cả** khai `/organelle=plastid:chloroplast` · **tất cả qua cổng IR ở 100,00 %** · CDS 74–87, tRNA 37, rRNA 8 |
| 50 bản ghi barcode (`PX676094–143`) | **543–625 bp, 0 gen, không IR** — đúng loại dữ liệu cho module `barcodes`, và **cổng IR chặn đúng** nếu đem chúng vào mô-đun plastome (đó là cổng đang làm việc, không phải lỗi dữ liệu) |
| **Table 1 của bài** (18 plastome mới: genome · LSC · SSC · IR · GC từng vùng) | **khớp từng dòng** — vd *C. altaica* 159.664 · 87.665 · 19.245 · 26.377 với GC 34,37/30,42/42,65; *C. phaenopyrum* IR **26.331** |
| Dải số trong **text** của bài (IR 26.311–26.396) | **26.311 không tái hiện được**: nó không có trong 18 plastome mới, không có trong 26 plastome công bố mình đo, và **không có trong chính Table 1 của bài** (cận dưới 26.331). Ghi lại làm phép đo, không phải lời buộc tội — quote bảng thì đúng, quote text thì sai 20 bp |
| **Table S6 (SSR)** trên *C. altaica* | bảng của họ = **57 dòng**; MISA của mình ở **mặc định (10,6,5,5,5,5)** cho 73, trong đó **57/57 dòng của họ khớp toạ độ**, 16 dòng còn lại là **thành phần của 8 cụm compound**. Chạy đúng ngưỡng **Methods của bài ghi** (10,5,4,3,3,3) cho **97** và chỉ khớp **55/57** ⇒ **Methods và Table S6 của bài không khớp nhau** (D-l: báo cáo cả hai) |
| **Table S5 (RSCU) · S7 (Ka/Ks) · S2 (barcode)** | có đủ cấu trúc để so; S7 thiếu **thiết kế so sánh**, S2 thiếu **mẫu thứ 51** — cả hai vào bảng "Waiting on data/people" của board |

**Kết luận của lần kiểm:** bộ dữ liệu **đủ** để bắt tay vào so — bốn mặt có đáp án số (cấu trúc/GC ·
RSCU · SSR · Ka/Ks) cộng barcode và bộ 49 genome cho π — và nó **đã trả tiền ngay** bằng một lỗi
thật của mình (§8.65). Đó là tiêu chuẩn đúng của một corpus: nó không chỉ để xác nhận.

### 8.67. Cấu trúc & GC so với Table 1 — 18/18 khớp kích thước, và một ô GC lệch 3 điểm (P8B-03, M8B.2, 10/10/2026)

**Việc.** Chạy `structure` + `genome_stats` trên **44 plastome** (18 mới + 26 đã công bố) và đặt
cạnh **Table 1** của bài — bảng in genome · LSC · SSC · IR · số gen · GC toàn hệ · **GC từng vùng**.

**Kết quả — bốn kích thước.**

| | |
|---|---|
| Dòng Table 1 máy đọc được từ text layer | 17 — và **17/17 khớp tuyệt đối** cả bốn số (genome · LSC · SSC · IR) |
| Dòng thứ 18 (*Crataegus sp.*, ô tên bị ngắt dòng nên regex bỏ sót) | đọc bằng tay ở dòng 388: **159.814 · 87.752 · 19.294 · 26.384** — và bản đo của mình cho `PX413297` là **đúng bốn số ấy** |
| **Tổng** | **18/18** |

Ví dụ nguyên văn: *C. altaica* 159.664 · 87.665 · 19.245 · 26.377 — bản đo của mình y hệt; và
*C. phaenopyrum* IR **26.331**, con số mà **text của bài ghi sai** (26.311, §8.66).

**Kết quả — GC từng vùng.** 17/18 dòng khớp tới mức in của bài (bài in 1–2 chữ số thập phân, mình
đo chính xác hơn: 34,4 so với 34,37 · 30,5 so với 30,46 — làm tròn, không phải lệch). **Một ô lệch
thật:**

| Dòng | Ô | Bài in | Mình đo |
|---|---|---|---|
| *C. jozana* (159.712) | GC **IR** | **45,62** | **42,64** |

Ba ô GC còn lại của **cùng dòng ấy** khớp tuyệt đối (36,64 · 34,37 · 30,55), và kích thước IR của
dòng ấy cũng khớp (26.384) — nên đây không phải khác định nghĩa vùng. Mọi plastome khác trong bảng
— kể cả 17 dòng còn lại — có GC IR nằm trong **42,62–42,67**; 45,62 nằm ngoài toàn bộ dải ấy. Ghi
lại **đúng như đo được**: một ô của bảng, lệch ~3 điểm, phần còn lại của dòng không lệch. Không kết
luận thay họ — chỉ nói nó không tái hiện được.

**26 plastome đã công bố không có bảng nào để so** — bài chỉ in Table 1 cho 18 mẫu mới. Bản đo của
mình trên nhóm ấy (IR 15.915–26.396 · LSC 87.464–88.307 · SSC 19.139–39.871) vì thế là **số của
mình, chưa có đối chứng**, và được ghi đúng như thế.

### 8.68. SSR so với Table S6 — 14/49 section khớp trọn, và hai nửa của bảng mã hoá cụm compound khác nhau (P8B-04, M8B.3, 10/10/2026)

**Việc.** Chạy `plastid run -m ssr` (MISA 2.1, **mặc định** 10,6,5,5,5,5 — bản đã pin, §7.8) trên
**50 plastome** của bộ, rồi đặt cạnh **Table S6** (49 section · **3.546 dòng**). Genome được khớp
**bằng chính toạ độ SSR của họ**, không bằng tên loài — vì bộ có năm bản ghi cùng tên
*Crataegus sp.* và nhiều cặp trùng kích thước, khớp theo tên sẽ lặng lẽ ghép sai cặp.

**Kết quả.**

| | |
|---|---|
| Section mà **mọi dòng của họ đều có trong bảng của mình** | **14/49** — đều là plastome **giải mới** trong nghiên cứu |
| Số dòng **thừa ra** ở những section ấy | 12–21 dòng, và chúng là **thành phần của các cụm compound** (đã kiểm ở §8.66 trên *C. altaica*: 73 − 16 = 57 = đúng bảng của họ) |
| Section **không** khớp | **35/49** — 4 plastome mới còn lại + **toàn bộ 26 plastome đã công bố** + 5 *Amelanchier* |

**Nguyên nhân của 35 section kia, đọc được từ ba dòng đầu tiên.** Ví dụ `NC085587`:

| | Ba SSR đầu |
|---|---|
| Bảng của họ | (10,19) · **(211,239)** · (1653,1662) |
| Bảng của mình | — · **(211,224) + (228,239)** · (1653,1662) |

Hai dòng `211–224` và `228–239` của mình cách nhau 4 bp ⇒ **một cụm compound**, và bảng của họ ghi
nó thành **một dòng duy nhất trải cả cụm**. Cùng dạng ở `OP964008` (241–273 so với 241–255 +
259–273) và `MK920297`. Tức **hai nửa của Table S6 mã hoá cụm compound theo hai cách khác nhau**:
nửa 18 plastome mới **bỏ hẳn** cụm, nửa còn lại **gộp cụm thành một dòng** — và nửa sau còn có
dòng mà MISA của mình **không hề tìm ra** (`NC085587` vị trí 10–19).

**Kết luận, nói đúng mức đã đo.** Ở **mặc định MISA**, bảng của mình **tái hiện trọn vẹn** phần SSR
đơn của 14 section plastome mới (mọi dòng của họ đều có, phần thừa là thành phần compound). Với 35
section còn lại thì **chưa so được**, và lý do không phải một tham số: phải biết **mỗi section do
lần chạy nào sinh ra** — nghi vấn đo được là nửa đã công bố lấy danh sách SSR từ **bài gốc của từng
bản ghi**, nơi công cụ và ngưỡng khác nhau. Đây là câu hỏi cho anh Điềm (bảng Waiting của board),
và **P8B-04 chỉ đóng được khi có câu trả lời** — không đóng bằng một con số trông hợp lý.

**Điều phép so này vẫn nói được ngay:** ngưỡng Methods của bài (10,5,4,3,3,3) cho 97 SSR trên
*C. altaica* và chỉ khớp 55/57 dòng của họ, còn **mặc định** cho 73 và khớp **57/57** — nên D-l có
câu trả lời đo được: **so bằng mặc định, và nói rõ thế**, chứ chạy theo ngưỡng họ *ghi* sẽ đi lệch
khỏi chính bảng họ *in*.

### 8.69. Ka/Ks so với Table S7 — bảng tái hiện từng chữ số, nhưng không phải so với *C. kansuensis* như text bài nói (P8B-08, M8B.5, 10/10/2026)

**Việc.** Table S7 = **42 PCG × 18 mẫu**. Bản đồ cột ↔ accession lấy từ Table S1 (13 cột tên loài
+ 5 cột mã `GSSZ · JRY · RR2H · RR3H · ZWSZ` = `PX413293–297`), không lấy từ tên trong bản ghi —
bộ có 5 bản ghi cùng tên *Crataegus sp.* và `PX413286` mang tên *C. chlorosarca* trong GenBank
nhưng *C. jozana* trong Table S1. Đây là phép so **chặt nhất** của giai đoạn: bảng in 6–8 chữ số,
và mục Waiting của board ghi thiết kế đã đóng vì Results nói *"relative to those of C. kansuensis"*
(`MF784433`, có trong bộ).

**Kết quả: bảng tái hiện được — và thiết kế đo được thì khác điều text nói.**

| | |
|---|---|
| Tham chiếu đo được | **`MN068257` — *Amelanchier cusickii*** (một bản ghi **ngoài nhóm** của chính bài) |
| Method đo được | **MLWL** (Tzeng, Pan & Li 2004) — **không** phải mặc định `MA` của công cụ |
| Công cụ, argv | **KaKs_Calculator 2.0.1** `-i <axt> -o <out> -m MLWL -c 11`, trên alignment **MAFFT** của từng cặp CDS đã cắt stop codon |
| Khớp | **629/714 ô so được — sai số 0 tới chữ số cuối cùng bảng in** |

**Đường tới thiết kế ấy: năm phép đo, không phải một lần thử.**

1. **`MF784433` cho số hoàn toàn khác**: rpoC2 của *C. monogyna* là **0,6667** (yn00 của mình) và
   **0,661981** (KaKs_Calculator `-m YN`) khi so với nó, còn bảng in **2,10854**.
2. **Cấu trúc bảng tự khai thiết kế**: 41/42 hàng của S7 **không đổi trong mỗi nhóm mẫu có CDS
   giống hệt nhau** (test md5 trên 18 mẫu × 42 gen; hàng phá vỡ duy nhất là `ycf1`) ⇒ giá trị là
   **hàm của trình tự mẫu** — đúng dạng của một phép so với **một tham chiếu cố định**, không phải
   trung bình hay so từng cặp một.
3. **`MF784433` bị chính bảng loại**: rpoB của *C. altaica* · GSSZ · ZWSZ **byte-identical** với
   rpoB của nó (md5 `eb99ab2c6eba` trên 3.207 bp), mà bảng cho ba mẫu ấy **0,1124499**; so với một
   tham chiếu giống hệt thì công cụ in **NA** (đo được: cặp giống hệt → Ka = Ks = 0 → NA). Cùng
   chuyện ở `matK`, nơi GSSZ và ZWSZ dùng chung chuỗi với `MF784433` và vẫn có giá trị.
4. **Quét 99 bản ghi** của bộ làm tham chiếu (method YN, ba mẫu thử): **không bản ghi nào** tái hiện.
5. **5 bản *Amelanchier* × 8 nhóm rpoC2 × 5 method**: **MLWL với `MN068257` khớp cả 8 nhóm**;
   `MN068255` khớp tới `ndhC`/`ccsA` thì lệch — hai bản Amelanchier khác nhau ở **6 gen** (mỗi gen
   1 base) và ở `ycf1`, và chính chỗ khác nhau ấy chỉ ra `MN068257`.

**42 hàng × 18 cột = 756 ô, phân loại đủ.**

| Số ô | Loại |
|---|---|
| **629** | khớp tuyệt đối, sai số 0 — trong đó **347 ô hai bên cùng bằng 0** và **282 ô giá trị thật** |
| **77** | bảng ghi `0` khi phép tính của mình **không xác định** (`Ks = 0` → công cụ in `NA`). Bảng của họ **không có ô NA nào** ⇒ NA được viết thành 0. Đo được **cả hai** quy ước in của công cụ: `Ka = 0` → `-0` (Excel đọc thành 0), `Ks = 0` → `NA` |
| **3** | `ycf1` lệch ~2×10⁻⁴ — và `ycf1` đúng là hàng duy nhất phá cấu trúc nhóm ở (2) |
| **2** | lệch thật: `clpP` của *M. germanica* **0,474917 vs 0,472932** (0,4 %; CDS của nó 591 bp so với 588 bp của 17 mẫu kia) · `ycf1` của *M. germanica* **0,381941 vs 0,344616** (10,8 %) |
| **42** | không so được: **18** vì bản ghi `MN068257` **không có `ndhI`** · **24** vì bản ghi mẫu **không có gen** (`ndhB` 11 · `ycf1` 11 · `matK` 1 · `infA` 1) |

**Điều dòng cuối vừa lộ ra — một phát hiện về chính bộ dữ liệu, không phải về bài.** S7 có giá trị
`ndhB` cho **11 mẫu mà bản ghi đã nộp không hề có `ndhB`**: **12/18 bản ghi mới thiếu trọn nhóm gen
IR** (`ycf1` · `ycf2` · `ndhB` · `rpl2` · `rpl23` · `rps7` — 11 bản ghi ở 74–75 CDS, 7 bản ghi đủ
83–85), một số thiếu thêm `matK` · `psbL` · `infA`. Tức **annotation đem phân tích không phải
annotation đã nộp**. Ghi lại đúng thế; nó là **đầu vào của P8B-07** (so annotation), không phải
một kết luận thay bài.

**Phiên bản — ghi lại, không làm mượt.** Methods của bài ghi **KaKs_Calculator v3.0**; bản ở đây là
**2.0.1** và nó tái hiện **đúng từng chữ số** bảng của họ: với method MLWL hai bản không khác nhau
ở mức in ra. Đây là chỗ hiếm hoi phiên bản công cụ *không* tạo khác biệt — ngược với §9.B-14
(MAFFT, nơi hai bản cho hai số Ts/Tv).

**Bên mình nói gì trên cùng cặp ấy.** Module `kaks` dùng **PAML yn00** (D4-01) — một estimator
khác, nên số không trùng và đó **không phải lỗi**: rpoC2 *C. monogyna* vs Amelanchier: yn00
**1,5263** so với MLWL **2,10854**; rpoB *C. altaica*: **0,1231** so với **0,112450**. Hai lần chạy
đầy đủ (80 gen × 18 mẫu, mỗi lần ~4 phút) nằm ở `.runs/phase8b/p8b08-vs-amelanchier/` và
`.runs/phase8b/p8b08-vs-kansuensis/`.

**Và bốn loài mà text bài gọi tên thì tái hiện ngay cả bằng estimator của mình.** Bài nói rpoC2
vượt 1 ở *C. monogyna* · *C. laevigata* · *C. songarica* · *M. germanica*; so với **Amelanchier**,
yn00 của mình cho **đúng bốn loài ấy** vượt 1 (1,5263 · 1,3684 · 1,3684 · 1,3684) và **không loài
nào khác**; còn so với ***C. kansuensis*** — đúng như text nói — thì **không mẫu nào** vượt 1
(pooled Ka/Ks cả 18 mẫu: **0,15–0,27**). Nghĩa là: **thứ hạng tái hiện, độ lớn thì không** — độ lớn
phụ thuộc estimator (yn00 vs MLWL) và tham chiếu (trong chi vs khác chi). Ghi lại như phép đo;
việc bài chọn tham chiếu nào để gọi tên "chọn lọc dương" là việc của bài.

### 8.70. RSCU so với Table S5 — hai bên lệch nhau ~1–2 %, và đó là *quy tắc đếm* chứ không phải dữ liệu (P8B-05, M8B.4, 10/10/2026)

**Việc.** 18 plastome chạy `plastid run -m codon_usage` (`.runs/phase8b/p8b05-px41328*…299/`), đặt
từng codon cạnh **Table S5** (64 dòng codon — 61 sense + 3 stop; bảng của mình in 61 codon sense,
stop bị loại theo scope). Hai bên in **luật tập gen** của mình: của mình là *mọi CDS, pseudogene
loại, **cả hai bản IR** được đếm, stop và codon cụt không đếm*; của họ, Methods chỉ nói
"PhyloSuite v1.2.2 … screen the PCGs" rồi vẽ bằng một web tool (`RSCU-Plot`) — **không nói tập gen,
không nói quy tắc đếm**.

**Kết quả.**

| Luật đếm | max &#124;ΔRSCU&#124; | trung bình &#124;Δ&#124; | trung bình tương đối |
|---|---|---|---|
| **cả hai bản IR** (luật của mình) | **0,1039** | 0,01475 | 1,98 % |
| **mỗi tên gen một lần** | **0,1018** | 0,01055 | 1,40 % |

Lệch **có hệ thống**, không phải nhiễu: phần lớn 61 codon lệch **cùng dấu ở cả 18 mẫu** — mạnh
nhất là `TTA` −0,046 · `AGT` −0,041 · `TCA` +0,041 · `CGT` −0,033 · `GTC` +0,031 · `GTA` −0,030 ·
`CTG` +0,029 · `CGA` +0,029 (dấu = mình − họ). Và nó **phân lớp theo bản ghi**: 11 bản ghi ở
74–75 CDS lệch trung bình **1,55 %**, 7 bản ghi đủ 83–85 CDS lệch **2,6 %** — đếm mỗi tên gen một
lần kéo nhóm đủ gen xuống **1,10 %** (C. crus-galli 0,0201 → 0,0083). Tức **quy tắc đếm gen lặp
của IR là một phần nguyên nhân, không phải toàn bộ**.

**Đã thử và loại, mỗi giả thuyết một số đo** (max &#124;Δ&#124; trên cả 18 mẫu):

| Giả thuyết | max &#124;Δ&#124; |
|---|---|
| tập CDS của chính bản ghi (luật của mình) | 0,1039 |
| bỏ `ycf15` (và bỏ luôn bản `rps12` thứ hai) | 0,1039 / 0,1092 |
| chỉ 42 gen của Table S7 | 0,1253 |
| tập gen **chung** của 18 mẫu (70 gen) | 0,0929 (trung bình của max từng mẫu) |
| cộng 8 gen mà bản ghi thiếu (`ndhB` · `rpl2` · `rpl23` · `rps7` · `ycf2` · `ycf1` · `matK` · `psbL`) | 0,2511 |
| cộng tRNA / cộng rRNA / cộng cả hai | 0,1070 / 0,2290 / 0,2686 |
| span CDS **kèm intron** | 0,8551 |
| **cả bộ gen** (frame 1) | 0,8493 |
| RSCU trung bình **không trọng số** theo gen | 0,3713 |
| bỏ **một** gen bất kỳ (74 phép thử) | tốt nhất 0,0331 (`rpl33`) — không gen nào giải thích được |
| **cộng** một gen bất kỳ của cả 99 bản ghi (300–3.000 bp) | tốt nhất 0,0350 — không gen nào giải thích được |

**Vì sao biết khác biệt nằm ở quy tắc đếm chứ không ở dữ liệu.** Chính **những CDS đã nộp ấy**
tái hiện Table S7 (Ka/Ks) **đúng từng chữ số** (§8.69) ⇒ trình tự hai bên là một. Và phía mình đã
được kiểm độc lập: trên bản vàng, cpstools và mình **đếm trùng từng codon** (21.861 codon sense,
0 lệch; RSCU lệch tối đa 0,0049 — đúng mức làm tròn 2 chữ số của họ, §8.61).

**Kết luận, nói đúng mức đo được.** Xu hướng mà bài đọc ra từ bảng — codon kết thúc A/T chiếm ưu
thế — **tái hiện** ở cả 18 mẫu (lệch 1–2 %, cùng dấu); nhưng bảng **không tái hiện tới chữ số in**
như ba mặt kia đã làm (Table 1 từng dòng · S6 ở mặc định MISA 57/57 · S7 629/714 ô sai số 0). Câu
hỏi còn lại là **quy tắc đếm của họ** (tập gen nào, IR đếm mấy lần), và nó nằm trong bảng
"Waiting on data / people" của board — **P8B-05 chỉ đóng bằng một câu trả lời, không bằng một con
số trông hợp lý**.

### 8.71. Annotation so với bản công bố của GeSeq — 341 khác biệt kết luận trên 18 bản ghi, và chúng đến từ đâu (P8B-07, M8B.7, 10/10/2026)

**Việc.** 18 bản ghi của bộ được annotate **leave-one-out** bằng chính `plastid annotate` (17 bản ghi
còn lại làm bộ tham chiếu, bản ghi đích bị loại bằng tên), rồi phân loại bằng **đúng luật của
`tools/corpus_bar.py`** (import `compare`/`summarise`/`by_gene` từ chính công cụ ấy — một luật, một
chỗ). 18/18 chạy xong, exit 0.

**Một chi tiết cơ chế, đo được khi chạy.** Manifest leave-one-out đã ghi trong cây là của bộ
**Iridaceae**; chạy với `--reference .panel/crataegus-18` mà không ghi lại manifest thì **mọi bản ghi
bị loại** và run dừng bằng một câu nêu tên (`every record … is excluded … leave-one-out leaves
nothing to search`) — tức cây chỉ giữ **một** thiết kế mỗi lúc. Đã ghi manifest Crataegus
(`plastid reference loo --family Rosaceae`), chạy, rồi **khôi phục** manifest Iridaceae
(`git status` sạch sau đó). Bản ghi đích bị loại **bằng tên** trong từng run, và tên ấy nằm trong
manifest của run.

**Kết quả — 18 bản ghi, 4.448 feature công bố so với 2.339 feature được gọi.**

| | |
|---|---|
| **Kết luận** | **341** — thiếu **21** · thừa **136** · extent **184** |
| Biên (trong 1 bp) | **288** — báo cáo, không tính vào bar |
| Từ vựng | **305** |
| Chặt (mọi đầu mút xê dịch đều tính là kết luận) | **629** |
| Mỗi bản ghi | 12–23 khác biệt kết luận (trung bình 18,9) |

**Ba thành phần, ba nguyên nhân — và hai trong ba nguyên nhân không phải lỗi của mình.**

1. **Thừa 136 — 110 là gen IR mà bản ghi đã nộp KHÔNG CÓ.** `rpl2` · `rpl23` · `ycf2` · `ndhB` ·
   `rps7` mỗi gen 22 lượt (= 11 bản ghi × 2 bản IR); cộng `ycf15` 14 · `ycf1` 11 · `infA` 1.
   12/18 bản ghi thiếu trọn nhóm gen IR (§8.69 đã đo từ phía S7), và bộ tham chiếu 17 bản ghi thì
   có — nên bộ **gọi đúng** những gen mà bản công bố không có. Đây là phát hiện về **bộ dữ liệu**,
   đã nêu ở §8.69.
2. **Extent 184 — 105 là SÁU gen của MỘT bản ghi tham chiếu.** `psbK` · `psbM` · `clpP` · `rpl32` ·
   `ndhG` · `psbJ` xuất hiện 15–18 lần mỗi gen; đo trên cả 18 bản ghi thì **`PX413299`
   (*M. germanica*) annotate cả sáu gen dài hơn 17 bản ghi kia** (`psbM` **270 bp so với 105**,
   `rpl32` 201 vs 159, `ndhG` 573 vs 531, `psbJ` 144 vs 123, `psbK` 195 vs 186, `clpP` 591 vs 588;
   nguyên văn trong bản ghi: `CDS complement(31219..31488) /gene="psbM"`). Engine **chuyển trung
   thành** cấu trúc exon của tham chiếu — đúng như thiết kế — nên một bản ghi tham chiếu annotate
   lệch sẽ lan sang **cả 18** run. Đây là **giới hạn đã biết** (P3B-10 ghi ở Residual risks: "the
   transfer is faithful to the chosen reference"), nay có số: **105/184 extent**.
3. **Thiếu 21 — 18 là `trnV-UAC`, và đây là một khiếm khuyết CỦA MÌNH, đo được.** Không run nào gọi
   được `trnV-UAC`, trong khi bản công bố có nó ở **cả 18** bản ghi. Truy nguyên: ARAGORN **có** gọi
   một tRNA tại đúng locus ấy nhưng **đọc sai danh tính** — trong bộ Crataegus nó gọi
   `tRNA-Stop(tta)` (91 base, mạch bổ sung, `c[54467,55134]`), và **trên bản vàng nó gọi
   `tRNA-Cys(aca)`** tại `c[51598,52275]` (chính là `trnV-UAC` của bản công bố). Bộ từ vựng không
   có tên cho "Stop", nên call bị bỏ và locus thành *không gọi được*: summary của run nói đúng
   **"Vocabulary loci not called: 1"**. Đây là **H-24** — ghi tên, chưa sửa (xem dưới).
4. **Từ vựng 305 — toàn bộ là tầng đặt tên, không phải tầng vị trí**: cùng toạ độ, khác tên.
   `trnI-CAU`→`trnM-CAU` 36 · `trnI-GAU`→`trnE-UUC` 36 · bốn rRNA (`rrn16S`→`rrn16` …) 36 mỗi cái ·
   `ndhI`→`ndhL` 18 (đúng cặp alias §8.31) · `trnfM-CAU`→`trnM-CAU` 18 · họ Ser/Gly
   (`trnG-GCC`→`trnS-CGA` 17 · `trnG-UCC`→`trnG-GCC` 17 · `trnS-GCU`↔`trnS-GGA` 16/2) — đúng họ
   H-04 mà §8.41 đo: ARAGORN đọc anticodon theo cách của nó, GeSeq theo cách của họ.

**Đặt cạnh bar của họ Iridaceae (§8.55), nói đúng chỗ khác nhau.** Bên ấy: **34 bản ghi**, mình
**418** khác biệt kết luận so với **588 của PGA**, chặt 1.070 so với 699. Bên này: **18 bản ghi**,
**341** khác biệt kết luận, chặt 629 — và **không có PGA** để so, vì bản công bố ở đây *là* GeSeq,
tức chính cái đáp án. Nên câu "hơn/kém PGA" không tồn tại ở bộ này; con số so được là **mật độ**:
**18,9 khác biệt kết luận mỗi bản ghi (Crataegus) so với 12,3 (Iridaceae)** — và mật độ ấy cao hơn
vì hai lý do đã gọi tên ở trên (bản ghi thiếu gen · một tham chiếu annotate lệch), không phải vì
engine đọc kém hơn.

**H-24, khiếm khuyết mở, ghi kèm bằng chứng.** `trnV-UAC` — một tRNA **có intron** — không được gọi
ở **18/18** bản ghi này và **cũng không** ở bản vàng (ARAGORN gọi `tRNA-Cys(aca)` ở locus ấy);
run báo đúng ở dòng "Vocabulary loci not called", nên nó không im lặng, nhưng một locus từ vựng
không gọi được là một khoảng trống thật. Ba hướng sửa, chưa chọn: chuyển tRNA từ tham chiếu như
CDS/rRNA đang làm · gọi bằng engine thứ hai · hay để nguyên và ghi. Chưa đủ dữ liệu để chọn trong
giai đoạn này.

### 8.72. π, vùng hotspot và 50 bản ghi barcode (P8B-06, M8B.6, 10/10/2026)

**Việc.** Hai nửa, hai câu hỏi khác nhau. **Nửa π**: chạy `diversity` (cửa sổ **600/200** — đúng
tham số DnaSP v6 của bài, và cũng là mặc định của mình) trên **44 plastome** (18 mới + 26 công bố);
bài **không in số π**, chỉ **gọi tên ba vùng đa dạng nhất** — nên thứ so được là **thứ hạng**.
**Nửa barcode**: 50 bản ghi `PX676094–143` so với **Table S2** (độ dài · GC) rồi đo π trong locus.

**Nửa π — 44 plastome, một alignment 166.154 cột, 157.073 cột dùng được.**

| | |
|---|---|
| π toàn alignment | **0,001633** |
| LSC · SSC · IRa · IRb | 0,002044 · 0,003533 · 0,000289 · 0,000249 |
| Spacer có giá trị | **91** |

**Thứ hạng, so với ba vùng bài gọi tên.** Bài nói ba vùng đa dạng nhất là `ndhC~trnV-UAC` ·
`psbZ~trnG-UCC` · `trnR-UCU~atpA`. Bản đo của mình trên 91 spacer:

| Vùng của bài | Hạng của mình | π |
|---|---|---|
| `trnR-UCU-atpA` | **4**/91 | 0,00839 |
| `psbZ-trnG-UCC` | **5**/91 | 0,00776 |
| `ndhC-trnV-UAC` | **9**/91 | 0,00742 |

Cả ba **đều nằm trong top 10** của mình, nhưng **không vùng nào đứng đầu**: vị trí số 1 là
**`petG-trnW-CCA`** với π **0,02042** — hơn vùng nhì (`rpl32-trnL-UAG`, 0,00874) **hơn hai lần**.
Khác biệt phạm vi phải nói kèm: bài lấy thứ hạng từ **cửa sổ trượt trên cả bộ gen** (DnaSP, 600/200),
mình lấy từ **π của từng spacer** trên cùng 44 bộ gen — hai cách đo gần nhau nhưng không phải một.
Ghi đúng thế: **ba vùng của bài tái hiện ở mức "trong top 10", không tái hiện ở mức "đứng đầu"**.

**Nửa barcode — 50/50 tái hiện tuyệt đối, và mẫu thứ 51 nói đúng tên.**

| | |
|---|---|
| **Độ dài** so với Table S2 | **50/50 khớp** (543–652 bp) |
| **GC** so với Table S2 | **50/50 khớp tới chữ số in** — lệch lớn nhất **0,0000** (bảng in 4–14 chữ số) |
| Dòng thiếu | `seq45` *Crataegus scabrifolia* isolate YNSZ3 — GenBank ID `-`, độ dài `----`: **mẫu thứ 51 chưa nộp**, đúng như Data availability statement nói (`PX676094–PX676143`) |
| π trong locus (50 bản ghi) | **0,018759** trên **455 cột dùng được** của alignment **806 cột** |
| Cùng spacer ấy trên 44 plastome | **0,00742** — tức mức lấy mẫu dày hơn (nhiều cá thể mỗi loài) cho π **cao hơn ~2,5 lần**, đúng tiền đề của bài |

**Một chỗ cơ chế, đo được, và nó sửa lại một dòng trong bảng Pitfalls của board.** Bảng ấy viết
"the 50 barcode records fed to a plastome module… use `barcodes`/`diversity`, never `structure`".
Đo lại: **cả `diversity` lẫn `barcodes` đều khai `requires_ir`**, và CLI bắt **mọi** genome trong run
phải qua cổng IR — nên 50 amplicon không IR bị **cả hai** chặn (`IR gate: PX676094.1: no valid IR
pair found`). Đó là cổng đang làm việc đúng, nhưng lời khuyên "dùng `diversity`" thì sai với đầu vào
amplicon; π ở đây vì thế được đo bằng **đúng hàm lõi mà mô-đun gọi**
(`core.align.nucleotide_diversity`) trên alignment của **đúng adapter MAFFT** (`external.mafft`) —
không phải một đường tính thứ hai. Dòng Pitfalls đã sửa tại chỗ.

### 8.73. SSR so với Table S6, lần hai: **cả hai ngưỡng trên cả 49 section** (P8B-04 tiếp, 10/10/2026)

**Việc.** §8.68 so ở **một** cài đặt và chỉ trên các section so được. Lần này chạy MISA (bản đã
vendored, dùng **đúng `_parse_gff` của mô-đun** — cùng cách đọc, chỉ đổi ngưỡng) trên **cả 50
plastome** ở **hai** cài đặt: **mặc định** `(1/10)(2/6)(3/5)(4/5)(5/5)(6/5)` và **ngưỡng Methods
của bài ghi** `(1/10)(2/5)(3/4)(4/3)(5/3)(6/3)`; mỗi section khớp với **một** genome theo toạ độ SSR
của họ (không theo tên — 5 bản ghi cùng tên *Crataegus sp.*), và **cùng genome ấy được đọc cho cả
hai cài đặt**, nên hai cột dưới đây so được với nhau.

| Cài đặt | Dòng của họ có mặt **đúng toạ độ** | Section **trọn vẹn** |
|---|---|---|
| **Mặc định** (10,6,5,5,5,5) | **2.581/3.546 — 72,8 %** | **14/49** |
| **Ngưỡng họ ghi** (10,5,4,3,3,3) | **3.070/3.546 — 86,6 %** | **4/49** |

**Hai con số này không mâu thuẫn — chúng trả lời hai câu hỏi khác nhau.** Ngưỡng thấp hơn ⇒ bảng
của mình **nhiều dòng hơn** (ZWSZ: 75 → 103 dòng) ⇒ *tìm thấy* nhiều dòng của họ hơn, **và** thừa
ra nhiều hơn ⇒ ít section trọn vẹn hơn. Nên: **14 section tái hiện trọn vẹn dưới mặc định** — và
**cả 14 đều thuộc nửa 18 plastome mới** — còn **tỉ lệ dòng phủ** thì cao hơn ở ngưỡng họ ghi.

**Cơ chế của phần không khớp, nay đo được ở mức từng dòng.** Section `Crataegus sp. ZWSZ`: 65 dòng
của họ, **56 khớp đúng toạ độ**, **65/65 khớp ở toạ độ *đầu*** — và **9/9 dòng lệch là một cụm
compound viết thành một dòng trải cả cụm**:

| Bảng của họ | Bảng của mình |
|---|---|
| `(242, 267)` | `(242,252)` + `(256,267)` |
| `(54421, 54488)` | `(54421,54430)` + `(54439,54449)` + `(54477,54488)` |
| `(52070, 52122)` | `(52070,52080)` + `(52113,52122)` |

Tức khi họ gộp cụm, **đầu trùng, cuối lệch** — đúng như §8.68 suy từ ba dòng đầu, nay đúng cho cả
9 dòng. Không dòng nào của họ là một SSR mình không có.

**Kết luận của P8B-04, ba điều Done-when hỏi.** (1) **Cài đặt nào tái hiện S6**: **mặc định** —
đó là cài đặt dưới đó các *section* tái hiện trọn vẹn, và nó khớp với điều §8.68 đã đo trên
*C. altaica* (57/57 so với 55/57). (2) **Mức khớp trên toàn bộ genome**: 72,8 % dòng đúng toạ độ ở
mặc định, 86,6 % ở ngưỡng họ ghi — cả hai trên **49/49 section**, không bỏ section nào. (3)
**Tên các cụm compound**: bảng của mình liệt **mọi thành phần**, bảng của họ (ở nửa đã công bố)
**gộp cụm thành một dòng trải**, còn ở nửa 18 plastome mới thì **không liệt thành phần** (73 dòng
của mình so với 57 của họ trên *C. altaica*: 16 dòng thừa là thành phần của 8 cụm).

**Còn lại, và nó không còn chặn phép so.** Câu hỏi "mỗi section của S6 do lần chạy nào sinh ra" vẫn
là câu hỏi cho anh Điềm — nhưng nay nó **không còn là điều kiện để so**: cả 49 section đã được so ở
cả hai cài đặt, và phần không khớp đã có nguyên nhân đo được (mã hoá cụm) chứ không phải một khoảng
trắng. Giữ nó lại trong bảng Waiting như một dữ kiện có thể làm phép so **chặt hơn**, không phải
như một việc phải xong.
