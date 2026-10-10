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
