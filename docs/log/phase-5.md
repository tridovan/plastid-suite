# Nhật ký đo lường — Giai đoạn 5 — hình & báo cáo, và các số đo cho D5-08

**Phần này là `PROJECT.md` §8** (**§8.34–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.34` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.34. Ba số đo cho annotator lục lạp — gen chồng lấn, chiều SSC, và độ sạch của bộ tham chiếu (D5-08, 05/10/2026)

**Bối cảnh.** Trong lúc chốt thiết kế D5-08 (annotator lục lạp tự viết) có ba câu hỏi quyết định
luật của công cụ, và cả ba đo được ngay trên dữ liệu đang có — không gọi công cụ ngoài, không cần
mạng.

**1. Genome vàng có 4 cặp gen chồng lấn — "chồng lấn" KHÔNG phải lỗi.** Đo bằng reader của mình
(`core/genbank.py` + `core/model.py`) trên `tests/data/NC_050833.1.gb`: 132 gen được đếm, giao
từng phần của từng cặp gen:

| Cặp | Chồng lấn |
|---|---|
| `psbD` × `psbC` | 53 bp |
| `rps12` × `rps12` | 114 bp (hai feature cùng tên — gen trans-spliced) |
| `ndhK` × `ndhC` | 10 bp |
| `atpE` × `atpB` | 4 bp |

Cả bốn cặp nằm trong hợp đồng **132 = 86 + 38 + 8**. ⇒ Luật trọng tài của annotator mới **giữ cả
hai + báo**; luật "điểm cao hơn thắng, bỏ cái kia" sẽ xoá ba gen thật và một bản sao `rps12` ngay
trên bản ghi chuẩn.

**2. Chiều SSC trên 34 bản công bố: 32 cùng chiều, 2 ngược — và hai ca ngược là flip-flop thật,
không phải lỗi.** Đo trên bộ 01 (`.data/thanhdiem_data/01_NCBI_tai_truc_tiep/gb`), từng bản:
`read_genbank` → `detect_ir` → `enforce_ir_gate` → `normalize`, rồi so **chiều tương đối bên trong
chính bản ghi** — strand đa số của các CDS có midpoint nằm trong vùng LSC, so với vùng SSC.
`normalize` chỉ **xoay gốc**, không lật chiều — nên so toạ độ thô giữa hai bản là vô nghĩa, phải so
quan hệ trong-bản. Kết quả: **32/34 "same" · 2/34 "opposite"** — `Gladiolus_x_gandavensis_NC_062351.1`
và `Iris_missouriensis_NC_042827.1`. Cả 34 bản đều qua cổng IR ở lần đo này (khớp §8.31: bộ công
bố **0** rớt cổng).

Kiểm tiếp câu hỏi "hai ca ngược có phải lỗi chú giải không": chạy `core/orf.check_genome` trên CDS
vùng SSC của cả hai, đối chứng là genome vàng. Kết quả: mỗi bản **12 CDS SSC**, trạng thái chỉ có
`ok`/`alternative_start` — **y hệt bản đối chứng, 0 lỗi**. Nếu cờ chiều bị ghi sai thì chuỗi đọc ra
là reverse-complement và dịch thành rác; sạch như vậy nghĩa là chiều *đúng theo chính bản ghi* —
SSC **thật sự đảo** so với phần còn lại của nó. Đây là **flip-flop SSC** đã biết của plastome. ⇒
Audit reference của D5-08 **báo, không chặn** ca này — chặn sẽ loại hai plastome công bố thật, một
trong đó là một bản *Iris* (chi bản thảo).

**3. 34 bản công bố làm chuẩn thì sạch — trừ đúng một bản.** Đo bằng BioPython `SeqIO` đọc thẳng
từng file của bộ 01: CDS thiếu `/translation`: **0/34** — nên `blastx` trên translation công bố là
khả thi cho cả bộ; bản ghi có ít nhất một feature gen (CDS/tRNA/rRNA) mang toạ độ **partial**
(`<`, `>`): **1/34** — `NC_024936.1` (*Iris gatesii*), đúng bản ghi mà §8.15 bẫy 1 đo PGA ra **0
gen, exit 0** khi dùng làm reference. ⇒ Luật "partial ⇒ lỗi, dừng" của audit chặn đúng **một** bản
công bố thật; người dùng sẽ phải loại nó khỏi bộ tham chiếu (hoặc thay bằng bản không partial).

**Ba số này là đầu vào cho các luật trong `docs/D5-08_plastid_annotator.md`** — mục "The
reference-set audit" và "Arbitration". Chúng không được dùng để chọn số; chúng chọn luật.
