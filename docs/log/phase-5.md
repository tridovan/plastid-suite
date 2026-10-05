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

---

### 8.35. Kiểm thử truy vấn số trên report **đã sinh** — và lớp va chạm buộc phải thu hẹp luật (P5-05, 05/10/2026)

**Bối cảnh.** §13.5 yêu cầu: mỗi số trong report phải trỏ về **đúng một ô** của bảng gốc. Hôm nay
việc đó đúng **theo cấu trúc** — template chỉ trỏ `{{bảng.cột}}`, và `find_bare_digits` đã chặn số
gõ tay trong template — nhưng chưa có gì kiểm **report đã sinh**. `report/audit.py` (P5-05) làm việc
đó: bỏ ra ba khoảng không do bảng nào trả lời (khối trích dẫn · nhãn dữ liệu của run · banner nháp),
rồi đòi mỗi số còn lại phải nằm trong **giá trị đã render** của một ô.

**1. Luật nghiêm (không miễn trừ gì thêm) hỏng trên chính bộ giao phẩm Pha 4 — 10 phát hiện, 6 lớp.**
Chạy trên **12 report** của hai run `.runs/phase4/standin-oldenlandia` và `standin-coleus`:

| Số | Bảng | Ô thứ hai cùng chứa số đó |
|---|---|---|
| `10` | `synteny_blocks` | id khối `NC_057983.1#10` (cột khoá `block`) |
| `1` | `barcodes_summary` | chuỗi luật `mean + 1 sd of this run's windows` |
| `6` · `5` | `popgen_*` | chuỗi cỡ nhóm `Coleus: 6, Oldenlandia: 5` |
| `1` | `popgen_*` | cùng chuỗi trên + một câu ghi chú |
| `80` | `kaks_*` | ghi chú `row 'reference': a gap starts after 80 base(s) — not a codon boundary` |

**Cả 10 phát hiện đều thuộc đúng một lớp: một lượng (ô số) va vào một số in bên trong một ô
*chữ* được report trích nguyên văn.** Không phát hiện nào là "số không có bảng nào mang" — tức
kiểm tra *orphan* chưa hề báo động giả.

**2. Luật sau khi đo: ô chữ có chữ cái là văn xuôi, không phải lượng.** Một giá trị của cột `str`
mà chứa chữ cái (câu luật, ghi chú, id khối) bị **bỏ khỏi thân report đang quét** và **không được
tính là ô mang số**; một ô `str` render đúng một con số trần (ví dụ `"80"`) thì vẫn là ô mang số.
Sau luật này: **12/12 report sạch, 0 phát hiện**, và cả bốn lý do miễn trừ đều thực sự có mặt
(12 report × 4 lý do).

**3. Bộ đếm của lần đo đó.** 12 report · **44 số phân biệt** được kiểm. Con số nhỏ vì report của
module là văn xuôi ngắn (bảng mới là nơi chứa số) — nên phép kiểm này canh **report**, không thay
`find_bare_digits`.

**4. Cắm thử một số mồ côi vào report thật thì bắt được.** Chèn `About 99.5 % of each genome
aligned` vào `03_synteny/report.md` của run `standin-oldenlandia` (đọc lại `result.json` bằng
`read_result`): kết quả **`orphan number '99.5' — no table cell carries it`**.

**5. `read_result` đọc lại được `result.json` thành `Result`.** Đo: `export_json` → `read_result`
→ `export_json` cho ra **chuỗi y hệt** trên fixture genome vàng, và ô trống (`null`) trở về đúng
`None`, không thành 0. Đúng **một** thứ không round-trip, có chủ ý: `Result.artifacts` — đầu ra
dạng file nằm **cạnh** `result.json` và manifest liệt kê, chép vào JSON sẽ nhân đôi mọi file.
`Table.emphasis` trước đây thiếu trong JSON và đã được bổ sung trong chính task này.

**Hệ quả.** Bản ghi của một run nay **đọc lại được** — hình, gói Supplementary và QC report
(P5-04) không phải chạy lại phân tích nào. Và luật "đúng một ô" có hiệu lực thật: nó bắt số lạ,
còn ca "hai dòng cùng render một giá trị" **không** bị coi là mơ hồ — hai ô không thể mâu thuẫn,
nên số vẫn có đúng một câu trả lời (đo trên chính genome vàng: hai bản IR đều 26.214 bp).

---

### 8.36. Hình vẽ: công thức xuất vector tất định, và dải vùng không được lấy từ mép cửa sổ (P5-01, 05/10/2026)

**Bối cảnh.** P5-01 dựng hệ thống hình: `[figures]` (matplotlib 3.11.2, extra tuỳ chọn, import
trong hàm) · `result/figure.py` (theme + `Figure` + driver) · `result/figures/` (một file mỗi hình,
registry) · `run/figures.py` (render ở **mức run**, đọc `result.json` đã ghi) · CLI `plastid run
--figures` và `plastid figures <run-dir>`.

**1. Cái gì trong file vector làm hai lần render khác nhau — đo bằng cách render cách nhau 1,2 giây.**
matplotlib đóng dấu thời gian thật vào cả ba định dạng, và **SVG còn thêm một salt ngẫu nhiên theo
tiến trình** cho id của các phần tử:

| Định dạng | Dấu bị đóng | Cách ghim | Còn khác không? |
|---|---|---|---|
| SVG | `<dc:date>2026-10-05T21:11:09.509605</dc:date>` **và** id kiểu `#m3b5f8bc5b3` | `rcParams["svg.hashsalt"]="plastid"` **+** `savefig(metadata={"Date": None})` | không |
| PDF | `/CreationDate (D:20261005211110+07'00')` | `savefig(metadata={"CreationDate": None})` | không |
| EPS | `%%CreationDate: Mon Oct  5 21:11:12 2026` | backend PostScript **không nhận metadata** ⇒ thay dòng đó bằng mốc cố định | không |

Đo riêng: SVG với `metadata={"Date": None}` mà **không** có `hashsalt` vẫn khác nhau (id phần tử
đổi theo tiến trình) — nên phải cả hai. `metadata` cũng **không** dùng chung một dict cho ba backend
được. Kiểm chứng end-to-end: chạy `plastid run --figures` hai lần trên cặp Iris đã pin ⇒
**svg · pdf · eps đều byte-identical**, và lần thứ hai là **cache hit** — hình vẫn ra, vì nó render
từ `result.json` đã ghi chứ không từ đối tượng `Result` trong bộ nhớ.

**2. Dải vùng LSC/SSC/IR không được lấy từ mép cửa sổ — cửa sổ chồng nhau 400 bp.** Đo trên chính
hình π của cặp Iris (600 bp/cửa sổ, bước 200 bp): lấy mép cửa sổ làm mép dải thì
`LSC 1–84000` và `IRa 83601–110200` **chồng nhau 399 bp**, SSC và IRb cũng vậy — tức là tô vùng này
đè lên vùng kia. Sửa: ranh giới đặt ở **trung điểm giữa cửa sổ cuối của vùng cũ và cửa sổ đầu của
vùng mới**, được `LSC 1–83800 · IRa 83800–110000 · SSC 110000–128600 · IRb 128600–154827`.
So với mốc thật của khung chuẩn hoá (`83127 · 109341 · 127510`): lệch **670 · 660 · 1090 bp** —
đúng độ phân giải của cách cửa sổ hoá, và hình **không** giả vờ chính xác hơn.

**3. Số in trên hình là chính ô của bảng, không có lần làm tròn thứ hai.** Dòng tiêu đề hình π in
`whole-alignment π 0.007181` — lấy từ `diversity_summary.format("pi")`, đúng chuỗi mà CSV và report
in. Test canh cả chữ lẫn số: chữ phải nằm trong SVG dưới dạng `<text>` (đặt `svg.fonttype: none` —
SVG là bản **để sửa**), còn PDF nhúng TrueType (`pdf.fonttype: 42` — bản **để nộp**).

**4. Hình phải được NHÌN, không chỉ assert.** Bài học §9.C-34 áp lại đúng lần này: hình render ra
hợp lệ, test xanh, nhưng **mắt mới thấy** dải vùng chồng nhau ở mục 2 — và thấy π tụt hẳn trong IR,
vọt lên ở SSC, khớp các số vàng (`LSC 0.009185 · IRa 0.001682 · SSC 0.014066`).

**Quyết định kiến trúc kèm theo (ghi để lần sau khỏi hỏi lại).** `TASKS.md` viết `core/figure.py`,
nhưng `tests/unit/test_architecture.py` chỉ cho `core/` import `plastid.core` — mà một hình là **hàm
của `Result`**, và `Result` nằm ở tầng `result`. Hơn nữa hai hình đọc bảng của **hai module khác
nhau**, việc mà một callback theo-từng-module không làm được. Nên: hình ở tầng `result`
(`docx_export.py` đã là tiền lệ cho một exporter phụ thuộc tuỳ chọn nằm ở đó), và render ở **mức
run** sau khi mọi module chạy xong — một hàm dùng chung cho `plastid run --figures` lẫn
`plastid figures <run-dir>`, nên cache không bao giờ chen vào được.

---

### 8.37. Bảy hình đầu tiên, và quy tắc mở rộng theo N đo bằng ba lần vẽ (P5-02, 05/10/2026)

**Bối cảnh.** P5-02 dựng sáu hình còn lại — `circular_map` · `ir_junction` · `rscu_heatmap` ·
`ssr_distribution` · `identity_synteny` · `barcode_assessment` — và **đo** quy tắc mở rộng D5-04
thay vì đoán. Cả bảy hình nay render từ **một thư mục run** (`plastid figures <run-dir>`), mỗi
hình khai báo bảng nguồn theo **(module, bảng)** và hình nào thiếu bảng thì **bị bỏ qua kèm lý do
có tên**, không vẽ rỗng.

**1. Quy tắc N (D5-04): đo ở N = 2 · 8 · 34 trên bộ dữ liệu anh Điềm, cùng một hình
(`identity_synteny`).**

| N | thời gian run `synteny` + vẽ | đỉnh RSS | hình đọc được? |
|---|---|---|---|
| 2 | 2,35 s | 124 MB | có |
| 8 | 6,16 s | 131 MB | có — 7 series, legend 4 cột |
| 34 | 22,07 s | 165 MB | **không** — legend phủ kín panel, 33 đường nhoè thành một khối xám |

Ngưỡng lấy từ chính ba lần đo: **trần 8 series**. Trên trần, hình đổi sang **dạng tổng hợp**: dải
phủ giữa đường thấp nhất và cao nhất của mọi bộ gen tại từng cửa sổ — **hai mép dải đều là giá trị
ô**, không sinh ra số mới, và bảng vẫn giữ nguyên từng series. Canvas và panel khối cũng nở theo N
(33 nhãn bộ gen trên chiều cao cố định là một vệt, không phải nhãn). Sau khi sửa, N = 34 **đọc
được**: dải phủ cho thấy rõ vùng tụt sâu quanh 110–120 kb — đúng vùng `ycf1`/SSC.

**2. Một lỗi vẽ chỉ mắt thấy, lần thứ hai trong cùng giai đoạn.** `barcode_assessment` vẽ dải ứng
viên theo toàn trục: hình hợp lệ, test xanh, nhưng **nhìn thì sai**. Nguyên nhân là một giả định
sai *trong bản khảo sát viết ra cho task này* — "locus vắt qua gốc thì `reference_start >
reference_end`". Đo lại trên hai run giao phẩm Pha 4: **mỗi run đúng một** ứng viên có
`start > end`, đó là **`rps16`** — gen hai phần (intron) chạy ngược — và `note` của nó **rỗng**.
Vắt qua gốc là ca *khác*, và thứ duy nhất nói ra điều đó là **`note` do module tự viết**
(`spans the origin: 153273–153724 and 1–71`). Hình nay: `start ≤ end` ⇒ một dải; `start > end` mà
note nói vắt gốc ⇒ hai dải; còn lại ⇒ **footprint** `[end, start]`.

**3. Đo lại bằng mắt, không chỉ bằng assert.** Nhìn bốn hình đầu tiên cũng bắt được: dải GC-skew
trên bản đồ vòng bị **khe hở 50%** (bề rộng bar tính bằng nửa cung), và nhãn gen **đè lên nhau** —
sửa bằng ngưỡng 4° cộng khoảng cách tối thiểu 8° giữa hai nhãn liền kề trên cùng một mạch.

**4. `read_run_results` + `Frame`.** Bản đồ vòng và sơ đồ IR junction đọc **hai module cùng lúc**
(`genome_stats` ở khung GenBank, `structure` ở khung chuẩn hoá) nên cần đổi khung: `Frame` nhận
đúng **hai ô** — chiều dài và `rotation_offset`. Genome vàng có offset 0 nên hai khung trùng nhau;
đó chính là lý do một hình bỏ qua phép đổi khung vẫn "đúng" trên ca kiểm chuẩn và **sai** trên một
bản ghi bắt đầu giữa IR. Test đặt offset = 250 để bắt buộc đi qua đường đó.
