# Nhật ký đo lường — Giai đoạn 3A — annotator mito tự viết (kèm §8.25, việc của giai đoạn 2)

**Phần này là `PROJECT.md` §8** (§8.20–8.25), tách ra ngày 04/10/2026 để `PROJECT.md` đọc được trọn
(§8 khi đó chiếm 1.531/2.379 dòng, 118/216 KB). **Số mục giữ nguyên** như khi còn nằm trong
`PROJECT.md` — mọi tham chiếu kiểu `§8.27` trong code, test, `docs/` và board vẫn trỏ đúng chỗ —
và **không bao giờ đánh số lại**. Mục mới của giai đoạn này ghi tiếp vào **cuối file này**, rồi
thêm một dòng vào bảng mục lục ở `PROJECT.md` §8.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "Nhật ký & tài liệu".

---

### 8.20. Writer GenBank — mảnh còn thiếu để một annotator tự viết có đầu ra dùng được (P3A-01, 03/10/2026)

**Vì sao cần.** Tới giờ cả toolkit chỉ **đọc** GenBank (`core/genbank.py`); chưa có đường nào
**ghi** ra. Đường `.fsa` + `.tbl` của P3-05 đủ cho NCBI nộp bài, nhưng `panel pin --file` và mọi
bước sau đó đều đọc bằng **reader của mình** — nên một annotator tự viết (Phase 3A) mà không có
writer thì sản phẩm của nó không vào được pipeline.

**Đo hình dạng bản ghi vàng trước khi viết** (không đoán):
* **265/265** feature gen đều có `/gene`, và với cả CDS · tRNA · rRNA thì **tên khác `/product`**;
  feature `gene` chỉ có `/gene`.
* **0 partial** trên toàn bản ghi · 2 `repeat_region` · 3 `misc_feature` (model không mang).
* Qualifier model **không giữ**, nên writer không ghi lại: `locus_tag` (268) · `db_xref` (266) ·
  `protein_id` (86) · `translation` (86) · `note` (51) · `rpt_type` (2).

**Cái bẫy đo được — BioPython viết qualifier trần SAI.** `/pseudo` và `/trans_splicing` trong mọi
bản ghi thật **không có giá trị**; BioPython gặp `{"pseudo": [""]}` thì ghi `/pseudo=""` — một file
khác — còn gặp `None` mới ghi `/pseudo` trần. Đối chiếu chính file vàng: nó viết `/trans_splicing`
trần, `/codon_start=1` và `/transl_table=11` (**số, không trích dẫn**) — đúng những gì BioPython
sinh ra với `None` và với `str`. Từ đó thành quy ước trong model: **`None` nghĩa là qualifier trần**
(§9.B-25).

**Một hàm qualifier dùng chung.** Đường `.tbl` (P3-05) đã có `_qualifiers_of`; writer GenBank cần
đúng thứ đó ⇒ theo luật DRY-sau-lần-trùng-thứ-hai nó thành `model.gene_qualifiers()`, **cả hai
đường gọi chung một hàm** — một feature không thể ra hai kiểu ở hai định dạng. Nhân đó chốt thêm
một quy tắc: **`/gene` chỉ ghi khi nó nói điều `/product` chưa nói** — đo trên `OX459128.1` (không
có `/gene` nào, tên feature là cả câu `/product`), ghi cả hai sẽ gán cho bản ghi một *ký hiệu gen*
là một câu văn.

**Kết quả đo:** đọc → ghi → đọc lại trên genome vàng cho **265/265** feature trùng khớp
(kind · tên · các part **theo đúng thứ tự transcript** · mạch · cờ partial · product ·
`transl_table` · `codon_start` · `exception` · `pseudo` · `trans_splicing`), trình tự và `source`
nguyên vẹn. Thêm 8 test cho writer (join mạch trừ, partial, tên-trùng-product, qualifier trần,
`source`, và ca từ chối ghi feature vượt quá độ dài trình tự). **324 test unit vẫn xanh** — đường
`.tbl` của P3-05 (đối chiếu byte-với-byte với bảng thật của NCBI) **không đổi hành vi**.

**Writer này KHÔNG phải trình sửa GenBank bất tổn thất** — nó ghi *model*, không ghi lại file gốc:
các qualifier kể trên rơi mất, `misc_feature` không có trong model, và phần LOCUS (division `PLN`,
ngày hôm nay) được sinh lại. Điều đó đúng cho thứ nó phục vụ — đầu ra của annotator mình — và được
nói thẳng trong docstring thay vì để người đọc tự phát hiện.

---

### 8.21. Bộ tham chiếu mito — dựng từ 8 bản ghi RefSeq, và cái bẫy nằm trong chính bộ tham chiếu (P3A-02, 03/10/2026)

**Nguồn:** 8 mitogenome RefSeq ghim bằng `panel pin` (5,2 MB, sha256 trong `PINNED.json`), chọn
theo **đúng luật taxonomy của Reference Panel** — cùng họ trước: 3 họ Lamiaceae (*Salvia
miltiorrhiza* · *Ajuga reptans* · *Rotheca serrata*) cho *Coleus*, 3 họ Rubiaceae (*Psychotria
viridis* · *Scyphiphora hydrophyllacea*) cho *Oldenlandia*, cộng 3 loài mẫu (Arabidopsis ·
Nicotiana · Oryza) làm outgroup. Kết quả: `mito_reference_proteins.faa` (**287 protein, 100 KB**)
+ `mito_reference_rrna.fa` (**25 rRNA, 49 KB**).

**Đo được — tự dịch protein là KHÔNG THỂ, không phải chỉ kém chính xác.** Cả **530/530** CDS tham
chiếu đều có `/translation` và **không bản nào khai `/transl_table`**; và CDS mito thực vật mở đầu
bằng codon như `ACG` mà **không bảng mã nào nhận làm codon khởi đầu** — dịch ở chế độ CDS thì
BioPython **ném lỗi trước khi ra được protein nào**. Nên bộ tham chiếu lấy **chính `/translation`
của người chú giải**: nó cũng xử lý luôn trans-splicing và các đầu do RNA editing, những thứ module
này không có việc phải suy diễn lại.

**Đo được — chính bộ tham chiếu bị ô nhiễm, và đây là cái bẫy lớn nhất của task.** Bản ghi mito
thực vật chú giải **ORF vô danh như CDS**: *Nicotiana* **121 cái**, *Salvia* **97**, *Oryza* **22**
— tổng **236**. Thêm 7 feature tên là cả câu (`hypothetical protein`, `LAGLIDADG endonuclease`).
Thêm gen **nguồn gốc lục lạp** nằm trong mito: `rbcL` · `petG` · `petL` · `atpB` · `atpE` (cả 5
trong bản ghi *Salvia*), và gen chuyển ngang `citB` · `epd` · `ltrA` · `rnhA` · `tuf` (*Psychotria*).
Thêm rRNA **sai kích thước**: *Psychotria* gọi một feature **722 bp** là "5S" và một **1112 bp** là
"18S" (khoảng thật: 5S 107–122 · 18S 1692–1935 · 26S 2907–3507). Một build ngây thơ sẽ nuốt hết.
Nên builder **trả về và in ra** mọi thứ bị loại, không bỏ im lặng.

**Kết quả:** phủ **31/31** locus của vocabulary, và artifact được **giữ đúng** bằng một test dựng
lại từ panel đã ghim rồi so với file đã commit — tức FASTA là **bằng chứng**, không phải file sót
lại. Đường pin (`panel pin` → `PINNED.json`) giữ sha256 của cả 8 nguồn.

**Lỗ hổng đo được, để mở có tên (H-06):** tham chiếu còn **23 gen mà toolkit chưa có tên**, trong
đó **10 gen mito chuẩn** — `rpl2` · `rpl10` · `rpl16` · `rpl23` · `rps2` · `rps7` · `rps12` ·
`rps14` · `rps19` · `sdh3`. `mito_gene_names.yaml` chỉ có 31 locus vì nó được dựng từ **một** bản
ghi; đây đúng là phép kiểm chéo mà chính file đó ghi là còn thiếu. Chúng **được giữ lại trong bộ
tham chiếu** với dấu `-` ở cột vocabulary, để một lần chạy có thể **nói ra** là nó tìm thấy `rps12`
thay vì im lặng về một gen nó không có tên.

### 8.22. Bộ gọi gen mito — chuyển chú giải từ bộ tham chiếu, và sáu cái bẫy đo được (P3A-03, 03/10/2026)

**Vì sao.** P3-04 chọn PMGA rồi không chạy nổi (11 GB + MAKER so với trần 8 GB, web service từ
chối kết nối). Luật của CLAUDE.md cho ca đó là *thay — hoặc tự viết phần đó*, và các mảnh đã có
sẵn: BLAST+ 2.17.0, bộ tham chiếu 287 protein + 25 rRNA (P3A-02), vocabulary 31 locus, các phép
kiểm ORF, reader/writer của mình. P3A-03 là **bộ gọi gen**: `blastx` cho CDS, `blastn` cho rRNA,
rồi `core/mito_caller.py` biến các HSP thành `Gene`.

**Kết quả đo trên bản ghi đáp án `OX459128.1`** (chạy thật, BLAST+ 2.17.0, bộ tham chiếu đã commit):

| Số đo | Giá trị |
|---|---|
| Gen gọi được | **31 CDS + 3 rRNA**, `not_called` **rỗng** |
| Span trùng khít bản công bố | **18/31**; 13 ca còn lại lệch **≤ 75 bp** ở một đầu |
| Mạch | khớp **31/31** (kể cả 2 gen mixed-strand: nad1, nad5) |
| Hồ sơ ORF của các call | **21 ok · 3 no_start · 7 no_stop** (bản ghi gốc: 21 ok · 4 pseudo · 3 no_stop · 3 no_start) |
| Exon | 11 gen một exon khớp; nad1/nad5 ra **4 exon so với 5** của bản ghi — hai exon 19 và 7 codon của bản ghi **không alignment nào tìm được**, và bộ gọi **không bịa** chúng |
| Gen tìm thấy mà chưa có tên (H-06) | **15 hit** được báo cáo: rps19 · rps14 · rpl16 · rps2 · ccmFN2 · rbcL · atpB · orfX · orfB |
| Hit không dùng | **68**, mỗi cái một lý do |
| Thời gian | vài giây cho cả hai lượt tìm (bộ tham chiếu 287 + 25 chuỗi) |

**Sáu cái bẫy, tất cả đều đo được và đều đã thành test:**

1. **BLAST cắt `sseqid` tại khoảng trắng đầu tiên.** Header tham chiếu kết thúc bằng tên loài
   (`…|Nicotiana tabacum`) nên nó về tới nơi là `…|Nicotiana`. Khoá chủ thể **phải** là ba field
   đầu; parser **từ chối** một id ngắn hơn thay vì đoán.
2. **Một header có thể là hai chuỗi khác nhau.** FASTA đã commit có **8 header trùng**: 6 là bản
   sao y hệt (gen nằm trong repeat), nhưng *Arabidopsis* `atp6` là **385 và 349 aa**, *Psychotria*
   `cox2` là **257 và 390 aa**. Nên một HSP thuộc về **độ dài nhỏ nhất chứa được khoảng subject
   của nó**, và HSP của hai chuỗi cùng khoá **không bao giờ được nối chung** — đó chính là cách
   một gen khảm được dựng lên.
3. **`-seg no`.** SEG che các đoạn ít phức tạp của *query* — hợp lý khi đi tìm protein lạ, sai khi
   chuyển chú giải cho một mitogenome thực vật vốn đầy protein màng kỵ nước. Đo được: `atp9` mất
   **137 → 67 bits**, 21 locus khác cũng mất điểm. Ngưỡng identity mới là thứ giữ rác ở ngoài.
4. **`makeblastdb -parse_seqids` là flag phá hoại.** Nó **lỗi thẳng**:
   `Error: Duplicate seq_ids are found: LCL|NC_006581.1|SDH3|-|NICOTIANA` (bộ tham chiếu có bản
   sao repeat một cách hợp pháp), và với id nó parse được thì nó **viết hoa** — phá mọi ký hiệu
   gen viết thường trong khoá.
5. **Mạch của rRNA đến từ phía subject, không phải query.** Tham chiếu lưu theo chiều 5'→3' của
   gen, nên `sstart > send` nghĩa là gen nằm trên **mạch âm** của query — ngược với trực giác và
   ngược với quy ước của blastx (nơi frame quyết định). Đo lại đúng cả ba rRNA của bản ghi đáp án.
6. **`blastx` chỉ nhận `-version`, không nhận `--version`.** `run_tool` chỉ thử `--version`/`-v`/`-V`
   ⇒ công cụ vừa trở thành mấu chốt của cả một module bị ghi là **"unknown"** trong khi mọi lần dò
   đều "thành công" (thất bại). Đã thêm `-version` vào danh sách dò, kèm test.

**Một phát hiện về chính bản ghi đáp án:** ba CDS `no_stop` của nó (ccmFc · ccmC · nad6) **kết thúc
đúng ở codon cuối của protein và bỏ codon stop nằm ngay sau** — đo được: bộ gọi lấy codon đó khi
chuỗi xác nhận, nên nad6/ccmFc của mình dài hơn bản ghi **3 bp**. Đây là khác biệt *có thật* giữa
hai cách chú giải, không phải lỗi của bên nào.

**Đã khai báo công cụ:** `environment.yml` ghim `blast=2.17.0` (trước đó là `blast>=2.8.1` chỉ vì
PGA cần), `docs/tools/blast.md` ghi nguồn, phiên bản, đường cài, ba quyết định đo được và lý do
không vendor (52 MB nhị phân).

---

### 8.23. Bộ gọi tRNA — ARAGORN 1.2.41, và 17 không phải là một tập hợp (P3A-04, 03/10/2026)

**Vì sao.** CDS đến từ BLAST+ (P3A-03); tRNA thì không thể chuyển từ bộ tham chiếu protein, nên
cần một bộ gọi riêng. Chọn **ARAGORN 1.2.41** (D3A-01): 147 kB trên bioconda, GPLv3, và chính là
công cụ mà PMGA chạy bên trong. Một bộ gọi tự viết rồi chỉnh cho ra đúng 17 sẽ là overfitting
đúng nghĩa.

**Kết quả đo trên `OX459128.1`:** ARAGORN gọi **17 tRNA**, bản ghi công bố có **17** — **nhưng
không phải cùng một tập**:

| Khác biệt | Chi tiết |
|---|---|
| Thiếu | `trnA-UGC 13244-13348` — bản ghi chú giải, ARAGORN không gọi (feature dài **105 bp**, hơn một tRNA thường) |
| Thừa | `trnL-AAG 36355-36427` — ARAGORN gọi, bản ghi không có gì ở đó |
| Lệch anticodon | tRNA-His ở 11741-11819: ARAGORN đọc **`trnH-AUG`**, bản ghi ghi **`trnH-GUG`** |
| Lệch toạ độ | **9** feature lệch 1–6 bp; lớn nhất là tRNA **vắt qua origin** — ARAGORN in `[258274,73]`, bản ghi viết `1-79` |
| Trùng khít | 15 feature còn lại |

Con số **17 = 17** vì vậy không phải bằng chứng khớp: nó che một ca thiếu, một ca thừa và một ca
lệch anticodon. Đây đúng là lý do bài kiểm tra so **từng feature**, không đếm tên: 17 tRNA của bản
ghi chỉ mang **14 tên khác nhau** (`trnH-GUG` · `trnK-UUU` · `trnM-CAU` mỗi tên hai lần) — gộp theo
tên sẽ ra 14 và sai.

**Ba cái bẫy đo được, đã thành test:**

1. **`-V` của ARAGORN không phải cờ phiên bản.** Nó thoát với mã **0** và in
   "No sequence file specified, type aragorn -h for help" — `run_tool` ghi nguyên câu đó làm
   **phiên bản** của công cụ gọi tRNA. Một chuỗi trông hợp lý trong hồ sơ provenance còn tệ hơn
   "unknown": không ai phát hiện được. Bộ dò giờ chỉ nhận dòng **có chữ số**, và adapter tự xác
   minh bản phát hành từ banner của `-h`.
2. **ARAGORN in anticodon dưới dạng DNA** (`(cat)`) còn tên công bố dùng RNA (`trnM-CAU`), và in
   axit amin bằng chữ (`tRNA-Met`). Adapter đổi Met→M, `cat`→`CAU` rồi đưa qua
   `core/gene_names.normalize_name` — **cùng một bộ chuẩn hoá** mà một bản ghi tải về phải đi qua.
3. **Gene vắt qua origin**: ARAGORN in `[258274,73]` (đầu > cuối). Đó không phải toạ độ GenBank;
   adapter tách thành hai part theo đúng thứ tự phiên mã, tức một `join()`.

**Đã khai báo:** `environment.yml` ghim `aragorn=1.2.41`, `docs/tools/aragorn.md` ghi nguồn
(bioconda — upstream `ansikte.se` trả HTTP 455), licence GPLv3, đường cài và ba quyết định trên.

---

### 8.24. Chạy trên mẫu chưa ai chú giải — `AP045469.1`, và tám locus không gọi được (P3A-07, 03/10/2026)

**Input:** `AP045469.1` (*Coleus hadiensis*, **313.585 bp**, **0 feature** — đo được), đưa vào qua
`plastid mito annotate` dưới dạng FASTA trần. Đây là hình dạng mọi lần chạy mito sau này sẽ đi:
chuỗi vào → GenBank ra → `panel pin --file` ghim → đọc lại offline.

**Kết quả:** **45 gen = 23 CDS + 19 tRNA + 3 rRNA**, ghim với sha256 `3e7c7114ae5d…`
(402.318 B, `parse_error: null`). `mito_genes` đọc thẳng từ bản ghim, **không cần công cụ, không
cần bộ tham chiếu**.

**Tám locus vocabulary KHÔNG gọi được — và đây là phát hiện, không phải lỗi cần giấu.** Bảng
`not_called` ghi rõ từng locus; đo thêm bằng cách nới E-value lên **10** (so với ngưỡng 1e-5 của
lần chạy) để trả lời câu "ngưỡng có quá chặt không?":

| Locus | Ở E-value ≤ 10 |
|---|---|
| `mttB` · `rps1` · `rps13` · `rps4` | **không có alignment nào cả** |
| `cox2` | một mảnh 40 aa, 67,5 % — quá ngắn để phủ |
| `ccmFn` · `matR` · `rps3` | mảnh 21–52 aa ở E-value 0,05–2,8 — không phải bằng chứng tương đồng |

Nên **không phải ngưỡng** đang giữ chúng ở ngoài. Câu hỏi "các gen đó vắng trong genome hay vắng
trong bản lắp ráp này" **không phải câu một bộ chuyển chú giải trả lời được** — nó cần người chú
giải, hoặc một phép kiểm tra từ read. Ghi lại đúng như vậy thay vì hạ ngưỡng cho đủ 31.

**Hồ sơ ORF của 23 call:** **10 ok · 9 no_start · 3 no_stop · 1 internal_stop**. Chín `no_start`
phần lớn là các call bắt đầu giữa protein (alignment không tới đầu N của tham chiếu) — module đánh
dấu `partial_start` chứ không dịch chuyển toạ độ. **Đây là prototype, không phải bản nộp** — dòng
scope trong mọi bảng và ghi chú của lần chạy đều nói vậy.

**Nhắc lại giới hạn (D3A-05):** không tự tìm gen mới, không mô hình intron, không dự đoán RNA
editing. Bốn gen không tìm thấy **không** được suy đoán là "không có".

---

### 8.25. Chạy thật hai module so sánh — và một lỗi im lặng sống sót vì chưa từng chạy dữ liệu thật (P2-07/P2-08, 03/10/2026)

**Bối cảnh.** P2-07 (π) và P2-08 (Ts/Tv) được viết và kiểm bằng unit test từ Phase 2, nhưng chưa
bao giờ chạy trên dữ liệu thật vì thiếu MAFFT 7.526. Hôm nay MAFFT đã cài và pin đúng §7.8, nên
đây là lần chạy thật đầu tiên — và nó tìm ra một lỗi.

**Lỗi đo được: MAFFT 7.526 in căn thẳng ra CHỮ THƯỜNG.** Bộ đếm phân loại cột theo `_DNA =
{"A","C","G","T"}` (chữ hoa), nên **mọi cột của một alignment chữ thường bị tính là "skipped"**:
lần chạy Ts/Tv thật đầu tiên báo **0 cột so sánh được và Ts/Tv để trống** — không phải một con số
sai, mà là *không có con số nào*, trong khi mọi bảng vẫn ra bình thường. Unit test không thấy được
vì chúng tự dựng `Alignment` bằng chữ hoa; adapter MAFFT **chưa từng có test nào**. Đã sửa hai đầu:
adapter **upper-case** khi đọc output (kèm comment ghi rõ số đo), và `Alignment` **từ chối** một
hàng có chữ thường — fail loud thay vì đếm ra số 0.

**Ts/Tv — tái lập đúng số vàng §6, bằng chính đường của toolkit:** thêm scope `cds` (các gen hai
genome dùng chung, mỗi gen một bản — IR tính một lần, align từng gen):

| | Lần chạy này | §6 |
|---|---|---|
| Gen dùng chung | **76** | 76 |
| Cột căn | **66.066** | 66.066 |
| Ts / Tv | **317 / 118** | 317 / 118 |
| Ts/Tv | **2,686** | 2,686 |

Scope mặc định (cả bộ gen) cho **1,8416** trên 152.338 cột — một con số khác, và mọi bảng đều ghi
rõ scope nào. Đây là lý do module có hẳn hai scope: §6 chỉ đúng trong scope hẹp.

**π — lần chạy thật đầu tiên:** π toàn alignment **0,007181**; theo vùng **LSC 0,009185 · IRa
0,001682 · IRb 0,001681 · SSC 0,014066**; kèm π theo gen và cửa sổ trượt 600/200. Tín hiệu đúng
sinh học: **hai bản IR gần như trùng nhau và thấp hơn hẳn** hai vùng đơn bản — một phép kiểm hình
dạng mà một lỗi "đếm ra 0" sẽ xoá sạch. **Còn thiếu:** π theo **intron/IGS** (H-07) — cần tách phép
suy intron ra tầng dùng chung trước.

