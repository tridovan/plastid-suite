# Nhật ký đo lường — Giai đoạn 4a — bộ dữ liệu anh Điềm làm bộ test

**Phần này là `PROJECT.md` §8** (§8.31–8.33), tách ra ngày 04/10/2026 để `PROJECT.md` đọc được trọn
(§8 khi đó chiếm 1.531/2.379 dòng, 118/216 KB). **Số mục giữ nguyên** như khi còn nằm trong
`PROJECT.md` — mọi tham chiếu kiểu `§8.27` trong code, test, `docs/` và board vẫn trỏ đúng chỗ —
và **không bao giờ đánh số lại**. Mục mới của giai đoạn này ghi tiếp vào **cuối file này**, rồi
thêm một dòng vào bảng mục lục ở `PROJECT.md` §8.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.31. Bộ dữ liệu anh Điềm gửi — bộ test đầu tiên không do dự án tự chọn (P4a-01, 04/10/2026)

**Bối cảnh.** Giai đoạn 1–4 được dev và test trên dữ liệu **tự chọn**. Anh Điềm gửi hai bộ — một
bộ tải trực tiếp từ database (chưa xử lý), một bộ đã qua xử lý sơ bộ — mỗi bộ ở hai dạng input
(`.gb` đã chú thích, `.fasta` chỉ có chuỗi) vì mỗi công cụ cần một loại input khác nhau. Đo ra:
**34 plastome họ Iridaceae × 2 bản chú giải** = 68 `.gb` + 68 `.fasta` · 28 MB · 136 file.

**Chỗ đặt và tốc độ.** Đặt trong cây làm việc tại `.data/thanhdiem_data/`, thêm `.data/` vào
`.gitignore` (28 MB không commit); **136/136 sha256 khớp bản gốc** trên `D:\`. Lý do đặt trong WSL
là tốc độ, và con số đo được nói rõ nó đến từ đâu:

| Việc | `/mnt/d` | trong WSL |
|---|---|---|
| đọc một lượt 28 MB | 0,55 s | 0,38 s (**1,4×**) |
| 300 file nhỏ (ghi / đọc) | 700 / 614 ms | 7 / 3 ms (**~100×**) |

Đọc thô nhanh hơn **không đáng kể** — thứ đáng kể là **thao tác file nhỏ**, và một run chính là
hàng nghìn file nhỏ (run Phase 4 có 616 lượt MAFFT ⇒ 1.232 file log). Đây là lý do thật, không
phải "copy cho nhanh".

**Hai bộ chú giải là hai cách đọc CÙNG MỘT chuỗi.** Đo: **34/34 chuỗi giống nhau từng byte** giữa
hai bộ. Nghĩa là mọi khác biệt giữa hai bộ **chỉ có thể là khác biệt chú giải** — một phép thử tự
nhiên mà dự án chưa từng có: cùng đầu vào, hai đường ống, so thẳng kết luận.

**Bản ghi vàng nằm trong bộ dữ liệu.** `01_…/Iris_domestica_NC_050833.1.gb` có sha256
`a946f23e6e5f6a886eea93c52d72497209531aaf855ec08853dbef7997638493` — **trùng khít** file fixture
đang vendored ở `tests/data/NC_050833.1.gb`. Bộ test chứa luôn chính cái hợp đồng nó kiểm.

**Bộ 02 = PGA chú giải lại genome vàng:** 153.724 bp · GC 37,85 % · IR 26.214 bp @ 100,00 % — nhưng
**85 CDS, không phải 86**, thiếu đúng **bản sao thứ hai của `rps12`**. Trên bộ 01, *E. bulbosa*
`PX516262.1` (bản ghi của chính nhóm) đọc ra **133 = 87 + 38 + 8**; bộ 02 ra **131 = 85 + 38 + 8**
— cùng một kiểu mất bản sao.

**Gate IR bắt được 3/34 bản ghi của PGA:**

| Bản ghi | IR theo annotation (PGA) | de-novo trên **cùng chuỗi đó** |
|---|---|---|
| *Chasmanthe floribunda* | 26.440 bp @ **29,24 %** | 26.440 bp @ **100,00 %** |
| *Iris hermona* | 26.180 bp @ **26,36 %** | 26.132 bp @ **100,00 %** |
| *Trimezia steyermarkii* | 26.169 bp @ **27,78 %** | 26.169 bp @ **100,00 %** |

*Chasmanthe* là ca rõ nhất: de-novo trả về **đúng span PGA đã ghi** nhưng ở 100 % ⇒ PGA chọn đúng
đoạn, **ghép cặp sai**. Bộ 01 của cả ba bản ghi đều qua gate ở 100 % ⇒ đây là lỗi của **bản chú
giải**, không phải của chuỗi. Tức **3/34 = 8,8 %** số bản ghi sẽ hỏng im lặng nếu ai đó tin IR
trong annotation — đúng lý do gate tồn tại (§9.A).

**34/34 bản ghi viết tên gen khác nhau giữa hai chú giải.** Phân loại:

| Khác biệt | Số bản ghi | Loại |
|---|---|---|
| `clpP1` → `clpP` | 23/34 | từ vựng |
| `pbf1` → `psbN` | 21/34 | từ vựng |
| `ndhL` → `ndhI` | 8/34 | từ vựng (cùng locus — xem dưới) |
| `atpL`/`psaL`/`bspB` → `atpI`/`psaI`/`psbB` (bộ ba của §8.27) | 7/34 | từ vựng |
| `lhbA` → `psbZ` | 1/34 | từ vựng |
| tRNA viết tên trần (không anticodon) | 8/34 | từ vựng |
| rRNA kiểu `16S ribosomal RNA` vs `rrn16` | 1/34 | từ vựng |
| **`ycf1` có CDS ở bản công bố, KHÔNG có ở PGA** | **24/34** | **khác kết luận** |
| **`rps12` mất một bản sao** | 12/34 | **khác kết luận** |

Hai dòng cuối **không phải khác cách viết mà khác kết luận**: ở **71 %** bộ dữ liệu, hai chú giải
bất đồng về việc một gen **có tồn tại hay không**. Đây là bài học của §8.27 ở quy mô 34 bản ghi —
và là lý do mọi so sánh theo tên phải nói ra nó đang khoá theo tên ai.

**Bản chú giải tự khai công cụ — và đây là bằng chứng, không phải suy đoán.** Cả **34/34** bản ghi
bộ 02 mang `/note="Annotation Method :: PGA-Plastid Genome Annotator"` trên feature `source`;
bộ 01: **0** bản ghi. Việc gọi bộ 02 là "bản PGA" vì thế đọc thẳng từ file, không phải đoán theo
kiểu đặt tên.

**`ndhI` (PGA) và `ndhL` (công bố): alias đã có từ 03/10, hôm nay đo lại (D4a-05).** Trên chính
genome vàng, hai bản chú giải đặt locus ở **cùng toạ độ `118156-118698(-)` với chuỗi 543 bp giống
hệt** (`ATGTTCCCT…TAA`) — đúng như §8.17 đã đo khi Chloe chú giải genome này, và mục `ndhI` ←
`ndhL` đã nằm trong `data/gene_aliases.yaml` từ đó. Bộ dữ liệu mới chỉ **xác nhận lại ở quy mô**:
locus được viết khác nhau ở **8/34** cặp, và kiểm bằng code thì **cả hai cách viết resolve về cùng
một khoá chuẩn `ndhI`** — nên một so sánh trải qua cả hai bộ vẫn ghép đúng locus này. Cặp genome
vàng **không đổi số** (cả hai bản NCBI đều viết `ndhL`). **Không có mục alias mới nào được thêm**:
một mục đã đúng thì việc phải làm là một test canh nó, không phải viết lại nó.

**Đã ghim — P4a-01, Done-when đạt:**

| Panel | Nội dung | Kiểm lại |
|---|---|---|
| `.panel/iridaceae-34/` | 34 bản công bố + `fasta.json` | 34 bản ghi · integrity OK (sha256 · md5 · bytes · số đo đọc lại) · **0** rớt gate · 2.924 CDS + 1.292 tRNA + 272 rRNA |
| `.panel/iridaceae-34-pga/` | 34 bản PGA | 34 bản ghi · integrity OK · **3** rớt gate, từng bản ghi được nêu tên · 2.887 CDS |
| `fasta.json` | 68 mục (2 × 34) | mỗi sha256 đọc lại đúng · **68/68 chuỗi FASTA trùng khít bản ghi đã pin** |

**Một khoảng trống đo được, thuộc P4a-04:** `plastid run -i` **chỉ đọc GenBank**
(`cli.py` → `read_genbank`); `read_fasta` mới chỉ có ở đường `mito annotate`. Muốn chạy Chloe trên
34 chuỗi trần thì `plastid run` cần một cửa FASTA — hoặc chạy trên bản ghi đã bóc chú giải (cách
P3A-07 đã làm), vì hai đường mang **cùng chuỗi** (đã kiểm 68/68).

**Bộ fixture đã vendor (P4a-02, Done-when đạt).** Bộ dữ liệu 28 MB không commit; cái được vendor là
**tập nhỏ nhất mà mỗi file mang sẵn một đáp án viết ra được** — 5 file · 1,3 MB ở
`tests/data/corpus/`, ba cặp:

| Cặp | Bản công bố | Bản PGA | Điều nó kiểm |
|---|---|---|---|
| genome vàng | `NC_050833.1.gb` (**đã có sẵn**, sha256 trùng khít file của bộ) | 85 CDS | hợp đồng §6 không phải thuộc tính của chuỗi, mà của chú giải |
| *Iris hermona* | 26.132 bp @ 100,00 % — **qua gate** | 26.180 bp @ **26,36 %** — **rớt gate** | một chuỗi, hai phán quyết; fixture hồi quy cho §9.A |
| *E. bulbosa* (bản ghi của nhóm) | **133 = 87 + 38 + 8** | **131 = 85 + 38 + 8** | cùng kiểu mất bản sao, trên chính bản ghi nhóm quan tâm |

12 test tier-A trong `tests/golden/test_corpus.py` khoá từng con số trên, **và khoá cả điều làm
cho cặp fixture có nghĩa: hai file của một cặp mang chuỗi giống nhau** (nếu điều đó thôi đúng thì
mọi kết luận "hai chú giải khác nhau" đang đo một thứ khác). Kèm một test mới canh **mọi file dưới
`tests/data` đều có trong `SHA256SUMS` và khớp** — một file đã vendor mà không ai kiểm được thì
không phải là file đã ghim. Cặp hermona còn được chạy qua CLI: **exit 0** và **exit 2**, không có
output module nào sau khi gate chặn. `make check` xanh: **545 fast + 42 golden**.

**Một đính chính trong chính phiên này:** bảng Phase 4a ban đầu ghi D4a-05 là "thêm alias
`ndhI`↔`ndhL`". Kiểm lại thì **mục đó đã nằm trong `data/gene_aliases.yaml` từ 03/10** (§8.17, đo
khi Chloe chú giải genome vàng — cùng toạ độ, cùng 543 bp, và phức hệ ndh chỉ có tiểu đơn vị A–K).
Hôm nay chỉ là đo lại ở quy mô 34 bản ghi, **không có mục alias mới nào được thêm** — việc phải
làm là một test canh nó, và test đó đã có.

### 8.32. Quy mô: 34 bản ghi tốn bao nhiêu (P4a-03, 04/10/2026)

**Máy đo — phải nói ra, vì một con số quy mô không có máy thì không phải là con số:** Intel
i7-1165G7 @ 2,80 GHz · 6 vCPU trong WSL2 · **8.135.488 kB RAM** (~7,8 GB), ~4 GB trống. Đúng lớp
máy phổ thông mà dự án nhắm tới, và đúng trần bộ nhớ đã ghi ở §9.28.

**Cách đo:** mỗi lần một module, `plastid run -m <module> -i <genome…>`, tham chiếu đứng đầu, và
**`--force` mọi lượt** — cache của run layer khoá theo `module_version` chứ không theo mã nguồn
(H-08), nên một lượt ăn cache là một lượt đo … việc tra cache. Số là `Elapsed` và `Maximum
resident set size` của `/usr/bin/time -v`; đã kiểm GNU time **có** tính RSS đỉnh của tiến trình
con (đo bằng một tiến trình con ngốn 300 MB: báo 317 MB).

**Chín module một-genome, cả 34 bản ghi:**

| Module | 34 lượt | RAM đỉnh |
|---|---|---|
| `genome_stats` | 20,5 s | 87 MB |
| `structure` | 21,6 s | 87 MB |
| `introns` | 21,9 s | 87 MB |
| `ssr` (MISA kèm theo) | 25,5 s | 87 MB |
| `codon_usage` | 22,8 s | 87 MB |
| `amino_acid` | 23,9 s | 87 MB |
| `orf` | 12,9 s | 51 MB |
| `gene_table` | 15,3 s | 51 MB |
| `submission` | 12,7 s | 52 MB |

**Cả một panel chỉ là chuyện vài phút, không cần chia lô.** 87 MB không phải ngẫu nhiên: đó là
giá của một lần BioPython đọc một bản ghi ~300 kB, và nó giống nhau ở mọi module có đọc gen.

**Hình dạng đắt tiền: một phép align cả bộ gen.** Ba module align **tham chiếu + mọi genome, cả
chuỗi**: `substitutions` (scope mặc định `whole`), `diversity` và `barcodes` — cùng một phép toán
trên cùng đầu vào, nên một đường cong mô tả cả ba:

| Số genome | Thời gian | RAM đỉnh | Kết cục |
|---|---|---|---|
| 4 | 1 ph 12 s | 91 MB | xong |
| 8 | 2 ph 44 s | 97 MB | xong |
| 16 | — | — | **vượt 7 phút**, dừng |
| 34 | — | — | **vượt 11 phút**, dừng |

**Giới hạn ở đây là sự kiên nhẫn, không phải RAM.** Align cả bộ gen của 34 plastome đỉnh dưới
100 MB trên máy 8 GB — trần bộ nhớ mà dự án lo (§9.28) **không phải** thứ chặn nó. Thứ chặn nó là
thời gian: MAFFT căn tiến dần trên 5,2 Mbp ở độ phân kỳ này tốn hàng chục phút một lõi, và chi phí
tăng **nhanh hơn** số genome — 4 → 8 genome gấp đôi thời gian, 8 → 16 gấp hơn hai lần nữa.

**Đường rẻ là đường mà hai bản thảo vốn đã dùng.** `substitutions --set scope=cds` thoát **trước**
khi dựng alignment cả bộ gen (mã module: `if scope == "cds": return _cds_scope(...)`) và align theo
từng gen. Đó đúng là scope **chỉ CDS · gen dùng chung · không tính IR** của số vàng §6, và là con số
hai bản thảo trích dẫn. Đo ở N = 34: **5 ph 24 s · 179 MB · xong** — trong khi scope cả bộ gen chưa
xong sau 11 phút trên đúng 34 chuỗi đó. Scope CDS tốn RAM hơn (179 MB so với 97 MB, vì giữ 79
alignment theo gen trong bộ nhớ) và là cái duy nhất chạy xong: trên máy 8 GB, **RAM là thứ rẻ hơn**.

**Hai module theo từng cặp, ở N = 34** — và đây là chỗ một dự đoán đọc từ code bị số đo cãi:

| Module | Lượt gọi tool | Thời gian | RAM đỉnh |
|---|---|---|---|
| `synteny` (`nucmer` → `delta-filter` → `show-coords` mỗi query) | 99 | **19,6 s** | 145 MB |
| `popgen` (MAFFT mỗi gen, thống kê tự viết) | 80 | **6 ph 13 s** | 179 MB |
| `kaks` (MAFFT + `yn00` mỗi gen dùng chung, mỗi cặp) | 5.202 MAFFT + 4.766 yn00 = **9.968** | **8 ph 26 s** | 98 MB |

**Số lượt gọi tool không phải thứ đắt — kích thước một alignment mới là.** 9.968 alignment hai
chuỗi ngắn xong trong 8 ph 26 s; **một** alignment 34 plastome nguyên vẹn thì 11 phút chưa xong.
Và 33 lượt `nucmer` — thứ được đoán là nửa chậm — hoá ra là phần nhanh nhất: **19,6 giây**. Đây
đúng là lý do phép đo này tồn tại: giá của một module **không đọc được từ vòng lặp của nó** (con số
viết tay ở bản nháp đầu nói "khoảng 2.600 mỗi loại"; đo được là 5.202 và 4.766).

**Quy tắc chia lô, viết ra để dùng:**

- Việc một-genome **không cần chia lô**: 34 genome × 9 module dưới 4 phút, dưới 100 MB.
- Align **cả bộ gen** thoải mái tới **8 genome** (2 ph 44 s), quá 16 là vượt 7 phút ⇒ trên mức đó
  dùng scope theo gen.
- Thiết kế **theo từng cặp** (nhiều lượt gọi nhỏ) chạy tốt ở N = 34: `kaks` 9.968 lượt gọi vẫn
  xong trong 8,5 phút.
- **Ngoại suy bị cấm**: mọi con số ở đây đo trên máy này (i7-1165G7, WSL2, trần 8 GB). Chúng chặn
  một *hình dạng* chi phí, không phải một lời hứa cho máy khác.

### 8.33. Ba cách đọc một chuỗi — công bố · PGA · Chloe trên 34 plastome (P4a-04, 04/10/2026)

**Bối cảnh.** Bộ dữ liệu có 34 chuỗi trần (`.fasta`) và hai bản chú giải của **chính 34 chuỗi đó**.
P4a-04 cho annotator của dự án (Chloe 0.1.18) chạy trên 34 chuỗi trần — lần đầu annotator gặp một
bộ dữ liệu thật ở n = 34 — rồi đặt ba cách đọc cạnh nhau. Chuỗi giống nhau từng byte ở cả ba
(đã kiểm: 34/34 cặp, và 68/68 FASTA khớp bản ghi đã pin, §8.31), nên **mọi khác biệt dưới đây là
khác biệt chú giải, không phải khác biệt dữ liệu.**

**Một điều kiện phải mở trước.** `plastid run -i` chỉ đọc GenBank, nên không trỏ được vào một
`.fasta` trần. Đã thêm cửa FASTA (chọn reader theo **dòng đầu của file**, không theo đuôi file —
`.fa`/`.fasta`/`.fna`/`.fas` đều có ngoài kia), và việc đó làm lộ ra lỗi thứ ba cùng họ: `orf`
**nói** trong comment rằng bản ghi không có CDS là input hợp lệ, nhưng bảng per-CDS được dựng
**trước** dòng bảo vệ nó, nên đường đó chưa từng chạy (§9.I.3).

**Kết quả: 34/34 chú giải được, không bản ghi nào lỗi.**

| | CDS | tRNA | rRNA |
|---|---|---|---|
| Bản công bố (bộ 01) | 2.924 | 1.292 | 272 |
| PGA (bộ 02) | 2.887 | 1.292 | 272 |
| **Chloe 0.1.18** | **2.922** | 1.292 | 272 |

Chloe: **0 CDS lỗi kiểm ORF** trên 2.922 CDS (201 `alternative_start`, còn lại `ok`), **0** gen
không xếp được nhóm, 340 tên in kèm synonym, ~35 s mỗi bản ghi. Mọi đầu ra đọc lại bằng chính
reader của dự án (đó là cách module `annotation` sinh ra bảng).

**Không bản ghi nào ba cách đọc đồng ý (0/34)** — kể cả ở mức đếm gen. Khác biệt chia hai loại, và
tỉ lệ giữa chúng mới là kết quả:

**Loại 1 — từ vựng: cùng locus, khác cách viết.** Ở 4 bản ghi, bản công bố viết `trmA-UGC` ·
`trnL-CAU` · `trnL` (trần) · `trnG-UCC` nơi **PGA và Chloe cùng** viết `trnA-UGC` · `trnI-CAU` ·
`trnL-UAG` · `trnG-GCC`. Không có anticodon leucine `CAU` — bộ ba Leu là UUA/UUG/CUU/CUC/CUA/CUG —
nên `trnL-CAU` là cách viết của `trnI-CAU`: đúng họ với alias `trnL-GUA → trnI-GAU` đã ghi ở §8.28.

**Loại 2 — kết luận khác nhau: gen có/không, hoặc bao nhiêu bản sao.**

| Khác biệt | Số bản ghi |
|---|---|
| `ycf1` có CDS ở bản công bố, **cả PGA và Chloe đều không gọi** | **24/34** |
| `rps12` (trans-spliced): Chloe gọi **2 bản sao**, hai bên kia gọi 1 | **22/34** |
| `pafII` (`ycf4`) chỉ bản công bố có | 2/34 |
| `rps19` chỉ bản công bố có | 1/34 |
| `rps16` chỉ bản công bố có (chính bản ghi *Libertia* ở §9.I.2) | 1/34 |

**Đọc kết quả này cho đúng — đây không phải bảng xếp hạng annotator.**

- Không bản ghi nào ba bên đồng ý, nên *"annotator nào đúng"* không phải câu bộ dữ liệu này trả
  lời được. Câu nó trả lời được: **ba đường ống đọc cùng một chuỗi ra ba tập gen khác nhau, và
  khác nhau ở những chỗ xác định** — không phải nhiễu rải rác.
- Năm khác biệt về gen đều xoay quanh hai locus khó nhất của plastome: **`ycf1`** (nằm ở ranh giới
  IR/SSC, gen dài nhất, nhiều bản ghi cụt hoặc chồng `ndhF` — đúng cái bẫy §8.11 đã ghi cho genome
  vàng, nơi ycf1 cụt và **không có CDS**) và **`rps12`** (trans-spliced, 3 exon, hai bản sao nằm ở
  hai IR — gọi 1 hay 2 bản sao là quyết định của annotator, không phải của chuỗi).
- **`132 = 86 + 38 + 8` không phải hằng số.** Trong bộ này bản công bố cho 84–87 CDS, Chloe 84–87,
  PGA 84–85. Đúng cái bẫy P2-10 đã đo (131 và 133 ở hai loài khác) — và là lý do không check nào
  được hard-code "86 CDS".

**§8.15 đo lại ở n = 34.** Kết luận cũ — PGA bị thay sau khi đo **một** genome (147 gen so với 132
của genome vàng) — nay có nền rộng hơn: trên 34 plastome, PGA **mất 37 CDS** so với bản công bố và
**3/34 bản ghi rớt gate IR** (§8.31). Chloe, trên cùng bộ dữ liệu: 34/34 chú giải, 0 CDS lỗi ORF.

**Còn lại của pha:** không — P4a-01…04 xong.
