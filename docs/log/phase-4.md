# Nhật ký đo lường — Giai đoạn 4 — so sánh đa genome

**Phần này là `PROJECT.md` §8** (§8.26–8.30), tách ra ngày 04/10/2026 để `PROJECT.md` đọc được trọn
(§8 khi đó chiếm 1.531/2.379 dòng, 118/216 KB). **Số mục giữ nguyên** như khi còn nằm trong
`PROJECT.md` — mọi tham chiếu kiểu `§8.27` trong code, test, `docs/` và board vẫn trỏ đúng chỗ —
và **không bao giờ đánh số lại**. Mục mới của giai đoạn này ghi tiếp vào **cuối file này**, rồi
thêm một dòng vào bảng mục lục ở `PROJECT.md` §8.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.26. Synteny chạy thật — "21 khối" không tái lập được, và con số thật của cặp Iris (P4-01, 04/10/2026)

**Bối cảnh.** P4-01 (M4.2, biểu đồ % tương đồng kiểu mVISTA) mang check của pha: *"biểu đồ phải ra
**21 khối**, %IDY 83,7–98,2 (mean 91,5) — kết quả của chính probe"*. D4-04 đã ghi trước: tham số
nucmer chính xác **không có** trong cây probe; phải dựng lại, và nếu không dựng lại được thì **đó là
một phát hiện cho §9, không phải một con số để chỉnh**.

**Đầu vào đã kiểm trước khi đo.** `/mnt/d/business/pts/probe/run_local/{target,compare}.fa` so
từng base với hai bản ghi vàng: `target.fa` **= NC_050833.1** (153.724 bp), `compare.fa`
**= NC_024936.1** (153.441 bp) — đúng cặp, đúng chiều. Vậy đầu vào không phải nghi vấn.

**Quét tham số (MUMmer 4.0.1, bioconda).** `--mum`/`--mumreference`/`--maxmatch`/mặc định ×
`delta-filter` {`-1`,`-g`,`-q`,`-r`,`-m`,`-h`,`-i`,`-o`} × `nucmer -l` ∈ {15…40} × `-c` ∈ {18…65}
× lọc độ dài tối thiểu {0…2.000 bp}: **không tổ hợp nào cho 21 khối.** Số khối dao động 7–43 và
%IDY luôn ở dải 97–99 (giá trị nhỏ nhất từng thấy: 72,34). Vài mốc:

| Tham số | Khối thô | Sau `delta-filter -1` | %IDY (min–max) |
|---|---|---|---|
| `--maxmatch` (mặc định `-l 20 -c 65`) | 12 | **7** | 89,66–99,11 |
| `--maxmatch -l 20 -c 40` | 22 | 7 | 72,34–100,00 |
| `--maxmatch -l 15 -c 15` | 435 | 7 | 59,76–100,00 |
| `--mum` | 10 | 7 | — |

**Kết luận đo được.** *I. domestica* × *I. gatesii* là **hai loài cùng chi**: vùng thẳng hàng
giống nhau **98–99 %**, không phải 91,5 %. Dải 83,7–98,2 (mean 91,5) mô tả một so sánh **khác
chi** — đo đối chứng ngay trong phiên này (*I. domestica* × *O. diffusa*, khác họ): 74 khối thô,
81,89–98,84, mean 89,54, cùng vùng giá trị. Nguồn gốc con số: dòng §7.1 *"mVISTA | sống | biểu đồ
% tương đồng | nucmer + matplotlib (21 khối, PNG 243.787 B)"* — đo **26/9** trong đợt khảo sát
công cụ, trên dữ liệu khác, rồi được chép vào §13.4/ROADMAP như thể là của cặp Iris. Ghi thành
phát hiện **§9.H**; check của pha thay bằng số đo được kèm tham số (bảng dưới).

**Số thật của module** — pipeline §7 ghi (`nucmer --maxmatch → delta-filter → show-coords`), tham
số truyền tường minh `-l 20 -c 65 -t 1`, lọc `-1`, target ở **hệ chuẩn hoá** (hai genome vàng vốn
đã bắt đầu tại LSC nên `rotation_offset = 0`):

| | |
|---|---|
| Khối thẳng hàng | **7** (đều xuôi; hai genome cùng chiều) |
| %IDY mỗi khối | 98,00 · 98,85 · 98,52 · 98,53 · 98,91 · **89,66** · 99,11 |
| %IDY min–max · mean (không trọng số) | **89,66–99,11** · **97,37** |
| Cặp base thẳng hàng | **152.332** / 153.724 = **99,09 %** coverage |
| Identity trên cặp base (bỏ gap) | **99,31 %** |
| Series | cửa sổ 600 bp / bước 200 bp, 769 cửa sổ |

**Định dạng delta của MUMmer 4 — đo và khoá bằng test.** Hai điều không có trong tài liệu nào của
pipeline và sai thì im lặng: (1) file delta ghi **một dòng `>` cho mỗi cặp trình tự**, mọi
alignment của cặp nằm dưới header đó, mỗi bản ghi kết thúc bằng dòng `0` — đọc theo kiểu "một
header một alignment" sẽ trả về **1 khối** cho một so sánh có 7; (2) trong danh sách gap, giá trị
`d` tiêu thụ `|d| − 1` cột **có cặp base**, rồi `d > 0` đẩy con trỏ target một mình, `d < 0` đẩy
query. Kiểm trên cả 12 alignment thô của cặp vàng: **mọi khối kết thúc đúng toạ độ header khai, và
`mismatch + base lẻ` = đúng số lỗi của header** (khối 1: 152 + 120 = 272). Parser ở
`core/delta.py`, fixture thật ở `tests/data/iris_pair_sample.delta`.

**Công cụ đã khai.** `mummer4=4.0.1` pin trong `environment.yml` (§7.8 trước đó ghi *không có
phiên bản* — giờ là bản đo được), `docs/tools/mummer4.md` ghi nguồn, dòng cài, argv và hai số đo
định dạng ở trên.

### 8.27. Ka/Ks chạy thật — và ba locus mà cách viết tên che mất (P4-02, 04/10/2026)

**Bối cảnh.** P4-02 dùng PAML `yn00` (D4-01) trên alignment codon của từng gen. Checklist có một
mục tưởng như phụ — *"gen mà một genome không mang phải được LIỆT KÊ, không được bỏ im lặng"* —
và chính nó tìm ra một lỗi thật.

**Cách chạy.** Mỗi gen dùng chung: MAFFT căn **chính chuỗi CDS** (không dịch protein rồi chiếu
ngược), cắt codon stop cuối, **kiểm khung ba-ngôi** (yn00 đọc 3 base một lần và *không* kiểm
khung — đo được: một gap 1 base làm nó báo stop codon giả rồi exit 255, tức là may mắn chứ không
phải kiểm tra), rồi mới gọi yn00. dN/dS đọc từ ma trận `2YN.dN`/`2YN.dS` và **đối chiếu chéo**
với dòng Yang & Nielsen của chính output — hai đường code trong PAML, lệch nhau thì không lấy
trung bình.

**Kết quả đo (cặp vàng, sau khi sửa tên — xem dưới):**

| | |
|---|---|
| Gen dùng chung | **79** |
| Chấm được | **78** |
| Không chấm được | **1** — `accD`: MAFFT đặt gap sau 1.445 base, không phải bội số của 3 ⇒ xếp `frameshift`, **không** đưa cho yn00 |
| Ks = 0 (tỉ số không xác định) | **25 gen** — kiểm lại: các gen này **giống hệt** hoặc lệch **1 base** giữa hai bản ghi, nên Ks = 0 là số thật |
| Pooled (trọng số theo site của yn00) | Ka **0,003586** · Ks **0,016550** · Ka/Ks **0,2167** |
| Ví dụ | `rbcL` Ka 0,0019 · Ks 0,0299 · **0,0635**; `matK` 0,0087/0,0174 = **0,5000**; `ndhF` 0,0094/0,0473 = 0,1987 |

**Lỗi tìm ra: ba locus bị cách viết tên che mất.** Lần chạy đầu báo `atpI`·`psaI`·`psbB`
"absent from this genome" và `atpL`·`psaL`·`bspB` "absent from the reference". Đo từng cặp:

| Cặp | Dài | Khác biệt | `/product` (giống nhau) |
|---|---|---|---|
| `atpI` (NC_050833.1) / `atpL` (NC_024936.1) | 744 bp | **0 base** | "ATP synthase CF0 A subunit IV" |
| `psaI` / `psaL` | 111 bp | 1 base | "photosystem I subunit VIII" |
| `psbB` / `bspB` | 1.527 bp | 7 base | "photosystem II p680 chlorophyll A apoprotein" |

Trên cả 13 bản ghi của panel: `atpL`/`psaL`/`bspB` chỉ xuất hiện ở **đúng một** bản ghi
(NC_024936.1); 9 bản ghi có chú giải dùng tên chuẩn. Ba alias vào `data/gene_aliases.yaml` kèm
bằng chứng, canonical theo tên chuẩn (luật D3-02).

**Hệ quả dây chuyền — số vàng Ts/Tv §6 đo lại.** Bộ gen dùng chung của cặp vàng từ **76 → 79**,
nên phép đo Ts/Tv theo scope CDS đổi theo: **2,686** (317 Ts / 118 Tv trên 66.066 cột) →
**2,6917** (323 Ts / 120 Tv trên 68.448 cột). Đây là **đo lại vì phép đo cũ đếm gen theo cách
viết tên**, không phải chỉnh số cho khớp công cụ; phạm vi (chỉ CDS · gen dùng chung · MAFFT ·
bỏ IR) không đổi. Số cũ giữ nguyên trong §8.25 như số đo của bộ hẹp hơn. Ba locus thêm vào chỉ
lệch ít nên tỉ số đổi rất nhỏ (+0,2 %) — nhưng nếu chúng lệch nhiều thì đây đã là một bảng
manuscript sai.

**Bài học.** *"Gen không được chia sẻ"* và *"gen bị viết tên khác"* là hai chuyện khác nhau, và
chỉ một mục checklist bắt **liệt kê** cái không chấm được mới làm lộ ra. Cùng họ với §8.9–8.10
(bốn tên chức năng) và §8.17 (`ndhI`/`ndhL`).

**Công cụ đã khai.** PAML 4.10.10 đã pin sẵn trong `environment.yml` (cho `mcmctree`); nay
`yn00` của cùng gói là estimator của pha này — `docs/tools/paml.md` ghi ctl, ba cái bẫy đo được
(không có flag phiên bản; đọc codon mà không kiểm khung; stop codon là thuốc độc) và cách gọi.

### 8.28. Di truyền quần thể viết tay — công thức, hai lần kiểm chứng, và cái bẫy Fst (P4-03, 04/10/2026)

**Bối cảnh (D4-03).** DnaSP là GUI Windows (không script được, không chạy trong Docker);
pixy/libsequence cần VCF hoặc dữ liệu read mà dự án không có. Nên các estimator được viết trong
`core/popgen.py` bằng Python thuần (chỉ stdlib) — bản cài tối thiểu chạy được, đúng luật
CLAUDE.md. Đổi lại, mỗi công thức phải có **nguồn** và mỗi hàm phải kiểm được bằng **vector giải
tích**.

**Kiểm chứng 1 — đối chiếu một implementation độc lập.** Tajima's D, π, θ(S) và S khớp
**scikit-allel 1.3.13** (gói Python độc lập, có bài bình duyệt) trên 4 vector ngẫu nhiên
(n = 4…12, S = 21…243), **chính xác đến 6 chữ số** — ví dụ n=12: D = −0.454318, π = 72.6970 ở cả
hai bên.

**Kiểm chứng 2 — công thức Fu & Li và một lần đo thất bại đáng ghi.** Công thức D\*/F\* có nguồn
Fu & Li (1993) và hằng số phương sai của Simonsen, Churchill & Aquadro (1995); chúng được đọc
trực tiếp trong mã nguồn **EggLib 3.6.1** (`Diversity1.cpp`) — một implementation đã công bố.
Nhưng **EggLib build từ sdist trên máy này trả về site view sai** (đo trên vector giải tích: site
thứ sáu ra `['A','A']` trong khi dữ liệu là A/C, vị trí 9.0 thay vì 5), nên nó **không dùng được
làm đối chứng số**. Vì vậy D\*/F\* được khoá bằng **vector tính tay từng bước** trong test
(n = 4, một singleton: D\* = −0,6123724356957945 · F\* = −0,4787135538781689). Đường khác cũng
đóng: `gestimator` của libsequence **không có trên bioconda** (package chỉ ship thư viện — đo:
`bin/` chỉ có `libsequenceConfig`), còn `Bio::PopGen` không còn trong BioPerl hiện hành.
Ghi lại vì "không kiểm được bằng công cụ" là một sự thật phải nói, không phải một chỗ để đoán.

**Kiểm chứng 3 — một định danh thú vị, và một cái bẫy thật.** Với **n = 4 và mọi đột biến đều là
singleton**, D\* ≡ D (đo: cả hai = −0,7545 ở S = 3…63; chứng minh được: tử số tỉ lệ 11 và phương
sai cũng tỉ lệ 11², nên hai thống kê trùng nhau). Không phải lỗi — nhưng là lý do phải in **cả
hai** cột chứ không suy cột này từ cột kia.

**Fst — viết Hudson trước, đo rồi bỏ.** Bản đầu dùng Hudson, Slatkin & Maddison (1992)
(`1 − H_w/H_b`). Đo trên vector giải tích: **hai nhóm có cùng thành phần haplotype (n = 2 mỗi
nhóm) cho Fst = −1**, vì H_w chỉ là một cặp — một artifact mẫu nhỏ sẽ vào bảng bản thảo như một
giá trị âm. Thay bằng **Nei's Gst (1973)**: hai nhóm cùng thành phần ⇒ **0** đúng, không haplotype
chung ⇒ **1** đúng, mất hết đa dạng ⇒ để trống. Gst đọc **tần số haplotype**: các gen mà mọi
trình tự đều khác nhau cho **cùng một giá trị Fst** bất kể chúng xa nhau bao nhiêu (đo trên panel:
33/80 gen cùng ra 0,1852) — mức độ khác biệt nằm ở π, và module ghi rõ điều đó.

**Kết quả đo trên panel stand-in (11 bản ghi, 2 nhóm Coleus 6 · Oldenlandia 5):** 80 gen được
căn; phần lớn mỗi nhóm n = 4 (gen có mặt ở cả 4 bản ghi có chú giải của nhóm). Ví dụ `rbcL`
Oldenlandia: n = 4 · S = 41 · π 0,01844071 · D 1,8434; `matK` Coleus: D −0,7284. **19 dòng để
trống** vì n < 4 hoặc không có site phân ly — mỗi dòng mang lý do trong cột Note ("n = 3: cần ít
nhất 4 trình tự"), và 3 bản ghi **không có CDS** (AP045470.1 · PQ325515.1 · PZ428230.1) được
liệt kê trong notes thay vì lặng lẽ góp zero. Fst giữa hai nhóm có đủ 80 dòng.

### 8.29. Barcode + vùng siêu biến chạy thật — π dùng lại, và hai lỗi tên lộ ra (P4-04, 04/10/2026)

**Bối cảnh.** P4-04 xếp hạng cửa sổ/vùng theo độ biến thiên, gọi tên các locus mã vạch ứng viên,
và **dùng lại** π của module diversity (cùng một hàm lõi, cùng window/step) chứ không tính lại.

**Đo được (cặp vàng, 600/200, 773 cửa sổ):** mean π **0,007332** · ngưỡng vùng **0,014157**
(mean + 1 sd của chính run) · **26 vùng siêu biến** · 10 locus ứng viên, **8 tìm thấy**:

| Locus | Kiểu | Vị trí (bp) | Dài | π | Hạng |
|---|---|---|---|---|---|
| atpF-atpH | intergenic | 12.889–13.263 | 375 | 0,024064 | 1 |
| psbK-psbI | intergenic | 7.284–7.677 | 394 | 0,015424 | 2 |
| ndhF | gene | 109.327–111.540 | 2.214 | 0,013134 | 3 |
| ycf1 | gene | 122.984–128.404 | 5.421 | 0,011852 | 4 |
| matK | gene | 1.663–3.231 | 1.569 | 0,010835 | 5 |
| trnL-trnF | intergenic | 47.944–48.324 | 381 | 0,010499 | 6 |
| rbcL | gene | 55.212–56.672 | 1.461 | 0,008898 | 7 |
| rps16 | gene | 4.697–5.808 | 252 (2 exon) | 0,000000 | 8 |
| trnH-psbA | intergenic | — | — | — | không tìm thấy |
| ITS | nuclear | — | — | — | không phải locus lục lạp |

**Bài kiểm "dùng lại": hai module, một con số.** Một tier-B test chạy `-m barcodes -m diversity`
trên cùng dữ liệu rồi so **từng cửa sổ**: 773/773 π trùng nhau. Đây là kiểm chứng *giữa hai
module*, không phải một bản cài lại.

**Phát hiện 1 — `trnH-psbA` không tồn tại như một spacer trong hai bản ghi vàng.** Bản đầu của
module báo locus này "not carried" — một lời khai **sai** về dữ liệu. Đo lại: `trnH-GUG` và `psbA`
**đều có**, nhưng **rps19 nằm giữa** chúng ở cả hai bản ghi (trình tự quanh điểm gốc:
`rpl2 → trnH-GUG → rps19 → [gốc] → psbA`). Dòng bảng giờ nói đúng: *"present but not adjacent in
this record's annotation"*. Kèm theo, code xử lý spacer **vắt qua origin** (bản ghi bắt đầu ở đâu
là quy ước của người nộp) được viết và có test riêng — đo được rằng một spacer hoàn toàn có thể
nằm hai bên vị trí 1.

**Phát hiện 2 — hai tên tRNA viết thường làm mất locus trong im lặng.** `NC_024936.1` viết
`trnk-UUU` và `trnl-GUA`; bản vàng viết `trnK-UUU` và `trnI-GAU`. Đo bằng trình tự: locus
`trnk-UUU` **giống hệt từng base** `trnK-UUU` (77 bp), và `trnl-GUA` **giống hệt từng base**
`trnI-GAU` (77 bp, cả hai bản IR). Hai sửa:
* **quy tắc cú pháp**: chữ cái axit amin trong `trnX-YYY` không phân biệt hoa/thường —
  `trnk-UUU` → `trnK-UUU`; `trnfM-CAU` **không** bị đụng (chữ `f` của nó là một phần của tên, ký
  tự kế tiếp là `M` chứ không phải `-`);
* **một alias mới**: `trnL-GUA` → `trnI-GAU` (GUA không phải anticodon của Leu; bản ghi không có
  `/product` lẫn `/anticodon` để đối chiếu, nhưng trình tự thì trùng khít).
Không sửa thì mỗi so sánh theo tên — gene table, thống kê tên, và bất kỳ module nào khoá theo
tên chuẩn — **âm thầm bỏ hai locus** khỏi bảng.

### 8.30. Bộ bảng so sánh cho hai bản thảo — ba run trên panel đã pin (P4-05, 04/10/2026)

**Bối cảnh.** P4-05 là chính sản phẩm của pha: nội dung so sánh cho hai bản thảo. Ba run, tất cả
trên panel đã pin (dữ liệu công khai, **không phải mẫu của nhóm**), không cần mạng:

| Run | Tham chiếu | Nội dung |
|---|---|---|
| `standin-oldenlandia` | `NC_057985.1` (*O. diffusa*, RefSeq) | substitutions (scope CDS) · kaks · synteny · barcodes · popgen |
| `standin-coleus` | `PX923121.1` (*C. amboinicus* — **bản ghi của chính nhóm**) | như trên |
| `standin-mito` | `OX459128.1` × `AP045469.1` | kaks · popgen — bản phác khác chi (D4-05) |

**Nhãn dữ liệu đi vào từng bảng.** Thêm `--label TEXT` cho `plastid run`: nhãn được ghép vào
`scope.label` của **mọi** bảng ở đúng một chỗ — runner, nơi duy nhất biến `Result` thành file —
nên nó có mặt ở header CSV, JSON và report mà không module nào phải nhớ. Panel chạy với
`--label "stand-in (public records, not the team's samples)"`.

**Một lỗi thật do panel làm lộ ra.** Scope CDS của `substitutions` yêu cầu *mọi* genome phải có
tên gen chung; panel có **3/11 bản ghi không có chú giải** (0 CDS) ⇒ bảng rỗng cho cả run
(`ValueError: table 'substitutions_by_gene' has no rows`, đo trên lần chạy đầu). Sửa: bản ghi
không có CDS **được nêu tên trong notes** và bị loại khỏi so sánh CDS, thay vì làm hỏng cả run;
không còn genome có chú giải nào ⇒ lỗi nói rõ lý do. Test riêng cho cả hai nhánh.

**Vài số đo được (hai run loài, 04/10/2026):**

**standin-oldenlandia** (tham chiếu NC_057985.1, 8 genome có chú giải): Ts/Tv scope CDS gộp **1.5157** = 16959 Ts / 11189 Tv trên **77 gen dùng chung** và 68743 cột; **bảng theo từng cặp** (`substitutions_by_genome.csv`) cho 0.9583–1.5173. Synteny: gần nhất `PP078674.1` 1 khối @ 98.39 %, xa nhất `NC_057679.1` 61 khối @ 89.07 %. Barcodes: 29 vùng siêu biến, 9/10 locus ứng viên, biến thiên nhất `trnH-psbA` (π 0.251894), `psbK-psbI` (π 0.212766). Ka/Ks pooled theo từng genome: 0.1585, 0.1007, 0.0986, 0.0853, 0.0857, 0.0959, 0.0868. Popgen: 80 gen, Coleus: 6, Oldenlandia: 5, 19/160 dòng không có kiểm định trung tính, 2/80 dòng Fst để trống.

**standin-coleus** (tham chiếu PX923121.1, 8 genome có chú giải): Ts/Tv scope CDS gộp **1.5395** = 15235 Ts / 9896 Tv trên **77 gen dùng chung** và 68743 cột; **bảng theo từng cặp** (`substitutions_by_genome.csv`) cho 1.3065–1.5424. Synteny: gần nhất `PZ428230.1` 1 khối @ 99.70 %, xa nhất `PP078674.1` 63 khối @ 89.10 %. Barcodes: 20 vùng siêu biến, 9/10 locus ứng viên, biến thiên nhất `trnL-trnF` (π 0.260571), `atpF-atpH` (π 0.250487). Ka/Ks pooled theo từng genome: 0.1663, 0.2468, 0.2067, 0.0853, 0.0831, 0.0856, 0.0922. Popgen: 80 gen, Coleus: 6, Oldenlandia: 5, 19/160 dòng không có kiểm định trung tính, 2/80 dòng Fst để trống.

**Một phát hiện về vận hành (ghi thành H-08, §9.C-35):** cache của run layer nằm ở
`.runs/.cache` **dùng chung giữa các run** và khoá cache **không gồm mã nguồn module** — sửa
`substitutions` xong chạy lại thì nó trả `cached` và phục vụ payload cũ (thiếu cả bảng mới), vì
`MODULE_VERSION` chưa đổi. Luật "đổi hành vi thì bump version" đã bị vi phạm trong chính phiên
này; bốn module đổi hành vi được bump lên **0.2.0**, cache cũ bị xoá, và việc mở H-08 ghi lại
cách sửa tận gốc.

**Một hạn chế đo được, phải nói ra:** căn CDS bằng MAFFT ở cặp **khác chi** đặt indel lệch khung
ở **19–21/80 gen** (cùng chi 5–11/80), nên các gen đó mang trạng thái `frameshift` chứ không có
số — MAFFT không phải công cụ căn theo codon, và đây là cái giá của việc căn thẳng trên chuỗi CDS
(đúng như yêu cầu "không dịch protein rồi chiếu ngược"). Con số nằm trong bảng chứ không bị giấu:
cột Status của `kaks_genes` in ra đúng lý do cho từng gen.

**Nửa mito — hai lỗi nữa, cả hai đều là "một cái tên che mất dữ liệu".** Lần chạy mito đầu
tiên chết với `table 'popgen_groups' has no rows`, và lý do không phải dữ liệu thiếu:
1. **Bản ghi tham chiếu đúng nhưng bản ghi truy vấn sai**: `AP045469.1` trong panel
   `oldenlandia-coleus-mito` là **bản ghi thô không có feature** (0 CDS — đúng như đo lúc pin);
   bản chú giải của chính dự án nằm ở `.panel/mito-prototype/` (P3A-07). Run dùng nhầm file thô ⇒
   không có gen nào để so. Sửa: dùng file prototype, và `popgen` nay **báo lỗi nói rõ** khi không
   có gen nào được ≥2 genome mang, thay vì để lỗi validation của bảng.
2. **Hai genome, không gen nào trùng tên**: `OX459128.1` **không có `/gene`** trên bất kỳ feature
   nào — "tên" của mỗi CDS thực ra là **câu `/product`** ("NADH dehydrogenase subunit 1") — trong
   khi bản chú giải của mình ghi ký hiệu chuẩn (`nad1`). Bảng từ vựng mito
   (`data/mito_gene_names.yaml`, có từ P3-04) map đúng những câu đó, nhưng các module so sánh
   **truyền bảng alias plastid một cách tường minh**, nên từ vựng mito không bao giờ được hỏi.
   Sửa: `core/gene_names.aliases_for(genome)` chọn bảng theo `/organelle` của chính bản ghi, và
   `cds_by_canonical_name` mặc định dùng nó; bốn module so sánh bỏ tham số bảng tường minh. Đo
   lại: OX 27 CDS có tên, prototype 23, **19 gen dùng chung** (trước đó: 0).

**Kết quả mito (`standin-mito`)** — 19 gen chung; kaks chấm **12/27**, pooled Ka 0,028443 · Ks
0,112924 · Ka/Ks 0,2519; phần còn lại mang lý do (6 frameshift · 8 không có ở genome này · 4
không có ở tham chiếu · 1 stop codon nội bộ). popgen: mỗi nhóm **1 trình tự** ⇒ mọi thống kê để
trống kèm lý do, **19/19 dòng Fst không tính được** — đúng như thiết kế sau khi thấy bẫy Fst.

**Điều chỉnh ở `popgen` sau khi thấy mito:** Fst cần **ít nhất 2 trình tự mỗi nhóm** — với nhóm
1 trình tự, đa dạng trong nhóm bằng 0 theo cấu trúc và Fst luôn ra đúng **1,0** cho hai nhóm khác
nhau (một con số trông như "khác biệt hoàn toàn" nhưng chỉ là artifact cỡ mẫu). Từ chối thay vì
đưa số, kèm lý do trong cột Note.

