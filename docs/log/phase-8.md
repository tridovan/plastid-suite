# Nhật ký đo lường — Giai đoạn 8 — kiểm thử chéo, hiệu năng, tự động hoá (M8.1 · M8.3 · M8.4)

**Phần này là `PROJECT.md` §8** (**§8.60–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.60` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.60. CPGview — khảo sát, đo độ truy cập, và quyết định KHÔNG đưa vào bộ công cụ (P8-01, 10/10/2026)

**Việc.** M8.1 và §13.8 nêu **CPGview** là một trong ba công cụ đối chiếu chéo, nhưng §7.8 chưa
từng có dòng nào về nó — luật phủ (`ROADMAP.md` → "How to use" 5) gọi đó là **lỗ hổng**, không
phải một ghi chú. P8-01 là quyết định: khai nó ra thành một dòng §7.8, hay ghi *khảo sát — không
đưa vào* kèm lý do đo được.

**Đo được (10/10/2026, từ máy dev).**

| Phép đo | Kết quả |
|---|---|
| Bài công bố | Liu và cs. 2023, *Molecular Ecology Resources* 23:694–704, doi:10.1111/1755-0998.13729 |
| Repo `github.com/cliu6688/CPGView` | **chỉ là trang giới thiệu** — không mã, không release, không tệp giấy phép, không CLI; README trỏ về web |
| Web service `http://www.1kmpg.cn/cpgview/` | trang ngoài **HTTP 200**, nhưng nội dung là `<frameset>` trỏ ứng dụng sang `http://47.96.249.172:16085/cpgview/home` |
| Ứng dụng (cổng 16085) | **không kết nối được** — `curl` trả `http_code=000`, TCP connect đóng |
| Phép kiểm chứng (control) | cùng IP đó **cổng 80 mở**, và một phép thử cổng cao độc lập **`portquiz.net:16085` mở, HTTP 200** ⇒ **không phải mạng phía mình chặn cổng cao**; ứng dụng thật sự không trả lời |
| Bản cài được duy nhất | container Figshare `cpgview_latest`: **2.871.624.635 B (~2,87 GB)**, md5 `4b371a5ce9410a4a433c9499a6b04bbd` (mục Figshare ghi CC BY 4.0) |
| Bản pip / conda | **không có** (tra PyPI và bioconda/conda-forge: 0 kết quả) |

**Quyết định (P8-01 · D-h, dev): khảo sát — KHÔNG đưa vào bộ công cụ.** Bốn lý do, ba cái đo được:

1. **Không có gì để pin.** Không pip, không conda, không CLI, repo rỗng ⇒ không có phiên bản để
   ghim và không có argv để ghi — đúng thứ mà luật "một công cụ chỉ nằm trên một máy thì không
   phải là một phần của bộ công cụ" (`CLAUDE.md`) cấm.
2. **Bản cài được là 2,87 GB cho một công cụ vẽ hình.** Bộ công cụ đã có hệ hình riêng (M5.1 —
   tám loại hình, SVG/PDF/EPS byte-identical), và bản cài tối thiểu pip-only là yêu cầu không
   thương lượng.
3. **Đường web đụng §10.12** (dữ liệu chưa công bố — còn treo), và hôm nay **web cũng không chạy**:
   phép kiểm một lần trên bản ghi đã công bố mà D-h đề nghị **không thực hiện được** — không phải
   vì chính sách, mà vì ứng dụng không trả lời (bảng trên).
4. **Nó không so một con số nào.** M8.1 là đối chiếu **số** (IR · π · RSCU · SSR · KaKs);
   CPGview vẽ hình và báo cấu trúc exon/intron — không có số nào của mình để đặt cạnh nó.

**Hệ quả cho M8.1.** Nhánh CPGview của M8.1 **đóng bằng phép đo này**, không bằng một lần chạy:
"đối chiếu với CPGview" từ đây có nghĩa là *đã khảo sát, không đưa vào, lý do đã đo* — và §13.8
thôi nêu tên một công cụ mà repo không thể giải trình. Nếu sau này service sống lại và nhóm muốn
một phép kiểm hình ảnh, việc mở lại có địa chỉ: dòng §7.8 này.

### 8.61. Đối chiếu chéo với CPStools 3.0 — năm mặt trùng nhau, một lỗi của mình, và bốn khiếm khuyết của công cụ kia (P8-02, 10/10/2026)

**Việc.** M8.1: đặt số của mình cạnh số của công cụ khác **trên cùng một bộ dữ liệu đã ghim**,
công cụ ghi tên + phiên bản + argv, và **mỗi khác biệt phải có nguyên nhân được gọi tên** — không
đếm. Vòng này so với **CPStools 3.0** (Xu Wenbo, MIT), mặt trùng nhau là IR · π · RSCU · SSR · KaKs.

**Hai bên, đo trên đúng một bộ vào.**

| | Bên mình | CPStools |
|---|---|---|
| Bộ vào | `tests/data/NC_050833.1.gb` sha256 `a946f23e…` (≡ `.panel/iridaceae-34/NC_050833.1.gb`) | cùng file |
| Cặp π/KaKs | `+ tests/data/NC_024936.1.gb` | cùng cặp |
| Bản chạy | `plastid run -m structure · ssr · codon_usage · diversity · kaks` | `cpstools 3.0` (`cpstools -v`), pip trong env riêng — `docs/tools/cpstools.md` |
| Aligner | MAFFT **7.526** (§7.8) | cùng bản, truyền qua `-m /root/toolbox/env/pts/bin/mafft` |
| Ka/Ks | PAML **yn00 4.10.10** (Yang & Nielsen 2000) | **KaKs_Calculator 2.0.1**, `-m YN` (cùng họ phương pháp) |

**Kết quả, năm mặt.**

| Mặt | Bên mình | CPStools | Kết luận |
|---|---|---|---|
| **IR** | LSC 1–83.127 · IRa 83.128–109.341 · SSC 109.342–127.510 · IRb 127.511–153.724 | `LSC:1-83127 IRb:83128-109341 SSC:109342-127510 IRa:127511-153724` | **Bốn khoảng trùng khít**; hai bản IR **đổi nhãn cho nhau** — quy ước, không phải khác biệt (cổng IR đo hai bản giống **100,00 %**, nên gọi bản nào là "a" là chuyện của người viết) |
| **SSR** | 37 = mono 31 · di 4 · tri 2 | 37 = mono 31 · di 4 · tri 2 | **37/37 dòng cùng toạ độ**; 35 dòng trùng cả lớp ngữ cảnh (Gene/Intron/IGS) — 2 dòng lệch là **cụm compound** đã ghi từ §8.1/§8.2, xem dưới |
| **RSCU** | 61 codon (bỏ stop) | 64 codon (**có** stop) | Trên **cùng một bộ gen**: đếm **trùng từng codon** (21.861 codon sense, 0 lệch); RSCU lệch tối đa **0,0049** = đúng mức làm tròn 2 chữ số của họ |
| **π** | theo gen · theo IGS · theo intron | theo gen · theo "IGS" (có cả intron) | **100/104 gen trùng tới 5 chữ số** (lệch trung bình 2×10⁻⁶); IGS: 77/83 trùng tới 5 chữ số |
| **Ka/Ks** | 79 gen dùng chung, 78 chấm điểm | 76 gen | Ka lệch trung bình **0,00014**; Ks 0,00105; Ka/Ks 0,0102 — nguồn lệch gọi tên ở dưới |

**Bốn khác biệt, bốn nguyên nhân.**

1. **Cụm compound (SSR, 2 dòng).** Golden có cụm `(TA)6 + (A)10` tại 48.659–48.679. Bên mình in
   **cả hai thành phần với span của cụm** (scope của bảng nói đúng thế); cpstools in **span riêng
   của từng thành phần** (48.659–48.670 · 48.670–48.679). Không bên nào sai — đây đúng là quy ước
   đã đo ở §8.1/§8.2 (khác biệt 1 bp ngày đó hoá ra là một cụm compound), nay gặp lại ở dạng đầy đủ.
2. **Bộ gen của RSCU (mọi codon).** cpstools khử trùng lặp theo **tên** (86 feature CDS → **79**
   locus; 7 tên có hai bản trong IR: `ycf2` · `ndhB` · `rpl2` · `rpl23` · `rps12` · `rps19` ·
   `rps7`), rồi mặc định `-l 300` bỏ tiếp 28 trình tự ngắn (**còn 51**), và bỏ trình tự không mở
   đầu bằng ATG (`rps19` · `rpl2` · `ndhD` — đúng 5 CDS mà §9.J-39 đã đo). Chạy `-l 0` còn **76**;
   trên **đúng 76 trình tự ấy**, phép đếm codon của hai bên **trùng từng dòng** — nên khác biệt
   không nằm ở công thức, mà ở **tập gen**, và tập gen ấy do tham số của họ quyết định.
3. **π theo IGS (6/83 dòng).** Sáu spacer ngắn có indel: bên mình lấy cột từ **alignment cả
   plastome**, họ align **từng đoạn tách rời**. Đo để chứng minh: lấy đúng hai trình tự họ tách ra
   (69 bp và **72 bp** — họ tách mỗi bản ghi một khác) đem align bằng MAFFT của mình → 74 cột, 67
   cột dùng được, **2** vị trí khác ⇒ π 0,0299. Cùng một locus, ba cách align cho **ba** đáp số:
   1 (mình) · 2 (cặp của họ qua aligner của mình) · 6 (họ) — ví dụ `rps3-rpl22`: π **0,0149** so với
   **0,0896**. Đây là khác biệt **phạm vi alignment**, không phải công thức π.
4. **Ka/Ks.** Cùng họ YN nhưng **hai bản cài khác nhau** (PAML yn00 4.10.10 vs KaKs_Calculator 2.0.1).
   Ka khớp tới ~4 chữ số; chỗ lệch lớn nhất là những gen **Ks rất nhỏ**, nơi mẫu số bé làm sai số
   Ks phóng đại thành sai số tỉ lệ — `psbK`: Ka 0,0080 vs 0,0073 nhưng Ks **0,0178 vs 0,0223** ⇒
   Ka/Ks 0,4494 vs 0,3279. Cả hai bên cùng gọi Ks là *không xác định* ở **cùng một họ gen** (mình
   25 dòng `ks = 0`, họ 23 dòng `NA`) — không có đột biến đồng nghĩa thì tỉ lệ không tồn tại, và
   **cả hai cùng nói thế** thay vì bịa một số.

**Một lỗi của mình, tìm ra bởi chính phép so này — và đã sửa.** Cột `Context` của module `ssr`
xếp **11/37 dòng sai**: bản ghi vàng chú giải `rps12` **trans-spliced** (exon 1 ở LSC, exon 2–3 ở
IR, **cả ba exon ghi trên mỗi feature CDS**) nên **span ngoài** của nó là 68.975–139.634 — gần cả
genome. Đọc span ngoài thì mọi SSR nằm giữa hai con số đó thành "Intron", **kể cả SSR nằm trong
`ycf1`, `ndhA`, `accD`**. Test hồi quy viết **trước** khi sửa (`tests/unit/test_ssr.py`, 5 ca),
luật sửa lấy đúng định nghĩa `core/spans.py` đã có: khoảng giữa hai part của một gen chỉ là intron
khi **hai đầu cùng một vùng**; cắt qua ranh giới vùng là **junction trans-splicing**, không phải
intron. Sau khi sửa: **37/37 dòng khớp cpstools** trừ đúng 2 dòng compound ở mục 1. Chi tiết:
§9.M-1.

**Bốn khiếm khuyết của cpstools, ghi lại vì chúng là dữ liệu để giải thích, không phải chuẩn**
(§7.9's rule) — đầy đủ ở `docs/tools/cpstools.md`: `Pi` chỉ chạy khi `-d` là đường dẫn **tuyệt
đối** (ghi file theo `abspath`, đọc theo tham số thô — đo được, `FileNotFoundError`); nó **từ chối
chạy lại** khi thư mục output còn (`The save directory has been existed…`); `KaKs` với `-m` không
hợp lệ in một lỗi mỗi gen rồi vẫn kết thúc bằng **"All KaKs calculations are completed"** — im lặng
báo thành công; và `RSCU` **tính cả ba codon stop** vào bảng RSCU.

**Điều phép so này KHÔNG kiểm.** CPStools không có mặt Ts/Tv — con số vàng ấy vẫn chỉ đối chiếu với
hai bài công bố (§6), không qua công cụ này. Nói ra để khoảng phủ không bị đọc thành "đã kiểm hết".

### 8.62. Đối chiếu chéo với CGAS 1.0.1 — bốn mặt, và một khoảng trống của chính mình (P8-03, 10/10/2026)

**Việc.** Mặt thứ hai của M8.1. CGAS là bộ 14 mô-đun (Abdullah, MIT) đã đọc mã từ 02/10 (§7.7 ·
§7.9); P8-03 hỏi đúng một câu: **chạy được gì ở đây, và cái không chạy được thì vì thiếu gì** —
thay vì báo một phép kiểm chưa từng chạy.

**Dựng và chạy.** `uv venv /root/toolbox/env/cgas` + `uv pip install -e .data/cgas-review` ⇒
`cgas --version` → **CGAS v1.0.1** (`docs/tools/cgas.md`). **R không có trên máy này**, và đây là
phép đo đáng giữ: các mô-đun 8 · 9 · 12 · 13 của CGAS **tính số bằng Python rồi mới gọi R để vẽ** —
thiếu R thì nó in `⚠ R is not installed … Skipping visualization` và **vẫn ghi đủ bảng số**. Nên
bốn mặt dưới đây chạy thật, chỉ phần hình là không. Không chạy được: mô-đun 1 (cần read FASTQ),
2 (chú giải — PGA, đã nghỉ hưu §8.56), 14 (cần IQ-TREE + MACSE).

**Bốn mặt, đo trên bản ghi vàng.**

| Mặt | Bên mình | CGAS | Kết luận |
|---|---|---|---|
| **Cấu trúc (mô-đun 7)** | LSC 83.127 · SSC 18.169 · IR 26.214 · GC 37,85/35,98/31,45/43,05 | y hệt, từng ô | **trùng khít cả bốn khoảng và bốn giá trị GC** |
| **RSCU (mô-đun 8)** | 61 codon sense, 26.263 lượt | 26.263 lượt sense | **0 codon lệch**; RSCU lệch tối đa 5×10⁻⁵ = đúng mức làm tròn 4 chữ số khi mình xuất bảng. CGAS in thêm **3 codon stop** (86 = 37+27+22 — đúng một stop cho mỗi feature CDS, tức nó đếm cả 86 CDS như mình, **khác cpstools**) |
| **SSR (mô-đun 12)** | 37 = LSC 27 · SSC 10 · IR 0 · mono 31 · di 4 · tri 2 | mặc định: **54** = LSC 39 · SSC 13 · **IR 2** · 31/12/3/4/3 | **Một tham số, không phải một bất đồng**: mặc định của CGAS là `{1:10, 2:5, 3:4, 4:3, 5:3, 6:3}`, MISA là `{1:10, 2:6, 3:5, 4:5, 5:5, 6:5}`. Chạy lại `-t "10,6,5,5,5,5"` ⇒ **37 = 27 · 10 · 0 · 31/4/2 — trùng khít bên mình**, kể cả số 0 ở IR |
| **Intron (mô-đun 11)** | **17** intron CDS | **17** intron CDS (in trong **15 dòng** — gen hai intron dùng chung một dòng, cột Intron I/II) | **Từng dòng trùng toạ độ và độ dài** (541–1.055 bp) |

**Phát hiện thứ nhất — CGAS im lặng ở đúng chỗ khó, và số 0 của nó là số sai.** Trên
`NC_024936.1` (*Iris gatesii*), mô-đun 7 in tên loài rồi **không tìm thấy vùng nào, không một cảnh
báo**, ghi vào bảng `LSC 0 · SSC 0 · IR 0`, và vẫn kết thúc bằng **"Successfully analyzed: 2/2
genomes"**. Bên mình trên chính bản ghi ấy: **LSC 82.702 · IR 26.184 · SSC 18.371**, hai bản IR
**100,00 %** giống nhau — và đường de-novo của mình cũng tìm ra đúng 26.184 bp, nên đây không phải
"bản ghi khó". Đây là họ lỗi §7.9 đã đo ở phía CGAS ("không thấy IR ⇒ im lặng"), nay thấy từ phía
bảng số: **bản ghi không sai, số 0 mới sai** — và là lý do cổng IR của dự án tồn tại (§9.A).

**Phát hiện thứ hai — một khoảng trống của chính mình: intron tRNA.** CGAS in thêm **8 intron
tRNA** trên bản vàng — `trnK-UUU` (1.403–4.019, **2.617 bp**) · `trnG-UCC` · `trnL-UAA` · `trnV-UAC`
· `trnI-GAU` ×2 · `trnA-UGC` ×2. Bản ghi vàng **có đúng 8 feature tRNA nhiều phần** (đo lại từ
chính bản ghi: 38 tRNA, 8 cái có intron) ⇒ **CGAS đúng, và mô-đun `introns` của mình không in
chúng**: nó suy intron từ **part của CDS** (`core/spans.py`), nên con số 17 của mình là *17 intron
CDS*, còn intron tRNA **không mô-đun nào in ra**. Đặt tên: **H-23** — không sửa trong task này
(mở rộng `introns` là việc của board, và phải đo lại cột π theo intron của `diversity` nếu span đổi).

**Điều phép so này KHÔNG kiểm.** Mô-đun 3 (so tên gen), 10 (SNP), 13 (π) và 14 (phylogeny) **không
chạy lần này** — không phải vì thiếu dữ liệu, mà vì phạm vi P8-03 chốt ở ba mặt có đầu ra tất định
(IR · codon · SSR) cộng mặt intron; ba mặt kia đã có phép so tương ứng ở §8.61 (π) và ở bộ vàng
(SNP), nên chúng không bị bỏ quên, chỉ là không lặp lại ở đây.
