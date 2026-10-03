# PlastidToolSuite — ngữ cảnh dự án

> Nguồn ngữ cảnh duy nhất của dự án. Đọc trước khi làm bất cứ việc gì.
> Cập nhật: 02/10/2026. Toàn bộ nội dung là **dữ kiện đã kiểm tra / đã chạy thật**.
> Mục lục nhanh: §1 bối cảnh · §2 mục tiêu & yêu cầu đầu ra · §3 quyết định · §4 vai trò · §5 nền tảng · §6 số vàng · §7 công cụ (nguồn tải: §7.8) · §8 nhật ký đã chạy · §9 lỗi đã gặp · §10 câu hỏi treo · §11 nguồn tài liệu · §12 lệnh cài · **§13 ROADMAP — mô-đun & thứ tự triển khai**.

---

## 1. Bối cảnh

Bộ công cụ phân tích **bộ gen lục lạp (plastome) và bộ gen ty thể** cho một nhóm nghiên cứu thực vật. **Công cụ học thuật, không thương mại.**

- Người dùng cuối: nhà sinh học, **không rành kỹ thuật**, dùng **máy tính phổ thông**.
- ~20 mô-đun, mỗi mô-đun là một việc phân tích: cấu trúc genome, nội dung gen, mã bộ ba (codon/RSCU), SSR, lặp lại phân tán/palindromic, chú giải, lắp ráp, so sánh nhiều genome, RNA editing, ước tính thời gian phân kỳ, chuyển gen giữa các bào quan (NUPT/MTPT)…
- Nguồn gốc: kế hoạch 22/9/2026 của cô Trang + 2 báo cáo kiểm tra công cụ của anh Điềm (cài thật — chạy thật).

## 2. Mục tiêu & sản phẩm đầu ra

| Sản phẩm | Nội dung |
|---|---|
| 2 bản thảo loài | Oldenlandia và Coleus — cần bảng số + hình + **bản nháp mục Results** |
| 1 bài báo công cụ | Eleutherine bulbosa (đang phản biện) — phần chính: NUPT/MTPT và cấu trúc đa hình ty thể |
| Phát hành | Bộ công cụ chạy trên máy phổ thông, có mã trích dẫn (DOI) |

### Yêu cầu bắt buộc cho đầu ra (đã chốt 02/10)

1. **Report ngắn gọn về các công cụ đã chạy kèm thông số.** Mỗi lần chạy phải xuất một bản ngắn: công cụ · phiên bản · tham số/ngưỡng đã dùng · thời gian chạy · dữ liệu vào.
2. **Output chia theo thư mục** để người dùng truy lại được — tách theo lần chạy / công việc / loại kết quả (bảng số, hình, dữ liệu trung gian, report), thay vì để chung một chỗ.
3. **Khắc phục nhược điểm của CGAS về đầu ra** (hiện trạng CGAS ghi ở §7.7): thư mục output tên cố định `ModuleN_...` nằm trong thư mục làm việc hiện tại ⇒ chạy lại cùng dữ liệu thì kết quả lẫn/đè lên lần trước; tên file phần lớn cố định; không có cấp "lần chạy" hay "dataset"; không có file mục lục và không có report công cụ + thông số.

## 3. Quyết định đã chốt

| Ngày | Quyết định |
|---|---|
| 22/9 | Xây bộ công cụ; mục tiêu máy phổ thông |
| 26/9 | Công cụ bị chặn → thay bằng bản local đã chạy thật (không còn phụ thuộc web) |
| 26/9 | PGA = GPL-3.0; Vmatch = chỉ miễn phí học thuật |
| 01/10 | **Làm phần phân tích trước, phần lắp ráp sau** (cô Trang đồng ý) |
| 01/10 | Đích: **Linux native + Windows thông qua Docker**; **macOS để sau** |
| 02/10 | Đầu ra phải có **report công cụ + thông số**; **output chia theo thư mục** để truy lại; **khắc phục nhược điểm đầu ra của CGAS** (§2 và §7.7) |
| 02/10 | ROADMAP chốt: **khung/hạ tầng là Giai đoạn 1** (làm trước tiên) · **phạm vi Ts/Tv = chỉ CDS, không tính IR** · `misa.ini` dùng **mặc định MISA** |
| 02/10 | Nguồn MISA đóng bằng md5 đối chiếu bản tải lại từ IPK (§7.8); MAFFT = conda-forge 7.526 (§7.8) |
| 02/10 | **MISA đóng gói kèm** bộ công cụ — pin đúng byte: zip 4.946 B · md5 `cf9d3579f7caacd02f99830b329c7d9c`; image phải có **perl** (bản anh Điềm dùng: 5.32.1). Giấy phép MISA: **không nêu rõ** (zip không kèm file giấy phép; trang IPK chỉ ghi trích dẫn Beier và cs. 2017) ⇒ đưa vào rà soát giấy phép Giai đoạn 9 |
| 02/10 | **Bộ số vàng tách 2 tầng**: tầng A — thuần Python, không cần công cụ ngoài, luôn chạy; tầng B — cần MISA/cpstools/MAFFT, thiếu công cụ thì **skip kèm lý do**, không im lặng bỏ qua |
| 02/10 | **Không dùng Snakemake** ở giai đoạn đầu (kế hoạch 22/9 mục 0.14 vốn chỉ là "cân nhắc"); chỉ xem lại khi pipeline vượt ~20 bước hoặc cần chạy lại từng phần |
| 02/10 | **Reference Panel**: chốt định nghĩa 4 tầng taxonomy (nguyên văn kế hoạch 22/9, mục 0.7 chế độ B) + 2 chế độ A/B — ghi ở §13.1 |
| 02/10 | **Hệ toạ độ** (§10.6): bảng theo từng genome ghi **toạ độ GenBank** (truy ngược về bản nộp NCBI); hình so sánh chéo nhiều genome dùng **hệ đã chuẩn hoá**; **mọi bảng ghi rõ đang ở hệ nào**; `rotation_offset` lưu trong `00_run_info` |
| 02/10 | **Report sinh từ `Result`** (§10.7): mọi số trong report trỏ về **đúng một ô bảng gốc**; template chỉ trỏ `{{bảng.cột}}` / `{{bảng[dòng].cột}}`, không gõ số; ref không giải được ⇒ **raise**, không in trống |
| 02/10 | **Reference Panel** (§10.8): engine **gợi ý** qua cache RefSeq có ngày (làm mới được); khi người dùng chốt danh sách thì **tải về + ghim checksum**; chạy lại **không cần mạng**. Mục (c) điền sau, không chặn code |
| 03/10 | **Backend (P1-12, dev quyết)**: giai đoạn 1–5 và 8 chạy **tại chỗ** trên máy phổ thông; giai đoạn 6 (thời gian phân kỳ) + 7 (lắp ráp) chạy **máy chủ từ xa** — máy 8 GB chỉ còn 1,5–2,5 GB trống mỗi việc (§9.14). Ghi lại theo luật phân quyết định trong CLAUDE.md |
| 03/10 | **Bảng gen xuất bản (P2-09, dev quyết)**: bảng in ra bài báo là bảng **phân nhóm theo chức năng**, không phải bản dump từng feature; nhóm gen là **dữ liệu** (`src/plastid/data/gene_groups.yaml`), không phải code. tRNA/rRNA khớp theo **tiền tố tên** (`trn`, `rrn`) chứ không theo danh sách anticodon — bản Oldenlandia/Coleus thật sẽ mang anticodon khác, danh sách cứng sẽ **rơi im lặng**. Gen không nhóm nào nhận ⇒ một dòng *unclassified* + ghi tên trong notes, **không bao giờ bị bỏ** |
| 03/10 | **Số đếm của bảng gen = số bản sao** (kể cả 2 bản IR), đúng quy ước 132 = 86 + 38 + 8, và bảng **tự đối chiếu với `count_genes()`** trong một bảng riêng — cùng kiểu đối chiếu đã bắt lỗi MISA ở P2-06. Một gen rơi giữa bộ đếm và bảng sẽ hiện ra thành *DISAGREES*, không phải thành một số nhỏ hơn không ai kiểm |
| 03/10 | **Xuất Word (D2-04 + P2-09, dev quyết)**: `python-docx` là **extra tuỳ chọn** (`pip install 'plastid-suite[docx]'`), import **bên trong hàm** nên bản cài tối thiểu vẫn thuần Python; thiếu nó ⇒ `ExportDependencyError` + exit 3, không phải traceback ImportError. `openpyxl` **chưa thêm** — chưa có gì dùng (YAGNI); sẽ vào cùng extra khi làm Excel ở Giai đoạn 5 |
| 03/10 | **Tài liệu Word phải tái lập được từng byte**: `python-docx` đóng dấu thời gian đồng hồ vào zip, nên gói được **dựng lại** với timestamp cố định (1980-01-01), entry sắp xếp, `create_system` ghim ⇒ cùng đầu vào ra cùng file. Đây là định nghĩa "Done when" của P2-09, không phải một lời hứa |
| 03/10 | **Gói .docx được tỉa**: template của `python-docx` mang theo cả catalogue style của Word (164 style, ~350 kB) + `stylesWithEffects.xml` (438 kB) + thumbnail + `customXml`; bản xuất chỉ giữ style tài liệu **thực sự dùng** cùng tổ tiên và default của chúng, mọi quan hệ trỏ tới phần đã bỏ đều bị cắt. Đo được: ~800 kB → **13 kB**, có test cấu trúc canh |
| 03/10 | **[ĐÃ THAY — xem dòng dưới] Tên gen: ánh xạ theo locus, giữ tên gốc (D3-02, dev quyết có bằng chứng)**: `pafI`·`pafII`·`pbf1`·`clpP1` trong `PX923121.1` là **cùng locus** với `ycf3`·`ycf4`·`psbN`·`clpP` — chứng minh bằng synteny + độ dài + `/product` trên 9 bản ghi (§8.9). Tên **chuẩn** giữ bản cổ điển (7/9 bản ghi + genome vàng); tên gốc đi kèm `source_name`, không bao giờ bị bỏ. Tên **in trong bản thảo** và **nộp NCBI** vẫn là quyết định của PI (D3-06). *Chiều chuẩn hoá ở dòng này bị đảo ở dòng dưới, sau khi biết hai bản ghi là của nhóm (§8.10)* |
| 03/10 | **Tên gen: khoá theo locus, in cả hai tên (D3-02, dev quyết có bằng chứng)**: bốn tên `pafI`·`pafII`·`pbf1`·`clpP1` là **cùng locus** với `ycf3`·`ycf4`·`psbN`·`clpP` (synteny + độ dài + `/product`, §8.9). **Khoá so sánh nội bộ = phương ngữ của nhóm** — vì `PX923121.1` là bản ghi **của nhóm** (Vu,H.-T. = cô Trang) và **đúng loài bản thảo**, và GeSeq sẽ sinh ra đúng bộ tên đó; tên cổ điển đi kèm `source_name`. **In ra thì in cả hai** (`pafI (ycf3)`) — cách các bảng gen công bố đang làm. `gene_groups.yaml` nhận thêm bốn tên vào đúng nhóm chức năng (§8.10) |
| 03/10 | **Chú giải: công cụ offline là pipeline, GeSeq là bước kiểm tra một lần rồi ghim (D3-01, dev quyết)**: bộ công cụ **phải chạy offline và đóng gói bằng Docker**, mỗi lần chạy phải ghi **tên · phiên bản · đường dẫn · argv** — công cụ web không ghi được cái nào trong số đó, nên đây là **ràng buộc, không phải sở thích**. Chọn **PGA** (GPL-3.0, batch, ra GenBank; ba cái bẫy của nó — thư mục ra phải rỗng, tên file có dấu cách, multi-FASTA — repo này đã tránh sẵn), plann.pl dự phòng. GeSeq chạy tay **một lần cho mỗi loài bản thảo** rồi **ghim** (`panel pin`) ⇒ bước thủ công chỉ xảy ra một lần, các lần sau offline và tái lập được. Kiểm 03/10: GeSeq không có bản chạy rời | Dev | `docs/D3-01_annotator.md` |
| 03/10 | **Đa genome (D2-01, dev quyết)**: giữ nguyên contract v0.1 `run(genome, params, ctx)`; thêm contract thứ hai `run(genomes, params, ctx)` chọn bằng `arity` khai báo lúc đăng ký — **không** nới rộng contract cũ, vì module một-genome sẽ phải nhận một dãy nó không dùng (seam 5). Runner cưỡng chế số genome trước khi chạy. Ghi ở `docs/API_v0.2.md` |
| 03/10 | **π: cửa sổ 600 bp · bước 200 bp (D2-02)** — quy ước kiểu mVISTA; window/step vào params, manifest và report |
| 03/10 | **ORF: đọc bảng mã của chính bản ghi, từ chối thay vì đoán (D3-03 + P3-01)**: `Gene` mang `/transl_table` và `/codon_start` của bản ghi; bảng mã ngoài {1, 11} ⇒ `OrfError` (hai bảng này chia sẻ toàn bộ axit amin, chỉ khác mã khởi đầu), `/codon_start` ≠ 1 ⇒ `OrfError` — dịch sai khung là kiểu "sai mà không báo". **ACG được nhận là mã khởi đầu nhưng xếp lớp riêng `non_canonical`**: nó *không* có trong danh sách của bảng 11, nhưng bản ghi đã công bố dùng nó cho rpl2/ndhD (RNA editing ACG→AUG); xếp vào `no_start` sẽ là báo động giả trên 2/86 CDS của genome vàng (§8.11) |
| 03/10 | **Lớp "annotated, not coding" khớp cả TÊN lẫn VỊ TRÍ (P3-01)**: ycf1 cụt chồng lấn 21 bp với CDS của `ndhF`, nên khớp theo vị trí không thôi thì ndhF "phủ" nó, còn khớp theo tên không thôi thì CDS ycf1 thật ở chỗ khác "phủ" nó — mỗi cách hỏng một kiểu im lặng (§8.11) |
| 03/10 | **File nộp NCBI: `.fsa` + `.tbl` thuần Python, luôn luôn; `.sqn` chỉ khi có `table2asn` *và* có template `.sbt` (D3-04, P3-05)**: web submission của NCBI nhận thẳng cặp `.fsa`+`.tbl`, nên bản cài tối thiểu vẫn nộp được. `.tbl` mang **tên gốc của bản ghi** (đổi tên là đổi nội dung nộp — D3-06 chỉ nói về *chữ in trong bản thảo*); toạ độ partial đi theo bản ghi, partial ở interval giữa ⇒ **từ chối**, không bỏ im lặng (§8.13) |
| 03/10 | **Chú giải: Chloe là annotator của pipeline; PGA còn lại cho genome ngoài angiosperm (D3-01 — sửa lại sau khi ĐO, dev quyết)**: ràng buộc cũ vẫn nguyên (offline, ghi được tên · phiên bản · đường dẫn · argv), nhưng *"công cụ offline nào"* chỉ là **giả thuyết cho tới khi chạy thật**. PGA trên genome vàng ra **147 gen** (7 reference) và **130** (1 reference) so với **132**, kèm ba kiểu hỏng im lặng (§8.15) ⇒ **thay bằng Chloe**: offline, có bộ reference chuẩn hoá của chính nó, và **đo qua reader của mình ra đúng 132 = 86 + 38 + 8, 0 CDS lỗi ORF** (§8.17). Giữ `--set annotator=pga` cho plastome không phải thực vật có hoa (Chloe chỉ làm angiosperm). GeSeq vẫn là kiểm tra chéo ghim một lần. **Bài học ghi lại:** công cụ thắng không nằm trong danh sách khảo sát ban đầu | Dev | `docs/D3-01_annotator.md`, `docs/tools/chloe.md` |
| 03/10 | **Nối công cụ chú giải mito (P3-04) — PMGA**, không phải MITOS2 (annotator của động vật có vú, bộ tham chiếu `refseq81m`). PMGA phân phối dạng **singularity container** (Figshare 27201798) + web server, không có repo git và không có tài liệu CLI công khai ⇒ **chưa viết adapter**: viết theo CLI tự nghĩ ra sẽ là một placeholder trông như thật. Cần: kéo container, ghi entry point vào `docs/` |
| 03/10 | **Bước nặng: chạy ở đâu, và kết quả chạy ở nơi khác phải ghim được (P1-12, làm nốt 03/10)**: luật backend đã có từ 03/10 (giai đoạn 1–5 và 8 chạy tại chỗ; 6–7 chạy máy chủ từ xa) nhưng **thiếu đường thực hiện** cho công cụ chạy bên ngoài: `panel pin` chỉ nhận accession. Bổ sung `plastid panel pin --file … --note …` (§8.19) để kết quả của web service / máy chủ khác vào được panel với đủ checksum + đọc lại + provenance — PMGA là ca đầu tiên dùng nó | Dev | `docs/tools/pmga.md`, §8.19 |
| 03/10 | **Chú giải mito: PMGA chạy trên web service rồi GHIM, không chạy local (dev quyết, xác nhận bởi anh Trí)**: máy này **không đủ mạnh để chạy annotate** — container PMGA 3,18 GB nén → **11 GB giải nén, 138.045 file**, cần Apptainer và bên trong chạy MAKER; trần WSL 8 GB, ổ C: thật chỉ còn 42 GB. Container đã được tải về, đọc lấy CLI thật, rồi **xoá** (14 GB = một phần ba chỗ trống). Chú giải mito đi theo **đường thứ hai của chính luật P1-12**: chạy PMGA trên web service một lần mỗi loài ⇒ `plastid panel pin --file … --note …` (§8.19) ghim kết quả với đủ checksum + đọc lại ⇒ `mito_genes` đọc offline. **`external/pmga.py` để mở có điều kiện:** viết khi có máy chạy nổi container, và số đầu tiên phải đo khi đó là RAM đỉnh + thời gian | Dev | `docs/tools/pmga.md`, TASKS.md P3-04 |
| 03/10 | **Dữ liệu thay thế cho dev (dev quyết, cô Trang chưa có mẫu)**: khi mẫu thật chưa có, **được tải bản công khai tương tự** (cùng chi, hoặc họ hàng gần nhất) để dev/test — chỉ qua `plastid panel pin` (accession · ngày · sha256), **mọi đầu ra phải ghi nhãn**, **không bao giờ** để số của bản thay thế vào bảng manuscript. Định nghĩa bộ dữ liệu commit ở `data/`; payload `.gb` không commit. Quy tắc ở CLAUDE.md → "Data — real samples vs development stand-ins"; bộ đầu tiên ở §8.6 |

## 4. Ai làm gì

| Người | Vai trò |
|---|---|
| **Vũ Thị Huyền Trang** (cô Trang) | Trưởng nhóm, không rành kỹ thuật; quyết định nội dung khoa học |
| **Nguyễn Thanh Điềm** (anh Điềm) | Sinh tin; đã cài thật và chạy thật 84 công cụ trên WSL2 Ubuntu 22.04 (RAM 3,7 GB); giữ dữ liệu đề tài và các env công cụ |
| **Người dùng (dev)** | Phát triển sản phẩm |
| Hermes | Điều phối, đo đạc, viết tài liệu |
| Claude Code | Viết code |

Xưng hô trong tài liệu tiếng Việt: **cô Trang**, **anh Điềm**.

## 5. Nền tảng kỹ thuật

- Đích phát hành: **Linux native** + **Windows thông qua Docker** (Docker chạy nhân Linux, cùng một codebase). **macOS để sau.**
- Không phát hành native Windows: đo được **win-64: 0/24 công cụ** biên dịch được; linux-64: 24/24; osx-64: 24/24; osx-arm64: 19/24 (thiếu bandage, sibelia, vmatch, oatk, mafft).
- GetOrganelle **không có** trên usegalaxy.org; có trên **usegalaxy.eu** (1.7.7.0+galaxy1) và usegalaxy.org.au (1.7.6.1 và 1.7.7.0+galaxy1); usegalaxy.fr không có.
- Máy dev: Windows 16 GB RAM, WSL2 Ubuntu 24.04; `.wslconfig`: memory=8GB, processors=6, swap=8GB, sparseVhd=true, autoMemoryReclaim=gradual. Docker Desktop + Antigravity đã kết nối; code nằm trong ext4 của WSL.
- Đo tốc độ ghi: `/mnt/d` **183 MB/s** · ext4 trong WSL **1,5 GB/s**.

## 6. Bộ SỐ VÀNG (dùng để đối chiếu)

Mẫu: **Iris domestica NC_050833.1** (bản nộp gốc: `MT001880.1`).

| Chỉ số | Số vàng | Nguồn |
|---|---|---|
| Độ dài | **153.724 bp** | GenBank |
| GC | **37,85 %** | GenBank |
| Số gen | **132** = 86 CDS + 38 tRNA + 8 rRNA | khớp bài công bố Ai và cs. 2019 (PubMed 33366388) |
| Feature `gene` trong GenBank | **133** (thừa 1 bản `ycf1` cụt, không có CDS) | tự đếm |
| Tên gen duy nhất / gen lặp trong IR | 112 / 21 | tự đếm |
| Vùng lặp đảo IR | **26.214 bp**, hai bản giống **100,00 %** | tự tính; chú giải GenBank (`repeat_region` 83.128–109.341 và 127.511–153.724) khớp |
| LSC / SSC | **83.127 / 18.169 bp** (tổng 4 vùng = 153.724) | tự tính |
| SSR | **37** = mono 31 · di 4 · tri 2 | MISA standalone 2.1 **và** cpstools 3.0 (2 công cụ độc lập, cùng ra 37) |
| Ts/Tv (so với Iris gatesii NC_024936.1) | **2,686** = 317 Ts / 118 Tv; **chỉ vùng CDS** (76 gen dùng chung, không tính vùng lặp IR); MAFFT; 66.066 cột căn — phạm vi chốt 02/10 | tự đo |
| Số công bố của loài này | *Iris domestica* 153.729 bp, 132 gen; *Belamcanda chinensis* 153.735 bp, 132 gen | 2 bài công bố |

---

## 7. Danh mục công cụ

### 7.1. Bảng tổng

| Việc | Công cụ | Trạng thái test |
|---|---|---|
| SSR | MISA standalone 2.1 (`misa.pl`) + cpstools 3.0 `SSRs` | đã chạy thật — 37 SSR (cả hai công cụ) |
| Repeat phân tán + palindromic | repeat-match (MUMmer3) + Vmatch 2.3.1 | đã chạy thật — 9 match ≥30 bp, có cặp đảo IR 26.214 bp |
| So sánh nhiều genome (khối LCB, synteny) | Sibelia 3.0.7 | đã chạy thật — 3 plastome |
| Biểu đồ % tương đồng theo toạ độ | nucmer (MUMmer4) `--maxmatch` → `delta-filter` → `show-coords` + matplotlib | đã chạy thật — 21 khối, %IDY 83,7–98,2 (mean 91,5) |
| Chú giải plastome | plann.pl (`github.com/daisieh/plann`) + PGA.pl | đã chạy thật, offline theo reference GenBank |
| Vẽ bản đồ vòng genome | pycirclize 1.10.1 | đã chạy thật — PNG 494.859 B @300 dpi |
| Lắp ráp từ reads ngắn | GetOrganelle 1.7.7.1 + SPAdes 4.3.0 | đã chạy với reads **mô phỏng** — 127.679 bp |
| Lắp ráp long-read (ít RAM) | raven-assembler 1.8.3 (+ OATK 1.0) | đã chạy thật — badread 150x ONT → 2 contig, dài nhất 107.394 bp |
| Lắp ráp long-read + DB HMM | OATK 1.0 + `embryophyta_pltd.fam` (OatkDB) | đã chạy thật — 126.660 bp, identity 100,00 %; file `.annot_pltd.txt` **rỗng** |
| Vẽ + thống kê đồ thị assembly (headless) | Bandage 0.9.0 | đã chạy thật — PNG 103.289 B trên GFA thật của OATK |
| Mạng haplotype | PyPopART | đã chạy thật — 13 trình tự matK → mafft → mjn/analyze/visualize, rc=0 |
| Ước tính thời gian phân kỳ | PAML mcmctree 4.10.10 · BEAST2 2.7.7 | xác nhận chương trình chạy; **chưa** chạy phân tích thật |
| RNA editing (đo từ RNA) | minimap2 + samtools mpileup (+ wgsim sinh dữ liệu) | đã chạy mô phỏng — 20 vị trí C→U: ngoài IR đúng 14/14 (recall 100 %), precision 74 %; 6 vị trí sót trong IR |
| Phân tích so sánh plastome (IR/Pi/RSCU/SSRs/KaKs) | cpstools 3.0 + CGAS | đã cài, đã dùng |

### 7.2. Không cài được / không chạy được (đã thử cả env nhóm và env cách ly)

| Công cụ | Kênh đã thử | Lệnh cài | Kết quả thô (log thật) |
|---|---|---|---|
| progressiveMauve | bioconda (env `~/toolbox/env`) | `conda install -c bioconda progressivemauve` | **KHONG-CHAY-DUOC** — segfault với cả input tối thiểu 40 bp; tái hiện 2 lần |
| PMAT2 | GitHub release (tar.gz 90 MB) | tải `pmat2.tar.gz` | **CAN-THU-CONG** — release là ảnh SIF (Apptainer/Singularity); máy không có Apptainer |
| PopART | web/GUI | ứng dụng GUI Java | **KHONG-CAI-DUOC** — chỉ có GUI, không có CLI |
| BandageNG | GitHub release (AppImage + `--appimage-extract`) | tải `.AppImage` | **KHONG-CHAY-DUOC** — binary đòi GLIBC 2.38 + libOpenGL; máy là Ubuntu 22.04, glibc 2.35 |
| r8s | nguồn (không có trên conda/conda-forge/PyPI/CRAN) | `curl` `r8s1.81.tar.gz` → `make -f Makefile.linux` | **KHONG-CAI-DUOC** — source 1.81 build lỗi (Makefile trỏ `/usr/include/sys/errno.h`); binary kèm theo là Mach-O của macOS |

### 7.3. Cài được nhưng chạy lỗi

| Công cụ | Kênh | Lệnh cài | Kết quả thô | Đã dùng thay bằng |
|---|---|---|---|---|
| flye 2.9.6 | conda | `conda install -p toolbox -c conda-forge -c bioconda flye` | **CAI-DUOC-TEST-LOI (RAM)** — assembly chết vì hết RAM (cần ≥8 GB; máy 3,7 GB) | raven-assembler 1.8.3 / OATK |
| Chloroplot 0.2.4 | conda (env R riêng `~/toolbox/r_bio`) + GitHub | `conda create -p ~/toolbox/r_bio -c conda-forge -c bioconda r-base bioconductor-biostrings bioconductor-genbankr bioconductor-cordon r-reutils r-circlize r-dplyr r-magrittr` rồi `remotes::install_github("shuyuzheng/Chloroplot")` | **CAI-DUOC-TEST-LOI** — env trộn hai dòng Bioconductor (Biostrings 2.70.3 của Bioc 3.19 cạnh GenomicFeatures 1.52.1 của Bioc 3.18) → `genbankr::readGenBank` lệch seqlevels; hạ cấp Biostrings 2.66 không được vì bản đó cần R < 4.3 | circlize trong chính env đó (PDF 6,1 KB) hoặc pycirclize |
| Deepred-Mt | GitHub clone | `git clone aedera/deepredmt` | **CAI-DUOC-THI-THIEU-DEP** — cần TensorFlow (~2 GB), RAM máy 3 GB | minimap2 + samtools mpileup (đo trực tiếp từ RNA-seq) |

### 7.4. Công cụ web (12: 7 còn sống, 5 đã chết) → bản local đã chạy thật

| Web | Trạng thái | Việc | Bản local dùng thay |
|---|---|---|---|
| GeSeq | sống | chú giải | plann.pl + PGA.pl |
| CPGAVAS2 | sống | chú giải + phân tích | plann.pl / PGA.pl |
| Chloe | sống | chú giải theo loài họ hàng | plann.pl |
| mVISTA | sống | biểu đồ % tương đồng | nucmer + matplotlib (21 khối, PNG 243.787 B) |
| MISA-web | sống | SSR | MISA standalone — 37 SSR |
| REputer | sống | repeat phân tán + palindromic | repeat-match (MUMmer3) + vmatch 2.3.1 |
| GetOrganelle (Galaxy) | sống | lắp ráp reads ngắn | GetOrganelle 1.7.7.1 local + SPAdes 4.3.0 |
| DOGMA | **chết** (DNS không phân giải) | chú giải | plann.pl |
| CpGAVAS bản 1 | **chết** | chú giải + phân tích | plann.pl |
| Verdant | **chết** | pipeline chú giải | plann.pl |
| CGAP | **chết** | phân tích so sánh | cpstools 3.0 + CGAS |
| PREP-Cp | **chết** | dự đoán RNA editing | không còn web tương đương (PREPACT, CURE-Chloroplast cũng `http=000`) → đo từ RNA-seq |

### 7.5. Ghi chú kỹ thuật đã gặp khi cài/chạy

- **Krait**: bản hiện tại là ứng dụng GUI (Python/PyQt), README khuyên dùng Krait2 ⇒ không chạy headless được.
- **`plann` trên PyPI** là ứng dụng lịch CalDAV (plann.no), không phải Plann chú giải plastome; Plann thật ở `github.com/daisieh/plann` (Perl).
- **OATK**: cần DB HMM ngoài (`embryophyta_pltd.fam`, 15,2 MB) và `nhmmscan` của HMMER ≥3.4. Lần chạy thật: assembly + GFA tốt nhưng `.annot_pltd.txt` rỗng.
- **treePL**: không có gói conda; phải tự build nlopt (nguồn trong `deps/`) và ADOL-C (không có trên conda). Đã build được cả hai.
- **REputer**: mreps và TRF tìm *tandem repeat* (khác lớp); REputer tìm *exact repeat* kể cả sợi đảo/palindromic ⇒ dùng repeat-match hoặc vmatch.
- **RNA editing**: PREP-Cp, PREPACT, CURE-Chloroplast đều không truy cập được (`http=000`) tại thời điểm kiểm 26/09.

### 7.6. Việc chưa chạy thử (giới hạn đã biết)

- NovoPlasty đã cài nhưng **chưa** chạy với seed thật.
- BEAST2/mcmctree mới xác nhận chương trình chạy, **chưa** chạy một phân tích phân kỳ thật (cây + hiệu chuẩn).
- GetOrganelle mới chạy với reads **mô phỏng**; chưa có reads Illumina thật của đề tài.
- OATK hmmannot chưa ra chú giải (file rỗng) — chưa kiểm `nhmmscan`/DB.
- Ngưỡng gọi RNA editing mới tối ưu trên dữ liệu mô phỏng.
- **macOS chưa kiểm lần nào**; mọi kết quả hiện tại là trên Ubuntu/WSL2.
- PMAT2 bản gốc chưa thử lại (cần Apptainer/Singularity — máy chưa có).
- Chưa có genome nhân (cần cho NUPT/MTPT).

---

### 7.7. Đầu ra của CGAS — hiện trạng (đọc trực tiếp mã nguồn 02/10/2026)

Nguồn: repo `github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS` — 19 file `.py`, 21.369 dòng, license MIT.

| Điều quan sát được trong code | Bằng chứng |
|---|---|
| Thư mục output **đặt cứng theo module** | module3 `module_3` · module5 `Module5_Gene_Comparative_Analysis` · module6 `Module6_Gene_Content_Tables` · module7 `Module7_Comparative_Analysis` · module8 `Module8_Codon_Usage_Analysis` · module9 `Module9_Amino_Acid_Analysis` · module10 `Module10_SNP_Analysis` · module11 `Module11_Intron_Analysis` · module12 `Module12_SSR_Analysis` |
| Input mặc định là **thư mục hiện tại**, output mặc định là tên cứng | module12: `parser.add_argument('-i','--input', default='.')` và `parser.add_argument('-o','--output', default=OUTPUT_FOLDER)`; `OUTPUT_FOLDER = "Module12_SSR_Analysis"` |
| CLI chỉ parse `-m/--modules/--all/--list`; `-i`/`-o` truyền xuống module dưới dạng tham số tự do | `cli.py` dùng `parse_known_args()`; dòng hướng dẫn in ra: *"Results are saved in Module<number>_ folders"* |
| Tên file **phần lớn cố định**, không mang tên dataset hay lần chạy | `ssr_summary.xlsx` · `ssr_by_species.xlsx` · `all_ssrs_detailed.xlsx` · `comparison_results.xlsx` · `nucleotide_diversity_results.txt` · `nucleotide_diversity_by_position.txt` · `complete_matrix.fasta` · `genes_only.fasta` · `introns_only.fasta` · `intergenic_only.fasta` · `genes_for_R_plot.txt` · `noncoding_for_R_plot.txt` |
| Chỉ một số ít file có timestamp | module5: `Chloroplast_Gene_Analysis_{TIMESTAMP}.xlsx` |
| **Module 1 là ngoại lệ**: có cây thư mục con đánh số | `01_QC` → `07_assembled_genomes` (mô tả trong `USAGE.md`); trong code: `self.reports_dir = self.output_dir / "06_reports"` |
| Không thấy **file mục lục / manifest** mô tả "file nào ở đâu" | tìm `manifest`/`index.txt`/`README.txt` trong `cgas/*.py`: không có; mỗi module tự in/ghi thống kê riêng (ví dụ module14 có `summary_file`) |
| Tài liệu có mô tả output theo từng module | `USAGE.md` có mục "Output Structure" ở nhiều chỗ (module 1 rõ nhất); README lưu ý "Keep the output directory outside of the input directory" |

Kết luận về câu hỏi "CGAS có phải output ra một đống không truy được?": **đúng phần cốt lõi** — không có cấp "lần chạy"/"dataset", tên file cố định nên chạy lại là lẫn/đè kết quả cũ, không có file mục lục và không có report về công cụ + thông số. **Không hoàn toàn đúng ở chỗ "không có gì"**: CGAS vẫn tách thư mục theo từng module và module 1 có cây thư mục con đánh số.

### 7.8. Nguồn tải / cách cài từng công cụ

Trạng thái: **ĐÃ GHI** = tài liệu của nhóm có ghi lệnh cài / địa chỉ; **CHƯA GHI** = tài liệu chỉ có tên công cụ, chưa ghi nguồn.

**Qua conda (conda-forge + bioconda)** — lệnh gộp đầy đủ ở §12:

| Công cụ | Kênh | Trạng thái |
|---|---|---|
| sibelia (3.0.7) · bandage (0.9.0) · oatk (1.0) · raven-assembler (1.8.3) · flye (2.9.6) · spades (4.3.0) · getorganelle (1.7.7.1) · novoplasty · mummer4 (nucmer) · mreps · trf · vmatch (2.3.1) · minimap2 · samtools · wgsim · blast · paml (mcmctree 4.10.10) · beast2 (2.7.7) · perl-bioperl · perl-xml-simple | bioconda | ĐÃ GHI |
| pycirclize (1.10.1) | conda-forge | ĐÃ GHI (kiểm 02/10: conda-forge 27 file, tới 1.10.1) |
| cpstools (3.0) | **bản đang chạy: pip/PyPI** (site-packages python3.12 của base env, tác giả Xu Wenbo, MIT) — trong `~/miniconda3/pkgs/` có tarball bioconda `cpstools-3.0-pyhdfd78af_0.tar.bz2` nhưng **không phải bản đang chạy** | Sửa lại: cpstools **không có trên conda-forge**; bioconda có (22 file, tới 3.0) |
| Vmatch — trang gốc của tác giả Kurtz | `http://www.vmatch.de` (bản đóng gói: `anaconda.org/bioconda/vmatch`) | ĐÃ GHI (kiểm 02/10) |
| GetOrganelle — cần tải bộ seed 1 lần (cần mạng) | `get_organelle_config.py -a embplant_pt,embplant_mt` | ĐÃ GHI |
| OATK — cần DB HMM ngoài | `github.com/c-zhou/OatkDB` → `v20230921/embryophyta_pltd.fam` (15,2 MB) | ĐÃ GHI |

**Qua git / pip:**

| Công cụ | Nguồn | Trạng thái |
|---|---|---|
| plann.pl (chú giải plastome, Perl) | `github.com/daisieh/plann` | ĐÃ GHI |
| PyPopART | `github.com/Adamtaranto/pypopart` → `pip install -e '.[viz]'` | ĐÃ GHI |
| Deepred-Mt | `github.com/aedera/deepredmt` (cần TensorFlow ~2 GB) | ĐÃ GHI |
| treePL | `github.com/blackrim/treePL` → phải tự build nlopt (nguồn trong `deps/`) + ADOL-C | ĐÃ GHI |
| PGA / PGA.pl | `github.com/quxiaojian/PGA` (bản 2.0: `github.com/quxiaojian/PGA2`); CGAS `INSTALL.md` cũng ghi PGA phải cài tay từ GitHub | ĐÃ GHI (kiểm 02/10) |
| Chloroplot + env R | `remotes::install_github("shuyuzheng/Chloroplot")`; env: r-base · bioconductor-biostrings · bioconductor-genbankr · bioconductor-cordon · r-reutils · r-circlize · r-dplyr · r-magrittr | ĐÃ GHI |
| CGAS | `github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS` (19 file `.py`, 21.369 dòng, MIT) | ĐÃ GHI (kiểm 02/10) |

**CHƯA GHI nguồn trong tài liệu của nhóm** (chỉ có tên hoặc chỉ có vị trí trên máy):

| Công cụ | Hiện biết | Nguồn chính thức |
|---|---|---|
| MISA standalone (`misa.pl` + `misa.ini`) | nằm ở `~/toolbox/src/misa/`, không có trong PATH → gọi `perl ~/toolbox/src/misa/misa.pl <fasta>` (perl 5.32.1); `misa.ini` phải cùng thư mục | **ĐÃ ĐÓNG (kiểm 02/10, đối md5 với bản tải lại từ nguồn):** https://webblast.ipk-gatersleben.de/misa/ → `misa_sourcecode_25082020.zip` · 4.946 B · md5 `cf9d3579f7caacd02f99830b329c7d9c` · bản **v2.1 (25/08/2020)** · `misa.pl` 20.148 B md5 `cfebc1d1a36fddfd6e0d897c96babd12` · `misa.ini` 158 B md5 `49a1c56048631647c0364532a6964cf9` · tác giả Thomas Thiel & Sebastian Beier (IPK) · trích dẫn Beier và cs. 2017, Bioinformatics 33:2583–2585 (doi:10.1093/bioinformatics/btx198) · Quyết định 02/10: **đóng gói kèm** (xem §3); image cần perl. |
| repeat-match | "có sẵn trong env toolbox" | thuộc **MUMmer3**: `https://mummer.sourceforge.net` (MUMmer3.23 tarball trên SourceForge); MUMmer4 (`nucmer`) thì qua bioconda |
| MAFFT | bản đang chạy: `/home/nthanhdiem/miniconda3/bin/mafft` → `v7.526 (2024/Apr/26)` | **ĐÃ ĐÓNG (kiểm 02/10):** cài từ **conda-forge** — `mafft 7.526` build `h4bc722e_0` (2.610.460 B, md5 `cf3b2acf649bc31e576ed43a2fa4aaa3`); conda-forge cũng có build `h4ab18f5_0` và `ha3fa977_0` của cùng 7.526; bioconda chỉ tới 7.525. Máy còn bản apt: **7.490-1 trên Ubuntu 22.04** (Ubuntu 24.04 là **7.505-1**) ⇒ report/log phải ghi **đường dẫn + phiên bản**, không ghi trần "MAFFT". Chạy lại hàm Ts/Tv bằng cả 2 bản mafft → cùng ra 2,686 ⇒ số vàng không phụ thuộc bản mafft |
| Krait | đã thử, là GUI nên không dùng | chưa ghi nguồn |
| PGA — giấy phép | PGA = GPL-3.0 | (đã ghi ở §3) |

### 7.8b. Cài đặt thật trên máy dev (03/10/2026 — gỡ chặn Giai đoạn 3)

Không có conda trên máy ⇒ dựng **micromamba 2.9.0** (`~/toolbox/bin/micromamba`), env dùng chung
`~/toolbox/env/pts`. **Không thêm gì vào repo**: dự án pin *nguồn gốc* công cụ, không pin *binary*.

| Công cụ | Bản | Kiểm chứng |
|---|---|---|
| **MAFFT** | **7.526** `h4bc722e_0` (conda-forge) | **khớp tuyệt đối bản đã pin**: tarball 2.610.460 B · md5 **`cf3b2acf649bc31e576ed43a2fa4aaa3`** · `mafft --version` → `v7.526 (2024/Apr/26)`. Gỡ chặn H-01 ⇒ P2-07/P2-08 chạy số thật được |
| **BLAST+** | 2.17.0+ (bioconda) | `blastn -version` → `2.17.0+`; PGA cần ≥ 2.8.1 |
| **Perl** | có sẵn | dùng cho MISA và PGA |
| **table2asn** | NCBI, `ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/` | **không có trên conda**; tải từ nguồn chính thức của NCBI |

**Phát hiện làm sai lệch kế hoạch — MITOS2 không dùng được cho ty thể thực vật.** `M3.1` ghi
"MITOS2/GeSeq (ty thể)". MITOS2 tìm gen bằng `blastx` với **protein ty thể của metazoan** (bộ tham
chiếu mặc định `refseq81m`) — sai giới, nên kết quả sẽ **sai im lặng** chứ không rỗng. Công cụ đúng
là **PMGA** (Plant Mitochondrial Genome Annotator, Li và cs. 2025, *Plant Communications*
6(3):101191): làm riêng cho ty thể thực vật, nhận FASTA một-nhiều, ra GFF3 + GenBank, và trong
benchmark của chính tác giả trên *Arabidopsis thaliana* / *Liriodendron tulipifera* cho **1 gen sai**
so với GeSeq 18 và AGORA 42. Module RNA editing của nó chạy **Deepred-mt** — công cụ §7.8 đã ghi.

*Giá phải trả:* PMGA bọc MAKER · CPGAVAS2 · BLASTn · ARGORN · tRNAscan-SE · MISA · TRF · Vmatch —
nặng hơn nhiều so với PGA (Perl + BLAST). Quyết ở P3-04, không phải bây giờ.
*Độ mạnh bằng chứng:* phạm vi của MITOS2 là chắc (trang phân phối và tài liệu nói thẳng); số
benchmark của PMGA đọc từ bản tóm tắt web, **chưa đọc bài gốc** — kiểm lại trước khi chọn.

### 7.8c. LỖI ĐÃ SỬA — công cụ chỉ nằm trên một máy thì không phải là một phần của bộ công cụ

**Chuyện đã xảy ra (03/10/2026, do chính tôi gây ra):** để gỡ chặn Giai đoạn 3, tôi cài MAFFT 7.526
và BLAST+ vào `~/toolbox/` trên máy dev — nhưng **không sửa gì trong repo**. Nếu dừng ở đó thì:

- mọi số tôi đo được **không ai tái lập được** — kể cả anh Điềm, kể cả chính tôi sáu tháng nữa;
- người cài `plastid-suite` theo hướng dẫn trong repo **không có công cụ nào**, và tier B sẽ skip
  sạch mà không rõ vì sao;
- Docker image ở Giai đoạn 9 sẽ không chạy được phần cần công cụ.

**Kiểm tra thì phát hiện lỗi nặng hơn:** `environment.yml` cũ liệt kê công cụ **không kèm phiên bản**
— trong đó có dòng `mafft` trần. Mà §9.B-13 **đã đo**: cùng một đầu vào, aligner khác nhau cho
**Ts/Tv 2,667** so với **2,686**. Tức một người cài theo `environment.yml` cũ sẽ ra **số vàng khác**
với cùng bộ gen — đúng loại "sai mà không báo" mà dự án đang săn, chỉ khác là nó nổ ở máy người khác.

**Đã sửa:** `environment.yml` pin theo đúng phiên bản §7.8 ghi; ba công cụ không có trên conda
(PGA, table2asn, công cụ ty thể thực vật) được ghi rõ nguồn và lý do ngay trong file; dry-run xác
nhận bộ giải chọn **`mafft 7.526 h4bc722e_0` từ conda-forge** — đúng bản đã pin, không rơi về 7.525
của bioconda (bioconda cũng có `mafft`, nên đây là bẫy thật, không phải giả định).

**Luật mới trong CLAUDE.md** (mục Invariants): công cụ phải được khai báo **trong repo** ngay khi
module dùng nó — `environment.yml` (pin phiên bản) · `src/plastid/vendor/` + `SOURCE.json` (script
nhỏ, như MISA) · hoặc `docs/` kèm nguồn và lệnh cài. Cài lên máy để **kiểm chứng** là việc nên làm;
dừng ở đó là việc không được phép — **sửa repo trong cùng phiên, không để "lát nữa"**.

### 7.9. Bài học từ mã nguồn CGAS (đọc trực tiếp 02/10/2026 — 18.353 dòng, 14 mô-đun + CLI)

Đọc để học cách tổ chức, và để không lặp lại lỗi của họ.

**Giữ lại được:** entry point theo tên việc (`cgas-ssr`, `cgas-codon`…) · module 1–2 đếm lỗi theo mẫu
và ghi file mẫu hỏng · class có tham số tiêm được (module 12) · `SUMMARY_REPORT.txt` của module 14
(có ngày, thư mục vào/ra, tham số — thiếu đúng phiên bản công cụ, dòng lệnh, thời gian) · header
docstring mỗi module · `test_files/` có dữ liệu thật.

**Phải tránh (kèm bằng chứng trong mã):**

| Bệnh | Bằng chứng | Bộ công cụ này làm khác |
|---|---|---|
| Sai mà không báo | CLI trả `False` nhưng `main()` bỏ qua → luôn exit 0; bare `except:` 5 chỗ ở module 8, 5 chỗ ở module 9; R lỗi vẫn in "ANALYSIS COMPLETE" | exit code thật · lần chạy hỏng vẫn để lại manifest ghi lý do · không bare except |
| Phân loại IR bằng substring | module12:223 `any(term in combined for term in ["inverted","ir",…])` — chuỗi `"ir"` khớp cả "pair", "circular", "direct" | nhận diện theo **cấu trúc** (độ dài, cặp, ranh giới) + cổng ≥ 99,9 % (§9A) |
| Không thấy IR ⇒ coi cả genome là LSC | module12:259-260 `return {'LSC': (1, seq_len)}` | "IR không xác định" là trạng thái riêng, chặn phân tích theo vị trí |
| SSR vắt ranh giới exon/intron ⇒ dán nhãn "intergenic" | `is_in_parts` đòi SSR nằm trọn trong 1 phần (module12:376-380) → rơi xuống nhánh IGS | có test cho ca vắt ranh giới |
| Không lưu phiên bản công cụ | có chạy `Rscript --version`, `blastn -version`… nhưng chỉ để kiểm tra rồi vứt | `run_tool()` ghi tên · phiên bản · nguồn · đường dẫn · argv · cwd · thời gian · exit code · log |
| Chạy lại là đè | `os.makedirs(exist_ok=True)` + tên file cố định; 12/14 module không có bảo vệ | mỗi lần chạy một thư mục riêng + `manifest.json` |
| Đường dẫn cứng | PGA dò `~/PGA`, `C:\PGA\…`, `f"/home/{username}/PGA/PGA.pl"`; `java -Xmx8g` cứng | đường dẫn công cụ ghi vào run info; không hardcode |
| Tài liệu nói dối | USAGE.md mô tả cây output của module 5/10/12: **0 tên file khớp code**; README dạy lệnh không tồn tại (`--module 5-14`, `--version`) | bộ vàng chạy trong CI — tài liệu trôi là build đỏ |
| Copy-paste | ~2.400 dòng keo R lặp ở 6 module; `find_genbank_files` 3 bản không tương thích (module 3 bỏ sót `.gbff`); 14 hàm wrapper CLI y hệt | engine chung + test chặn tầng import (`tests/unit/test_architecture.py`) |
| Không test được | 0 unit test; pyproject trỏ `tests/` không tồn tại; logic dính `os.getcwd()` lúc import, `os.chdir()`, `sys.exit` trong thư viện | hàm thuần `run(genome, params, ctx) -> Result`; test là hợp đồng |

**Hai chỗ trùng khớp đáng sợ với rủi ro mình đã đo:** (1) họ phân loại IR bằng substring và im lặng
coi cả genome là LSC — đúng loại lỗi IR âm thầm mình đo được 2/9 ca; (2) họ dán nhãn SSR vắt ranh giới
thành "intergenic" — mình đã có test cho ca này.

## 8. Nhật ký các lần đã chạy thật

### 8.1. Cài đặt công cụ — môi trường đã dùng

| | |
|---|---|
| Máy | WSL2 Ubuntu 22.04, glibc 2.35, 8 nhân, RAM 3,7 GB |
| Conda | 24.11.3, kênh conda-forge + bioconda |
| Env | `~/toolbox/env/toolbox` (có sẵn) · `~/toolbox/env/thaythe` (cài mới đợt 26/9) · `~/toolbox/pyenv` (Python cho docx) · `~/toolbox/r_bio` (R cho Chloroplot) |
| Dữ liệu & log | `ke_hoach_new/_thay_the/` |
| Kết quả cài | **84 công cụ đã thử · 64 chạy được · 2 cài được nhưng test lỗi · 17 không cài được/không chạy** (trong 17: 7 công cụ web còn sống, 5 web đã chết) |
| Vị trí không nằm trong PATH | `misa.pl` + `misa.ini` ở `~/toolbox/src/misa/` · `nucmer`, `repeat-match` trong env `mummer4` (phải activate) · `trf` không có |

### 8.2. Chạy chuỗi trên máy anh Điềm (26/09/2026 — 15 giây, 18/19 mốc đạt)

| Mốc | Kết quả thô |
|---|---|
| Độ dài / GC | 153.724 bp · 37,85 % |
| Gen | 132 = 86 CDS + 38 tRNA + 8 rRNA (feature `gene` = 133; tên duy nhất 112; gen lặp IR 21) |
| IR (tự tính) | 26.214 bp · hai bản giống 100,00 % |
| IR (chú giải GenBank) | 26.214 bp — `repeat_region` 83.128–109.341 và 127.511–153.724 |
| 4 vùng | IRa 26.214 + IRb 26.214 + LSC 83.127 + SSC 18.169 = 153.724 |
| Xoay chiều | điểm xoay đặt giữa vùng LSC; IR vẫn 26.214 bp |
| SSR — MISA standalone 2.1 | **37** (1 bp: 31 · 2 bp: 4 · 3 bp: 2; "SSRs present in compound formation: 3"; file `.misa` có 34 dòng) |
| SSR — cpstools 3.0 | **37**; `NC_050833_SSRs_loc_results.txt` có cột vị trí gen: `Gene`/`IGS`/`Intron` |
| Ts/Tv | 317 Ts / 118 Tv = **2,686** (MAFFT, 76 gen CDS dùng chung với Iris gatesii) |
| Công cụ dò được | `misa.pl`+`misa.ini` = có · `cpstools` 3.0 = có · `mafft` = có · `blastn` = có · `nucmer` = không (phải activate env mummer4) · `trf` = không |
| Spike RAM | **chưa chạy** — trong `D:\diem\bai bao tool` không có read thô; trên ổ D: chỉ có `D:\diem\Huy_genome\RO0172510A19_010.fastq.gz` (45 MB, nghi amplicon 16S). GetOrganelle có trong env `verify_plastid` và `cgas`; RAM trống 1,9 GB / 3,7 GB |
| Ts/Tv đã có trong bộ Eleutherine | `KQ_THONG_NHAT\08_SNP\Complete_Substitution_Analysis.csv` — dải **1,095–1,777** (33 taxon × *E. bulbosa*; × *I. domestica* = 1,471 · × *I. gatesii* = 1,548), tính trên **cả plastome** |

### 8.3. Đo trên máy Hermes (26/09/2026 — Python thuần, không thư viện ngoài)

| Việc | Lệnh / cách làm | Kết quả thô |
|---|---|---|
| Lấy danh sách genome | NCBI E-utilities `esearch` (`db=nuccore&term=Iridaceae[Organism] AND refseq[filter] AND chloroplast`) | 76 plastome RefSeq |
| Tải genome | `efetch` (`id=NC_050833.1&rettype=gbwithparts&retmode=text`) | 302 KB trong **2,3 giây**; gộp 20 genome = 5,9 MB |
| Nhận diện IR — cách viết đơn giản (k-mer + mở rộng khớp chính xác) | so với chú giải `repeat_region` của NCBI | **đúng 7/9, sai 2/9 (22 %) và không ca nào báo lỗi**: NC_063920 đo 22.998 bp (thật 26.162) · NC_063921 đo 7.727 bp (thật 26.183) |
| Nguyên nhân 2 ca sai | đo mức giống nhau của 2 bản IR | NC_063920: khớp tổng 92,63 % (1.929 base khác), lõi hoàn hảo 22.998 bp · NC_063921: khớp tổng 53,65 %, lõi 9.179 bp |
| Hằng số ghép cặp | c = (cuối gen rRNA bản 1) + (đầu gen rRNA bản 2); phép 0-based | **K = c − 2** → IR = **26.214 bp (lệch 0 bp)**; K = c − 1 → khớp 28,56 % (ra 1 bp) |
| Cổng kiểm tra | "hai bản IR giống ≥ 99,9 %" | bắt **2/2** ca sai · **báo oan 0/7** ca lành; ngưỡng 95 % chỉ bắt 1/2 |
| Điểm bắt đầu genome | đối chiếu 20 genome RefSeq | 19/20 bắt đầu ở LSC; 1/20 (NC_085783.1) bắt đầu trong IR. Tài liệu dẫn 27,67 % trên GenBank rộng (CPJSdraw, 8.229 genome) |
| Hiệu năng | chạy trên 20 genome | **4,7 giây**, RAM 331 MB (máy WSL 2,8 GB) |
| Chuỗi end-to-end bản Python thuần | `chain_e2e.py` | 8 mốc, **9,8 giây**: 153.724 bp · GC 37,85 % · IR 26.214 / 100,00 % · 4 vùng khớp · xoay chiều về 109.341 → IR không đổi · SSR tự viết **36** (1bp:30 · 2bp:4 · 3bp:2) · Ts/Tv difflib = **2,667** · máy này không có công cụ chuẩn |
| Chỗ lệch 1 SSR | so danh sách với cpstools | cụm compound **48.659–48.679 = `(TA)6` + `(A)10`**; bản tự viết cắt mất cái `(A)10` |
| Bảng độ nhạy theo ngưỡng mono | cùng dữ liệu | ngưỡng ≥8 → **133** SSR · ≥10 → **36** · ≥12 → **16** |
| Phân bố SSR theo vùng (cpstools) | `_SSRs_loc_results.txt` | IGS 18 · Intron 8 · Gene 11; motif: A 16 · T 14 · TA 3 · AT 1 · C 1 · CTT 1 · TCT 1 |
| So sánh cách căn | difflib vs MAFFT | **2,667** (304 Ts/114 Tv) vs **2,686** (317 Ts/118 Tv) |

### 8.4. Môi trường dev (01/10/2026)

| Hạng mục | Kết quả đo |
|---|---|
| Máy Windows | RAM **15,7 GB** (trống 6,2 GB) · C: 62,3 GB trống / 363,8 GB · D: 60,6 GB trống / 97,7 GB |
| WSL | Ubuntu 24.04.2 · WSL 2.7.14.0 · nhân 6.18.33.2 · WSLg 1.0.73.2. Trước khi sửa `.wslconfig`: RAM 2,8 GB, 8 core |
| `.wslconfig` | memory=8GB, processors=6, swap=8GB, sparseVhd=true, autoMemoryReclaim=gradual (bản cũ lưu `.wslconfig.bak`) |
| Docker | Có Docker Desktop trên Windows; `docker` gọi từ WSL báo *"The command 'docker' could not be found in this WSL 2 distro"* trước khi bật WSL Integration |
| Công cụ trong WSL (máy này, trước khi cài gì) | `misa.pl`, `cpstools`, `mafft`, `blastn`, `mummer`, `trf` = **đều không có** |
| Python | có sẵn `uv` 0.11.16 + CPython **3.11.15** (`/root/.local/bin/python3.11`) |
| Claude Code | bản Linux **2.1.286** tại `/root/.local/bin/claude` (installer `https://claude.ai/install.sh`); `which claude` trỏ bản Linux |
| Cấu hình Claude Code | khối `env` chuyển từ `C:\Users\ADMIN\.claude\settings.json` sang `/root/.claude/settings.json` (quyền 600): `ANTHROPIC_BASE_URL=https://opencode.ai/zen/go`, `ANTHROPIC_MODEL=deepseek-v4.1-flash`, `CLAUDE_CODE_EFFORT_LEVEL`, `ANTHROPIC_API_KEY` |
| Chạy thử Claude Code | tự gọi `make golden` → báo đúng 6/6 chỉ số của bộ số vàng; tự gọi `make test` → đọc đúng lý do các test bị skip |
| Cảnh báo đã gặp | *"Ignoring 17 permissions … this workspace has not been trusted"* (đã đặt `projects["/root/work/plastid-suite"].hasTrustDialogAccepted = true`) · `[claude-code:unrecognized_model] deepseek-v4.1-flash` (vẫn gọi được; API 2,2 giây/lần) |

### 8.5. Đo lại bằng implementation trong repo (02/10/2026 — P1-03)

Chạy đúng 20 genome Iridaceae của §8.3 qua `core/ir.py` của bộ công cụ:

| Việc | Kết quả đo |
|---|---|
| Detector de novo (port `strict_best` từ `ir_probe.py`) | khớp GenBank **7/9**, sai đúng 2 ca như §8.3: NC_063921 → không ra kết quả (**fail ồn ào**), NC_063920 → 22.998 bp (**im lặng** — lõi khớp chính xác) |
| Cổng IR trên cặp chú giải | bắt **2/2** ca sai · **báo oan 0/7** ca lành — tái lập đúng §8.3 |
| NC_085783.1 (không chú giải, bắt đầu trong IR) | de novo ra **21.148 bp @ 100,00 %**, cổng đạt; chuẩn hoá xoay được về đầu LSC |
| **Phát hiện mới — `detect_v2` nở quá**: `ir_probe_v2.py` (mở rộng chịu mismatch) không phải bản đã kiểm chứng | NC_086581: ra **27.601 bp** trong khi GenBank và `strict_best` là **26.217 bp**. Nguyên nhân: nó nhận từng khối 100 bp theo *tỉ lệ giống của toàn vùng đã mở rộng*, nên vượt ranh giới IR mà không dừng. **Bản đã kiểm chứng là `strict_best` trong `ir_probe.py`** (đường chéo mirror K=20 + run khớp chính xác + cửa sổ hợp lệ 10–35 kb, ≥ 95 %) — đã port đúng bản này |
| **Đo lại NC_063921** | với cách ghép cấu trúc C = a1 + b2 − 1: **98,87 %**, không phải 53,65 % như ghi ở §8.3. Kết luận cổng **không đổi** (98,87 < 99,9 ⇒ cắm cờ) |
| Đối chiếu chéo annotation ↔ de novo | có: NC_063920 ghi nhận lệch **3.163 bp** (26.161 vs 22.998) — không im lặng; golden lệch 0 |

### 8.6. Bộ dữ liệu công khai thay thế — Oldenlandia / Coleus (03/10/2026 — D-04)

**Đây KHÔNG phải mẫu của nhóm** — là bản GenBank công khai tải về để dev (quy tắc: CLAUDE.md →
"Data — real samples vs development stand-ins"). Lệnh tái lập (một lần; sau đó chạy offline):

```
plastid panel pin --groups data/oldenlandia-coleus.yaml      --out .panel/oldenlandia-coleus
plastid panel pin --groups data/oldenlandia-coleus-mito.yaml --out .panel/oldenlandia-coleus-mito
```

Payload `.gb` **không commit** (`.panel/` bị gitignore); `PINNED.json` giữ ngày tải + sha256 từng file.
Bảng dưới **sinh thẳng từ `PINNED.json`** (không gõ tay): `pin()` đọc lại từng file bằng
`core/genbank.py` rồi đo bằng `core/ir.py` — cổng IR chỉ **ghi lại**, không chặn ở bước pin
(seam 4: cổng chỉ có hiệu lực lúc chạy). GC% không nằm ở đây — nó là số của module `genome_stats`
lúc chạy, không phải thuộc tính provenance.

| Accession | Dài (bp) | CDS | tRNA | rRNA | Topology | IR (bp, identity) · cổng |
|---|---|---|---|---|---|---|
| `NC_057985.1` | 153.653 | 85 | 37 | 8 | circular | 26.286 (1,0000) PASS |
| `PP078674.1` | 153.676 | 86 | 37 | 8 | circular | 26.323 (1,0000) PASS |
| `PQ325515.1` | 153.640 | **0** | **0** | **0** | circular | 26.323 (1,0000) PASS |
| `NC_057983.1` | 152.327 | 85 | 37 | 8 | circular | 25.286 (1,0000) PASS |
| `PX977793.1` | 152.414 | 85 | 37 | 8 | circular | 25.308 (1,0000) PASS |
| `PX923121.1` | 152.544 | 88 | 37 | 8 | circular | 25.711 (1,0000) PASS |
| `PZ428230.1` | 152.510 | **0** | **0** | **0** | **linear** | **15.756** (1,0000) PASS |
| `NC_069583.1` | 152.484 | 87 | 36 | 8 | circular | 25.709 (1,0000) PASS |
| `AP045470.1` | 152.481 | **0** | **0** | **0** | circular | 25.709 (1,0000) PASS |
| `NC_058328.1` | 152.190 | 88 | 37 | 8 | circular | 25.792 (1,0000) PASS |
| `NC_057679.1` | 152.603 | 85 | 37 | 8 | circular | 25.716 (1,0000) PASS |
| `OX459128.1` | 258.274 | 31 | 17 | 3 | circular | IRDetectionError |
| `AP045469.1` | 313.585 | **0** | **0** | **0** | linear | IRDetectionError |

**sha256 ổn định:** pin lại cùng ngày 03/10 — cả **13/13 sha256 không đổi**, tức NCBI trả đúng
từng byte như lần tải đầu (bằng chứng bộ dữ liệu tái lập được).

Số đo mới, phải nhớ:

- **5 bản GenBank bị loại vì trùng khít RefSeq** (so chuỗi trực tiếp, không suy từ metadata):
  `MT767006.1`≡`NC_057983.1` · `MT767008.1`≡`NC_057985.1` · `MT473748.1`≡`NC_058328.1` ·
  `MW538954.1`≡`NC_057679.1` · `OP611428.1`≡`NC_069583.1`.
- **3/11 plastome tải về không có annotation** (`PQ325515.1`, `PZ428230.1`, `AP045470.1` — 0 CDS):
  dùng được cho phân tích theo vị trí (SNP · π · SSR), **không** dùng được cho phân tích theo gen.
  Bài học: `entry_from_summary` cũ đặt `annotated=True` vô điều kiện — đã sửa thành chỉ tin RefSeq;
  annotation của bản GenBank **đo sau khi tải**, không giả định.
- **`PZ428230.1` — ĐÃ KẾT LUẬN: lỗi ở dữ liệu, không phải ở detector.** IR de novo ra **15.756 bp**
  trong khi cả 5 loài Coleus còn lại đều ~25,7 kb, và bản ghi khai báo **`linear`**. Đo để phân biệt
  "detector cắt cụt IR thật" với "bản ghi hỏng": mở rộng vùng IR ra mỗi bên rồi đo lại identity
  — với `PZ428230.1`, identity **rơi ngay lập tức** (1,0000 → 0,9556 ở +500 bp → 0,9187 ở +1 kb),
  tức **không hề có** cặp repeat ~25,7 kb nào bị bỏ sót; bản ghi này đúng là thiếu ~10 kb IR.
  (Bản cùng loài `PX923121.1` cho cùng dạng đường rơi sau 25.711 bp — tức detector hành xử nhất quán.)
  **Kết luận: loại `PZ428230.1` khỏi mọi phân tích theo vị trí** (và nó vốn đã trống annotation).
  **Bài học:** khi bản ghi không có chú giải thì không có gì đối chiếu IR — cổng chỉ kiểm *identity*,
  không kiểm *độ dài*, nên một repeat 15,7 kb vẫn PASS. Xem §9.A.
- **IR của 12/13 plastome là de novo**, chỉ `NC_058328.1` có annotation để đối chiếu — nghĩa là cột IR
  ở bảng trên chưa có cross-check từ chú giải cho phần lớn bản ghi (đúng loại chỗ §9A cảnh báo).
- Đồng nghĩa taxonomy: *Oldenlandia diffusa* hiện mang tên **Scleromitrion diffusum** trong NCBI —
  cột "genus" của panel theo lineage sẽ là `Scleromitrion`, không phải `Oldenlandia`.
- Mito: cả hai đều `IRDetectionError: no valid IR pair found (annotation absent; de-novo scan failed)`
  — **đúng**, mito không có IR canonical; module đòi cổng IR sẽ không chạy được trên mito cho tới khi
  có đường riêng (Phase 3+).

### 8.7. Bảng gen xuất bản + xuất Word (03/10/2026 — P2-09, D2-04)

Chạy thật: `plastid run -m gene_table -i tests/data/NC_050833.1.gb --docx --dataset golden`.

**Bảng phân nhóm, đo trên genome vàng** (132 bản sao, 18 nhóm, 112 tên phân biệt):

| Nhóm | Bản sao | Ghi chú |
|---|---|---|
| Subunits of photosystem I / II | 5 · 15 | |
| Subunits of cytochrome b/f · ATP synthase | 6 · 6 | |
| Large subunit of rubisco | 1 | |
| Subunits of NADH dehydrogenase | 12 | gồm **ndhL** (mảnh) |
| Ribosomal RNA · Transfer RNA | 8 · 38 | khớp theo tiền tố `rrn` / `trn` |
| Small · Large subunit of ribosome | 15 · 11 | |
| DNA-dependent RNA polymerase | 4 | |
| Other genes (6 nhóm) | 6 | infA · matK · cemA · accD · ccsA · clpP |
| Genes of unknown function | 5 | ycf1, ycf2 (×2), ycf3, ycf4 |
| **Total** | **132** | = 86 CDS + 38 tRNA + 8 rRNA, khớp `count_genes()` |

Đối chiếu trong `gene_reconciliation`: CDS 86/86 · tRNA 38/38 · rRNA 8/8 · Total 132/132 — **`agrees` cả 4**.
Gen không nhóm nào nhận: **0**.

**Tài liệu Word — đo được, không phải cảm nhận:**

- **Tái lập từng byte**: hai lần dựng cách nhau 2,1 giây cho **cùng sha256** (mốc thời gian zip là
  2 giây/lần, nên khoảng cách đó thật sự vượt qua ngưỡng); chạy lại qua **cache** cũng ra **đúng byte đó**
  (`e4a65ddf…` cho gene_table). Test canh: mọi entry mang `date_time = 1980-01-01`, thứ tự sắp xếp,
  `create_system = 0`.
- **Kích thước gói**: ~800 kB (template thô) → **13 kB**; style giữ lại **< 25** trên **164** của template.
- **Render thật** bằng LibreOffice (`soffice --headless`) → PDF → PNG và **xem tận mắt**: tiêu đề theo
  style `Title`, mục theo `Heading 1`, lưới provenance (Dataset · Run · Input + sha256 · Coordinates ·
  Module · Generated by), bảng có caption ở trên + dòng scope, header tô nền **lặp lại ở mỗi trang**,
  hàng không bị cắt ngang trang, cột số căn phải, dòng **Total in đậm**, chân trang
  `PlastidToolSuite · <accession> · <module> · Page X of Y`.
- **Bảng 16 cột** (`genome_stats`) tự chuyển sang **khổ ngang** theo **bề rộng nội dung ước lượng**, không
  theo số cột — bản đầu chọn theo số cột nên đẩy một bảng 7 cột, 1 hàng sang trang ngang gần như trắng.

### 8.8. Bản nháp bản thảo trên bộ dữ liệu thay thế (03/10/2026 — P2-10)

Công thức ở `docs/P2-10_manuscript_draft.md`. Chạy thật **một lệnh cho mỗi loài**, 7 module:

```
plastid run -m genome_stats -m structure -m gene_table -m codon_usage -m amino_acid \
            -m introns -m ssr -i <genome>.gb --docx --dataset <nhãn>
```

**Phát hiện quan trọng — bộ "thay thế" chứa đúng hai loài bản thảo:**

- `PP078674.1` = ***Scleromitrion diffusum*** — **tên hiện hành của *Oldenlandia diffusa***
- `PX923121.1` = ***Coleus amboinicus*** — **đúng loài thứ hai**

Tức hai bản ghi này **không phải loài họ hàng** mà là bản công khai **của chính hai loài**.

⚠️ **ĐÍNH CHÍNH (xem §8.10):** kết luận trên còn thiếu một nửa. `PX923121.1` **là bản ghi do nhóm
mình nộp** (Vu,H.-T. = cô Trang), nên với *C. amboinicus* nhãn đúng là
**`own-record-coleus-amboinicus`**, không phải `standin-…`. `PP078674.1` (*O. diffusa*) là của
Shen,C. nên nhãn `standin-oldenlandia-diffusa` **vẫn đúng**. Bản chạy lại theo nhãn mới nằm ở
`.runs/standin/own-record-coleus-amboinicus/`; các số trong bảng dưới không đổi (cùng một bản ghi).

**Nhãn từng loài (đính chính §8.10):** *C. amboinicus* = **`own-record-coleus-amboinicus`** (bản ghi
của nhóm) · *O. diffusa* = **`standin-oldenlandia-diffusa`** (bản ghi của Shen,C.).

**Số đo được (mỗi số là một ô trong bảng của lần chạy, không gõ tay):**

| | *S. diffusum* `PP078674.1` | *C. amboinicus* `PX923121.1` |
|---|---|---|
| Chiều dài · GC | 153.676 bp · 37,61 % | 152.544 bp · 37,75 % |
| Gen | **131 = 86 CDS + 37 tRNA + 8 rRNA** | **133 = 88 CDS + 37 tRNA + 8 rRNA** |
| LSC · SSC · IR | 83.835 · 17.195 · 26.323 | 83.404 · 17.718 · 25.711 |
| Cổng IR | PASS 100,00 % (de novo) | PASS 100,00 % (de novo) |
| Bảng gen ↔ `count_genes()` | 86/86 · 37/37 · 8/8 · 131/131 — `agrees` | 88/88 · 37/37 · 8/8 · 133/133 — `agrees` |
| Codon | 26.085 codon · ENC 50,65 · CAI 0,5603 · CBI 0,2645 · Fop 0,5067 | 26.721 codon · ENC 49,80 · CAI 0,5343 · CBI 0,2794 · Fop 0,5170 |
| Amino acid | GRAVY −0,0856 · aromaticity 0,1116 | GRAVY −0,0962 · aromaticity 0,1105 |
| Intron | 16 intron / 11 gen · 12.287 bp · dài nhất 1.142 | 17 intron / 12 gen · 12.134 bp · dài nhất 1.022 |
| Trans-splicing | 2 điểm nối, đều rps12 — `agrees` với chú giải | 2 điểm nối, đều rps12 — `agrees` |
| SSR | 32 = mono 32 · di 0 · tri 0 (3 compound) | 32 = mono 31 · di 1 · tri 0 (2 compound) |
| Internal stop | 0 | 0 |

**Điểm đáng chú ý:** trans-splicing được nhận đúng trên **hai genome mới**, không chỉ genome vàng —
`annotation_check` ghi `agrees` cả hai, tức logic tách intron/trans-splicing của P2-05 không phải
may mắn trên một ca. **Số gen khác 132** (131 và 133) là **đúng**: 132 là con số của genome vàng
*Iris domestica*, không phải hằng số cho mọi loài — quy ước đếm mới là hằng số.

**Chưa có trong bản nháp:** mục π và Ts/Tv (cần MAFFT 7.526 — máy này không có), và **hình** (Giai
đoạn 5). Bản nháp ghi rõ phần thiếu là phần vắng, không phải phần ngắn im lặng.

### 8.9. Bốn tên gen lạ ở *C. amboinicus* — cùng locus, khác tên (03/10/2026 — D3-02)

`PX923121.1` (*C. amboinicus*) và `PX977793.1` (*O. corymbosa*) có **bốn tên gen không khớp** bộ
nhóm chức năng: `pafI` · `pafII` · `pbf1` · `clpP1`. Tra theo **vị trí** (không tra theo tên) trên
9 bản ghi có chú giải:

| Tên trong `PX923121.1` | Nằm giữa | Ở 7 bản ghi còn lại, đúng chỗ đó là | Độ dài |
|---|---|---|---|
| `pafI` | `psaA` … `trnS-GGA` | **`ycf3`** | 507 bp (510 bp ở cả hai bên) |
| `pafII` | `psaI` … `cemA` | **`ycf4`** | 531 bp (555 bp ở cả hai bên) |
| `pbf1` | `psbT` … `psbH` | **`psbN`** | **132 bp — khớp tuyệt đối** |
| `clpP1` | `rps12` … `psbB` | **`clpP`** | 591 bp |

Ba lớp bằng chứng độc lập: **synteny** (hàng xóm hai bên giống hệt, 7 bản ghi) · **độ dài** (khớp
từng cặp, khớp cả theo chi: 507 bp ở Oldenlandia / 510 bp ở Coleus) · **`/product`** ("photosystem I
assembly factor I" ↔ "photosystem I assembly protein Ycf3"). Ngoài ra `pbf1` = **132 bp** trùng khít
`orf43/psbN` (43 codon + stop).

**Vì sao tên đổi:** `ycf` = *hypothetical chloroplast open reading frame* — tên đặt khi chưa biết gen
làm gì. Đã biết rồi thì đổi sang tên chức năng: Ycf3/Ycf4 là yếu tố lắp ráp PSI (PafI/PafII); `psbN`
→ `pbf1` do Krech và cs. 2013 đề xuất sau khi chứng minh sản phẩm gen tham gia sinh tổng hợp **cả**
PSI lẫn PSII. Một bảng gen công bố ghi thẳng `pafI(ycf3)`, `pafII(ycf4)`.

**Quyết định (dev, có bằng chứng):** tên chuẩn giữ **tên cổ điển** — `ycf3` · `ycf4` · `psbN` ·
`clpP` — vì **7/9 bản ghi và genome vàng** dùng chúng, mà genome vàng là hợp đồng của dự án.
Tên gốc **không bị bỏ**: đi kèm là `source_name`, truy ngược được về bản ghi. Lý do phải ánh xạ:
so sánh chéo genome khoá theo **tên gen**, để nguyên thì `pafI` và `ycf3` thành hai gen khác nhau và
bảng π theo gen / bảng gen xuất bản **tách đôi một gen** — đúng loại "sai mà không báo".

**Còn lại của PI:** tên **in trong bản thảo** và tên **nộp NCBI** (D3-06) — hai việc đó mới đổi nội
dung công bố. Ghi chú: hai bản ghi dùng tên mới là **Direct Submission** nộp 24/01/2026 và
07/02/2026 bởi cùng một nhóm (`Vu,H.-T.` và cs.) — tức đây là **quy ước của một nhóm**, không phải
chuẩn cộng đồng; nhóm mình cần xác nhận có phải bản ghi của mình không.

**Bài học phương pháp:** lần đầu tôi xếp việc này vào cột PI vì *chưa có bằng chứng*. Đo xong thì nó
thành quyết định của dev. Escalate khi thiếu bằng chứng là đúng; **giữ nguyên** mức escalate sau khi
đã có bằng chứng là sai — ở đây nó chỉ nằm trong board một lúc, chưa gây thiệt hại, nhưng cơ chế thì
sai: một việc dev làm được mà bị xếp cho PI sẽ **chặn im lặng** cho tới khi có người hỏi.

### 8.10. Hai bản ghi là CỦA NHÓM — và cả một phương ngữ đặt tên (03/10/2026)

**Xác nhận của anh Trí:** `Vu,H.-T.` = **Vũ Thị Huyền Trang** = **cô Trang, PI**. (Tên này cũng chính
là tên file kế hoạch gốc ở §11 — lẽ ra phải nhận ra sớm hơn.)

**Ai nộp bản ghi nào** (đo từ khối REFERENCE của cả 11 plastome):

| Bản ghi | Loài | Người nộp |
|---|---|---|
| `PX923121.1` | ***Coleus amboinicus*** | **Vu,H.-T. and Nguyen,T.-V.H.** — 24/01/2026 |
| `PX977793.1` | *Oldenlandia corymbosa* | **Vu,H.-T. and Duong,T.-H.** — 07/02/2026 |
| `PP078674.1` | *Scleromitrion diffusum* | Shen,C. |
| `PQ325515.1` | *Scleromitrion diffusum* | huang,z. |
| `AP045470.1` · `PZ428230.1` · `NC_057679.1` · `NC_058328.1` · `NC_069583.1` | — | Furuta,Y. và cs. · Rahmat,E. và cs. · NCBI Genome Project · Zhao,F. và cs. |

**Hệ quả 1 — nhãn `standin-coleus-amboinicus` là SAI.** `PX923121.1` **đúng loài bản thảo *C. amboinicus***
và là **bản ghi do nhóm mình nộp**. Bản nháp P2-10 cho loài này dựng trên **dữ liệu của chính nhóm**, không
phải dữ liệu thay thế. Còn `PP078674.1` (*O. diffusa* = *S. diffusum*) là của **Shen,C.** — bản đó **đúng
là** dữ liệu thay thế. Ghi chú: mẫu đã công bố **không chắc** là mẫu bản thảo sẽ dùng; vẫn phải hỏi anh Điềm.

**Hệ quả 2 — tên gen trong hai bản ghi đó là do nhóm đặt, qua GeSeq.** Khối COMMENT ghi
`Assembly Method :: GetOrganelle v.1.7+` và bản ghi thiếu `##Annotation-Data`; tài liệu cho thấy GeSeq
chú giải theo tương đồng với bộ tham chiếu dùng **tên chức năng mới** ⇒ đầu ra là `pafI`·`pafII`·`pbf1`.
Nghĩa là **Giai đoạn 3 chạy GeSeq sẽ tự sinh ra đúng bộ tên đó** — việc chuẩn hoá tên không phải việc
phụ, mà là việc bắt buộc.

**Hệ quả 3 — khác biệt không phải 4 gen, mà là CẢ MỘT PHƯƠNG NGỮ.** Đo trên 11 bản ghi:

| Khía cạnh | Dạng | Số bản ghi |
|---|---|---|
| tRNA anticodon | hyphen `trnA-UGC` | **8** |
| | ngoặc `trnA(UGC)` | 1 (`NC_058328.1`) |
| | **trần, không anticodon** `trnK` · `trnI` | 2 (`PX923121.1` · `NC_069583.1`) |
| rRNA | trơn `rrn16` | 6 |
| | có S `rrn16S` | 3 |
| protein | tên cổ điển `ycf3`·`ycf4`·`psbN`·`clpP` | 7 + genome vàng |
| | tên chức năng `pafI`·`pafII`·`pbf1`·`clpP1` | 2 (của nhóm) |

**Tên trần giải được bằng TRÌNH TỰ, không bằng phỏng đoán** — so với các tRNA đã có tên ở bản ghi khác:

- `trnK` ở `PX923121.1` ≡ **`trnK-UUU`** — **100,0 %** trùng với `NC_069583.1` và `NC_057679.1`
- `trnI` ở `NC_069583.1` ≡ **`trnI-GAU`** — **100,0 %** trùng với `PX977793.1`

*(Thử đọc anticodon theo offset cố định đã KHÔNG hội tụ: offset 34 chỉ khớp 106/293 — anticodon không
nằm ở một vị trí cố định trong bộ này. So trình tự là cách đúng.)*

**Quyết định (dev, có bằng chứng):**
- **Khoá so sánh nội bộ = phương ngữ của nhóm** (`pafI`·`pafII`·`pbf1`·`clpP1`), tên cổ điển đi kèm ở
  `source_name`. Lý do không phải "đa số": bản ghi *C. amboinicus* **là của nhóm**, là **loài bản thảo**,
  và GeSeq — công cụ nhóm dùng — sẽ sinh ra đúng bộ tên đó. Chuẩn hoá về tên cổ điển sẽ làm **bài báo
  của nhóm in tên khác với bản ghi của chính nhóm**.
- **In ra thì in CẢ HAI**: `pafI (ycf3)`. Đây là cách các bảng gen công bố đang làm, và nó gỡ được câu
  hỏi "tên nào" cho người đọc — mỗi tên vẫn truy được về bản ghi tương ứng.
- **`gene_groups.yaml` phải nhận thêm bốn tên này** vào đúng nhóm chức năng (pafI/pafII → photosystem I;
  pbf1 → photosystem II; clpP1 → protease). Hiện chúng rơi vào *unclassified* — đúng như thiết kế
  fail-loud, và chính nó đã chỉ ra vấn đề này.
- **D3-01 (chọn công cụ) giờ có bằng chứng**: quy trình thật của nhóm là **GetOrganelle → GeSeq**.

### 8.11. Kiểm tra ORF — P3-01 (03/10/2026)

**Đã làm:** `core/orf.py` (phân loại) + `modules/orf.py` (báo cáo) + `Gene` nay mang `/transl_table`
và `/codon_start` của chính bản ghi. Không có CDS nào bị viết lại để "cho qua".

**Đo trên genome vàng (NC_050833.1) — đúng hợp đồng P3-01:**

| Chỉ số | Số đo |
|---|---|
| CDS kiểm tra | **86** (khớp `count_genes().cds`, có bảng đối chiếu riêng) |
| `ok` | **81** |
| `alternative_start` | **5** = rps19 **GTG** ×2 · rpl2 **ACG** ×2 · ndhD **ATC** |
| `no_start` · `no_stop` · `internal_stop` · `not_multiple_of_three` | **0 · 0 · 0 · 0** |
| `gene` feature không có CDS nào cùng tên phủ | **1** — ycf1 cụt 108.448–109.347 (`pseudo`) |

Mọi CDS khai `/transl_table=11` và `/codon_start=1`; module đọc đúng khai báo đó (D3-03).

**Hai phát hiện mới, mỗi cái thành một luật:**

1. **ACG không có trong danh sách mã khởi đầu của bảng 11 (NCBI)** — nhưng bản ghi đã công bố dùng nó:
   đo trên panel được **4 ACG** (đều ở `ndhD`) cộng **2** ở genome vàng (`rpl2`); tài liệu cho thấy
   ACG→AUG do RNA editing ở thực vật. Nếu xếp ACG vào `no_start` thì đúng hai CDS của genome vàng bị
   báo động giả. Quyết định: ACG **được nhận là start** nhưng mang lớp riêng `non_canonical`
   (khác `alternative` là codon có trong bảng) — vừa không báo động giả, vừa không im lặng.
2. **ycf1 cụt chồng lấn 21 bp với CDS của `ndhF`** (109.327–109.347). Nên lớp "annotated, not coding"
   phải khớp **cả tên lẫn vị trí**: khớp theo vị trí không thôi thì `ndhF` "phủ" nó và nó biến mất;
   khớp theo tên không thôi thì CDS ycf1 thật ở chỗ khác (122.984–128.404) sẽ "phủ" nó. Cả hai cách
   đều là im lặng — mỗi cách hỏng một kiểu.

**Từ chối thay vì đoán:** bảng mã không nằm trong {1, 11} (hai bảng này chia sẻ toàn bộ axit amin,
chỉ khác mã khởi đầu — D3-03) → `OrfError`; `/codon_start` ≠ 1 → `OrfError`. Dịch sai khung là kiểu
"sai mà không báo", thà dừng.

**Kiểm chứng:** `make check` xanh — **229 fast + 30 golden**; chạy thật `plastid run -m orf` trên
genome vàng sinh đủ 5 bảng + report, 0 CDS lỗi.

### 8.12. Chuẩn hoá tên gen — P3-02 (03/10/2026)

**Đã làm:** `core/gene_names.py` (**một** chỗ duy nhất đổi tên gen) · `data/gene_aliases.yaml`
(4 locus ngữ nghĩa, mỗi entry mang bằng chứng) · `modules/gene_names.py` (báo cáo trên nhiều bản
ghi) · `gene_groups.yaml` nhận 4 tên mới · normalizer được nối vào `gene_table`, `genome_stats`,
`orf`.

**Ba loại khác biệt, xử lý và báo cáo riêng:**

| Loại | Ví dụ | Nguồn tri thức |
|---|---|---|
| cú pháp | `trnL(UAA)` → `trnL-UAA` · `rrn16S` → `rrn16` | quy tắc, không cần dữ liệu ngoài |
| ngữ nghĩa | `ycf3` → **`pafI`** (in ra `pafI (ycf3)`) | `gene_aliases.yaml` + D3-02 |
| trình tự | `trnK` → `trnK-UUU` | so 100 % với bản sao có tên ở bản ghi khác |

**Đo trên panel 11 bản ghi (chạy thật qua CLI):**

| Chỉ số | Số đo |
|---|---|
| Gen được đếm | 1.048 |
| Tên khác nhau *như bản ghi viết* | **155** |
| Tên chuẩn khác nhau sau chuẩn hoá | **115** |
| Viết khác cú pháp | 34 (30 tên ngoặc ở `NC_058328.1` + `rrn*S`) |
| Một locus hai tên | 8 (4 locus × 2 chi) |
| Giải bằng trình tự | **2** — `trnK` → `trnK-UUU`, `trnI` → `trnI-GAU` (đúng như §8.10 đo) |
| Không giải được | **0** |
| *Unclassified* (gene_groups) | **0** trên cả 11 bản ghi |

Genome vàng: 112 tên đến **nguyên vẹn** (normalizer không viết lại tên nào) và cho **112** tên
chuẩn — không mất, không gộp. Panel ∪ genome vàng = **156** tên (khớp board).

**Ba việc sửa được trong lúc làm:**

1. **Feature `gene` không mang trình tự.** Bản ghi `PX923121.1` có cả feature `gene` (trùm cả
   intron) lẫn feature `tRNA` cùng tên `trnK`. Giải tên theo *từng feature* thì so intron với
   tRNA → báo unresolved oan. Sửa: chỉ giải từ feature `tRNA`, rồi áp tên đã giải cho mọi feature
   cùng tên. (Test panel bắt được, không phải review.)
2. **CLI chặn panel.** `cli.py` so `len(genomes) != arity` — tài liệu API v0.2 ghi arity ≥ 2 là
   "**at least**", runner cũng đã làm đúng, chỉ CLI sai. Hệ quả thật: module `arity=2` không chạy
   được trên panel 11 bản ghi. Đã sửa + test hồi quy.
3. **Dạng in ra:** mọi output in **`tên chuẩn (tên đồng nghĩa)`** — `pafI (ycf3)` — kể cả khi bản
   ghi vốn đã dùng tên chuẩn, để người đọc không phân biệt được bản ghi nào dùng phương ngữ nào.
   Tên gốc vẫn nằm trong `name_mappings` (cột `name`) và trong chính bản ghi.

**Ràng buộc đã biết (chưa xử lý, cần team biết):** test quét template cấm **mọi** chữ số ngoài
`{{…}}`, nên template không thể viết tên gen có số (`ycf3`, `rps19`) hay mã quyết định (`D3-02`).
Hai template mới phải né bằng cách để dữ liệu trong bảng. Guard này vẫn giữ nguyên (không nới);
nếu sau này prose cần tên gen, cách đúng là tham chiếu một ô bảng.

**Kiểm chứng:** `make check` xanh — **255 fast + 30 golden**; chạy thật `plastid run -m gene_names`
trên 11 bản ghi panel.

### 8.13. Chú giải & file nộp — P3-03 · P3-04 · P3-05 (03/10/2026)

**P3-03 · PGA (adapter + module).** `external/pga.py` chạy PGA qua `run_tool()` (argv lấy từ
README chính thức: `perl PGA.pl -r <refdir> -t <targetdir> -o <outdir>`), trong thư mục tạm mới
mỗi lần, tên file đã slugify, **một bản ghi mỗi lần chạy** — ba thói quen xấu của PGA (xoá thư mục
ra, vỡ ở tên có dấu cách, multi-FASTA) đã bị chặn từ trước. Đầu ra GenBank được **đọc lại bằng
reader của mình** rồi đi qua đúng bộ kiểm tra ORF + chuẩn hoá tên như một bản ghi tải về.

*Đo được ở máy này:* thiếu `PGA.pl` ⇒ `ToolError` có tên + exit 3. Đường đi từ "annotator ghi ra
GenBank" tới "bộ gen đi qua kiểm tra" chạy với một **stub** thay PGA (test ghi rõ điều đó): ra
**132 = 86 + 38 + 8**, 0 CDS lỗi ORF, 5 alternative start, 4 tên có đồng nghĩa. **Đây là kiểm tra
dây nối của mình, không phải bằng chứng PGA chú giải đúng** — stub chép lại chính chú giải của
genome vàng. Test thật (tier B) skip kèm lý do; `PTS_REQUIRE_TOOLS=1` sẽ biến skip thành fail.

**P3-04 · mito — quyết định đã chốt, adapter chưa viết.** Chọn **PMGA** (không phải MITOS2). PMGA
không nằm trên conda, không có repo git: bài báo phân phối nó dưới dạng **singularity container**
(Figshare 27201798) và web server, và **không có tài liệu CLI công khai** nào tra được. Viết
`external/pmga.py` theo một CLI tự nghĩ ra sẽ là placeholder trông như thật ⇒ dừng lại, ghi yêu cầu
cụ thể vào TASKS.md. Hai số đo để phiên sau khỏi đo lại: bản ghi mito thay thế `OX459128.1`
**không có qualifier `/gene`** (mọi "tên" trong model là chuỗi `/product`, ví dụ "NADH
dehydrogenase subunit 1") ⇒ muốn báo "tên gen mito" thì phải có từ vựng product→symbol trước; và
CDS của nó có start `CCT`, `CGC` — không bảng mã nào nhận, đúng như bộ kiểm tra ORF sẽ báo.

**P3-05 · file nộp NCBI.** `.fsa` + `.tbl` viết thuần Python; `.sqn` là tuỳ chọn (`--set
sqn=true`) và cần cả `table2asn` **lẫn** template `.sbt` — template chứa tên người và cơ quan, bộ
công cụ không tự bịa. Thiếu công cụ ⇒ `ToolError` tên rõ ràng + exit 3 (`table2asn` không có ở
máy này; bản tải dở 12 MB nằm ở `~/toolbox/bin/table2asn.gz`).

*Quy tắc định dạng lấy từ tài liệu NCBI, không lấy từ trí nhớ* — mỗi lỗi trong số này chỉ lộ ra ở
bước upload, sau khi đã làm xong hết:

| Quy tắc | Ghi chú |
|---|---|
| 5 cột, tab: start · stop · feature key · **qualifier name (cột 4)** · value (cột 5) | định nghĩa cột theo tài liệu NCBI |
| `join()` **không tồn tại** trong định dạng này | interval **xếp tầng**, key chỉ ở dòng đầu |
| mạch âm | interval **đảo** (toạ độ lớn trước), từng interval một — nên gen trans-spliced khác mạch vẫn viết được |
| `<` luôn ở cột 1, `>` luôn ở cột 2, **bất kể mạch** | và cả hai đều là **tiền tố**: `>1050`, không phải `1050>` |
| partial ở interval **giữa** | NCBI từ chối ("internal partials") ⇒ bộ ghi **raise**, không bỏ dấu |

*Đo được:* vòng tròn khép kín trên genome vàng — **265/265** feature (133 `gene` + 86 CDS + 38 tRNA
+ 8 rRNA) đọc lại đúng kind, tên, interval (kể cả partial), product và bảng mã; 0 vấn đề first-pass.
Chạy thật qua CLI: `submission.fsa` (156.300 B) + `submission.tbl` (15.503 B) nằm trong thư mục
run, có trong `manifest.json`, và lần chạy thứ hai (cache hit) phục hồi **byte-identical**.

*Đối chiếu với bảng THẬT (bổ sung cuối phiên).* Ba file `.tbl` do **chính pipeline của NCBI** sinh ra
(đầu ra VADR, vendor ở `tests/data/real_*.tbl` kèm nguồn + sha256 + giấy phép public domain) đã được
đọc lại và ghi lại **byte-identical từng dòng một** — tức là mọi quy ước "chết ở bước upload" đều đã
đối chiếu với bảng thật, không chỉ với tài liệu. Đọc chúng làm **sửa writer hai lần**: (1) dòng
interval chỉ có **3 trường**, không đệm 2 cột rỗng như bản đầu; (2) một file `.tbl` **có thể chứa
nhiều bản ghi** — một file VADR chứa **8 bản ghi**, và parser cũ **gộp im lặng** tất cả feature của 8
bản ghi vào một danh sách, tức là gán feature của bản ghi này cho sequence của bản ghi khác. Nay
`parse_tbl` trả về `TblRecord` (seq_id + features) và có test riêng cho ca 8 bản ghi.
*Hai khoảng trống còn lại:* (1) quy ước **join/xếp tầng** vẫn chỉ đối chiếu với ví dụ trong tài liệu
NCBI — không tìm được file thật nào có feature nhiều đoạn; (2) nửa `.sqn` chưa chạy lần nào (thiếu
template của nhóm).

**Cập nhật cuối phiên — `table2asn` đã cài và đã khai báo.** Tải trọn 20.213.344 B (nén) →
56.847.504 B, `table2asn -h` → **1.29.324**; md5 nhị phân `8b5b6cbede9d38b47c7c4cae29a25a0a`.
Khai báo theo invariant trong `docs/tools/table2asn.md` (nguồn + dòng tải nguyên văn + lý do không
vendored). Đo được hai lỗi có tên, đúng như thiết kế: thiếu công cụ ⇒ `ToolError` chỉ luôn
`--set table2asn=`; có công cụ nhưng thiếu template ⇒ `ToolError` chỉ trang sinh `.sbt` của NCBI —
**exit 3**, manifest ghi `status: failed` kèm thông báo, run vẫn truy lại được.
Thử hai cách viết `.sbt` bằng tay: cả hai bị `Error loading template file` — nên **không đoán**
định dạng đó nữa; `.sbt` là file của nhóm, sinh từ form NCBI. Nửa `.sqn` vì vậy vẫn chờ đúng một
thứ: template.

**Một ràng buộc mới phát hiện:** test quét template cấm **mọi** chữ số ngoài `{{…}}`, nên template
không thể viết tên gen có số (`ycf3`, `rps19`), mã quyết định (`D3-02`), hay `ASN.1`. Guard vẫn
giữ nguyên (không nới); hai template mới né bằng cách để dữ liệu trong bảng. Nếu sau này prose cần
tên gen, cách đúng là tham chiếu một ô bảng — và nếu team muốn nới guard thì phải là một quyết định
riêng, có test âm đi kèm.

### 8.14. Đường mito — P3-04 phần chạy được (03/10/2026)

**Đã làm:** `modules/mito_genes.py` + `data/mito_gene_names.yaml` (31 locus) + `load_mito_aliases()`.
Module **không bao giờ** gọi `ctx.require_ir()` — có test chạy nó với `ctx.ir=None`, và một test
đối chứng cho thấy module có gate vẫn dừng trên đúng context đó (gate vẫn thật, chỉ là nó nằm trước
phân tích theo vị trí, không nằm trước việc đếm gen).

**Đo trên mitogenome thay thế `OX459128.1`:** **51 = 31 CDS + 17 tRNA + 3 rRNA**; **31/31 CDS được
gọi bằng ký hiệu chuẩn** (`nad1`·`nad2`·`cox1`·`atp6`·`ccmFn`·`matR`·`sdh4` …); **0 tên không map
được**; run không đụng gate.

**Vì sao phải có từ vựng riêng:** bản ghi này **không có qualifier `/gene`** trên bất kỳ feature
nào, nên "tên" mà reader trả về thực chất là **câu `/product`** ("NADH dehydrogenase subunit 1").
Và **vì sao file phải tách khỏi `gene_aliases.yaml`**: đo được chuỗi product **trùng giữa hai bào
quan** — "ribosomal protein S3" là product của cả CDS plastid (genome vàng) lẫn CDS mito. Dùng
chung một bảng thì một gen plastid có thể bị đổi sang ký hiệu mito. Có test canh đúng ca này.

**Chất lượng bằng chứng, nói thẳng:** chuỗi synonym thì **đo được** (trích nguyên văn từ
`OX459128.1`); còn ký hiệu là **danh pháp chuẩn của mito thực vật** và **chưa** đối chiếu với một
chú giải tham chiếu đã công bố. Đây là file có bằng chứng yếu nhất trong repo, và chính nó nói ra
điều đó; việc đối chiếu thuộc Giai đoạn 8 (kiểm thử chéo). Tên nào không map được thì **báo nguyên
dạng**, không đoán — có bảng `unmapped_names` riêng.

**Còn lại của P3-04:** adapter PMGA. Không viết, vì PMGA chỉ phát hành dạng singularity container
(Figshare 27201798) và không có tài liệu CLI công khai — viết theo CLI tự nghĩ ra là placeholder.

### 8.15. PGA chạy thật trên genome vàng — và ba cái bẫy (03/10/2026)

**Bối cảnh:** mục "chờ anh Trí chạy PGA một lần" trên board hoá ra **không phải chờ người**: cùng
một cái máy, việc bị chặn là *quyền* của agent, không phải khả năng. Anh Trí cấp quyền trong
`.claude/settings.local.json` (file này đã vào `.gitignore`), và PGA chạy được ngay sau đó.

**Đo được** (genome vàng đã gỡ chú giải làm target, reference là các plastome chi *Iris*):

| Cấu hình | Kết quả | So với chuẩn 132 = 86 + 38 + 8 |
|---|---|---|
| 7 reference (6 RefSeq + *I. gatesii*) | **147 = 95 CDS + 44 tRNA + 8 rRNA** | thừa 15 gen |
| 1 reference (*I. gatesii*) | **130 = 87 CDS + 35 tRNA + 8 rRNA** | thiếu 2, lệch 4 |

→ **Done-when của P3-03 ("tái tạo đúng số gen của genome vàng") KHÔNG đạt** với PGA. Đây là số đo
về annotator, không phải lỗi module: module chạy đúng, đọc lại output, và báo cáo trung thực.

**Chênh lệch đến từ đâu (đo được, không suy đoán):** 4 CDS trượt kiểm tra ORF (2 `no_start`,
1 `no_stop`, 1 `internal_stop`); **4 locus bị gọi hai lần dưới hai tên** — `clpP`/`clpP1`,
`ycf3`/`pafI`, `ycf4`/`pafII`, `psbN`/`pbf1` — vì PGA lấy tên từ **từng reference**, mà bộ reference
có cả hai phương ngữ; 4 gen không nhóm chức năng nào nhận.

**Ba cái bẫy, mỗi cái thành một lỗi có tên hoặc một lần thử lại (thay vì "sai mà không báo"):**

1. **Reference có toạ độ partial ⇒ PGA chú giải RỖNG và vẫn exit 0.** Đo được: reference
   `NC_024936.1` (có `<153124`, `>83020`) → 0 gen, "Finish annotating", exit 0. Adapter nay **từ
   chối theo tên**, liệt kê file nào có partial.
2. **Tên file tạm của PGA chứa `$`** (bảng ký tự random của chính nó), rồi
   `system("blastn … $reference1_random …")` đưa chuỗi đó qua shell lần nữa ⇒ shell nuốt `$…` ⇒
   blastn nhận đường dẫn không tồn tại ⇒ **0 gen, exit 0**. Xác suất 1 − (63/64)¹⁴ ≈ **20 %/lần**.
   Adapter **thử lại một lần** và ghi vào report khi phải thử lại.
3. **Dò version suýt xoá dữ liệu.** `run_tool` hỏi version bằng cách chạy `perl <tool> -e …`; PGA.pl
   nhận cờ lạ, chạy với tham số **mặc định**, và việc đầu tiên nó làm là `rm -rf gb` — tức xoá thư
   mục tên `gb` trong thư mục làm việc hiện tại. Nay `tool_version` **không bao giờ chạy script
   Perl** (hỏi version của chính `perl`), và binary được dò trong thư mục tạm. Regression test:
   `tests/unit/test_run.py::test_the_version_probe_never_runs_a_perl_tool`.

**Một phát hiện phụ, ảnh hưởng D3-02:** ba bản ghi RefSeq trong bộ reference
(`NC_085782.1`, `NC_085783.1`, `NC_086581.1`) **cũng dùng phương ngữ chức năng**
(`clpP1`·`pafI`·`pafII`·`pbf1`). Nghĩa là đây không chỉ là "quy ước của nhóm" như §8.10 viết — nó
đã xuất hiện trong các bản nộp RefSeq khác cùng chi. Quyết định D3-02 (khoá = tên chức năng, in
cả hai) vì thế **đứng vững hơn**, không yếu đi.

### 8.16. Đổi công cụ chú giải: Chloe thay PGA (03/10/2026)

**Vì sao đổi — theo nguyên tắc "sản phẩm quyết định, không phải bộ công cụ" (CLAUDE.md).** PGA đã
được đo (§8.15) và không tái tạo được số vàng: 147 gen (7 reference) và 130 gen (1 reference) so
với **132 = 86 + 38 + 8**. Nguyên nhân nằm ở chính thiết kế của nó: PGA bắt người dùng tự đưa bộ
reference, và lấy tên gen từ từng reference ⇒ gọi cùng một locus hai lần dưới hai tên.

**Chloe** (Ian Small, ĐH Tây Úc; `github.com/ian-small/Chloe.jl`, Julia, CLI
`julia --project=. chloe.jl annotate <fasta>`) có sẵn **bộ reference chuẩn hoá của chính nó**
(10 plastome tham chiếu, tải qua Julia artifact `Chloe_refs 1.0.1`) — đúng chỗ PGA thiếu.

**Đo trên genome vàng đã gỡ chú giải (cùng một input như §8.15):**

| Công cụ | Tổng | CDS | tRNA | rRNA | Thời gian |
|---|---|---|---|---|---|
| **Chloe** | **134** | **88** | **38** ✓ | **8** ✓ | 11,1 giây |
| PGA 7 reference | 147 | 95 | 44 | 8 | ~13 giây |
| PGA 1 reference | 130 | 87 | 35 | 8 | ~7 giây |
| *genome vàng (sự thật)* | *132* | *86* | *38* | *8* | — |

**tRNA và rRNA khớp tuyệt đối (38/38, 8/8).** CDS lệch đúng 2, và lệch đó **giải thích được**:
Chloe ghi `rps12` thành **4 dòng gen** (ghép theo exon 5'/3' của cả hai bản IR), genome vàng ghi 2
(mỗi bản IR một dòng, trans-spliced) ⇒ +2. Đây là khác biệt **cách biểu diễn**, không phải lỗi chú
giải — và bộ chuẩn hoá của mình gộp được nó.

**Còn phải làm trước khi chốt (ghi vào TASKS.md):** (1) viết `external/chloe.py` + đọc GFF bằng
reader của mình (GFF, không phải GenBank — Biopython đọc được), đi qua đúng bộ kiểm tra ORF +
chuẩn hoá tên; (2) gộp rps12 theo exon rồi đếm lại ⇒ kỳ vọng **132 = 86 + 38 + 8**; (3) **giấy
phép: Chloe không có file LICENSE** — cùng loại với MISA, phải vào đợt rà soát M9.3 trước khi phát
hành; (4) Chloe chỉ dành cho **thực vật có hoa** (angiosperm) — đúng phạm vi bộ công cụ, nhưng phải
nói rõ trong tài liệu; (5) phụ thuộc **Julia** (conda-forge) — nặng hơn perl, chấp nhận được vì
annotator là bước tuỳ chọn và Docker image chứa được.

**Chưa chốt vội:** cần đọc lại output bằng reader của mình và đếm lại sau khi gộp rps12 — nếu ra
đúng 132 thì D3-01 được sửa thành **Chloe là pipeline offline, GeSeq ghim là kiểm tra chéo**.

**Đính chính (03/10/2026, sau khi đọc lại GFF bằng reader của mình — §8.17):** hai điều ở trên
phải sửa. (1) **134/88 là do cách đếm, không phải do Chloe**: gộp các dòng exon theo `parent` mà
chúng chia sẻ thì file ra **86 CDS** — đúng số vàng, và **không cần "gộp rps12"** như dự kiến.
(2) Bảng trên đo trên GFF **mặc định**, tức GFF có **biến đổi trình tự** — số gen vẫn đúng nhưng
toạ độ nằm ở hệ khác; xem §8.17.

### 8.17. Chloe chạy thật qua pipeline của mình — P3-03 chốt (03/10/2026)

**Kết quả: 132 = 86 + 38 + 8, và 0 CDS hỏng ORF** (80 `ok` + 6 `alternative_start`). Đây là
Done-when của P3-03, đo bằng **reader của mình** (`core/gff.py`) trên **file GFF thật của Chloe**,
rồi chạy lại lần nữa end-to-end qua CLI với Chloe thật.

| | |
|---|---|
| Gen chú giải được | **132** = 86 CDS + 38 tRNA + 8 rRNA (khớp số vàng) |
| Kiểm tra ORF | **80 ok + 6 alternative_start, 0 lỗi** (bản ghi vàng: 81 + 5, 0 lỗi) |
| Tên gen | 5 tên có tên đồng nghĩa, **0 tên không xếp được nhóm** |
| Chạy end-to-end | 22 giây (10 giây chú giải + ~12 giây Julia khởi động/biên dịch) |
| Provenance ghi được | `chloe 0.1.18` (pinned) · `julia version 1.13.0` · argv đầy đủ · rc 0 |

**File thật được lưu vào repo:** `tests/data/chloe_NC_050833.1.gff` (sha256 `c8b660b7…`, ghi ở
`tests/data/PROVENANCE.md`). Nhờ nó, ba con số trên kiểm lại được **trên máy không có Julia** —
đúng tinh thần "công cụ chỉ nằm trên một máy thì không phải một phần của bộ công cụ" (§7.8c).

**Phát hiện lớn nhất — Chloe mặc định BIẾN ĐỔI trình tự, và đếm gen không thấy điều đó.**
`--no-transform` được chính Chloe mô tả là *"do not flip and orient sequence to standard
configuration"*. Đo trên genome vàng, hai lần chạy chỉ khác đúng cờ đó:

| | Mặc định (có biến đổi) | Có `--no-transform` |
|---|---|---|
| Số gen | **132** — y hệt | **132** |
| Tên gen | y hệt | y hệt |
| Kiểm tra ORF | **12/86 CDS hỏng** (8 `no_start`, 4 `no_stop`) | **0 hỏng** |
| Trình tự viết kèm | `.chloe.fa` khác input ở **13.001 vị trí**, đúng dải **109.342–127.510** | không viết file nào |

Dải đó **chính là SSC** (18.169 bp): Chloe **đảo ngược chiều vùng sao chép đơn** rồi chú giải hệ
toạ độ đó. Ví dụ cụ thể: `ndhH` trong file mặc định nằm ở **114.734–115.915 (+)** — dịch xuôi vùng
đó ra `IPTPITGRKIQIETNSSRGPESIKLEFWTLNSLYP…`, không phải ndhH; còn ndhH thật ở 120.937–122.118 (−),
dịch ngược ra `MTVPATRKDFMIVNMGPHHPSMHGVLRLILTLDGE…`. Vì **số gen và tên gen không đổi**, một
pipeline chỉ đếm gen sẽ báo "khớp số vàng" trong khi một phần ba bộ gen nằm sai khung — đúng loại
"sai mà không báo" mà dự án này tồn tại để chặn.

⇒ Adapter **luôn truyền `--no-transform`**, và **từ chối mọi lần chạy có viết file FASTA biến đổi**
(so từng base với input đã đưa vào). Có test cho cả hai đường: stub (không cần công cụ) và chạy
thật với Chloe (tier B, đo lại đúng hiện tượng trên).

**Hai khác biệt tên gen lộ ra khi đối chiếu với bản ghi vàng:**

1. **`ndhI` (Chloe) vs `ndhL` (bản ghi vàng) — cùng locus, đã chốt.** Cùng vị trí, **cùng
   protein** (543 bp, dịch ra một chuỗi y hệt); `/product` của chính bản ghi vàng ghi *"NADH-
   plastoquinone oxidoreductase subunit I"*; họ ndh của plastome có tiểu đơn vị A–K, **không có L**
   ⇒ thêm cặp alias thứ năm vào `data/gene_aliases.yaml` (canonical `ndhI`, synonym `ndhL`).
2. **`trnG-GCC` (Chloe) vs `trnG-UCC` (bản ghi vàng) — CHƯA chốt, và không tự chốt.** Bản ghi vàng
   dùng tên `trnG-UCC` cho **cả hai** bản tRNA glycine, dù hai bản **khác trình tự**; Chloe gọi bản
   không intron là `trnG-GCC`. Bộ tRNA chuẩn của plastome thực vật có hoa có **cả hai** gen
   (`trnG-UCC` có intron và `trnG-GCC`), nên khả năng cao bản ghi vàng đặt trùng tên — nhưng **đo
   không đủ để kết luận**: đọc anticodon theo offset cố định chỉ khớp **11/38** tRNA của bản ghi
   vàng (đúng như ghi chú của bộ chuẩn hoá), nên không có cách đọc máy móc nào đủ tin ở đây. Đây là
   **câu hỏi cho cô Trang** (§10): tên nào in trong bảng bản thảo.

**Phạm vi & giấy phép (cho M9.3 và cho phần Methods):** Chloe chỉ chú giải **thực vật có hoa**
(angiosperm) — đúng phạm vi ba loài bản thảo, nhưng phải nói rõ trong tài liệu; plastome ngoài
nhóm này thì dùng `--set annotator=pga`. Checkout **không có file LICENSE**; bản công bố trên
Zenodo (doi:10.5281/zenodo.18194040, v0.1.18) ghi **MIT** — hai chuyện khác nhau, đã ghi cả hai
vào `docs/tools/chloe.md` để đợt rà soát M9.3 quyết.

**Một lỗi nhỏ của chính mình, ghi lại vì cùng loại:** helper `_length()` trong test tier B đọc
độ dài ở `LOCUS.split()[1]` — đúng với bản ghi tự chế mà nó được viết ra để phục vụ, **sai** với
bản ghi thật (`LOCUS NC_050833 153724 bp`, vị trí 1 là tên bản ghi). Nó chỉ lộ ra khi test Chloe
chạy trên genome vàng, vì máy này không có PGA nên nhánh test kia chưa từng chạy. Một helper chỉ
được chạy với đúng fixture nó sinh ra thì chưa phải là helper đã được kiểm.

---

### 8.18. Chạy thật `table2asn` — bản `.tbl` của mình bị NCBI chê, và ba lỗi đã sửa (03/10/2026)

**Vì sao chạy:** nửa `.sqn` của P3-05 bị xếp là "chờ file `.sbt` của nhóm". Form của NCBI cần đăng
nhập, nhưng **`.sbt` chỉ là một file ASN.1 văn bản theo một cấu trúc xác định** — tìm được các mẫu
công khai (Cenote-Taker3, enviro-lab/ncbi-submit) và viết một bản **stand-in** cho repo
(`tests/data/standin_template.sbt`, tên/địa chỉ là chỗ trống, dòng tiêu đề ghi rõ "NOT FOR
SUBMISSION"). Có nó thì chạy được cả bước.

**Kết quả lần chạy đầu — bản `.tbl` của mình KHÔNG đạt:** NCBI's validator trả về **13 lỗi
(ERROR) và 27 cảnh báo** cho genome vàng. Ba nguyên nhân, cả ba đều là của mình:

| Lỗi NCBI báo | Nguyên nhân | Sửa |
|---|---|---|
| `NoSourceDescriptor` + `BioSourceMissing` (2) | `.tbl` **không có feature `source`** — người nộp phải nói trình tự là gì | ghi feature `source` từ **qualifier của chính bản ghi** (organism · organelle · cultivar · specimen_voucher…) |
| `SEQ_FEAT.StartCodon` (2) | `rpl2` mở đầu bằng **ACG** — không phải mã khởi đầu của bảng 11, vì RNA editing C→U biến nó thành AUG | `/exception` của bản ghi **đi theo gene** vào bảng (trước đây bị bỏ) |
| `SEQ_INST.BadProteinStart` (2) | cùng CDS đó, không có mã khởi đầu nên dịch ra protein có gap ở đầu | ghi **start partial** (`<`) — giữ nguyên toạ độ, chỉ thêm dấu |

**Thêm một phát hiện nhỏ nhưng chặn đường:** `mol_type` và `db_xref` **không phải qualifier của
định dạng bảng** — table2asn gọi chúng là "Unrecognized qualifier name" và **exit 2** (dù vẫn ghi
`.sqn`). Hai giá trị đó NCBI tự suy ra, nên bảng để chúng ở lại trên genome và **không ghi vào
`.tbl`**; module nói rõ điều này trong notes. Bảng nộp thật của *Dioscorea alata* có `mol_type` và
cũng exit 2 vì nó.

**Đối chứng — và đây là chỗ chứng minh phần còn lại không phải lỗi của mình.** Lấy **một bảng nộp
thật đã được chấp nhận**: `Dalata-cp.tbl` + `Dalata-cp.fasta` của *Dioscorea alata*
(bredeson/Dioscorea-alata-genomics, bảng dùng cho bản nộp GenBank thật), chạy qua **đúng
validator đó**:

| | Lỗi | Trong đó `GeneXrefStrandProblem` | Cảnh báo |
|---|---|---|---|
| Bảng nộp thật (Dioscorea) | **9** | **7** | 19 |
| Bảng của mình — trước khi sửa | 13 | 7 | 27 |
| **Bảng của mình — sau khi sửa** | **7** | **7** | 98 |

**7 lỗi còn lại là `GeneXrefStrandProblem`** — lớp lỗi mà định dạng này sinh ra cho **mọi** gene
nằm trong IR (có hai bản, một + một −) khi có feature `gene` riêng: bảng nộp thật của Dioscorea
cũng đúng 7 lỗi cùng lớp, còn có thêm 2 lỗi vị trí `rps12` mà bảng của mình không có. Nói cách
khác: **bảng của mình giờ ngang một bản nộp đã được chấp nhận**, và lỗi còn lại là thuộc tính của
định dạng + genome có IR, không phải của writer.

**Cảnh báo:** 98, trong đó **86 là `SEQ_FEAT.GenCodeMismatch`** với nội dung tự mâu thuẫn —
"conflict between CDS (code 11) and BioSource… (uses code 11)" — tức bảng mã hai bên *giống nhau*;
nó xuất hiện khi bản ghi khai `organelle`, và bản ghi vàng (RefSeq, đã được chấp nhận) cũng ghi
`/transl_table=11` trên mọi CDS. **Không hành động được.**

**Module giờ đọc file `.val` chứ không chỉ giữ nó:** số lỗi/cảnh báo vào bảng tổng, phân loại theo
mã vào bảng `validator_findings`, `.val` đi kèm như artifact, và báo cáo in ra đúng câu "đọc
`.val` trước khi nộp" — một file không ai đọc thì không phải là một phép kiểm.

**Ý nghĩa cho sản phẩm:** trước lần chạy này, bộ công cụ có thể tạo ra một cặp `.fsa` + `.tbl`
*trông* đúng và **bị NCBI từ chối ở khâu upload** — sau khi đã làm xong mọi thứ. Đúng loại lỗi mà
P3-05 tồn tại để chặn, và nó chỉ lộ ra khi có công cụ của chính NCBI để đối chiếu.

---

### 8.19. Nửa còn lại của luật backend: ghim kết quả công cụ chạy ở nơi khác (03/10/2026)

**Luật đã có, không phải quyết định mới.** §3 (P1-12) và §13.10 đã chốt: bước nặng thì **chạy ở
đâu cho phù hợp** — máy phổ thông khi chạy được, web/máy chủ khác khi không. Nhưng kiểm lại thì
thiếu đúng một mảnh: `plastid panel pin` **chỉ nhận accession NCBI**, tức là kết quả của một công
cụ chạy ở nơi khác **không ghim được**. Luật có mà đường thực hiện thì không.

**Đã bổ sung:** `plastid panel pin --file <path> [--note "PMGA web, 03/10/2026"]` — ghim một bản
ghi do nơi khác sinh ra, với **đúng những gì một lần tải NCBI có**: bytes · sha256 · md5 · bản ghi
được **đọc lại bằng reader của mình** (kèm đo IR như mọi pin khác) · và `source`/`note` nói rõ nó
từ đâu tới. Hai khác biệt so với ghim theo accession: khoá là **id của chính bản ghi** (không có
accession thì lấy tên file — đúng ca một mẫu chưa có accession), và `PINNED.json` được **thêm vào**
chứ không bị thay, vì ghim một file là ghim *vào* một panel đã có.

Đo trên mitogenome thay thế: `pinned OX459128.1: 381931 bytes · sha256 e31da1d3… (258274 bp ·
31 CDS + 17 tRNA + 3 rRNA · circular)` — khớp đúng số §8.14. Bản ghi không đọc được vẫn được ghim
kèm `parse_error` (từ chối ghim sẽ để lại một con số mà không ai giữ nguồn).

**Chưa chạy, và đo ra thì KHÔNG chạy được — cả hai đường (03/10).** Tôi đã khép task này trên một
đường chưa từng chạy; kiểm lại thì:
* trang chủ `www.1kmpg.cn/pmga/` trả **HTTP 200 nhưng chỉ là redirect** sang
  `47.96.249.172:16084/home`, và host đó **từ chối kết nối** (`Connection refused`, ~0,2–2 giây)
  trên **cả ba cổng 16084 · 16081 · 16086** — tức không có gì đang lắng nghe; `www.1kmpg.cn/mgavas/`
  cũng chỉ là frame trỏ vào `47.96.249.172:16081`;
* nên **cả hai đường đều không chạy được hôm nay**: local vì máy yếu, web vì **service của tác giả
  đang chết**. Phép kiểm số gen mito vì thế **vẫn mở**, và cách đóng nó là *chạy lại*, không phải
  một task mới: `curl -sS -o /dev/null -w '%{http_code}' http://47.96.249.172:16084/home` — khi nào
  có đáp án thì một lần chạy web + `panel pin --file` là xong.
* bước `panel pin --file` thì **đã test**, nhưng bằng **một file stand-in** (bản ghi NCBI copy ra
  đĩa), **không phải output PMGA thật** — hai chuyện khác nhau, không được gộp. Và một điểm
không phải việc của dev: đẩy genome **chưa công bố** của nhóm lên server bên thứ ba là **quyết định
của cô Trang** (dữ liệu rời khỏi lab) — chạy kiểm bằng genome công khai thì không vướng gì.

**Hệ quả cho P3-04:** PMGA giờ có **đủ hai đường** đúng như luật: (1) container local — đã có trên
máy này (11 GB, cần Apptainer; **không phải đường cho máy 4–8 GB**), (2) chạy web service một lần
rồi **ghim file** — chạy được ngay, không cần container. Adapter `external/pmga.py` vì thế là tiện
ích về sau, **không còn là thứ chặn task**; điều kiện viết nó là có máy chạy nổi container, và số
đầu tiên phải đo khi đó là **RAM đỉnh + thời gian**.

---


**Đo được thêm, cùng ngày — dung lượng ĐĨA THẬT.** Container 3,18 GB giải nén ra 11 GB, và nó bị
xoá ngay sau khi đọc xong CLI, vì: `df` **bên trong WSL báo còn ~905 GB — đó là trần đĩa ảo 1 TB,
không phải dung lượng thật**. Ổ **C: thật chỉ còn 42 GB (đã dùng 89 %)**, và file `ext4.vhdx`
(54 GB) **không co lại khi xoá file bên trong** — WSL từ chối chế độ sparse của chính nó
("currently disabled due to potential data corruption"), nên dòng `sparseVhd=true` trong
`.wslconfig` **không có tác dụng**; muốn lấy lại phải compact thủ công từ Windows. Bài học: mọi
thứ tải/cài trong WSL đều nằm trong một file không tự nhỏ đi, nên **phải đo `df -h /mnt/c`**, và
cân nhắc đường không-cần-cài (web service + `panel pin --file`) trước.

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

---

## 9. Lỗi đã gặp & số đo (tổng hợp)

### A. Lỗi "sai mà không báo" ở khâu nhận diện IR

1. Cách viết đơn giản khớp chú giải NCBI **7/9** ca, **sai 2/9** và **không báo lỗi** (số ở §8.3).
2. **2/9 IR có hai bản đã phân hoá** (khớp tổng 92,63 % và 53,65 %) ⇒ phép so khớp **không cho phép gap** sai ở nhóm này.
3. Ranh giới IR trong chú giải GenBank là **quy ước của người chú giải**, không phải kết quả so khớp.
4. **Hằng số ghép cặp phải là c − 2**; lệch 1 base → IR đo ra **1 bp** mà chương trình **không báo lỗi**.
5. **Điểm bắt đầu genome**: 1/20 genome RefSeq bắt đầu ngoài LSC; tài liệu dẫn 27,67 % trên GenBank rộng.
6. IR sai → phân vùng LSC/SSC sai → các phân tích theo vị trí (ranh giới, SNP, π, biểu đồ % tương đồng) sai theo, output không có dấu hiệu.
7. Hiệu năng: 20 genome trong **4,7 giây**, RAM 331 MB.
8. **Mới 03/10 (§8.6, `PZ428230.1`) — cổng PASS trên một IR sai độ dài.** Bản ghi không có chú giải,
   detector de novo tìm được một cặp lặp **giống nhau 100 %** nhưng chỉ dài **15.756 bp**, trong khi cả
   5 loài Coleus còn lại đều ~25,7 kb → cổng **PASS** (cổng chỉ kiểm *identity*, không kiểm *độ dài*).
   Đo mở rộng vùng ra hai bên cho thấy identity rơi ngay lập tức ⇒ **dữ liệu sai**, không phải detector
   bỏ sót. Hệ quả: **bản ghi không chú giải thì không có gì đối chiếu IR** — cần kiểm tra chéo độ dài IR
   với các bản cùng chi trước khi dùng cho phân tích theo vị trí (ghi thành việc mở, chưa làm).

### B. Đếm gen, SSR, Ts/Tv

8. **"135 gen" không tái lập được** từ bất kỳ cách đếm nào: 132 (CDS+tRNA+rRNA) · 133 (feature `gene`) · 112 (tên duy nhất) · 21 (gen lặp IR).
9. **1 bản `ycf1` cụt** gây lệch cách đếm: tại 108.448–109.347 chỉ có `gene` + `misc_feature` (không CDS); bản thật ở 122.984–128.404 mới có CDS.
10. Bản tự viết đếm SSR ra **36** (thiếu 1); MISA và cpstools đều ra **37**.
11. **"mono 10 · di 6 · tri 5 · tetra 4 · penta 4 · hexa 4" là ngưỡng MISA, không phải số lượng.** Phân loại thật: mono 31 · di 4 · tri 2.
12. **Số SSR phụ thuộc ngưỡng**: ≥8 → 133 · ≥10 → 36 · ≥12 → 16. `misa.ini` của anh Điềm dùng `(1/10)(2/6)(3/5)(4/5)(5/5)(6/5)`; tài liệu MISA ghi mặc định 4 cho unit 5–6 bp.
13. **Ts/Tv phụ thuộc cách căn**: difflib 2,667 vs MAFFT 2,686.
14. **Mới 03/10/2026 — parser GFF của MISA đọc thiếu một nửa từ vựng (im lặng).** MISA ghi SSR
    monomeric là `monomeric_repeat`, nhưng **mọi** unit dài hơn là `microsatellite` (lớp nằm
    trong cột Note), và một cụm compound là **một** dòng `repeat_region` chứa nhiều thành phần.
    Parser cũ chỉ nhận `monomeric_repeat` ⇒ **mất 7/37 SSR** (bảng per-SSR còn 30 dòng) trong khi
    bảng tổng vẫn ghi 37 vì lấy từ file `.statistics` — hai bảng âm thầm lệch nhau.
    Phát hiện bằng cột đối chiếu của P2-06. **Bài học:** fixture test cũ tự bịa type
    `dinucleotide_repeat` mà MISA không bao giờ ghi ⇒ test xanh nhưng sai; fixture phải lấy
    **hình dạng dòng thật**.
15. **Mới 03/10/2026 — thứ tự các phần của CDS trans-spliced.** Sắp xếp lại các part, hoặc đọc
    tất cả theo mạch của gen, đều cho ORF sai. Đo trên **cả hai** bản rps12 của genome vàng:
    chỉ cách "đọc từng part theo mạch của chính nó, giữ nguyên thứ tự annotation" mới ra
    `ATG…TAA` (sắp xếp lại → `AAT…TAT`; đọc xuôi hết → `ATA…ATT`).
16. **Mới 03/10/2026 — khoá bảng trùng thì tham chiếu rơi vào một dòng bất kỳ.** Gen nằm trong IR
    có 2 bản nên `intron_id`/tên gen lặp; `Table` đã được vá để **từ chối** khoá trùng.
17. **Mới 03/10/2026 — `yaml.YAMLError` KHÔNG phải `ValueError`.** File YAML hỏng làm `load_groups`
    ném traceback thay vì thông báo; đã dịch lỗi tại biên.
18. **Mới 03/10/2026 — `annotated=True` là một lời khai chưa ai kiểm.** `entry_from_summary` cũ
    đặt `annotated=True` cho mọi bản ghi; thực đo trên bộ dữ liệu thay thế: **3/11 plastome có
    0 CDS** (GenBank không mặc nhiên có chú giải). Đã sửa: chỉ tin RefSeq, bản GenBank **đo sau
    khi tải**.
19. **Mới 03/10/2026 — dòng `Total` của bảng đối chiếu tự tố cáo chính nó.** Lần chạy đầu của
    `gene_table`: CDS 86/86 · tRNA 38/38 · rRNA 8/8 đều `agrees`, nhưng dòng **Total ra
    `in_table = 0`** — `by_kind.get("Total")` trên một `Counter` chỉ có ba khoá CDS/tRNA/rRNA.
    Đúng loại lỗi mà cột đối chiếu sinh ra để bắt: nếu bảng chỉ in số mà không tự so, một dòng
    tổng bằng 0 sẽ nằm giữa bảng và **không ai kiểm**. Đã sửa (`_in_table`), thành test.

20. **Mới 03/10/2026 — biến đổi trình tự làm đúng SỐ GEN mà sai TOÀN BỘ KHUNG (§8.17).**
    Chloe mặc định đảo chiều vùng sao chép đơn rồi chú giải hệ toạ độ đó. Hệ quả đo được: **132 gen
    — y hệt, tên y hệt**, nhưng **12/86 CDS hỏng ORF** khi đọc với trình tự người dùng đưa vào, và
    vị trí các gen ở SSC lệch hẳn (ndhH ở 114.734–115.915 (+) thay vì 120.937–122.118 (−)).
    **Bài học:** một con số đếm được không chứng minh được hệ toạ độ; phải có một phép kiểm *theo
    vị trí* (ở đây: dịch ORF) mới thấy. Đây là lý do P3-03 bắt buộc "đọc lại bằng reader của mình"
    chứ không tin output của annotator.
21. **Mới 03/10/2026 — annotator và bản ghi vàng gọi CÙNG một locus bằng hai tên khác nhau.**
    `ndhI` (Chloe) = `ndhL` (bản ghi vàng): cùng vị trí, **cùng protein**, và `/product` của bản
    ghi vàng ghi "subunit I" ⇒ thành cặp alias thứ năm. Ca thứ hai **chưa chốt**: `trnG-GCC` vs
    `trnG-UCC` — bản ghi vàng dùng **một tên cho hai bản tRNA khác trình tự**, mà bộ gen chuẩn có
    cả hai gen. Không tự sửa: đo anticodon bằng offset cố định chỉ khớp 11/38 tRNA ⇒ **không đủ
    bằng chứng**, đưa vào §10 cho PI quyết.

22. **Mới 03/10/2026 — ba thứ định dạng nộp bắt buộc mà model gene không mang (§8.18).**
    Chạy `table2asn` thật trên genome vàng: **13 lỗi**, sửa ba chỗ còn **7**. Hai chỗ đầu là
    *passthrough* (feature `source`, qualifier `/exception`) — bản ghi có, mình bỏ. Chỗ thứ ba là
    một quy ước: **CDS không mở đầu bằng mã khởi đầu hợp lệ thì phải ghi start partial**, nếu
    không NCBI báo `BadProteinStart`. Bài học: một model chỉ giữ những gì module *hiện tại* cần
    thì mỗi định dạng xuất ra sẽ phát hiện lại từ đầu những gì nó bỏ.
23. **Mới 03/10/2026 — "exit 0" không phải là tiêu chí của `table2asn`.** Nó **exit 2** khi chỉ
    *cảnh báo* về một qualifier lạ (`mol_type`, `db_xref`) trong khi `.sqn` vẫn đúng; bảng nộp
    thật của *Dioscorea alata* cũng exit 2 vì lý do đó. Nghĩa là: đọc `.val`/`.stats`, đừng tin
    mã thoát. Và ngược lại — `.tbl` thiếu feature `source` thì **vẫn exit 0** nhưng bị đánh 2 lỗi.
24. **Mới 03/10/2026 — phép kiểm mới chỉ có giá trị khi có ĐỐI CHỨNG.** 7 lỗi
    `GeneXrefStrandProblem` trông như lỗi của mình cho tới khi chạy **cùng validator** trên một
    bảng nộp đã được chấp nhận (Dioscorea) và thấy nó cũng đúng 7 lỗi cùng lớp. Không có đối
    chứng thì hoặc là sửa những thứ không hỏng, hoặc là bỏ qua những thứ hỏng thật.
25. **Mới 03/10/2026 — thư viện viết GenBank ghi ra cú pháp không ai dùng, và round-trip bằng
    chính nó thì vẫn xanh.** BioPython gặp qualifier trần (`/pseudo`, `/trans_splicing`) dưới dạng
    **chuỗi rỗng** thì ghi `=""`; gặp `None` mới ghi `/pseudo` trần. File `=""` vẫn đọc lại được
    bằng chính BioPython, nên một test kiểu "ghi rồi đọc bằng thư viện đã ghi" **xanh trong khi
    file khác cú pháp NCBI**. Bắt được chỉ vì có **bản ghi vàng thật** để đối chiếu cú pháp và vì
    **reader của mình** mới là thứ đọc lại. Bài học: round-trip bằng chính thư viện đã ghi không
    phải bằng chứng — phải đối chiếu với một file do NCBI viết (§8.20).

### C. Môi trường & vận hành

20. **`chloroplast[filter]` của NCBI bỏ sót bản ghi ghi "plastid"** (đo 03/10/2026, khi chạy
    Reference Panel thật): `Iridaceae[Organism] AND refseq[filter] AND chloroplast[filter]` → **71**;
    đổi sang `plastid[filter]` (hoặc `(chloroplast OR plastid)`) → **76**, khớp con số §8.3. Ca cụ thể:
    **NC_024936.1** (*Iris gatesii* — chính genome so sánh của bộ số vàng) có title
    *"Iris gatesii voucher T06-52 plastid, complete genome"*: khớp `refseq[filter]` và
    `Iridaceae[Organism]` nhưng **không** khớp `chloroplast[filter]`. Bài học: lọc theo `[filter]`
    của NCBI phải kiểm bằng số đếm, không tin tên filter.
    *Xác nhận kèm (cùng lần chạy):* `plastid panel pin` tải lại NC_050833.1 từ NCBI ra sha256
    `a946f23e…` — **khớp đúng byte** với dữ liệu vàng đang vendor trong repo; NC_024936.1 ra
    `fd441727…` cũng khớp. Nguồn dữ liệu vàng được xác nhận lại end-to-end.

14. **Máy 8 GB không đủ cho phần lắp ráp:** WSL2/Docker mặc định ~50 % RAM host → còn **1,5–2,5 GB mỗi việc**; 120.000 cặp read là hết RAM (phải hạ còn 20.000); flye cần ≥8 GB.
15. **Dò công cụ theo PATH là không đủ**: `misa.pl` ở `~/toolbox/src/misa/`; `nucmer`/`repeat-match` phải activate env `mummer4`; `trf` không có.
16. **Đường dẫn có dấu tiếng Việt** gây lệch đường dẫn giữa hai máy (thư mục "chạy Trí").
17. Script gửi kèm **không tới được máy anh Điềm** — anh làm ở thư mục khác (`D:\diem\bai bao tool\ke_hoach_new\chạy Trí\probe`).
18. `docker` không gọi được từ WSL trước khi bật WSL Integration (thông báo lỗi ở §8.4).
19. Chuỗi phân tích trên máy 2,8 GB dùng 331 MB và 9,8 giây.
20. **Mới 03/10/2026 — bộ tỉa style đọc sai chỗ, và bản render vẫn "trông ổn".** Style trong
    OOXML được tham chiếu bằng **phần tử con** (`<w:pStyle w:val="Heading1"/>`), không phải
    thuộc tính. Bản đầu đọc `element.get("w:pStyle")` ⇒ không thấy style nào được dùng ⇒ tỉa
    gần hết, kể cả `Title` và `Heading 1`. **LibreOffice vẫn render ra hình như bình thường**
    vì nó lùi về mặc định của chính nó — lỗi chỉ lộ ra khi mở PDF cạnh bản đúng, và chỉ bị bắt
    vì có **test cấu trúc** ("mọi style tài liệu dùng phải còn trong styles.xml"). Bài học:
    tỉa gói OOXML thì phải canh bằng cấu trúc, không bằng mắt.

---

## 10. Câu hỏi còn treo

1. ~~Phạm vi con số Ts/Tv~~ — **ĐÃ CHỐT 02/10**: cổng dùng **chỉ vùng CDS · 76 gen dùng chung · căn bằng MAFFT**, **không tính vùng lặp IR** ⇒ mốc **2,686** (317 Ts / 118 Tv). Số Ts/Tv của bộ Eleutherine (1,095–1,777) là phạm vi khác (SNP cả plastome) — chỉ dùng cho nội dung riêng của bộ đó và **phải dán nhãn phạm vi** ở mọi chỗ xuất ra.
2. ~~Ngưỡng `misa.ini`~~ — **ĐÃ ĐÓNG 02/10**: file `misa.ini` đang dùng trên máy anh Điềm **chính là file gốc trong zip của IPK** (md5 `49a1c56048631647c0364532a6964cf9` khớp bản tải lại): `1-10 2-6 3-5 4-5 5-5 6-5`, interruptions 100, GFF true ⇒ **dùng đúng mặc định của MISA**, không phải tuỳ chỉnh.
3. Ai có máy Mac để kiểm.
4. **Dữ liệu còn thiếu:** có lấy được reads Illumina thật của 1 mẫu không (để đo ngưỡng bộ nhớ — hiện mới có read mô phỏng).
5. **Dữ liệu còn thiếu:** có lấy được **genome nhân** của mẫu không — NUPT/MTPT cần genome nhân để chạy (hiện chưa có).
6. ~~Hệ toạ độ sau chuẩn hoá~~ — **ĐÃ CHỐT 02/10**: bảng theo từng genome ghi **toạ độ GenBank**; hình so sánh chéo nhiều genome dùng **hệ đã chuẩn hoá** (các vùng thẳng hàng nhau); mọi bảng ghi rõ hệ đang dùng; `rotation_offset` lưu trong `00_run_info`.
7. ~~Cách sinh report~~ — **ĐÃ CHỐT 02/10**: sinh từ đối tượng `Result`. Mỗi số trong report trỏ về **đúng một ô** bảng gốc; template chỉ trỏ (`{{bảng.cột}}`, `{{bảng[dòng].cột}}`), không gõ số; ref không giải được ⇒ **raise**, không in trống. Yêu cầu "kiểm thử truy vấn số" (§13.5) được thoả bằng cấu trúc, không bằng kỷ luật.
8. ~~Reference Panel: snapshot hay truy vấn sống~~ — **ĐÃ CHỐT 02/10**: engine gợi ý qua **cache RefSeq có ngày** (lệnh làm mới tường minh); người dùng chốt danh sách thì **tải + ghim checksum**; chạy lại offline, tái lập được.
9. Mục **(c)** của Reference Panel còn thiếu
10. **Tên gen `trnG` in trong bảng bản thảo là tên nào (§8.17).** Bản ghi vàng NC_050833.1
    dùng `trnG-UCC` cho **cả hai** bản tRNA glycine (hai bản khác trình tự), Chloe gọi bản không
    intron là `trnG-GCC`, và bộ tRNA chuẩn của plastome có **cả hai** gen. Đây là chuyện *chữ in
    trong bản thảo / bản nộp*, tức việc của cô Trang (D3-06), không phải việc dev tự quyết — dev
    đã ghi lại đủ bằng chứng để cô Trang quyết trong một lần xem.
 — đoạn văn bản nguồn bị cụt.

## 11. Nguồn tài liệu (ngoài file này)

Nằm trong `D:\business\pts\`:

| File | Nội dung |
|---|---|
| `Ke_hoach_xay_dung_bo_cong_cu_phan_tich_bo_gen_luc_lap_22.9.2026_VuThiHuyenTrang.docx` | Kế hoạch gốc (91.729 ký tự) |
| `Bao_cao_tool_thay_the_NguyenThanhDiem.docx` | 5 + 3 + 12 công cụ bị chặn đều có thay thế đã chạy thật |
| `Danh_sach_tool_khong_cai_duoc_NguyenThanhDiem.docx` | 84 công cụ đã thử · 64 chạy · 17 không chạy |
| `REVIEW_2026-09-26_3_tai_lieu_va_muc_tieu_may_pho_thong.docx` | Bản kỹ thuật đầy đủ: mục 9 (khả thi theo khối) · mục 10–12 (đo rủi ro đường nối, chạy thử liền mạch, đóng vòng bằng công cụ chuẩn) |
| `BAO_CAO_4_CAU_HOI_PlastidToolSuite_2026-09-26.docx` | Bản cho người không chuyên (khả thi / vấn đề / giải quyết được không / giải quyết thế nào) + mục 8 chiến lược tích hợp |
| `KE_HOACH_TRIEN_KHAI_PlastidToolSuite_2026-09-26.docx` | Kế hoạch triển khai (môi trường dev, kiểm thử, đóng gói, lộ trình 8 tuần) |
| `HUONG_DAN_SETUP_MOI_TRUONG_CODE_2026-10-01.docx` | Hướng dẫn setup WSL / Claude Code / Docker |
| `thanhdiem_result/` | Trả lời phiếu việc của anh Điềm + kết quả thô (`.misa`, `.statistics`, cpstools, MAFFT per-gene) |
| `probe/` | Script + dữ liệu dùng để đo các lỗi ở §8–9 (`chain_e2e.py`, `chain_tools.py`, `ir_gate.py`, `iridaceae20.gb`) |

## 12. Lệnh cài gọn cho máy mới (nguyên văn từ báo cáo anh Điềm)

```bash
conda create -y -p ~/toolbox/env/thaythe -c conda-forge -c bioconda \
  sibelia bandage getorganelle spades oatk mummer4 mreps trf vmatch \
  raven-assembler novoplasty minimap2 samtools wgsim blast paml beast2 \
  mafft perl perl-bioperl perl-xml-simple
conda install -y -c conda-forge pycirclize
conda install -y -c bioconda cpstools        # SỬA 02/10: cpstools KHÔNG có trên conda-forge (§7.8)
git clone https://github.com/daisieh/plann ~/toolbox/src/plann
git clone https://github.com/Adamtaranto/pypopart ~/toolbox/src/pypopart \
  && pip install -e '~/toolbox/src/pypopart[viz]'
curl -L -o embryophyta_pltd.fam \
  https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam
get_organelle_config.py -a embplant_pt,embplant_mt   # 1 lần, cần mạng
# MISA standalone (misa.pl + misa.ini) — nguồn IPK, kiểm md5 trước khi dùng:
curl -L -o misa_sourcecode_25082020.zip \
  https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip
md5sum misa_sourcecode_25082020.zip    # phải ra cf9d3579f7caacd02f99830b329c7d9c
unzip misa_sourcecode_25082020.zip -d ~/toolbox/src/misa
```

**Sửa 02/10/2026:** thêm `mafft` và `perl` (thiếu thì không ra được số vàng Ts/Tv 2,686 và SSR 37);
tách `cpstools` sang kênh bioconda; thêm lệnh tải MISA kèm md5. Lệnh cũ trong báo cáo anh Điềm giữ
nguyên ở §8.1 để đối chiếu.

---

## 13. ROADMAP — thứ tự phát triển

**Thứ tự đã đổi theo quyết định 01/10: phân tích trước, lắp ráp sau.** Làm từ trên xuống; mỗi giai đoạn ra một thứ dùng được ngay.

**Điều kiện chung cho mọi giai đoạn** (đọc trước khi làm bất kỳ giai đoạn nào):

- **Áp dụng cho cả lục lạp và ty thể** — trừ phần ước tính thời gian phân kỳ (chỉ cần cho bài báo công cụ).
- **Cổng bắt buộc ở mọi phân tích theo vị trí: hai bản IR phải giống ≥ 99,9 %.** Không đạt ⇒ cờ đỏ, **dừng**, không cho phân tích theo vị trí chạy tiếp. Lý do đo được: cách nhận diện IR thông thường sai âm thầm **2/9 ca (22 %)**, cổng này bắt đúng 2/2 ca sai và báo oan 0/7 ca lành (chi tiết §9 mục A).
- **Cách đếm gen (không ghi trần "132 gen"):** **132 = 86 CDS + 38 tRNA + 8 rRNA**, không dùng feature `gene` (feature `gene` = 133, thừa 1 bản `ycf1` cụt; tên gen duy nhất = 112; gen lặp trong IR = 21). Khớp bài công bố Ai và cs. 2019.
- **Phạm vi con số Ts/Tv dùng làm cổng:** **chỉ vùng CDS · 76 gen dùng chung · căn bằng MAFFT** ⇒ 2,686 (317 Ts / 118 Tv, 66.066 cột căn). Số của bộ Eleutherine (1,095–1,777) là **khác phạm vi** (SNP của cả plastome) — không được lấy ra đối chiếu trực tiếp. **Phạm vi này đã chốt 02/10** (chỉ CDS, không tính IR).
- Mỗi lần chạy phải xuất **report công cụ + thông số** (tên công cụ · phiên bản · đường dẫn · tham số/ngưỡng · thời gian) và **chia output theo thư mục** để người dùng truy lại — yêu cầu đã chốt 02/10 (§2). Đây là điểm CGAS đang thiếu (§7.7).
- **Chuẩn hoá bộ gen dùng chung**: xác định điểm bắt đầu/chiều và tách 4 vùng LSC/IRa/SSC/IRb — làm ở giai đoạn 1 rồi mọi giai đoạn dùng lại.
- **Chạy ở đâu**: các giai đoạn 1–5 và 8 chạy trên **máy phổ thông**; giai đoạn 6 (thời gian phân kỳ) và giai đoạn 7 (lắp ráp) cần **máy chủ từ xa** (máy 8 GB không đủ — §9 mục 14).
- **Ca kiểm thử chuẩn**: bộ dữ liệu **Eleutherine bulbosa** + **bộ số vàng** ở §6.
- **Kiểm thử ngay trong từng giai đoạn, không dồn về cuối**: mỗi giai đoạn phải có mục *“Kiểm tra trong giai đoạn”* và chỉ được coi là xong khi số kiểm tra khớp; sai thì **dừng lại sửa ngay trong giai đoạn đó**, không để lỗi trôi sang giai đoạn sau. Giai đoạn 8 chỉ làm phần **kiểm thử chéo**.
- **Chưa chốt (dev tự quyết khi lập kế hoạch với Claude)**: cách tổ chức mô-đun/plugin, định dạng dữ liệu trung gian giữa các bước, và quy trình chạy tự động nhiều bước.

### 13.1. Giai đoạn 1 — Khung chạy & dữ liệu dùng chung

Làm trước tiên: đây là thứ mọi giai đoạn sau phụ thuộc, nếu không làm thì giai đoạn 2 không có cổng nào để thi hành.

**Ba thứ bắt buộc phải có trong khung (chốt khi duyệt plan Giai đoạn 1, 02/10)** — thiếu một thứ là Giai đoạn 2 không có cổng để chạy:

1. Bộ số vàng chạy **một lệnh** ra đủ số §6 (153.724 bp · 132 gen · IR 26.214 · LSC/SSC · 37 SSR) — theo 2 tầng ở §3.
2. **Cổng IR ≥ 99,9 % nằm trong khung ngay từ đầu**, không phải thêm sau.
3. **Quy ước output theo yêu cầu 02/10**: mỗi lần chạy có report công cụ + thông số, thư mục tách theo lần chạy/dataset (điểm CGAS thiếu — §7.7).

- Mô hình dữ liệu lõi (bộ gen · gen · vùng) dùng chung cho mọi phần
- Cơ chế cắm thêm phần (plugin) cho từng mô-đun
- Schema kết quả chuẩn hoá (một kiểu kết quả chung, xuất được nhiều định dạng)
- **Engine sinh report diễn giải dạng câu**: điền số thật vào template câu soạn sẵn (không dùng mô hình ngôn ngữ tự suy diễn), luôn gắn nhãn *“bản nháp — cần người kiểm tra”*
- Log & provenance: thư mục `00_run_info` ghi phiên bản công cụ, tham số, checksum đầu vào (dùng cho Materials & Methods)
- File cấu hình mặc định tập trung (`defaults.yaml`)
- **Citation registry dùng chung**: tên công cụ + bài trích dẫn, để report tự lấy (không gõ tay từng report)
- **Bộ dữ liệu mẫu + bộ kiểm thử vàng**: chạy một lệnh ra đúng bộ số vàng §6 (2 tầng A/B — §3)
- **Reference Panel & tìm bộ gen so sánh** (4 tầng taxonomy) + cơ chế tải/đệm dữ liệu NCBI — engine này **chỉ đề xuất**, người dùng vẫn là bên chốt danh sách genome đưa vào chạy (định nghĩa 4 tầng + 2 chế độ: xem khối ngay dưới)
- Metadata nhóm mẫu dùng chung (nhãn nhóm cho mọi hình/bảng so sánh)
- **Chuẩn hoá điểm bắt đầu & chiều** cho mọi bộ gen (mới lắp ráp lẫn tải về) trước khi vào bất kỳ so sánh nào
- Xử lý gia tăng: chỉ chạy phần mới, không chạy lại từ đầu
- Cách đóng gói dự án (`pyproject.toml` + conda)
- Execution backend cho bước nặng RAM (chạy tại chỗ hay đẩy ra máy chủ từ xa)
- Cân nhắc dùng engine điều phối có sẵn (Snakemake) thay vì tự viết

#### Reference Panel — định nghĩa 4 tầng taxonomy (nguyên văn kế hoạch 22/9, mục 0.7 chế độ B)

| Tầng | Là gì | Dùng cho |
|---|---|---|
| Tầng 0 | các accession khác của cùng loài | C.8 di truyền quần thể · phần "trong loài" của C.6 |
| Tầng 1 | cùng chi | nhánh B · cặp Ka/Ks gần (C.2) · lớp mVISTA gần (C.4) |
| Tầng 2 | cùng họ, khác chi | lớp mVISTA xa (C.4) · ancestral state (C.9) |
| Tầng 3 | cùng bộ, khác họ | đối chứng Ka/Ks tránh bão hoà site đồng nghĩa — **không dùng để gợi ý outgroup** |

Engine có **2 chế độ**, không phải 1:

- **Chế độ A — "lấy hết trong 1 cấp phân loại"** (mặc định cấp họ): lọc RefSeq trước, rồi lọc tay theo 3 tiêu chí — toàn vẹn genome · đầy đủ chú giải · đúng vị trí phân loại. Dùng cho phát sinh loài C.5 và tổng quan cấp họ/chi.
- **Chế độ B — "phân tầng đại diện"**: 4 tầng ở bảng trên.

Điều dễ hiểu sai:
(a) engine có **2 chế độ**, không phải 1 — A và B ở trên;
(b) **outgroup phát sinh loài do người dùng tự chọn hoàn toàn** — tầng 3 không được dùng để gợi ý;
(c) *chưa ghi — đoạn văn bản nguồn bị cụt, cần bổ sung.*

**Kiểm tra trong giai đoạn:** `make test` + `make golden` chạy **một lệnh** ra đủ bộ số vàng §6 · report của lần chạy thử có đủ tên công cụ + phiên bản + tham số · mở thư mục output là biết file nào của việc gì.

### 13.2. Giai đoạn 2 — Đọc và thống kê bộ gen

Dữ liệu vào: file GenBank (kể cả tải từ NCBI). Chưa cần dữ liệu đọc thô.

- Nội dung gen & thống kê (số gen theo nhóm, GC%, độ dài)
- Bảng gen xuất bản (Word)
- Cấu trúc bộ gen: ranh giới LSC/SSC/IR, %GC theo vùng, GC-skew
- Codon usage bias (RSCU, ENC, CAI, CBI, Fop…)
- Thành phần amino acid (kèm GRAVY/Aromaticity)
- SNP / thay thế theo từng vị trí
- Intron: số lượng, độ dài, vị trí
- SSR (vi vệ tinh): phân loại theo vùng và motif
- Đa dạng nucleotide (π) theo gen/intron/vùng liên gen, sliding-window

**Áp dụng cho:** cả **lục lạp** và **ty thể**.

**Ra được:** bảng số + hình cho **2 bản thảo Oldenlandia / Coleus** (bản nháp mục Results).

**Kiểm tra trong giai đoạn:** bộ số vàng §6 — độ dài **153.724 bp** · **132 gen = 86 CDS + 38 tRNA + 8 rRNA** (cách đếm ở điều kiện chung) · IR **26.214 bp** (2 bản giống ≥ 99,9 %) · LSC **83.127** / SSC **18.169** · **37 SSR** (mono 31 · di 4 · tri 2) khi có MISA/cpstools (tầng B — §3). Lệch ⇒ dừng, sửa ngay tại đây.

### 13.3. Giai đoạn 3 — Chú giải bộ gen

Chạy được trên bộ gen đã có sẵn.

- Chú giải: PGA (lục lạp), MITOS2/GeSeq (ty thể)
- Kiểm tra & chuẩn hoá tên gen
- Kiểm tra khung đọc (start/stop/internal codon)
- Chuẩn bị hồ sơ nộp NCBI (FASTA/TBL/SQN)

**Kiểm tra trong giai đoạn:** đọc lại file kết quả bằng chính bộ đọc GenBank (không hỏng định dạng) · kiểm tra khung đọc · đếm lại gen theo cách đếm đã chốt (132). *Riêng bộ gen tự lắp ráp* (chỉ có ở giai đoạn 7): phải đối chiếu họ hàng gần hoặc read bắc qua ranh giới LSC–IR / SSC–IR, nếu không thì dán nhãn *chưa xác định* — ở giai đoạn 3 phần này bỏ qua.

### 13.4. Giai đoạn 4 — So sánh nhiều bộ gen

- Áp lực chọn lọc Ka/Ks
- Synteny/collinearity kiểu mVISTA: % tương đồng theo toạ độ
- Ứng viên DNA barcode & vùng biến động cao
- Di truyền quần thể kiểu DnaSP (Hd, θ, Tajima's D, Fu & Li, Fst…)
- Dự đoán RNA editing (cần dữ liệu RNA)

**Áp dụng cho:** cả **lục lạp** và **ty thể**.

**Ra được:** nội dung so sánh cho 2 bản thảo.

**Kiểm tra trong giai đoạn:** Ts/Tv = **2,686** cho cặp *Iris domestica* × *Iris gatesii* (**chỉ CDS · 76 gen · MAFFT** — đúng phạm vi ở điều kiện chung) · biểu đồ tương đồng phải ra **21 khối**, %IDY 83,7–98,2 (mean 91,5). Số không khớp ⇒ dừng, kiểm lại dữ liệu vào và tham số.

### 13.5. Giai đoạn 5 — Hình và báo cáo

- Thống nhất theme màu/font; các loại hình: bản đồ plastome vòng · RSCU heatmap · π sliding window (có tô nền LSC/SSC/IR) · sơ đồ IR junction · mVISTA-style %IDY + synteny tuyến tính · phân bố SSR/repeat theo vùng & motif · cây có annotation · NUPT/MTPT theo vị trí gen · mạng haplotype · cấu hình đa hình ty thể · **đánh giá DNA barcode ứng viên**; xuất vector (SVG/PDF/EPS); **quy tắc mở rộng hình theo số bộ gen**
- Thống kê: mô tả, kiểm định phi tham số, tương quan, hiệu chỉnh đa kiểm định (FDR/Bonferroni), diễn giải dạng câu
- Xuất bản: Excel tổng hợp đa sheet · Word (bảng gen/SSR/thống kê) · gói Supplementary Materials · báo cáo PDF 1 trang/bộ gen · báo cáo QC kỹ thuật khách quan
- **Yêu cầu 02/10:** mỗi lần chạy xuất **report công cụ + thông số**, và **output chia theo thư mục** để người dùng truy lại (§2)

**Kiểm tra trong giai đoạn:** sinh lại được bảng/hình từ dữ liệu trung gian (không làm tay) · **kiểm thử truy vấn số**: mỗi số trong report phải trỏ được về **đúng một ô** trong bảng gốc (nếu không, lỗi “số trông hợp lý nhưng sai” vẫn lọt) · report công cụ + thông số có đủ tên, phiên bản, đường dẫn, ngưỡng.

### 13.6. Giai đoạn 6 — Nội dung bài báo công cụ

- NUPT/MTPT — chuyển gen giữa nhân, lục lạp, ty thể (trọng tâm)
- Cấu trúc đa hình ty thể + tái tổ hợp qua trình tự lặp (Chimeric ORF/CMS)
- So sánh song song lục lạp – ty thể trong cùng loài
- Ancestral state reconstruction
- Ước tính thời gian phân kỳ (rất nặng; cần máy chủ)

**Kiểm tra trong giai đoạn:** phần cấu hình ty thể **có mốc công bố để đối chiếu ngay tại đây** — các mitogenome có cấu trúc phức tạp đã công bố: *Damnacanthus indicus* (qua trình tự lặp trực tiếp), *Salvia* (S. officinalis 2023 · S. splendens 2025 · S. miltiorrhiza 2022), *Hemiboea*, *Coffea arabica* (2023) — đều là mốc trong kế hoạch 22/9. Riêng **NUPT/MTPT chưa có số công bố** ⇒ mốc kiểm duy nhất hiện có là đối chiếu với bài Oldenlandia/Coleus, thực hiện ở giai đoạn 8. Giai đoạn này **chạy trên dữ liệu đã công bố / GenBank**; phần cần đồ thị lắp ráp thì chuyển sang giai đoạn 7.

### 13.7. Giai đoạn 7 — Lắp ráp (đã chuyển xuống sau)

Cần dữ liệu đọc thô; nặng về bộ nhớ.

- Lắp ráp + QC hàng loạt (fastp → GetOrganelle)
- Kiểm tra đồ thị assembly (Bandage)
- Cấp dữ liệu cho phần cấu hình ty thể ở giai đoạn 6
- **Đường chạy từ xa** cho máy yếu: usegalaxy.eu (`get_organelle_from_reads`) — cắt bớt read (subsample) trước khi upload, **pin phiên bản công cụ + ghi lại history ID**
- **Ca read Illumina thật**: hiện **chưa có** dữ liệu read của đề tài (đang nợ) ⇒ mới chỉ có cổng dựa trên read mô phỏng

**Kiểm tra trong giai đoạn:** so với kết quả đã chạy thật — GetOrganelle + SPAdes cho **127.679 bp** từ 20.000 cặp read mô phỏng · OATK: 126.660 bp, identity **100,00 %** · Bandage ra PNG (103.289 B). Khi có read thật: bổ sung mục đo **đỉnh RAM**.

### 13.8. Giai đoạn 8 — Kiểm thử chéo, hiệu năng và tự động hoá

**Không phải nơi test lần đầu.** Mỗi giai đoạn đã tự kiểm bằng mục *“Kiểm tra trong giai đoạn”*. Giai đoạn này chỉ làm phần không thể làm sớm hơn:

- **Đối chiếu chéo với công cụ khác** trên cùng bộ dữ liệu: CGAS · CPStools · CPGview
- **Đối chiếu phần NUPT/MTPT** với bài Oldenlandia/Coleus (mốc duy nhất hiện có cho phần này)
- **Kiểm thử batch lớn** (10–50+ bộ gen) để biết giới hạn bộ nhớ và thời gian
- **Tự động hoá**: CI chạy test mỗi khi sửa code; bộ vàng §6 chạy lại tự động

### 13.9. Giai đoạn 9 — Đóng gói, giấy phép & công bố

- README / hướng dẫn dùng nhanh / tài liệu sử dụng
- Đóng gói pip/conda + Dockerfile (Windows chạy qua Docker)
- **Giấy phép**: `LICENSE.md` + rà soát công khai các công cụ có ràng buộc (PGA = GPL-3.0; Vmatch chỉ miễn phí học thuật; **MISA = không nêu rõ** — zip không kèm giấy phép, chỉ có trích dẫn Beier và cs. 2017) và cách ly nếu cần — làm **trước khi phát hành**
- **Bảng mô-đun × nền tảng × RAM** + tuyên bố **3 mức hỗ trợ**: Linux native (chính thức) · Windows qua Docker (chính thức) · macOS (thực nghiệm, chưa kiểm)
- **Tuyên bố “thứ không làm”** (không Windows native; không lắp ráp dưới ngưỡng RAM; không dự đoán RNA editing nếu không có dữ liệu RNA-seq) — để chặn kỳ vọng sai
- Video hướng dẫn ngắn · bản thảo công bố công cụ · DOI (Zenodo) + CITATION.cff
- Giao diện web: để ngỏ, làm sau khi lõi ổn định

**Kiểm tra trong giai đoạn:** chạy lại bộ số vàng **trong Docker trên Windows** — cùng dữ liệu vào phải ra cùng số (3 môi trường: WSL · CI · Docker).


---

### 13.10. Đã chạy thử được gì (chỉ ghi chỗ có bằng chứng — nguồn §7.1, §8)

| Việc | Kết quả thô |
|---|---|
| Reading/thống kê gen | 132 gen = 86 CDS + 38 tRNA + 8 rRNA trên NC_050833.1 |
| Cấu trúc bộ gen | IR 26.214 bp (2 bản giống 100,00 %) · LSC 83.127 · SSC 18.169 · xoay chiều giữ nguyên IR |
| SNP/thay thế | Ts/Tv = 2,686 (76 gen CDS, MAFFT); bộ Eleutherine (cả plastome) có dải 1,095–1,777 |
| SSR | 37 SSR — MISA standalone 2.1 **và** cpstools 3.0 cùng ra 37; mono 31 · di 4 · tri 2; IGS 18 · Intron 8 · Gene 11 |
| Synteny/mVISTA | nucmer → 21 khối, %IDY 83,7–98,2 (mean 91,5); PNG 243.787 B |
| RNA editing | mô phỏng 20 vị trí C→U: recall 100 % ngoài IR, precision 74 %; 6 vị trí sót trong IR |
| Mạng haplotype | PyPopART trên 13 trình tự matK → rc=0 |
| Bản đồ vòng genome | pycirclize 1.10.1 → PNG 494.859 B @300 dpi |
| Chú giải | plann.pl + PGA.pl chạy offline theo reference GenBank |
| Lắp ráp (giai đoạn 6) | GetOrganelle + SPAdes: 20.000 cặp read mô phỏng → 127.679 bp · raven 1.8.3: 2 contig, dài nhất 107.394 bp · OATK 1.0: 126.660 bp, identity 100,00 % · flye: chết vì hết RAM (máy 3,7 GB) |
| Đồ thị assembly | Bandage 0.9.0 → PNG 103.289 B |
| Thời gian phân kỳ | mcmctree (PAML 4.10.10) + BEAST2 2.7.7 chạy được, **chưa** chạy phân tích thật; treePL đã build nlopt + ADOL-C |
| **Chưa chạy thử lần nào** | Bảng gen Word · codon usage · amino acid · intron · π · Ka/Ks · DNA barcode · di truyền quần thể · NUPT/MTPT · cấu trúc ty thể · ancestral state · thời gian phân kỳ (phân tích thật) · kiểm tra khung đọc · nộp NCBI · thống kê · xuất bản · kiểm thử · đóng gói |

### 13.11. Bảng tra để đối chiếu với kế hoạch 22/9 (chỉ dùng khi cần tra ngược tài liệu gốc)

| Giai đoạn mới | Mã trong kế hoạch 22/9 |
|---|---|
| 1 — Khung chạy & dữ liệu dùng chung | 0.1 · 0.2 · 0.3 · 0.3b · 0.4 · 0.5 · 0.6 · 0.7 · 0.8 · 0.9 · 0.10 · 0.11 · 0.12 · 0.13 · 0.14 |
| 2 — Đọc và thống kê bộ gen | B.1 · B.2 · B.3 · B.4 · B.5 · B.6 · B.7 · B.8 · B.9 |
| 3 — Chú giải bộ gen | A.2 · A.3 · A.4 · A.5 |
| 4 — So sánh nhiều bộ gen | C.2 · C.4 · C.6 · C.8 · C.3 |
| 5 — Hình và báo cáo | E.1–E.14 · D.1–D.5 · F.1–F.5 |
| 6 — Nội dung bài báo công cụ | C.1 · C.10 · C.7 · C.9 · C.5 · C.5b |
| 7 — Lắp ráp | A.1 · A.1b |
| 8 — Kiểm thử chéo & tự động hoá | G.1–G.5 |
| 9 — Đóng gói, giấy phép & công bố | H.1–H.6 |

*Ghi chú: kế hoạch 22/9 định nghĩa đủ các mã 0.1–0.14, E.1–E.14, D.1–D.5, G.1–G.5 (đã kiểm bằng cách đếm mã trong file kế hoạch). Bảng công cụ 17/09 chỉ định nghĩa tới E.13, D.1–D.3 và không có G.3 — đây là bảng cũ hơn, dùng bảng 22/9 làm chuẩn.*
