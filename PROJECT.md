# PlastidToolSuite — ngữ cảnh dự án

> Nguồn ngữ cảnh duy nhất của dự án. Đọc trước khi làm bất cứ việc gì.
> Cập nhật: 07/10/2026. Toàn bộ nội dung là **dữ kiện đã kiểm tra / đã chạy thật**.
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
| 10/10/2026 | **Vòng làm việc này CHỈ làm lục lạp — phần ty thể hoãn có tên, không phải "đang chờ dữ liệu" (dev quyết, nhắc lại lần thứ ba nên chốt thành luật phạm vi)**: D5-07 đọc như một *thứ tự* ("lục lạp trước") và vì thế các phiên sau vẫn quay lại bàn việc ty thể. Nay là **phạm vi**: mọi mô-đun · lần chạy · bản thảo · kiểm chứng chéo · bước phát hành của vòng này phục vụ **plastome**. **Hoãn có tên**: M6.1 · M6.2 · M6.3 (cả Pha 6) · M7.3 · M8.2 — đánh dấu `HOÃN` trong `ROADMAP.md`, không xuất hiện trên board nào, không nằm trong bảng "waiting" (chúng không chờ ai cả). **Giữ nguyên**: mã ty thể đã có (`mito_annotate` · `mito_genes` · core · panel · test) — không xoá, không mở rộng, `make check` vẫn chạy; và **cấu trúc** một-bào-quan-một-lần-chạy + ba nhóm trong registry (D5-06) chính là chỗ đã chừa để mở lại mà không phải đập lại | Dev | §13 "Hai bào quan" · CLAUDE.md → "Priority today" · ROADMAP (Pha 6 · M7.3 · M8.2). **Hệ quả nói thẳng:** bài báo công cụ (deliverable 2 — NUPT/MTPT + cấu trúc mito) **đi cùng Pha 6**, vì ba mô-đun ấy *chính là* nội dung của nó; vòng này còn lại hai bản thảo loài và bản phát hành toolkit |
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
| 03/10 | **Backend (P1-12, dev quyết)**: giai đoạn 1–5 và 8 chạy **tại chỗ** trên máy phổ thông; giai đoạn 6 (thời gian phân kỳ) + 7 (lắp ráp) chạy **máy chủ từ xa** — máy 8 GB chỉ còn 1,5–2,5 GB trống mỗi việc (§9.28). Ghi lại theo luật phân quyết định trong CLAUDE.md. **[ĐÃ THAY 08/10/2026 — riêng giai đoạn 7 (lắp ráp): nay tuỳ máy người dùng, xem dòng 08/10 bên dưới; giai đoạn 6 vẫn chạy máy chủ từ xa]** |
| 03/10 | **Bảng gen xuất bản (P2-09, dev quyết)**: bảng in ra bài báo là bảng **phân nhóm theo chức năng**, không phải bản dump từng feature; nhóm gen là **dữ liệu** (`src/plastid/data/gene_groups.yaml`), không phải code. tRNA/rRNA khớp theo **tiền tố tên** (`trn`, `rrn`) chứ không theo danh sách anticodon — bản Oldenlandia/Coleus thật sẽ mang anticodon khác, danh sách cứng sẽ **rơi im lặng**. Gen không nhóm nào nhận ⇒ một dòng *unclassified* + ghi tên trong notes, **không bao giờ bị bỏ** |
| 03/10 | **Số đếm của bảng gen = số bản sao** (kể cả 2 bản IR), đúng quy ước 132 = 86 + 38 + 8, và bảng **tự đối chiếu với `count_genes()`** trong một bảng riêng — cùng kiểu đối chiếu đã bắt lỗi MISA ở P2-06. Một gen rơi giữa bộ đếm và bảng sẽ hiện ra thành *DISAGREES*, không phải thành một số nhỏ hơn không ai kiểm |
| 03/10 | **Xuất Word (D2-04 + P2-09, dev quyết)**: `python-docx` là **extra tuỳ chọn** (`pip install 'plastid-suite[docx]'`), import **bên trong hàm** nên bản cài tối thiểu vẫn thuần Python; thiếu nó ⇒ `ExportDependencyError` + exit 3, không phải traceback ImportError. `openpyxl` **chưa thêm** — chưa có gì dùng (YAGNI); sẽ vào cùng extra khi làm Excel ở Giai đoạn 5 |
| 03/10 | **Tài liệu Word phải tái lập được từng byte**: `python-docx` đóng dấu thời gian đồng hồ vào zip, nên gói được **dựng lại** với timestamp cố định (1980-01-01), entry sắp xếp, `create_system` ghim ⇒ cùng đầu vào ra cùng file. Đây là định nghĩa "Done when" của P2-09, không phải một lời hứa |
| 03/10 | **Gói .docx được tỉa**: template của `python-docx` mang theo cả catalogue style của Word (164 style, ~350 kB) + `stylesWithEffects.xml` (438 kB) + thumbnail + `customXml`; bản xuất chỉ giữ style tài liệu **thực sự dùng** cùng tổ tiên và default của chúng, mọi quan hệ trỏ tới phần đã bỏ đều bị cắt. Đo được: ~800 kB → **13 kB**, có test cấu trúc canh |
| 03/10 | **[ĐÃ THAY — xem dòng dưới] Tên gen: ánh xạ theo locus, giữ tên gốc (D3-02, dev quyết có bằng chứng)**: `pafI`·`pafII`·`pbf1`·`clpP1` trong `PX923121.1` là **cùng locus** với `ycf3`·`ycf4`·`psbN`·`clpP` — chứng minh bằng synteny + độ dài + `/product` trên 9 bản ghi (§8.9). Tên **chuẩn** giữ bản cổ điển (7/9 bản ghi + genome vàng); tên gốc đi kèm `source_name`, không bao giờ bị bỏ. Tên **in trong bản thảo** và **nộp NCBI** vẫn là quyết định của PI (D3-06). *Chiều chuẩn hoá ở dòng này bị đảo ở dòng dưới, sau khi biết hai bản ghi là của nhóm (§8.10)* |
| 03/10 | **[ĐÃ THAY 07/10/2026 — xem dòng D5-09 bên dưới] Tên gen: khoá theo locus, in cả hai tên (D3-02, dev quyết có bằng chứng)**: bốn tên `pafI`·`pafII`·`pbf1`·`clpP1` là **cùng locus** với `ycf3`·`ycf4`·`psbN`·`clpP` (synteny + độ dài + `/product`, §8.9). **Khoá so sánh nội bộ = phương ngữ của nhóm** — vì `PX923121.1` là bản ghi **của nhóm** (Vu,H.-T. = cô Trang) và **đúng loài bản thảo**, và GeSeq sẽ sinh ra đúng bộ tên đó; tên cổ điển đi kèm `source_name`. **In ra thì in cả hai** (`pafI (ycf3)`) — cách các bảng gen công bố đang làm. `gene_groups.yaml` nhận thêm bốn tên vào đúng nhóm chức năng (§8.10) |
| 03/10 | **[ĐÃ THAY 03/10 — xem dòng "Chloe là annotator" bên dưới] Chú giải: công cụ offline là pipeline, GeSeq là bước kiểm tra một lần rồi ghim (D3-01, dev quyết)**: bộ công cụ **phải chạy offline và đóng gói bằng Docker**, mỗi lần chạy phải ghi **tên · phiên bản · đường dẫn · argv** — công cụ web không ghi được cái nào trong số đó, nên đây là **ràng buộc, không phải sở thích**. Chọn **PGA** (GPL-3.0, batch, ra GenBank; ba cái bẫy của nó — thư mục ra phải rỗng, tên file có dấu cách, multi-FASTA — repo này đã tránh sẵn), plann.pl dự phòng. GeSeq chạy tay **một lần cho mỗi loài bản thảo** rồi **ghim** (`panel pin`) ⇒ bước thủ công chỉ xảy ra một lần, các lần sau offline và tái lập được. Kiểm 03/10: GeSeq không có bản chạy rời | Dev | `docs/D3-01_annotator.md` |
| 03/10 | **Đa genome (D2-01, dev quyết)**: giữ nguyên contract v0.1 `run(genome, params, ctx)`; thêm contract thứ hai `run(genomes, params, ctx)` chọn bằng `arity` khai báo lúc đăng ký — **không** nới rộng contract cũ, vì module một-genome sẽ phải nhận một dãy nó không dùng (seam 5). Runner cưỡng chế số genome trước khi chạy. Ghi ở `docs/API_v0.2.md` |
| 03/10 | **π: cửa sổ 600 bp · bước 200 bp (D2-02)** — quy ước kiểu mVISTA; window/step vào params, manifest và report |
| 03/10 | **Codon usage & amino acid dùng NCBI bảng 1 (D2-03, dev quyết)** — bảng 1 và bảng 11 (bảng plastome khai) **chia sẻ toàn bộ axit amin**, chỉ khác mã khởi đầu (D3-03), nên một bảng mã trả lời cả hai. **Quy ước kèm theo: mã khởi đầu được đọc như axit amin của chính nó** — genome vàng có 5/86 CDS mở đầu bằng GTG/ACG/ATC (§9.J), nên bảng in ra **phải ghi rõ quy ước và bảng mà bản ghi khai**, không ngụ ý dùng bảng của bản ghi. Bảng mã biên dịch trong `core/codon.py`; con số khoá ở `src/plastid/data/defaults.yaml` | Dev | §9.J |
| 03/10 | **ORF: đọc bảng mã của chính bản ghi, từ chối thay vì đoán (D3-03 + P3-01)**: `Gene` mang `/transl_table` và `/codon_start` của bản ghi; bảng mã ngoài {1, 11} ⇒ `OrfError` (hai bảng này chia sẻ toàn bộ axit amin, chỉ khác mã khởi đầu), `/codon_start` ≠ 1 ⇒ `OrfError` — dịch sai khung là kiểu "sai mà không báo". **ACG được nhận là mã khởi đầu nhưng xếp lớp riêng `non_canonical`**: nó *không* có trong danh sách của bảng 11, nhưng bản ghi đã công bố dùng nó cho rpl2/ndhD (RNA editing ACG→AUG); xếp vào `no_start` sẽ là báo động giả trên 2/86 CDS của genome vàng (§8.11) |
| 03/10 | **Lớp "annotated, not coding" khớp cả TÊN lẫn VỊ TRÍ (P3-01)**: ycf1 cụt chồng lấn 21 bp với CDS của `ndhF`, nên khớp theo vị trí không thôi thì ndhF "phủ" nó, còn khớp theo tên không thôi thì CDS ycf1 thật ở chỗ khác "phủ" nó — mỗi cách hỏng một kiểu im lặng (§8.11) |
| 03/10 | **File nộp NCBI: `.fsa` + `.tbl` thuần Python, luôn luôn; `.sqn` chỉ khi có `table2asn` *và* có template `.sbt` (D3-04, P3-05)**: web submission của NCBI nhận thẳng cặp `.fsa`+`.tbl`, nên bản cài tối thiểu vẫn nộp được. `.tbl` mang **tên gốc của bản ghi** (đổi tên là đổi nội dung nộp — D3-06 chỉ nói về *chữ in trong bản thảo*); toạ độ partial đi theo bản ghi, partial ở interval giữa ⇒ **từ chối**, không bỏ im lặng (§8.13) |
| 03/10 | **[ĐÃ THAY 05/10 — xem dòng D5-08; và ĐÃ THAY LẦN NỮA 10/10/2026 — xem dòng D5-08 · §8.56: cả hai annotator ngoài đã nghỉ hưu theo P3B-11] Chú giải: Chloe là annotator của pipeline; PGA còn lại cho genome ngoài angiosperm (D3-01 — sửa lại sau khi ĐO, dev quyết)**: ràng buộc cũ vẫn nguyên (offline, ghi được tên · phiên bản · đường dẫn · argv), nhưng *"công cụ offline nào"* chỉ là **giả thuyết cho tới khi chạy thật**. PGA trên genome vàng ra **147 gen** (7 reference) và **130** (1 reference) so với **132**, kèm ba kiểu hỏng im lặng (§8.15) ⇒ **thay bằng Chloe**: offline, có bộ reference chuẩn hoá của chính nó, và **đo qua reader của mình ra đúng 132 = 86 + 38 + 8, 0 CDS lỗi ORF** (§8.17). Giữ `--set annotator=pga` cho plastome không phải thực vật có hoa (Chloe chỉ làm angiosperm). GeSeq vẫn là kiểm tra chéo ghim một lần. **Bài học ghi lại:** công cụ thắng không nằm trong danh sách khảo sát ban đầu | Dev | `docs/D3-01_annotator.md`, `docs/tools/chloe.md` |
| 03/10 | **[ĐÃ THAY 03/10 (cùng ngày) — xem dòng "Bộ gọi gen mito tự viết" bên dưới] Nối công cụ chú giải mito (P3-04) — PMGA**, không phải MITOS2 (annotator của động vật có vú, bộ tham chiếu `refseq81m`). PMGA phân phối dạng **singularity container** (Figshare 27201798) + web server, không có repo git và không có tài liệu CLI công khai ⇒ **chưa viết adapter**: viết theo CLI tự nghĩ ra sẽ là một placeholder trông như thật. Cần: kéo container, ghi entry point vào `docs/` |
| 03/10 | **Bước nặng: chạy ở đâu, và kết quả chạy ở nơi khác phải ghim được (P1-12, làm nốt 03/10)**: luật backend đã có từ 03/10 (giai đoạn 1–5 và 8 chạy tại chỗ; 6–7 chạy máy chủ từ xa) nhưng **thiếu đường thực hiện** cho công cụ chạy bên ngoài: `panel pin` chỉ nhận accession. Bổ sung `plastid panel pin --file … --note …` (§8.19) để kết quả của web service / máy chủ khác vào được panel với đủ checksum + đọc lại + provenance — PMGA là ca đầu tiên dùng nó | Dev | `docs/tools/pmga.md`, §8.19 |
| 03/10 | **[ĐÃ THAY 03/10 (cùng ngày) — xem dòng "Bộ gọi gen mito tự viết" bên dưới: PMGA không chạy nổi ở đây, đường chú giải mito là bộ gọi của dự án; dòng này giữ làm bằng chứng đo] Chú giải mito: PMGA chạy trên web service rồi GHIM, không chạy local (dev quyết, xác nhận bởi anh Trí)**: máy này **không đủ mạnh để chạy annotate** — container PMGA 3,18 GB nén → **11 GB giải nén, 138.045 file**, cần Apptainer và bên trong chạy MAKER; trần WSL 8 GB, ổ C: thật chỉ còn 42 GB. Container đã được tải về, đọc lấy CLI thật, rồi **xoá** (14 GB = một phần ba chỗ trống). Chú giải mito đi theo **đường thứ hai của chính luật P1-12**: chạy PMGA trên web service một lần mỗi loài ⇒ `plastid panel pin --file … --note …` (§8.19) ghim kết quả với đủ checksum + đọc lại ⇒ `mito_genes` đọc offline. **`external/pmga.py` để mở có điều kiện:** viết khi có máy chạy nổi container, và số đầu tiên phải đo khi đó là RAM đỉnh + thời gian | Dev | `docs/tools/pmga.md` (P3-04; nay H-05) |
| 03/10 | **Dữ liệu thay thế cho dev (dev quyết, cô Trang chưa có mẫu)**: khi mẫu thật chưa có, **được tải bản công khai tương tự** (cùng chi, hoặc họ hàng gần nhất) để dev/test — chỉ qua `plastid panel pin` (accession · ngày · sha256), **mọi đầu ra phải ghi nhãn**, **không bao giờ** để số của bản thay thế vào bảng manuscript. Định nghĩa bộ dữ liệu commit ở `data/`; payload `.gb` không commit. Quy tắc ở CLAUDE.md → "Data — real samples vs development stand-ins"; bộ đầu tiên ở §8.6 |
| 03/10 | **Bộ gọi gen mito tự viết: ranh giới lấy từ alignment, xếp hạng theo identity, và mọi thứ bị loại đều được kể ra (D3A-02 + P3A-03, dev quyết có đo)** — PMGA không chạy nổi ở đây nên Phase 3A tự viết phần chuyển chú giải: `blastx` (protein tham chiếu) + `blastn` (rRNA) qua `run_tool`, rồi `core/mito_caller.py` biến HSP thành gen. Bốn luật, mỗi luật có số đo đỡ lưng (§8.22): (1) **ranh giới từ alignment**, không bao giờ từ ORF-walking — codon stop chỉ được lấy khi chính chuỗi xác nhận, đọc **đúng một** codon, không đi tìm; (2) **chọn tham chiếu theo identity cao nhất** trong số các ứng viên qua được ngưỡng — đo trên `atp6`: xếp theo bitscore kéo gen lệch **381 bp** về thượng nguồn vì một alignment dài 418 aa chỉ khớp 75,6 %, trong khi bản khớp 93,4 % nằm đúng chỗ; tổng sai lệch toạ độ 310 bp so với 655 bp; (3) **ngưỡng 50 % identity / 50 % coverage**, in ra trong mọi bảng và trong `defaults.yaml`; (4) **mọi HSP không dùng đều được liệt kê kèm lý do** — dưới ngưỡng, trùng khung đọc, hay là **bản sao thứ hai** của một phần gen. Kèm theo: `-seg no` (mask làm atp9 mất một nửa điểm), **không** dùng `-parse_seqids` (nó **lỗi** trên chính bộ tham chiếu), và gen tìm thấy mà vocabulary chưa có tên thì **báo cáo chứ không đặt tên** (H-06) | Dev | §8.22, `docs/tools/blast.md` |
| 03/10 | **Ba quyết định phụ của Giai đoạn 3A — ghi lại đây để các mã được trích trong code có chỗ tra (D3A-01 · D3A-04 · D3A-05)**: **D3A-01** = bộ gọi tRNA dùng **ARAGORN 1.2.41** (bioconda, GPLv3; chính là công cụ PMGA chạy bên trong) thay vì tự viết rồi chỉnh cho ra đúng 17 — đó sẽ là overfitting (§8.23); **D3A-04** = call đa exon chỉ báo **các exon mà alignment chứng minh được**, không bao giờ tự nhận một join (gen trans-spliced hiện ra đúng các phần của nó); **D3A-05** = phạm vi prototype: **chỉ chú giải + nộp** — không tìm gen mới, không mô hình intron, không dự đoán RNA editing | Dev | docstring `core/mito_caller.py` · `modules/mito_annotate.py`; §8.23–8.24 |
| 04/10 | **Synteny/mVISTA: nucmer (MUMmer4) pin 4.0.1, tham số truyền tường minh (D4-02 + D4-04, dev quyết có đo)**: đường ống đúng như §7 ghi (`nucmer --maxmatch → delta-filter -1 → show-coords -rcl`) qua `run_tool`, đủ name · version · argv trong manifest. `-l 20 -c 65` là mặc định của MUMmer nhưng **truyền tường minh** và để trong `src/plastid/data/defaults.yaml` — số khối đổi theo chúng (12 khối thô ở mặc định, 22 ở `-c 40`), một mặc định là con số không ai viết ra; `-t 1` để số luồng không đổi được kết quả. Check của pha *"21 khối, %IDY 83,7–98,2 (mean 91,5)"* **không tái lập được** từ đầu vào của probe (quét đầy đủ ở §8.26): cặp cùng chi thẳng hàng 98–99 %, dải số đó thuộc một so sánh khác chi và đến từ dòng khảo sát §7.1 đo trên dữ liệu khác ⇒ check thay bằng số **đo được** cho cặp Iris: **7 khối · 89,66–99,11 · mean 97,37 · coverage 99,09 %**, kèm hệ chuẩn hoá. Phát hiện ở §9.H | Dev | §8.26, `docs/tools/mummer4.md` |
| 04/10 | **Ka/Ks: PAML `yn00` (D4-01, dev quyết)**: estimator chuẩn pairwise (Yang & Nielsen 2000), **đã pin sẵn** trong `environment.yml` (`paml=4.10.10`, dùng cho `mcmctree` ở giai đoạn 6 — cùng gói, khác chương trình). Loại: NG86 tự viết (chệch, và một bản cài sai sẽ là "sai mà không báo" không ai kiểm được) và KaKs_Calculator2 (công cụ thứ hai cho cùng một việc). Alignment **theo từng gen trên chính chuỗi CDS** (MAFFT), cắt stop cuối, **kiểm khung ba-ngôi trước khi gọi tool** vì yn00 không kiểm; `dN`/`dS` đọc từ ma trận và đối chiếu chéo với dòng Yang & Nielsen của output. Pooled = trung bình có trọng số theo số site của chính yn00, không phải trung bình các tỉ số | Dev | §8.27, `docs/tools/paml.md` |
| 04/10 | **Di truyền quần thể: tự viết bằng Python thuần (D4-03, dev quyết có đo)**: DnaSP là GUI Windows (không script được, không vào Docker), pixy/libsequence đòi VCF hoặc dữ liệu read; các estimator là công thức đóng đã công bố, kiểm được bằng vector giải tích, và bản cài tối thiểu vẫn thuần Python (CLAUDE.md). Bằng chứng: Tajima's D · π · θ(S) khớp **scikit-allel 1.3.13** đến 6 chữ số trên 4 vector; Fu & Li D*/F* khoá bằng vector tính tay (EggLib build trên máy này hỏng, `gestimator` không có trên bioconda — đo và ghi lại). **Fst = Nei's Gst**: bản Hudson viết trước cho **−1** với hai nhóm cùng thành phần ở n nhỏ ⇒ bỏ, có ghi lý do | Dev | §8.28 |
| 04/10 | **Nhãn dữ liệu thay thế đi vào MỌI bảng (`--label`, dev quyết)**: một bộ dữ liệu thay thế không bao giờ được đi mà không có nhãn (CLAUDE.md). Nhãn được áp **một lần** trong runner — nơi duy nhất biến `Result` thành file — bằng cách ghép tiền tố vào `scope.label` của từng bảng, nên nó có mặt ở header CSV, JSON và report mà không module nào phải nhớ. Không có nhãn thì không đổi gì | Dev | §8.30, `docs/P4-05_comparison_tables.md` |
| 04/10 | **Nửa mito của pha 4 là bản phác khác chi (D4-05, dev quyết)**: hai mitogenome từ hai chi khác nhau, một trong hai là **chú giải nguyên mẫu của chính dự án** (§8.24) ⇒ chạy `kaks` + `popgen`, dán nhãn *cross-genus sketch*, đọc như minh hoạ máy móc chứ không phải kết quả về loài. Một so sánh mito đáng vào bản thảo phải đợi bản ghi thứ hai cùng chi | Dev | §8.30 |
| 04/10 | **Bộ dữ liệu anh Điềm gửi là BỘ TEST — không phải stand-in, cũng không phải mẫu bản thảo (D4a-01, dev quyết)**: 34 plastome họ Iridaceae × 2 bản chú giải (công bố NCBI · PGA) + 34 chuỗi trần, đặt tại `.data/thanhdiem_data/` (**trong cây làm việc, đã ignore git**; 136 file · 28 MB · sha256 khớp từng file với `D:\business\pts\thanhdiem_data`). Không loài nào là loài bản thảo ⇒ **không đầu ra nào từ đây vào bảng bản thảo**; nhưng `PX516262.1` *E. bulbosa* **là bản ghi của chính nhóm** ⇒ nhãn theo **từng bản ghi**, không theo cả bộ. Đây là bộ dữ liệu đầu tiên không do dự án tự chọn, và nó mang sẵn ba dạng đáp án: bản công bố · PGA · chính chuỗi | Dev | §8.31 |
| 04/10 | **Ghim bộ test: mỗi bộ chú giải một thư mục panel riêng, `.fasta` đi đường manifest (D4a-03, dev quyết có đo)**: `panel pin --file` khoá bản ghi theo **id của chính bản ghi**, mà PGA đặt id theo voucher (`Aristea_ecklonii_voucher_KK2022_1_NC_072561`) còn NCBI đặt theo accession (`NC_072561.1`) ⇒ ghim chung một thư mục thì **một genome nằm hai lần dưới hai tên, không có gì ghép cặp**; tách `.panel/iridaceae-34/` và `.panel/iridaceae-34-pga/`, cặp được nhận ra bằng **tên file dùng chung**. Chuỗi trần **không** đi được đường GenBank: đo trên một `.fasta` ⇒ pin theo tên file, `parse_error: no GenBank records found`, CLI in UNVERIFIED và **exit 1**; thay bằng `fasta.json` (68 mục: sha256 · bytes · độ dài · sha256 của bản ghi đã pin) — và mỗi chuỗi FASTA **được kiểm là trùng khít** chuỗi của bản ghi đã pin, nên bộ FASTA không mang thêm thông tin nào ngoài cái panel đã có | Dev | §8.31 |
| 04/10 | **`ndhI`/`ndhL`: alias ĐÃ CÓ, hôm nay chỉ đo lại trên 34 bản ghi (D4a-05)**: `src/plastid/data/gene_aliases.yaml` đã mang mục `ndhI` ← `ndhL` từ **03/10** (§8.17, đo khi Chloe chú giải genome vàng: cùng 543 bp, cùng vị trí, và phức hệ ndh của lục lạp chỉ có tiểu đơn vị A–K). Đo lại trên bộ dữ liệu anh Điềm: locus này được viết khác nhau ở **8/34** cặp, và alias resolve **cả hai cách viết về cùng một khoá chuẩn** (`ndhI`) — kiểm bằng code, không phải bằng mắt. **Không thêm mục mới, không có quyết định mới**; việc của P4a-02 là một test canh nó. **[ĐÃ THAY 07/10/2026 — D5-09: chiều `canonical` nay theo bản ghi ★, khoá chuẩn không cố định là `ndhI`; phép đo 04/10 ở dòng này không đổi]** | — | §8.31, `src/plastid/data/gene_aliases.yaml` (đã có) |
| 04/10 | **Bản công bố (bộ 01) là ĐÁP ÁN cho bài kiểm annotator; PGA là ý kiến thứ ba (D4a-04, dev quyết)**: bộ 01 là bản NCBI đã công bố, dùng làm chuẩn để đối chiếu; PGA đã bị thay ở §8.15/§8.16 sau khi đo **147 gen** trên genome vàng (so với 132) ⇒ chỉ dùng làm đối tượng so sánh, **không bao giờ làm chuẩn** | Dev | §8.31 |
| 04/10 | **Nhật ký §8 tách ra `docs/log/`, số mục giữ nguyên (dev quyết, có đo)**: §8 chiếm **1.531/2.379 dòng · 118/216 KB** — hơn nửa PROJECT.md — nên đọc một mục (§8.27, 3,8 KB) phải mở cả file 216 KB, và file đã vượt hạn mức Read mặc định 2.000 dòng. Tách thành **5 file theo giai đoạn** (`phase-1` · `phase-2-3` · `phase-3a` · `phase-4` · `phase-4a`; `phase-5` tiếp theo sau), `PROJECT.md` §8 còn **bảng mục lục** (33 dòng lúc tách). **Số mục tuyệt đối không đổi** vì **271 tham chiếu** (đếm lúc tách; nay đã tăng) trỏ theo số (27 trong docstring code · 36 trong test · 23 trong `docs/` · 115 trong bốn file phiên · 70 trong PROJECT.md) — đánh số lại là phá chúng. Là **phép di chuyển thuần**: nối 5 file lại **giống hệt từng byte** §8 trong git (đã kiểm), phần §1–§7 và §9–§13 không đổi. Luật mới ở CLAUDE.md → "The run log", và `tests/unit/test_project_docs.py` fail build nếu một `§8.N` được trích mà không có mục, nếu index thiếu mục, hoặc nếu một số nằm ở hai file | Dev | `docs/log/`, CLAUDE.md |
| 04/10 | **Thư viện vẽ hình: matplotlib, để trong extra `[figures]` (D5-01, dev quyết)**: đúng khuôn mẫu D2-04 đã dùng cho `python-docx` — import **bên trong hàm**, thiếu thì `ExportDependencyError` + exit 3 chứ không phải traceback, nên **bản cài tối thiểu vẫn thuần Python** (biopython + PyYAML, CLAUDE.md). Đây là thư viện **offline, cài bằng pip, ghi được SVG · PDF · EPS** và vẽ được cả bản đồ vòng, heatmap lẫn đồ thị thống kê — không có lựa chọn thứ hai nào làm đủ cả ba. Hệ quả phải xử lý: matplotlib **đóng dấu thời gian + chuỗi creator** vào SVG/PDF ⇒ phải ghim như bản xuất Word đã ghim (timestamp cố định), nếu không thì "sinh lại được từ dữ liệu trung gian" là câu sai | Dev | §3 (P5-01) |
| 04/10 | **Thống kê: scipy trong extra `[stats]`, và MỌI thống kê phải có vector đối chiếu (D5-02, dev quyết)**: lý do D4-03 viết `popgen` bằng tay là *công cụ* (DnaSP = GUI Windows, pixy/libsequence đòi VCF) — lý do đó **không áp dụng cho scipy**, vốn cài được bằng pip và script được. Viết tay Mann–Whitney có hiệu chỉnh ties chính là cách sinh ra "số trông hợp lý nhưng sai". Bù lại, mỗi hàm thống kê đều bị **ghim bằng vector tính tay hoặc nguồn tham chiếu** như `popgen` đã làm | Dev | §3 (P5-03) |
| 04/10 | **Hình là HÀM THUẦN của `Result`, không bao giờ tính lại (D5-03, dev quyết)**: hình chỉ đọc ô của bảng, nên không thể có chuyện hình và bảng cho hai số khác nhau cho cùng một đại lượng — đúng lỗi §8.25 trong bộ áo mới. Kiểm bằng test: dữ liệu hình đọc phải bằng đúng các ô của bảng nguồn | Dev | §3 (P5-01) |
| 04/10 | **Quy tắc mở rộng hình theo N: ĐO, không đoán (D5-04, dev quyết)**: vẽ cùng một hình ở N = 2 · 8 · 34 trên bộ dữ liệu anh Điềm (bộ test đầu tiên đủ lớn để một hình 34 series lộ chuyện), rồi viết quy tắc từ chính ba phép đo đó — trần số series, hay đổi sang dạng tổng hợp khi vượt. Một hình 34 series không đọc được là một hình **hỏng**, không phải một hình dày | Dev | §8 (P5-02) |
| 05/10 | **Hai bào quan — mỗi lần chạy đúng một bào quan, mô-đun khai báo bào quan nó phục vụ (D5-06, dev quyết)**: bộ công cụ phục vụ **lục lạp và ty thể**, nhưng **một lần chạy chỉ thuộc một bào quan** — xác định theo qualifier `/organelle` của chính bản ghi; chuỗi trần thì bộc lộ bằng **đường vào** (`plastid mito annotate` là đường ty thể). **Bộ dữ liệu đầu vào cũng vậy: một lần chạy, một panel — một bào quan, không trộn plastome với mitogenome**; ngoại lệ duy nhất có tên: phân tích **chuyển đoạn giữa các ngăn** — nhân ↔ lục lạp ↔ ty thể (NUPT/MTPT, M6.1). Mô-đun chia **ba nhóm — dùng chung · chỉ lục lạp · chỉ ty thể** — cộng một danh sách **chuyển tiếp**: bảy phân tích dùng chung về bản chất nhưng hôm nay còn gắn vùng lục lạp (bảng ở §13, mục "Hai bào quan"). Tín hiệu duy nhất hôm nay là cờ `requires_ir`, mà nó trả lời "có cần cổng IR không", **không** trả lời "phục vụ bào quan nào" (`gene_table` không cần IR nhưng nội dung là bảng plastome; `codon_usage` không theo vị trí nhưng vẫn bị cổng chặn) ⇒ chuẩn hoá: mô-đun **khai báo nhóm của nó** trong registry, **artifact tự viết phải khai `/organelle`** (hôm nay bản chú giải mito không khai — H-11), và chạy sai bào quan — hoặc một bộ dữ liệu trộn hai bào quan — thì **dừng có tên**; hôm nay cổng IR mới chặn được một chiều, bộ trộn thì chưa (H-09) | Dev | §13 (H-09 · H-10 · H-11) |
| 05/10 | **Annotate lục lạp sẽ do dự án tự viết — theo khuôn CGAS nhưng sửa đúng các lỗ đã đo (D5-08, dev quyết)**: bước annotate trong pipeline (assembly → **annotate** → thống nhất tên gene → phân tích — **[thứ tự đã đổi 07/10/2026: đồng nhất tên chuyển lên TRƯỚC chú giải, xem dòng D5-09]**) là công cụ của mình, lấy **mô hình của CGAS** — bộ tham chiếu · chú giải từng mẫu · report, để phần thống nhất tên và phân tích là các bước sau (toolkit đã có) — nhưng **thay lõi tìm–quyết**: nhiều nguồn BLAST theo loại vùng (`blastx` trên `/translation` công bố cho CDS — không tự dịch, §8.21; `blastn` cho rRNA · IR · vùng không mã hoá), trọng tài **cấu trúc trước, điểm sau** (điểm chỉ so trong cùng chương trình), **chồng lấn giữ cả hai + báo** (đo trên genome vàng: 4 cặp gen chồng lấn thật), nhiều reference ⇒ **một call mỗi locus + liệt kê ứng viên thua** (PGA: 147 gen vì nhân bản, §8.15), mọi ca bị loại **kèm lý do**. **Chloe bị bỏ khỏi kế hoạch** (dev quyết; code hôm nay vẫn mặc định `chloe` tới khi engine mới đạt Done-when — xem CLAUDE.md). Done-when đo trên **bản công bố**: genome vàng ra đúng **132 = 86 + 38 + 8** với mọi ca lệch liệt kê; corpus 34 bản ghi **leave-one-out**. Thiết kế + 4 điểm còn mở (D-a tRNA · D-b bộ tham chiếu · D-c phase · D-e lõi chung): `docs/D5-08_plastid_annotator.md` | Dev | `docs/D5-08_plastid_annotator.md` |
| 05/10 | **Hình nằm ở tầng `result/`, render ở MỨC RUN từ `result.json` — không phải `core/figure.py` như board viết (dev quyết, có đo)**: (a) `tests/unit/test_architecture.py` chỉ cho `core/` import `plastid.core`, mà một hình là **hàm của `Result`** — `Result` ở tầng `result`; (b) **hai hình đọc bảng của hai module khác nhau** (`circular_map`, `ir_junction`: `genome_stats` ở khung GenBank + `structure` ở khung chuẩn hoá), việc mà một callback theo từng module không làm được. Hệ quả: hình render **sau** vòng lặp module, đọc `result.json` đã ghi — nên **một run được phục vụ hoàn toàn từ cache vẫn có hình**, và `plastid run --figures` với `plastid figures <run-dir>` là **cùng một hàm**. Đo: chạy hai lần trên cặp Iris đã pin ⇒ svg · pdf · eps **byte y hệt**; lần thứ hai là cache hit | Dev | §8.36 |
| 05/10 | **Luật "đúng một ô" cho report đã sinh: ô CHỮ có chữ cái là văn xuôi, không phải lượng (dev quyết, có đo)**: luật nghiêm (mọi số phải nằm trong giá trị đã render của một ô) hỏng trên chính bộ giao phẩm Pha 4 với **10 phát hiện, 6 lớp** — và **cả 10** là một lượng va vào một số in bên trong một ô *chữ* được report trích nguyên văn (id khối `…_072561.1#10`, chuỗi luật `mean + 1 sd …`, cỡ nhóm `Coleus: 6, Oldenlandia: 5`, một câu ghi chú). Luật nay: giá trị của cột `str` mà chứa chữ cái bị **bỏ khỏi thân report đang quét** và **không tính là ô mang số**; ô `str` render đúng một số trần thì vẫn là ô mang số. Ca *hai dòng render cùng một giá trị* **không** bị coi là mơ hồ — hai ô không thể mâu thuẫn (genome vàng: hai bản IR đều 26.214 bp) | Dev | §8.35 |
| 05/10 | **Ưu tiên hiện tại: LỤC LẠP trước — ty thể giữ nguyên và chừa hướng tích hợp sau (D5-07, cả nhóm chốt)**: dữ liệu **ty thể thực vật còn rất ít**, nên mọi phát triển mô-đun mới tập trung vào **lục lạp**; đường ty thể đã có (`mito_annotate` · `mito_genes` · bộ tham chiếu mito) **giữ nguyên, không mở rộng**. "Chừa hướng" là các điểm mở của D5-06, giữ nguyên hiệu lực: từ vựng tên theo bào quan · nhóm mô-đun khai báo · bộ dữ liệu một bào quan · bản ghi tự khai `/organelle` — nhờ đó sau này thêm phân tích ty thể **không phải đập lại cấu trúc**. Hệ quả xếp việc: **H-10** (mở 7 phân tích chuyển tiếp cho ty thể) và **H-06** (mở rộng từ vựng mito lên 10 locus chuẩn) **tạm hoãn** tới khi có dữ liệu ty thể; M6.2/M6.3 vẫn chờ dữ liệu như cũ | Cả nhóm (cô Trang) | §13 "Hai bào quan" (H-06 · H-10) |
| 07/10 | **Đồng nhất tên gen là một bước riêng, chạy TRƯỚC chú giải và TRÊN bộ tham chiếu (D5-09, dev quyết — đảo chiều D3-02)**: pipeline nay là *(nguồn tham chiếu → audit → người dùng chọn một bản ghi ★ → `gene_normalize` chuẩn hoá DB)* → lắp ráp → chú giải → *`gene_normalize` chế độ kiểm* → phân tích; mẫu **thừa hưởng** từ vựng A qua bước chú giải, và chạy lại ở chế độ kiểm mà có gì đổi ⇒ **đó là phát hiện**, không phải việc thường. **Từ vựng A = tên của bản ghi ★** — **loại tham chiếu thứ 6** (không chọn lại theo taxonomy mỗi lần chạy, không phải của công cụ), do **người dùng cuối chọn từ báo cáo audit**; tên in trong bản thảo đổi theo ⇒ PI ký nhận. **Chiều chuẩn hoá đảo** — 5/9 locus đã đo đổi cách viết (`ycf3`·`ycf4`·`psbN`·`clpP`·`ndhL` thay `pafI`·`pafII`·`pbf1`·`clpP1`·`ndhI`); ba module đọc khoá chuẩn (`kaks`·`substitutions`·`popgen`) ⇒ **Ts/Tv 2,6917 là dự đoán cho tới khi đo lại — [ĐÃ ĐO 08/10/2026: giữ nguyên 2,6917 (323 Ts / 120 Tv), và phép đo cho thấy nó *bất biến theo chiều* vì việc lật chiều là đổi nhãn chứ không chọn lại tập gen dùng chung — §8.45]**, không giả định. Mô-đun mới `gene_normalize` (**không nới `gene_names`**); slot của bước mới còn mở — ROADMAP → Open decisions | Dev · người dùng cuối (chọn ★) | `docs/D5-09_reference_naming.md` |
| 08/10 | **Nơi chạy phần lắp ráp: tuỳ máy người dùng (xác nhận 08/10/2026)** — máy **đủ RAM** thì chạy **tại chỗ**; **không đủ** thì đề xuất **đường từ xa Galaxy** (usegalaxy.eu, §13.7 / M7.4). Luật 03/10 ("giai đoạn 6–7 chạy máy chủ từ xa") nay chỉ còn nguyên với **giai đoạn 6 — thời gian phân kỳ**. Chính sách đưa dữ liệu chưa công bố lên máy chủ công khai: **chưa quyết định** (§10.12) | anh Điềm | §13 (điều kiện chung) · ROADMAP → "Cross-phase rules" |
| 08/10 | **Người dùng cuối, và cách dùng (định nghĩa, xác nhận 08/10/2026)**: người dùng cuối là **nhà nghiên cứu công nghệ sinh học / tin sinh học**; giao diện là **dòng lệnh trên Ubuntu** (Windows qua Docker — như quyết định 01/10). Định nghĩa này là thứ D5-09 §7.1 dùng cho việc chọn ★, và là người phải đọc được báo cáo audit | anh Điềm | CLAUDE.md · `docs/D5-09_reference_naming.md` §7.1 |
| 10/10 | **Thời gian phân kỳ và ancestral state KHÔNG phải mô-đun của bộ công cụ — xoá khỏi kế hoạch (dev quyết, 10/10/2026)**: hai phân tích đó cần **cây · hiệu chuẩn hoá thạch · một câu hỏi khoa học**, không cần bộ gen — chúng là phân tích của **cô Trang** cho bài báo, chạy bằng công cụ của cô (`mcmctree` · `BEAST2` ghi ở §7 là *đã khảo sát, không đưa vào*). Bộ công cụ làm **code**; *phân tích gì · khẳng định gì · đưa gì vào bài báo* là của PI (CLAUDE.md → who decides). Hệ quả: M6.4 và M6.5 ra khỏi ROADMAP, phase 6 còn ba mô-đun (NUPT/MTPT · cấu trúc mito · so sánh lục lạp–ty thể), và **không còn phần nào cần máy chủ từ xa** | Dev | §13.6 · ROADMAP → Phase 6 |
| 09/10 | **Cổng audit chặn run, và "bỏ file" là một hành động được khai ra (P3B-09, dev quyết, có đo)**: D5-08 viết "mọi ca audit hỏng thì dừng run — exit 2, không chú giải gì", và vòng lặp là *"fix or drop the bad file → run again"*. Cổng nay nằm trong `plastid_annotate._reference_set` — **đúng tập bản ghi run sắp tìm, trước HSP đầu tiên** — nên mọi đường vào (lệnh `annotate` lẫn `plastid run -m plastid_annotate`) đều bị chặn như nhau. Nửa "drop" của vòng lặp là **`--exclude <accession>`**: khai trên dòng lệnh, ghi vào manifest (nguyên argv), nêu trong nhãn bộ tham chiếu và trong note của `result.json`; một cái tên không khớp bản ghi nào thì **dừng run**. Lý do phải có nửa ấy, đo được: bộ đã pin chặn đúng một bản ghi (`NC_024936.1`, `rps19` partial), nên nếu chỉ có "dừng" thì annotator **không chạy được trên chính bộ tham chiếu của dự án** cho tới khi người dùng tự dựng một thư mục khác. Một pin **không** waive được BLOCK — `--exclude` không phải lối vòng qua cổng, nó là "file ấy không được dùng" | Dev | §8.49 · §9.L-5 · §9.L-6 |
| 09/10 | **`plastid annotate` = một loạt run, mỗi mẫu một thư mục (P3B-09, dev quyết)**: D5-08 viết "mỗi mẫu một artifact, run một report"; chỗ dựa là §2 ("output tách theo lần chạy/dataset, không ghi đè tại chỗ") và seam 3 ("lệnh thứ hai không được thành đường code thứ hai") — nên **mỗi mẫu đi qua đúng thân `_execute_run`** mà `plastid run` và `plastid mito annotate` dùng, ra một thư mục run với manifest · report · artifact riêng, dưới `<out>/<accession>/`. Mẫu hỏng: **được nêu tên, thư mục run của nó ghi `status: failed` + lý do, batch đi tiếp**, exit code khác 0 (mã của lỗi đầu tiên). Manifest thêm khoá `dataset` (tuỳ chọn — schema 2 không đổi) để một run hỏng trước khi ghi gì vẫn nói được nó là mẫu nào | Dev | §8.49 |

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
| Ts/Tv (so với Iris gatesii NC_024936.1) | **2,6917** = 323 Ts / 120 Tv; **chỉ vùng CDS** (79 gen dùng chung, không tính vùng lặp IR); MAFFT; 68.448 cột căn — phạm vi chốt 02/10, **bộ gen đo lại 04/10**, **chiều tên đo lại 08/10 (không đổi)** | tự đo — xem hai ghi chú dưới bảng |
| Số công bố của loài này | *Iris domestica* 153.729 bp, 132 gen; *Belamcanda chinensis* 153.735 bp, 132 gen | 2 bài công bố |

> **Ghi chú 04/10/2026 — dòng Ts/Tv đo lại.** Phạm vi không đổi (chỉ CDS · gen dùng chung ·
> MAFFT · bỏ IR), nhưng **bộ gen dùng chung** trước đây đếm theo *cách viết tên*: ba locus
> `atpI`·`psaI`·`psbB` mà NC_024936.1 viết là `atpL`·`psaL`·`bspB` bị coi là gen nó không có.
> Đo lại bằng ba tín hiệu (trình tự, vị trí, `/product`) ⇒ cùng locus ⇒ bộ dùng chung là **79
> gen**, và số là **2,6917** (323 Ts / 120 Tv, 68.448 cột). Số cũ **2,686** (317/118, 76 gen,
> 66.066 cột) vẫn nằm nguyên trong nhật ký §8.25 như số đo của bộ hẹp hơn. Chi tiết: §8.27.

> **Ghi chú 08/10/2026 — dòng Ts/Tv đo lại lần hai, dưới chiều tên mới.** D5-09 đảo chiều chuẩn
> hoá tên (từ vựng in ra nay theo bản ghi ★) và **ba mô-đun đọc khoá chuẩn** (`kaks` ·
> `substitutions` · `popgen`) ⇒ dòng này là **dự đoán** cho tới khi đo lại. Đo lại **08/10/2026**
> theo đúng thứ tự pipeline (`gene_normalize` ghi mẫu về từ vựng A → `substitutions --set
> scope=cds`): **2,6917 = 323 Ts / 120 Tv · 79 gen · 68.448 cột** — không một đơn vị nào lệch.
> Phép đo còn cho thấy **vì sao nó bất biến theo chiều**: bảng alias ánh xạ mọi thành viên của
> một locus về locus đó, nên lật chiều là *đổi nhãn* tập gen dùng chung chứ không *chọn lại* nó
> (đo: cả hai bản ghi Iris đều 86 CDS → 79 tên chuẩn, 0 locus có hai cách viết). Nhãn in ra thì
> đổi — `substitutions_by_gene` nay in `ycf3`·`ycf4`·`psbN`·`clpP`·`ndhL`. Chi tiết: §8.45.

---

## 7. Danh mục công cụ

### 7.1. Bảng tổng

| Việc | Công cụ | Trạng thái test |
|---|---|---|
| SSR | MISA standalone 2.1 (`misa.pl`) + cpstools 3.0 `SSRs` | đã chạy thật — 37 SSR (cả hai công cụ) |
| Repeat phân tán + palindromic | repeat-match (MUMmer3) + Vmatch 2.3.1 | đã chạy thật — 9 match ≥30 bp, có cặp đảo IR 26.214 bp **[ĐÃ ĐO LẠI 10/10/2026 — §8.58: bản vàng cho **9** match ≥30 bp (4 dispersed · 5 palindromic) và cặp IR 26.214 bp, **hai engine độc lập trùng nhau từng dòng**; số in ra là số của bản vàng, đo lại, **không** kế thừa dòng khảo sát này — trùng số là trùng hợp, đã ghi]** |
| So sánh nhiều genome (khối LCB, synteny) | Sibelia 3.0.7 | đã chạy thật — 3 plastome · **[không đưa vào bộ công cụ — 07/10/2026: tuyến synteny của dự án là `nucmer`, D4-02/M4.2; ghi chú có ngày ở §13.0]** |
| Biểu đồ % tương đồng theo toạ độ | nucmer (MUMmer4) `--maxmatch` → `delta-filter` → `show-coords` + matplotlib | đã chạy thật — 21 khối, %IDY 83,7–98,2 (mean 91,5) *(con số này đo 26/9 trên **dữ liệu khác**; cặp Iris = **7 khối · 97,37** — §8.26, §9.H)* |
| Chú giải plastome | **công cụ tự viết** (D5-08) — hiện tạm: Chloe 0.1.18; PGA đã đo và bị loại (§8.15) | **[ĐÃ THAY 10/10/2026: động cơ tự viết đã đạt Done-when của nó — §8.55: bản vàng 132 = 86 + 38 + 8 với 35 khác biệt được liệt kê; corpus 418 so với PGA 588 (30/34 bản ghi tốt hơn), đọc strict là 1070 so với 699 — nêu kèm. Tập chloe · pga trong `annotation` nay đến hạn nghỉ: **P3B-11**]** |
| Vẽ bản đồ vòng genome | pycirclize 1.10.1 | đã chạy thật — PNG 494.859 B @300 dpi |
| Lắp ráp từ reads ngắn | GetOrganelle 1.7.7.1 + SPAdes 4.3.0 | đã chạy với reads **mô phỏng** — 127.679 bp |
| Lắp ráp long-read (ít RAM) | raven-assembler 1.8.3 (+ OATK 1.0) | đã chạy thật — badread 150x ONT → 2 contig, dài nhất 107.394 bp |
| Lắp ráp long-read + DB HMM | OATK 1.0 + `embryophyta_pltd.fam` (OatkDB) | đã chạy thật — 126.660 bp, identity 100,00 %; file `.annot_pltd.txt` **rỗng** |
| Vẽ + thống kê đồ thị assembly (headless) | Bandage 0.9.0 | đã chạy thật — PNG 103.289 B trên GFA thật của OATK |
| Mạng haplotype | PyPopART | đã chạy thật — 13 trình tự matK → mafft → mjn/analyze/visualize, rc=0 |
| Ước tính thời gian phân kỳ **[ĐÃ THAY 10/10/2026 — không đưa vào bộ công cụ: không phải mô-đun, mà là phân tích của cô Trang cho bài báo; xem §13.6]** | PAML mcmctree 4.10.10 · BEAST2 2.7.7 | xác nhận chương trình chạy; **chưa** chạy phân tích thật |
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
| Chloe | sống | chú giải theo loài họ hàng | plann.pl **[ĐÃ THAY 10/10/2026 — P3B-11, §8.56: Chloe nghỉ hưu cùng `pga`; phép đo của nó ở §8.17, đường install giữ ở dòng cập nhật D3-01]** |
| mVISTA | sống | biểu đồ % tương đồng | nucmer + matplotlib (21 khối — probe 26/9, dữ liệu khác; cặp Iris = 7 khối, §8.26; PNG 243.787 B) |
| MISA-web | sống | SSR | MISA standalone — 37 SSR |
| REputer | sống | repeat phân tán + palindromic | repeat-match (MUMmer3) + vmatch 2.3.1 |
| GetOrganelle (Galaxy) | sống | lắp ráp reads ngắn | GetOrganelle 1.7.7.1 local + SPAdes 4.3.0 |
| DOGMA | **chết** (DNS không phân giải) | chú giải | plann.pl |
| CpGAVAS bản 1 | **chết** | chú giải + phân tích | plann.pl |
| Verdant | **chết** | pipeline chú giải | plann.pl |
| CGAP | **chết** | phân tích so sánh | cpstools 3.0 + CGAS |
| PREP-Cp | **chết** | dự đoán RNA editing | không còn web tương đương (PREPACT, CURE-Chloroplast cũng `http=000`) → đo từ RNA-seq |

*Cột "bản local dùng thay" là ảnh chụp 26/9. Riêng các dòng **chú giải**, đường hiện tại khác: từ
03/10 lục lạp dùng `annotation` (Chloe mặc định — §8.15–8.17; **[ĐÃ THAY 10/10/2026: động cơ tự viết theo D5-08 đã đạt Done-when (§8.55) và là đường chạy của `plastid annotate`; tập chloe · pga đến hạn nghỉ — P3B-11]**),
ty thể dùng `mito_annotate` (P3A); `plann.pl`/`PGA.pl` chỉ còn là lịch sử khảo sát.*

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

**Qua conda (conda-forge + bioconda)** — lệnh cài đầy đủ và có pin: [`environment.yml`](environment.yml):

| Công cụ | Kênh | Trạng thái |
|---|---|---|
| sibelia (3.0.7) · bandage (0.9.0) · oatk (1.0) · raven-assembler (1.8.3) · flye (2.9.6) · spades (4.3.0) · getorganelle (1.7.7.1) · novoplasty · mummer4 (nucmer) · mreps · trf · vmatch (2.3.1) · minimap2 · samtools · wgsim · blast · paml (mcmctree 4.10.10) · beast2 (2.7.7) · perl-bioperl · perl-xml-simple | bioconda | ĐÃ GHI |
| pycirclize (1.10.1) | conda-forge | ĐÃ GHI (kiểm 02/10: conda-forge 27 file, tới 1.10.1) |
| cpstools (3.0) | **bản đang chạy: pip/PyPI** (site-packages python3.12 của base env, tác giả Xu Wenbo, MIT) — trong `~/miniconda3/pkgs/` có tarball bioconda `cpstools-3.0-pyhdfd78af_0.tar.bz2` nhưng **không phải bản đang chạy** | Sửa lại: cpstools **không có trên conda-forge**; bioconda có (22 file, tới 3.0). **[ĐO 10/10/2026 — P8-02, §8.61: trên máy dev nay là **pip 3.0 trong env riêng** `/root/toolbox/env/cpstools` (Python 3.12), cùng `KaKs_Calculator` **2.0.1** (bioconda `kakscalculator2`) ở `/root/toolbox/env/kakscalculator`; đây là **công cụ đối chiếu, không phải phụ thuộc mô-đun** (D-j) — `docs/tools/cpstools.md` là dòng khai của nó]** |
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
| PGA / PGA.pl **ĐÃ THAY 10/10/2026 — P3B-11, §8.56: toolkit không còn chạy PGA; đầu ra của nó sống tiếp trong `tests/data/corpus/*.pga.gb` như đối chứng đã đo** | `github.com/quxiaojian/PGA` (bản 2.0: `github.com/quxiaojian/PGA2`); CGAS `INSTALL.md` cũng ghi PGA phải cài tay từ GitHub. **Nguồn đã đọc mã 05/10/2026: commit `0683a0ec2b2061c49c60ec6c5a2293cd08623253` (29/10/2020) · `PGA.pl` 16.320 dòng · sha256 `0d55b70e53de5d677522cb8fed646898c2ea28bc3f56275a9190128bfd65c3c`** | ĐÃ GHI (kiểm 02/10) |
| Chloroplot + env R | `remotes::install_github("shuyuzheng/Chloroplot")`; env: r-base · bioconductor-biostrings · bioconductor-genbankr · bioconductor-cordon · r-reutils · r-circlize · r-dplyr · r-magrittr | ĐÃ GHI |
| CGAS | `github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS` (19 file `.py`, 21.369 dòng, MIT) | ĐÃ GHI (kiểm 02/10). **[ĐO 10/10/2026 — P8-03, §8.62: chạy thật tại đây (`uv venv /root/toolbox/env/cgas` + `pip install -e .data/cgas-review`, `cgas --version` → **v1.0.1**), dùng làm **công cụ đối chiếu** — không phải phụ thuộc mô-đun; dòng khai đầy đủ ở `docs/tools/cgas.md`]** |

**CHƯA GHI nguồn trong tài liệu của nhóm** (chỉ có tên hoặc chỉ có vị trí trên máy):

| Công cụ | Hiện biết | Nguồn chính thức |
|---|---|---|
| MISA standalone (`misa.pl` + `misa.ini`) | nằm ở `~/toolbox/src/misa/`, không có trong PATH → gọi `perl ~/toolbox/src/misa/misa.pl <fasta>` (perl 5.32.1); `misa.ini` phải cùng thư mục | **ĐÃ ĐÓNG (kiểm 02/10, đối md5 với bản tải lại từ nguồn):** https://webblast.ipk-gatersleben.de/misa/ → `misa_sourcecode_25082020.zip` · 4.946 B · md5 `cf9d3579f7caacd02f99830b329c7d9c` · bản **v2.1 (25/08/2020)** · `misa.pl` 20.148 B md5 `cfebc1d1a36fddfd6e0d897c96babd12` · `misa.ini` 158 B md5 `49a1c56048631647c0364532a6964cf9` · tác giả Thomas Thiel & Sebastian Beier (IPK) · trích dẫn Beier và cs. 2017, Bioinformatics 33:2583–2585 (doi:10.1093/bioinformatics/btx198) · Quyết định 02/10: **đóng gói kèm** (xem §3); image cần perl. |
| repeat-match | "có sẵn trong env toolbox" | **[ĐÃ THAY 10/10/2026 — P2B-02, §8.58, đo được: nhị phân đang chạy thuộc gói `mummer4-4.0.1` (bioconda), không phải MUMmer3 như dòng cũ ghi; nguồn MUMmer3 (`mummer.sourceforge.net`, tarball 3.23) **không phải bản đã đo**]** `repeat-match` **không có cờ phiên bản** (`--version` là option không hợp lệ) nên bản phát hành đọc từ `mummer` cùng thư mục — lệch pin thì dừng; toạ độ **1-based**, hậu tố `r` nằm ở **đầu của bản sao đọc trên sợi đảo** (khung reverse-complement) — chi tiết ở `docs/tools/mummer4.md` |
| MAFFT | bản đang chạy: `/home/nthanhdiem/miniconda3/bin/mafft` → `v7.526 (2024/Apr/26)` | **ĐÃ ĐÓNG (kiểm 02/10):** cài từ **conda-forge** — `mafft 7.526` build `h4bc722e_0` (2.610.460 B, md5 `cf3b2acf649bc31e576ed43a2fa4aaa3`); conda-forge cũng có build `h4ab18f5_0` và `ha3fa977_0` của cùng 7.526; bioconda chỉ tới 7.525. Máy còn bản apt: **7.490-1 trên Ubuntu 22.04** (Ubuntu 24.04 là **7.505-1**) ⇒ report/log phải ghi **đường dẫn + phiên bản**, không ghi trần "MAFFT". Chạy lại hàm Ts/Tv bằng cả 2 bản mafft → cùng một kết quả ⇒ số vàng không phụ thuộc bản mafft (số vàng sau đó đo lại trên bộ dùng chung 79 gen — §6, §8.27) |
| Krait | đã thử, là GUI nên không dùng | chưa ghi nguồn |
| PGA — giấy phép | PGA = GPL-3.0 | (đã ghi ở §3) |

**Bổ sung 10/10/2026 — CPGview (P8-01, §8.60): KHẢO SÁT, KHÔNG ĐƯA VÀO.** Nó được M8.1/§13.8 nêu
tên mà §7.8 chưa từng có dòng nào — luật phủ gọi đó là lỗ hổng; đây là dòng còn thiếu.

| Công cụ | Nguồn | Trạng thái |
|---|---|---|
| CPGview (Liu và cs. 2023, *Mol Ecol Resour* 23:694–704, doi:10.1111/1755-0998.13729) | web `http://www.1kmpg.cn/cpgview/` · repo `github.com/cliu6688/CPGView` · container Figshare `cpgview_latest` **2.871.624.635 B**, md5 `4b371a5ce9410a4a433c9499a6b04bbd` | **KHẢO SÁT — KHÔNG ĐƯA VÀO** (đo 10/10/2026, §8.60): repo GitHub **chỉ là trang giới thiệu** (0 mã, 0 release, 0 giấy phép); không pip/conda/CLI ⇒ **không có phiên bản để pin, không có argv để ghi**; bản cài duy nhất là container **2,87 GB** cho một công cụ **vẽ hình** (bộ công cụ có hệ hình riêng — M5.1) và bản cài tối thiểu phải pip-only; đường web vướng §10.12, và **hôm nay ứng dụng không trả lời** (cổng 16085 — cùng IP đó cổng 80 mở, control `portquiz.net:16085` mở, nên không phải mạng phía mình) |

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

*Giá phải trả:* PMGA bọc MAKER · CPGAVAS2 · BLASTn · ARAGORN · tRNAscan-SE · MISA · TRF · Vmatch —
nặng hơn nhiều so với PGA (Perl + BLAST). **[ĐÃ THAY 03/10/2026 — P3-04 đã quyết: PMGA không chạy
nổi ở đây (container 3,18 GB nén → 11 GB giải nén, cần Apptainer + MAKER, trần WSL 8 GB; web
service từ chối kết nối) ⇒ đường chú giải mito là bộ gọi gen tự viết, §8.22.]**
*Độ mạnh bằng chứng:* phạm vi của MITOS2 là chắc (trang phân phối và tài liệu nói thẳng); số
benchmark của PMGA đọc từ bản tóm tắt web, **chưa đọc bài gốc** — giữ làm kỷ lục; PMGA đã bị bỏ ở
P3-04 nên không còn phải "kiểm lại trước khi chọn".

### 7.8c. LỖI ĐÃ SỬA — công cụ chỉ nằm trên một máy thì không phải là một phần của bộ công cụ

**Chuyện đã xảy ra (03/10/2026, do chính tôi gây ra):** để gỡ chặn Giai đoạn 3, tôi cài MAFFT 7.526
và BLAST+ vào `~/toolbox/` trên máy dev — nhưng **không sửa gì trong repo**. Nếu dừng ở đó thì:

- mọi số tôi đo được **không ai tái lập được** — kể cả anh Điềm, kể cả chính tôi sáu tháng nữa;
- người cài `plastid-suite` theo hướng dẫn trong repo **không có công cụ nào**, và tier B sẽ skip
  sạch mà không rõ vì sao;
- Docker image ở Giai đoạn 9 sẽ không chạy được phần cần công cụ.

**Kiểm tra thì phát hiện lỗi nặng hơn:** `environment.yml` cũ liệt kê công cụ **không kèm phiên bản**
— trong đó có dòng `mafft` trần. Mà §9.B-14 **đã đo**: cùng một đầu vào, aligner khác nhau cho
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

**Nhật ký đã tách ra `docs/log/`** (04/10/2026) để đọc được một mục mà không phải mở cả
file: §8 khi đó chiếm 1.531/2.379 dòng và 118/216 KB. **Số mục giữ nguyên** — `§8.27` vẫn
là `§8.27` — nên mọi tham chiếu trong repo (đếm lúc tách: 27 trong code, 36 trong test, 23 trong
`docs/`, 115 trong **bốn** file phiên — nay năm, thêm `docs/handoff-archive.md`) vẫn trỏ đúng. Mục mới: ghi vào file của giai đoạn đang làm, rồi
thêm **một dòng** vào bảng dưới. Luật ở `CLAUDE.md` → "The run log".

| § | Nội dung | File |
|---|---|---|
| 8.1 | Cài đặt công cụ — môi trường đã dùng | [`docs/log/phase-1.md`](docs/log/phase-1.md) |
| 8.2 | Chạy chuỗi trên máy anh Điềm (26/09/2026 — 15 giây, 18/19 mốc đạt) | [`docs/log/phase-1.md`](docs/log/phase-1.md) |
| 8.3 | Đo trên máy Hermes (26/09/2026 — Python thuần, không thư viện ngoài) | [`docs/log/phase-1.md`](docs/log/phase-1.md) |
| 8.4 | Môi trường dev (01/10/2026) | [`docs/log/phase-1.md`](docs/log/phase-1.md) |
| 8.5 | Đo lại bằng implementation trong repo (02/10/2026 — P1-03) | [`docs/log/phase-1.md`](docs/log/phase-1.md) |
| 8.6 | Bộ dữ liệu công khai thay thế — Oldenlandia / Coleus (03/10/2026 — D-04) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.7 | Bảng gen xuất bản + xuất Word (03/10/2026 — P2-09, D2-04) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.8 | Bản nháp bản thảo trên bộ dữ liệu thay thế (03/10/2026 — P2-10) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.9 | Bốn tên gen lạ ở *C. amboinicus* — cùng locus, khác tên (03/10/2026 — D3-02) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.10 | Hai bản ghi là CỦA NHÓM — và cả một phương ngữ đặt tên (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.11 | Kiểm tra ORF — P3-01 (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.12 | Chuẩn hoá tên gen — P3-02 (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.13 | Chú giải & file nộp — P3-03 · P3-04 · P3-05 (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.14 | Đường mito — P3-04 phần chạy được (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.15 | PGA chạy thật trên genome vàng — và ba cái bẫy (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.16 | Đổi công cụ chú giải: Chloe thay PGA (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.17 | Chloe chạy thật qua pipeline của mình — P3-03 chốt (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.18 | Chạy thật `table2asn` — bản `.tbl` của mình bị NCBI chê, và ba lỗi đã sửa (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.19 | Nửa còn lại của luật backend: ghim kết quả công cụ chạy ở nơi khác (03/10/2026) | [`docs/log/phase-2-3.md`](docs/log/phase-2-3.md) |
| 8.20 | Writer GenBank — mảnh còn thiếu để một annotator tự viết có đầu ra dùng được (P3A-01, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.21 | Bộ tham chiếu mito — dựng từ 8 bản ghi RefSeq, và cái bẫy nằm trong chính bộ tham chiếu (P3A-02, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.22 | Bộ gọi gen mito — chuyển chú giải từ bộ tham chiếu, và sáu cái bẫy đo được (P3A-03, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.23 | Bộ gọi tRNA — ARAGORN 1.2.41, và 17 không phải là một tập hợp (P3A-04, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.24 | Chạy trên mẫu chưa ai chú giải — `AP045469.1`, và tám locus không gọi được (P3A-07, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.25 | Chạy thật hai module so sánh — và một lỗi im lặng sống sót vì chưa từng chạy dữ liệu thật (P2-07/P2-08, 03/10/2026) | [`docs/log/phase-3a.md`](docs/log/phase-3a.md) |
| 8.26 | Synteny chạy thật — "21 khối" không tái lập được, và con số thật của cặp Iris (P4-01, 04/10/2026) | [`docs/log/phase-4.md`](docs/log/phase-4.md) |
| 8.27 | Ka/Ks chạy thật — và ba locus mà cách viết tên che mất (P4-02, 04/10/2026) | [`docs/log/phase-4.md`](docs/log/phase-4.md) |
| 8.28 | Di truyền quần thể viết tay — công thức, hai lần kiểm chứng, và cái bẫy Fst (P4-03, 04/10/2026) | [`docs/log/phase-4.md`](docs/log/phase-4.md) |
| 8.29 | Barcode + vùng siêu biến chạy thật — π dùng lại, và hai lỗi tên lộ ra (P4-04, 04/10/2026) | [`docs/log/phase-4.md`](docs/log/phase-4.md) |
| 8.30 | Bộ bảng so sánh cho hai bản thảo — ba run trên panel đã pin (P4-05, 04/10/2026) | [`docs/log/phase-4.md`](docs/log/phase-4.md) |
| 8.31 | Bộ dữ liệu anh Điềm gửi — bộ test đầu tiên không do dự án tự chọn (P4a-01, 04/10/2026) | [`docs/log/phase-4a.md`](docs/log/phase-4a.md) |
| 8.32 | Quy mô: 34 bản ghi tốn bao nhiêu (P4a-03, 04/10/2026) | [`docs/log/phase-4a.md`](docs/log/phase-4a.md) |
| 8.33 | Ba cách đọc một chuỗi — công bố · PGA · Chloe trên 34 plastome (P4a-04, 04/10/2026) | [`docs/log/phase-4a.md`](docs/log/phase-4a.md) |
| 8.34 | Ba số đo cho annotator lục lạp — gen chồng lấn · chiều SSC · độ sạch bộ tham chiếu (D5-08, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.35 | Kiểm thử truy vấn số trên report đã sinh — 10 va chạm buộc thu hẹp luật (P5-05, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.36 | Hình vẽ — công thức xuất vector tất định, và dải vùng không lấy từ mép cửa sổ (P5-01, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.37 | Bảy hình đầu tiên, và quy tắc mở rộng theo N đo bằng ba lần vẽ (P5-02, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.38 | Thống kê — vector tính tay, và một so sánh thật cho ra 0 khác biệt có ý nghĩa (P5-03, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.39 | Xuất bản phẩm — Excel tất định, gói Supplementary, và 34 trang PDF (P5-04, 05/10/2026) | [`docs/log/phase-5.md`](docs/log/phase-5.md) |
| 8.40 | Bộ tham chiếu lục lạp — 34 bản ghi Iridaceae pin lại, và luật leave-one-out thành dữ liệu (P3B-01, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.41 | Ba con số của D5-09 §9.3 — biên, phép thử hàng xóm, `identity_floor` (P3B-02, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.42 | Khoá cache mang thêm bốn thứ — chú giải, thứ tự, mã nguồn, thứ được yêu cầu (P3B-05, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.43 | Cổng kiểm bộ tham chiếu — lớp chặn, lớp báo, lần đọc thứ hai bằng MAFFT, ngoại lệ buộc vào sha256 (P3B-03, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.44 | `gene_normalize` — bốn phép thử, ba trạng thái, tham chiếu chuẩn có ghi lại; và chiều alias đảo theo bản ghi ★ (P3B-04, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.45 | Ts/Tv đo lại dưới chiều tên mới — 2,6917 giữ nguyên, và phép đo cho thấy nó *bất biến theo chiều* (P3B-06, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.46 | Lõi gọi gen dùng chung — tách khỏi đường mito, ngưỡng plastid đo trên bản ghi vàng, và ba lỗi im lặng (P3B-07, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.47 | Annotator plastid chạy thật trên genome vàng — 86 CDS (nhân bản IR), ba cặp chồng lấn §8.34, và hai khác biệt mức kết luận (P3B-08, 08/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.48 | Một locus bị chú giải hai lần: bản cụt không được tìm kiếm, và 0,5 nằm giữa một khoảng trống đo được — `ycf1` của bản ghi vàng nay khớp bản công bố từng base (H-19, 09/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.49 | `plastid annotate`: cổng audit trước run (exit 2, 13 giây, không artifact), batch một mẫu một thư mục, bào quan khai theo nhóm trong registry — và ba phát hiện đo được: `rps19` mất 54 bp khi bỏ bản ghi bị chặn (H-20) · luật BLOCK partial quá rộng (H-21) · `.docx` chưa từng chạy được vì thiếu citation (P3B-09, 09/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.59 | **Bản thảo hai loài chạy lại** (P2B-03, M2.10) — **hai run mỗi loài** (arity: CLI từ chối trộn, exit 2, đo được) · π cùng chi: Oldenlandia 0,023522 (95 % cột dùng được) · Coleus 0,008462 (98 %) · Ts/Tv scope CDS **1,4616** / **1,5134** (từng cặp 0,9583–1,5337 · 1,3387–1,6240) · repeats 4 + 4, mỗi loài dài nhất **đúng bằng IR của nó** (26.323 · 25.711) · **lỗi 1**: bộ so sánh vắt hai chi cho IRa/IRb bất đối xứng (3.008 vs 25.215 cột) — chuyển sang bộ **cùng chi**, nơi Coleus cho IRa = IRb = 0,001486 trên đúng 25.681 cột · **lỗi 2**: template `substitutions` in câu scope cứng ("whole aligned plastome") trên run `scope=cds` — nay in ô `scope` của chính run · bộ golden vẫn xanh: 2,6917 · 37 · 132 = 86 + 38 + 8 | [`docs/log/phase-2b.md`](docs/log/phase-2b.md) |
| 8.58 | **Module `repeats`** (P2B-02, D-f) — match chính xác cực đại ≥30 bp, hai lớp dispersed/palindromic: bản vàng **9 match** (4 dispersed · 5 palindromic), dài nhất **26.214 bp = chính cặp IR** (IRa ↔ IRb) · **hai engine độc lập trùng nhau từng dòng** (`repeat-match` + Vmatch 2.3.1) · **nguồn thật của `repeat-match` là gói `mummer4-4.0.1`**, không phải MUMmer3 như §7.8 ghi — sửa tại chỗ, và nó **không có cờ phiên bản** nên bản phát hành đọc từ `mummer` cùng thư mục · hậu tố `r` = đầu của bản sao **trong khung reverse-complement** (cặp IR là phép kiểm) · Vmatch chỉ miễn phí học thuật nên **không** là engine của mô-đun (bản cài tối thiểu giữ nguyên), nửa tandem của M2.12 vẫn mở | [`docs/log/phase-2b.md`](docs/log/phase-2b.md) |
| 8.57 | **Span intron & IGS vào `core/spans.py`** (P2B-01) — một phép suy, hai mô-đun dùng: `introns` in span, `diversity` đo π trên đúng span ấy · bản vàng: **17 intron** khớp từng dòng với bảng `introns`, **91 IGS** (79 LSC · 12 SSC · 32.134 bp) · π theo intron 0,004138–0,012658 · π theo IGS 0,000000–0,082353 · **ba cặp intron trong IR đo y hệt nhau** — lý do đo được để D-g loại IR khỏi tập spacer · định nghĩa IGS sửa sau khi ĐO: gen chiếm **exon** của nó, không chiếm khoảng giữa (span ngoài nuốt mất 3 spacer thật: `rps12-clpP` 138 · `trnK-UUU-matK` 260 · `matK-trnK-UUU` 789) · ba lỗi sửa kèm: tra vùng trộn hai hệ toạ độ · `diversity` không đi qua cổng IR chung · **`_check_span` không thể fail** (§9.H-6) | [`docs/log/phase-2b.md`](docs/log/phase-2b.md) |
| 8.56 | **Nghỉ hưu `chloe` · `pga`** (P3B-11, thi hành mặc định của bảng quyết định "retires when P3B-10 passes"): xoá 9 tệp / **2.406 dòng** — `annotation` · `external/chloe.py` · `external/pga.py` · `core/gff.py` (chỉ Chloe dùng) · 3 tệp test · fixture · `docs/tools/chloe.md`; **giữ** `tests/data/corpus/*.pga.gb` + `tools/corpus_bar.py` (đối chứng đã đo, nửa còn lại của thanh chuẩn) · **lỗ hổng lộ ra khi dọn**: vế "tìm không ra gì là run hỏng" của P3B-10 chỉ có test ở đường mito, câu "asserted on both paths" trong TASKS.md là **sai** — đã viết `test_plastid_annotate_module.py` (3 ca) **trước** khi xoá · registry 22 → 21 mô-đun, nhóm plastid 7 → 6 (tripwire `test_organelle` phải đọc bảng §13 mới sửa) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.55 | **Cấu trúc exon lấy từ chính bản tham chiếu** (`core/exon_structure.py`: junction = exon k + cả intron + exon k+1, tìm bằng blastn có chuỗi alignment nên ranh giới map chính xác kể cả qua gap; theo **từng cặp exon** nên `rps12` trans-spliced vẫn đặt được exon 26 bp; bốn điều nó từ chối làm) · stop codon đọc **cả khi alignment dừng sớm** (65 khác biệt còn lại đúng dạng này, `psaJ` 30) · **H-20**: cấu trúc quyết định giữa các ứng viên — loại ứng viên có protein là **suffix** của ứng viên khác (cờ chỉ đường plastid bật) · **THANH CHUẨN ĐẠT: 418 so với PGA 588** (30/34 tốt hơn), bản vàng 44 → 35 khác biệt — **kèm cảnh báo đo được: cách đọc strict là 1070 so với 699, KHÔNG đạt**, và 621/652 cái `boundary` là đầu tRNA của ARAGORN (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.54 | `lhbA` là **call trùng** tại đúng `psbZ` 35.895–36.083 (ba tín hiệu §8.9 khớp tuyệt đối; 33 bản ghi `psbZ` / 1 `lhbA`; bộ từ vựng 80 → 79 tên) · `rps12` vắt ranh giới IR nay được **viết hai lần** như bản công bố (chung exon 1; còn thiếu exon thứ ba 26 bp) · **con số 86 của bản vàng đang đúng vì một cái sai bù một cái sai** — sửa H-18 một mình làm nó rơi xuống 85 · corpus 692 → **681** (PGA 588): 22/34 bản ghi viết `rps12` một lần, nên bản sao thành extra — kể cả tính nó là khác biệt trình bày thì **659 so với 576, vẫn không đạt** (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.53 | Bản vàng: **45 khác biệt** so với bản công bố (18 kết luận · 19 boundary · 8 vocabulary) thành test so **tập hợp** — một khác biệt đổi, biến mất hay mới xuất hiện đều fail kèm hai danh sách; bốn khiếm khuyết đã biết nằm **trong** danh sách chứ không bị nó che (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.52 | ARAGORN `i(a,b)`: trường thứ hai là **độ dài** intron, không phải toạ độ kết thúc — đọc sai cho tRNA nối ra 104–134 nt (bất khả), đọc đúng ra 72–95 và khớp **chính xác** exon của bản công bố ở `trnL-UAA` (85) và `trnA-UGC` (73); corpus **811 → 692**, thanh chuẩn vẫn không đạt (PGA 588, nay đo đủ 34 bản ghi) · rủi ro mới **H-22**: trần intron 3.000 bp của `-i` im lặng, intron dài nhất đo được 2.671 bp (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.51 | Corpus chạy lại sau hai bản sửa (thiếu 167→47, thừa 161→41) · thước đo vào repo (`tools/corpus_bar.py`, 19 test cho chính luật, `make corpus-bar`) · ngưỡng `boundary` 1 bp khai ra kèm lý do đo được · **thanh chuẩn vẫn KHÔNG đạt**: 811 so với PGA 588, và 0/34 bản ghi tốt hơn · ba phát hiện: `lhbA` là call trùng (86 của bản vàng đúng nhờ nó bù cho `rps12` thiếu một bản sao) · `rps12` thiếu exon giữa 26 bp và cả bản sao IRa · exon đầu 6–40 bp của `petB`/`petD`/`rpl16`/`rps16` (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.50 | P3B-10 mở màn: ARAGORN bỏ qua tRNA có intron trong im lặng (30 → 38, bản vàng đạt 132 = 86 + 38 + 8) · bản sao IR bị đặt lệch đúng một chiều dài IR (số gen vẫn đúng nên không test nào bắt) · corpus 34 bản ghi chạy xong và **thanh chuẩn chưa đạt**: 1.579 khác biệt mức kết luận so với PGA 637 (P3B-10, 10/10/2026) | [`docs/log/phase-3b.md`](docs/log/phase-3b.md) |
| 8.66 | **Bộ Crataegus vào cây làm việc** (P8B-02) — 117 file copy vào `.data/`, **byte-identical** (sha256 từng file) · **100/100 bản ghi đọc được** bằng reader của mình · 50 plastome: tất cả khai `/organelle`, **tất cả qua cổng IR 100,00 %** · 50 barcode 543–625 bp không gen không IR (cổng chặn đúng nếu đem vào mô-đun plastome) · **Table 1 của bài khớp từng dòng** (vd *C. altaica* 159.664 · 87.665 · 19.245 · 26.377 · GC 34,37/30,42/42,65) · con số **26.311 trong text bài không tái hiện được** ở bất kỳ bản ghi nào kể cả Table 1 của chính bài (cận dưới 26.331) · **Table S6 (SSR) = MISA mặc định**, không phải ngưỡng Methods của bài ghi (57 so với 97; 57/57 khớp ở mặc định, 55/57 ở ngưỡng họ ghi) ⇒ D-l · thiếu: mẫu barcode **thứ 51** và **thiết kế Ka/Ks** (bảng Waiting của board) | [`docs/log/phase-8b.md`](docs/log/phase-8b.md) |
| 8.65 | **Cờ compound của MISA — một cụm, hai cách viết, bảng chỉ đọc được một** (P8B-01) — lỗi **thứ hai** trong hai ngày do một phép so với bên ngoài tìm ra: `_parse_gff` suy cờ từ "dòng này chứa >1 thành phần", nên mọi thành phần nằm **trên dòng riêng** với nhãn `compound_repeat` đều ra `no` ⇒ **PX413282: 0/16** thành phần được gắn cờ trong khi MISA báo 8 cụm (**bản vàng: 2/6** — nó chỉ có một cụm viết theo cách thứ nhất, nên che mất lỗi) · phép đối chiếu GFF-vs-statistics **không bắt được** vì so **số lượng** chứ không so **cờ** (cùng họ với §9.M-1) · sửa: nhận cả hai cách viết, test hồi quy viết trước, test cũ được **thu hẹp phạm vi** (câu "MISA ghi một cụm thành một dòng" của nó nay đo được là chưa đủ) · mô-đun **0.3.0**; bản vàng **không đổi số nào** (37 = 31·4·2), cờ đổi 2 → **6**, và ghi chú report nay nói rõ MISA đếm **mỗi cụm một lần** còn bảng đếm **mọi thành phần** | [`docs/log/phase-8b.md`](docs/log/phase-8b.md) |
| 8.64 | **CI trên mọi lần push** (P8-05, M8.4) — hai job: `fast` (ruff · test nhanh · **bộ vàng tier A**, thuần Python, không công cụ không mạng) và `golden-tier-b` (micromamba với **perl 5.32.1 · mafft 7.526 · mummer4 4.0.1 · paml 4.10.10 · blast 2.17.0 · aragorn 1.2.41 · vmatch 2.3.1**, `PTS_REQUIRE_TOOLS=1`, `table2asn` tải từ NCBI **kèm kiểm md5**) · **ba lỗi thật lộ ra ngay lần chạy đầu**: scipy 1.18 trả `pvalue=nan` (1.17 trả 1.0) và `kruskal` trả nan ở **cả hai** bản mà module `stats` **xuất thẳng số nan vào bảng FDR** → nay thành dòng *not computable* kèm lý do, hai test hồi quy · một test **hardcode `/root/toolbox/bin/table2asn`** (đường dẫn tuyệt đối bị cấm) ném `PermissionError` trên runner thay vì skip → dùng fixture `submission_tools` · env micromamba **không nằm trên PATH** của bước `run:` → đưa lên `GITHUB_PATH` + in từng phiên bản đã pin · **bốn lần chạy**: [đỏ vì ba lỗi thật](https://github.com/tridovan/plastid-suite/actions/runs/38024381329) → [xanh](https://github.com/tridovan/plastid-suite/actions/runs/38024744478) → **[đỏ vì cố ý phá số vàng](https://github.com/tridovan/plastid-suite/actions/runs/38024933018)** (`37.863555723244254 == 37.8536 ± 0.001`; `'37.85 %' not in report`) → [xanh](https://github.com/tridovan/plastid-suite/actions/runs/38025149827) | [`docs/log/phase-8.md`](docs/log/phase-8.md) |
| 8.63 | **`diversity` và `barcodes` ở N = 34** (P8-04, M8.3) — hai số đo cuối cùng, giao thức cũ (một mô-đun một lần · `--force` vì cache khoá theo `module_version` — H-08 · `/usr/bin/time -v`): `diversity` **13 ph 14 s · 138 MB**, π **0,040086** trên 120.964/191.707 cột (SSC 0,1211 · LSC 0,0439 · IR ≈ 0,0135) · `barcodes` **12 ph 27 s · 138 MB**, π trung bình 0,050020, **13 vùng siêu biến · 10 locus ứng viên · 8 tìm thấy** · **dự đoán cũ sai theo hướng đáng giữ**: *not measured* được viết vì `substitutions scope=whole` **bị dừng ở trần 11 phút**, mà **một cái trần không phải một phép đo** — hai mô-đun này **xong** trong ~13 phút · `docs/M8.3_scale.md` nay **không còn dòng "not measured"** nào | [`docs/log/phase-8.md`](docs/log/phase-8.md) |
| 8.62 | **Đối chiếu chéo với CGAS 1.0.1** (P8-03, M8.1) — bốn mặt trên bản ghi vàng, cùng bộ vào: **cấu trúc mô-đun 7 trùng khít** (LSC 83.127 · SSC 18.169 · IR 26.214 · GC 37,85/35,98/31,45/43,05) · **RSCU mô-đun 8: 26.263 lượt codon sense, 0 codon lệch** (RSCU lệch tối đa 5×10⁻⁵ = làm tròn khi mình xuất; CGAS đếm cả 86 CDS **và** in 3 codon stop 86 = 37+27+22) · **SSR mô-đun 12: mặc định 54 so với 37**, chạy lại `-t "10,6,5,5,5,5"` ⇒ **37 = 27·10·0 · 31/4/2 trùng khít** — khác biệt là **một tham số**, không phải bất đồng · **intron mô-đun 11: 17 intron CDS trùng từng dòng** (in trong 15 dòng) · **hai phát hiện**: (1) CGAS **im lặng** trên `NC_024936.1` — in `LSC 0 · SSC 0 · IR 0`, không cảnh báo, vẫn tính "2/2 genomes analysed", trong khi bên mình có LSC 82.702 · IR 26.184 (100,00 %) và **cả đường de-novo cũng tìm ra** ⇒ số 0 của nó là số sai, và đây là lý do cổng IR tồn tại; (2) **khoảng trống của mình**: CGAS in **8 intron tRNA** (trnK-UUU 2.617 bp…) mà **không mô-đun nào của mình in** — `introns` suy từ part của CDS ⇒ **H-23** · R không có trên máy, nhưng 8·9·12·13 **tính số bằng Python rồi mới gọi R để vẽ** (đo được: thiếu R thì skip hình, số vẫn ra) · không chạy: mô-đun 1 (cần read) · 2 (PGA, đã nghỉ hưu) · 14 (thiếu IQ-TREE/MACSE) | [`docs/log/phase-8.md`](docs/log/phase-8.md) |
| 8.61 | **Đối chiếu chéo với CPStools 3.0** (P8-02, M8.1) — năm mặt trên **cùng một bộ vào đã ghim** (`NC_050833.1` sha256 `a946f23e…`; cặp vàng cho π/KaKs), hai bên ghi tên + phiên bản + argv: **IR trùng khít** bốn khoảng (chỉ đổi nhãn IRa/IRb — quy ước) · **SSR 37/37 cùng toạ độ · cùng 31/4/2**, 2 dòng lệch là cụm compound `(TA)6+(A)10` (quy ước §8.1/§8.2) · **RSCU trên cùng bộ gen đếm trùng từng codon** (21.861 codon sense, 0 lệch; lệch RSCU tối đa 0,0049 = làm tròn 2 chữ số của họ) — khác biệt ở lần chạy mặc định là **tập gen** (họ khử trùng lặp theo tên 86→79, `-l 300` bỏ 28, bỏ trình tự không mở đầu ATG) + họ **tính cả codon stop** · **π 100/104 gen trùng tới 5 chữ số** (lệch trung bình 2×10⁻⁶), IGS 77/83, 6 dòng lệch là **phạm vi alignment** (đo: cùng một locus ba cách align cho 3 đáp số 1·2·6 vị trí khác) · **Ka/Ks**: Ka lệch 0,00014, chỗ lệch lớn là gen **Ks rất nhỏ** (psbK Ka 0,0080 vs 0,0073 nhưng Ks 0,0178 vs 0,0223 ⇒ tỉ lệ 0,4494 vs 0,3279), hai bên cùng gọi Ks không xác định ở cùng họ gen · **một lỗi của mình do chính phép so tìm ra và đã sửa**: cột `Context` của `ssr` sai 11/37 dòng vì đọc span ngoài của `rps12` trans-spliced (§9.M-1) · bốn khiếm khuyết của cpstools ghi ở `docs/tools/cpstools.md` · **Ts/Tv không nằm trong phép so này** (cpstools không có mặt ấy) | [`docs/log/phase-8.md`](docs/log/phase-8.md) |
| 8.60 | **CPGview — khảo sát, KHÔNG đưa vào bộ công cụ** (P8-01, D-h; lỗ hổng phủ của M8.1/§13.8 nay đã đóng): bài Liu và cs. 2023 · repo GitHub **chỉ là trang giới thiệu** (0 mã, 0 release, 0 giấy phép) · không pip/conda/CLI ⇒ không có gì để pin · web trả HTTP 200 ở trang ngoài nhưng **ứng dụng (cổng 16085) không kết nối được** — control `portquiz.net:16085` mở nên không phải mạng phía mình · bản cài duy nhất là container Figshare **2,87 GB** (md5 `4b371a5c…`) cho một công cụ **vẽ hình** · nhánh CPGview của M8.1 đóng bằng phép đo, không bằng một lần chạy | [`docs/log/phase-8.md`](docs/log/phase-8.md) |


---

## 9. Lỗi đã gặp & số đo (tổng hợp)

### A. Lỗi "sai mà không báo" ở khâu nhận diện IR

1. Cách viết đơn giản khớp chú giải NCBI **7/9** ca, **sai 2/9** và **không báo lỗi** (số ở §8.3).
2. **2/9 IR có hai bản đã phân hoá** (khớp tổng 92,63 %; ca thứ hai **đo lại ở §8.5 = 98,87 %**, không phải 53,65 % của bản đo đầu — cổng vẫn cắm cờ vì < 99,9 %) ⇒ phép so khớp **không cho phép gap** sai ở nhóm này.
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
   với các bản cùng chi trước khi dùng cho phân tích theo vị trí (việc mở: **H-13**).

### B. Đếm gen, SSR, Ts/Tv

9. **"135 gen" không tái lập được** từ bất kỳ cách đếm nào: 132 (CDS+tRNA+rRNA) · 133 (feature `gene`) · 112 (tên duy nhất) · 21 (gen lặp IR).
10. **1 bản `ycf1` cụt** gây lệch cách đếm: tại 108.448–109.347 chỉ có `gene` + `misc_feature` (không CDS); bản thật ở 122.984–128.404 mới có CDS.
11. Bản tự viết đếm SSR ra **36** (thiếu 1); MISA và cpstools đều ra **37**.
12. **"mono 10 · di 6 · tri 5 · tetra 4 · penta 4 · hexa 4" là ngưỡng MISA, không phải số lượng.** Phân loại thật: mono 31 · di 4 · tri 2.
13. **Số SSR phụ thuộc ngưỡng**: ≥8 → 133 · ≥10 → 36 · ≥12 → 16. `misa.ini` của anh Điềm dùng `(1/10)(2/6)(3/5)(4/5)(5/5)(6/5)`; tài liệu MISA ghi mặc định 4 cho unit 5–6 bp.
14. **Ts/Tv phụ thuộc cách căn**: difflib 2,667 vs MAFFT 2,686 — đo trên bộ 76 gen lúc đó;
    số vàng nay là **2,6917** trên 79 gen (§8.27). Phép đo này là về *aligner*, không về bộ gen.
15. **Mới 03/10/2026 — parser GFF của MISA đọc thiếu một nửa từ vựng (im lặng).** MISA ghi SSR
    monomeric là `monomeric_repeat`, nhưng **mọi** unit dài hơn là `microsatellite` (lớp nằm
    trong cột Note), và một cụm compound là **một** dòng `repeat_region` chứa nhiều thành phần.
    Parser cũ chỉ nhận `monomeric_repeat` ⇒ **mất 7/37 SSR** (bảng per-SSR còn 30 dòng) trong khi
    bảng tổng vẫn ghi 37 vì lấy từ file `.statistics` — hai bảng âm thầm lệch nhau.
    Phát hiện bằng cột đối chiếu của P2-06. **Bài học:** fixture test cũ tự bịa type
    `dinucleotide_repeat` mà MISA không bao giờ ghi ⇒ test xanh nhưng sai; fixture phải lấy
    **hình dạng dòng thật**.
16. **Mới 03/10/2026 — thứ tự các phần của CDS trans-spliced.** Sắp xếp lại các part, hoặc đọc
    tất cả theo mạch của gen, đều cho ORF sai. Đo trên **cả hai** bản rps12 của genome vàng:
    chỉ cách "đọc từng part theo mạch của chính nó, giữ nguyên thứ tự annotation" mới ra
    `ATG…TAA` (sắp xếp lại → `AAT…TAT`; đọc xuôi hết → `ATA…ATT`).
17. **Mới 03/10/2026 — khoá bảng trùng thì tham chiếu rơi vào một dòng bất kỳ.** Gen nằm trong IR
    có 2 bản nên `intron_id`/tên gen lặp; `Table` đã được vá để **từ chối** khoá trùng.
18. **Mới 03/10/2026 — `yaml.YAMLError` KHÔNG phải `ValueError`.** File YAML hỏng làm `load_groups`
    ném traceback thay vì thông báo; đã dịch lỗi tại biên.
19. **Mới 03/10/2026 — `annotated=True` là một lời khai chưa ai kiểm.** `entry_from_summary` cũ
    đặt `annotated=True` cho mọi bản ghi; thực đo trên bộ dữ liệu thay thế: **3/11 plastome có
    0 CDS** (GenBank không mặc nhiên có chú giải). Đã sửa: chỉ tin RefSeq, bản GenBank **đo sau
    khi tải**.
20. **Mới 03/10/2026 — dòng `Total` của bảng đối chiếu tự tố cáo chính nó.** Lần chạy đầu của
    `gene_table`: CDS 86/86 · tRNA 38/38 · rRNA 8/8 đều `agrees`, nhưng dòng **Total ra
    `in_table = 0`** — `by_kind.get("Total")` trên một `Counter` chỉ có ba khoá CDS/tRNA/rRNA.
    Đúng loại lỗi mà cột đối chiếu sinh ra để bắt: nếu bảng chỉ in số mà không tự so, một dòng
    tổng bằng 0 sẽ nằm giữa bảng và **không ai kiểm**. Đã sửa (`_in_table`), thành test.

21. **Mới 03/10/2026 — biến đổi trình tự làm đúng SỐ GEN mà sai TOÀN BỘ KHUNG (§8.17).**
    Chloe mặc định đảo chiều vùng sao chép đơn rồi chú giải hệ toạ độ đó. Hệ quả đo được: **132 gen
    — y hệt, tên y hệt**, nhưng **12/86 CDS hỏng ORF** khi đọc với trình tự người dùng đưa vào, và
    vị trí các gen ở SSC lệch hẳn (ndhH ở 114.734–115.915 (+) thay vì 120.937–122.118 (−)).
    **Bài học:** một con số đếm được không chứng minh được hệ toạ độ; phải có một phép kiểm *theo
    vị trí* (ở đây: dịch ORF) mới thấy. Đây là lý do P3-03 bắt buộc "đọc lại bằng reader của mình"
    chứ không tin output của annotator.
22. **Mới 03/10/2026 — annotator và bản ghi vàng gọi CÙNG một locus bằng hai tên khác nhau.**
    `ndhI` (Chloe) = `ndhL` (bản ghi vàng): cùng vị trí, **cùng protein**, và `/product` của bản
    ghi vàng ghi "subunit I" ⇒ thành cặp alias thứ năm. Ca thứ hai **chưa chốt**: `trnG-GCC` vs
    `trnG-UCC` — bản ghi vàng dùng **một tên cho hai bản tRNA khác trình tự**, mà bộ gen chuẩn có
    cả hai gen. Không tự sửa: đo anticodon bằng offset cố định chỉ khớp 11/38 tRNA ⇒ **không đủ
    bằng chứng**, đưa vào §10 cho PI quyết.

23. **Mới 03/10/2026 — ba thứ định dạng nộp bắt buộc mà model gene không mang (§8.18).**
    Chạy `table2asn` thật trên genome vàng: **13 lỗi**, sửa ba chỗ còn **7**. Hai chỗ đầu là
    *passthrough* (feature `source`, qualifier `/exception`) — bản ghi có, mình bỏ. Chỗ thứ ba là
    một quy ước: **CDS không mở đầu bằng mã khởi đầu hợp lệ thì phải ghi start partial**, nếu
    không NCBI báo `BadProteinStart`. Bài học: một model chỉ giữ những gì module *hiện tại* cần
    thì mỗi định dạng xuất ra sẽ phát hiện lại từ đầu những gì nó bỏ.
24. **Mới 03/10/2026 — "exit 0" không phải là tiêu chí của `table2asn`.** Nó **exit 2** khi chỉ
    *cảnh báo* về một qualifier lạ (`mol_type`, `db_xref`) trong khi `.sqn` vẫn đúng; bảng nộp
    thật của *Dioscorea alata* cũng exit 2 vì lý do đó. Nghĩa là: đọc `.val`/`.stats`, đừng tin
    mã thoát. Và ngược lại — `.tbl` thiếu feature `source` thì **vẫn exit 0** nhưng bị đánh 2 lỗi.
25. **Mới 03/10/2026 — phép kiểm mới chỉ có giá trị khi có ĐỐI CHỨNG.** 7 lỗi
    `GeneXrefStrandProblem` trông như lỗi của mình cho tới khi chạy **cùng validator** trên một
    bảng nộp đã được chấp nhận (Dioscorea) và thấy nó cũng đúng 7 lỗi cùng lớp. Không có đối
    chứng thì hoặc là sửa những thứ không hỏng, hoặc là bỏ qua những thứ hỏng thật.
26. **Mới 03/10/2026 — thư viện viết GenBank ghi ra cú pháp không ai dùng, và round-trip bằng
    chính nó thì vẫn xanh.** BioPython gặp qualifier trần (`/pseudo`, `/trans_splicing`) dưới dạng
    **chuỗi rỗng** thì ghi `=""`; gặp `None` mới ghi `/pseudo` trần. File `=""` vẫn đọc lại được
    bằng chính BioPython, nên một test kiểu "ghi rồi đọc bằng thư viện đã ghi" **xanh trong khi
    file khác cú pháp NCBI**. Bắt được chỉ vì có **bản ghi vàng thật** để đối chiếu cú pháp và vì
    **reader của mình** mới là thứ đọc lại. Bài học: round-trip bằng chính thư viện đã ghi không
    phải bằng chứng — phải đối chiếu với một file do NCBI viết (§8.20).

### C. Môi trường & vận hành

27. **`chloroplast[filter]` của NCBI bỏ sót bản ghi ghi "plastid"** (đo 03/10/2026, khi chạy
    Reference Panel thật): `Iridaceae[Organism] AND refseq[filter] AND chloroplast[filter]` → **71**;
    đổi sang `plastid[filter]` (hoặc `(chloroplast OR plastid)`) → **76**, khớp con số §8.3. Ca cụ thể:
    **NC_024936.1** (*Iris gatesii* — chính genome so sánh của bộ số vàng) có title
    *"Iris gatesii voucher T06-52 plastid, complete genome"*: khớp `refseq[filter]` và
    `Iridaceae[Organism]` nhưng **không** khớp `chloroplast[filter]`. Bài học: lọc theo `[filter]`
    của NCBI phải kiểm bằng số đếm, không tin tên filter.
    *Xác nhận kèm (cùng lần chạy):* `plastid panel pin` tải lại NC_050833.1 từ NCBI ra sha256
    `a946f23e…` — **khớp đúng byte** với dữ liệu vàng đang vendor trong repo; NC_024936.1 ra
    `fd441727…` cũng khớp. Nguồn dữ liệu vàng được xác nhận lại end-to-end.

28. **Máy 8 GB không đủ cho phần lắp ráp:** WSL2/Docker mặc định ~50 % RAM host → còn **1,5–2,5 GB mỗi việc**; 120.000 cặp read là hết RAM (phải hạ còn 20.000); flye cần ≥8 GB.
29. **Dò công cụ theo PATH là không đủ**: `misa.pl` ở `~/toolbox/src/misa/`; `nucmer`/`repeat-match` phải activate env `mummer4`; `trf` không có.
30. **Đường dẫn có dấu tiếng Việt** gây lệch đường dẫn giữa hai máy (thư mục "chạy Trí").
31. Script gửi kèm **không tới được máy anh Điềm** — anh làm ở thư mục khác (`D:\diem\bai bao tool\ke_hoach_new\chạy Trí\probe`).
32. `docker` không gọi được từ WSL trước khi bật WSL Integration (thông báo lỗi ở §8.4).
33. Chuỗi phân tích trên máy 2,8 GB dùng 331 MB và 9,8 giây.
34. **Mới 03/10/2026 — bộ tỉa style đọc sai chỗ, và bản render vẫn "trông ổn".** Style trong
    OOXML được tham chiếu bằng **phần tử con** (`<w:pStyle w:val="Heading1"/>`), không phải
    thuộc tính. Bản đầu đọc `element.get("w:pStyle")` ⇒ không thấy style nào được dùng ⇒ tỉa
    gần hết, kể cả `Title` và `Heading 1`. **LibreOffice vẫn render ra hình như bình thường**
    vì nó lùi về mặc định của chính nó — lỗi chỉ lộ ra khi mở PDF cạnh bản đúng, và chỉ bị bắt
    vì có **test cấu trúc** ("mọi style tài liệu dùng phải còn trong styles.xml"). Bài học:
    tỉa gói OOXML thì phải canh bằng cấu trúc, không bằng mắt.

---

35. **Mới 04/10/2026 — cache không biết code đã đổi, nếu version không đổi.** Cache của run layer
    nằm ở `.runs/.cache` — **dùng chung giữa các run** — và khoá cache gồm `module_id ·
    module_version · params · sha256(input) · IR fingerprint`, **không** gồm mã nguồn module. Đo
    được: sửa scope CDS của `substitutions` xong chạy lại, module trả về **`cached`** và phục vụ
    đúng payload cũ (thiếu cả bảng mới), vì `MODULE_VERSION` vẫn là 0.1.0. Luật của repo là
    *đổi hành vi thì bump version* — lần này luật đó bị vi phạm trong chính phiên, và cái bắt được
    là con số "cached" in ra ở dòng cuối. Đã bump 0.2.0 cho bốn module đổi hành vi và ghi thành
    việc mở (H-08): hoặc cache khoá thêm hash mã nguồn, hoặc CI kiểm rằng một module đổi hành vi
    mà không đổi version thì không được cache.

36. **Mới 04/10/2026 — log của công cụ bị GHI ĐÈ: mọi run chỉ giữ lần gọi cuối cùng.** Bộ đếm
    thứ tự log dùng `glob("*__{tên}__*.stdout.log")` — pattern này **không khớp** tên file nó sắp
    ghi (`001__mafft__stdout.log`, không có dấu chấm trước `stdout`), nên bộ đếm luôn bằng 0 và
    **mỗi lần gọi ghi đè lên lần trước**. Đo được trên một run thật: 616 lần gọi MAFFT ⇒ **1 file
    log**; cả run chỉ còn output của lần gọi cuối. Hệ quả: câu "output thô của công cụ đi kèm run"
    — nền tảng của "số tái lập được từ thư mục run" — **sai với mọi công cụ trừ cái cuối**, và sai
    từ Phase 1 mà không test nào bắt (test tier-B chỉ kiểm *một* file log tồn tại). Sửa pattern +
    test hồi quy (ba lần gọi ⇒ ba cặp file, nội dung lần đầu còn nguyên); đo lại: run mito từ 1
    lên **74 log MAFFT · 24 log yn00**. Ghi chú kèm: một module **cache hit** không chạy công cụ
    nào nên run đó **không có log** — trạng thái `cached` trong manifest là dấu hiệu để đọc.

37. **Mới 08/10/2026 — id của run mang mốc đồng hồ, nên TÊN thư mục run không phải thứ tự thời
    gian.** Đo được khi chạy `make check`: máy này (WSL2) **nhảy lùi ~1 giây** giữa hai lần chạy
    liên tiếp — manifest của run thứ hai ghi `started_at 09:54:26` và `finished_at 09:54:25`, còn
    run thứ ba nhận lại đúng mốc giây của run đầu nên được đặt hậu tố `-2`. Hệ quả: `sorted()`
    trên tên thư mục xếp `…-095425`, `…-095425-…-2`, `…-095426` — tức run **thứ ba** đứng trước
    run **thứ hai** — và test `test_second_identical_run_is_cached_and_identical` đọc nhầm thư mục
    của run thứ hai, đúng chỗ nó vừa xoá file để thử khôi phục từ cache ⇒ đỏ. **Bản thân sản phẩm
    không sai**: `RunLayout.create` vẫn bảo đảm "một run một thư mục, không ghi đè" nhờ hậu tố
    `-2` (đã có test riêng cho cả ca đồng hồ đứng yên lẫn ca đồng hồ lùi). Sai ở **giả định của
    test** — lấy thứ tự tên làm thứ tự thời gian; sửa bằng cách nhận diện thư mục của một run qua
    **hiệu tập hợp trước/sau**, không qua `sorted()`. Luật kèm theo: một thư mục run chỉ được
    nhận diện bằng chính nó, không bằng vị trí của nó trong danh sách.

### D. Chuyển chú giải mito bằng BLAST (P3A-03)

1. **`sseqid` bị cắt tại khoảng trắng** ⇒ một khoá chủ thể dựng từ bốn field sẽ **không bao giờ**
   khớp hit nào. Header tham chiếu kết thúc bằng tên loài nên chỗ cắt luôn nằm ở đó.
2. **`makeblastdb -parse_seqids` lỗi trên chính bộ tham chiếu đã commit** (`Duplicate seq_ids …
   LCL|NC_006581.1|SDH3|-|NICOTIANA`) và **viết hoa** id khi parse được. Flag trông "đúng chuẩn"
   nhất lại là flag phá pipeline.
3. **`blastx` chỉ nhận `-version`** — `--version`, `-v`, `-V` đều trả rc=1 kèm USAGE. `run_tool`
   chỉ thử ba cái sau ⇒ phiên bản bị ghi là "unknown" trong im lặng.
4. **SEG mask ăn mất điểm của gen thật**: atp9 137 → 67 bits, 21 locus khác cũng giảm. Mask được
   thiết kế cho query ít phức tạp; mitogenome thực vật thì toàn protein màng.
5. **Xếp hạng ứng viên theo bitscore tạo gen khảm.** `atp6`: alignment 418 aa khớp 75,6 % thắng
   một alignment khớp 93,4 % ⇒ call bị kéo **381 bp** về thượng nguồn, ôm luôn một ORF phân kỳ.
   Xếp theo identity: tổng sai lệch 13 gen còn lại **310 bp** thay vì **655 bp**.
6. **Hai chuỗi cùng một header khoá**: nếu nối HSP của cả hai vào một call, coverage bị chia sai
   mẫu số và hai phân tử khác nhau bị ghép thành một gen.
7. **Codon stop của tham chiếu không nằm trong alignment** (protein không có stop). Mở rộng +3 chỉ
   khi chuỗi xác nhận; mở rộng mà không kiểm là cách bịa ra một đầu gen, còn đi tìm stop là
   ORF-walking (D3A-02 cấm).
8. **Ba CDS `no_stop` của chính bản ghi đáp án** (ccmFc · ccmC · nad6) bỏ codon stop nằm ngay sau
   ⇒ cùng một gen, hai cách chú giải, lệch 3 bp. Không bên nào sai; phải **kể ra**, không chỉnh.

### E. Bộ gọi tRNA (ARAGORN, P3A-04)

1. **`aragorn -V` thoát mã 0 và in một câu lỗi** ⇒ `run_tool` ghi câu đó làm *phiên bản*. Cùng họ
   với ca BLAST+ chỉ nhận `-version`: bộ dò phiên bản tin vào quy ước `--version`/`-v`/`-V`, và
   hai công cụ liên tiếp đã cho thấy quy ước đó không tồn tại. Luật mới: một dòng chỉ là phiên
   bản nếu nó **có chữ số**.
2. **Đếm tRNA theo tên là sai**: 17 feature chỉ có 14 tên (3 tên lặp hai lần). Gộp theo tên ra 14.
3. **17 = 17 không có nghĩa là cùng một tập**: một locus thiếu, một locus thừa, một anticodon
   lệch. Một phép so chỉ đếm số lượng sẽ báo "khớp" cho cả ba sai khác đó.
4. **Gene vắt qua origin** bị ARAGORN in là `[258274,73]` — nếu ghi thẳng vào GenBank thì đó là
   toạ độ vượt quá độ dài bản ghi. Phải tách thành `join()`.

### F. Bộ chuyển chú giải trên một mẫu hoàn toàn mới (P3A-07)

1. **`AP045469.1` (313.585 bp, 0 feature): chỉ 23/31 CDS gọi được.** Bốn locus không có alignment
   nào ngay cả ở E-value 10, một locus chỉ có mảnh 40 aa, ba locus mảnh 21–52 aa ở E-value
   0,05–2,8. **Không phải ngưỡng** — hạ ngưỡng không tìm thêm được gì.
2. **Bài học về cách báo cáo:** một prototype trung thực phải in ra *cái nó không tìm thấy* kèm
   bằng chứng đã cố nới ngưỡng. "Không tìm thấy" khác "không tồn tại", và bộ chuyển chú giải chỉ
   trả lời được câu thứ nhất.
3. **23 call vẫn cho 13 CDS không qua phép kiểm ORF** (9 `no_start`, 3 `no_stop`, 1
   `internal_stop`) — trên một genome chưa ai chú giải thì đó là *thông tin*, không phải lỗi: nó
   nói chỗ nào alignment không xác lập được đầu gen.

### G. So sánh chuỗi — lỗi chỉ hiện ra khi chạy dữ liệu thật (P2-07/P2-08)

1. **MAFFT 7.526 in chữ thường** ⇒ bộ đếm chỉ nhận `ACGT` hoa ⇒ **mọi cột bị bỏ qua trong im
   lặng**: Ts/Tv ra *trống* (không phải sai), π ra 0. Unit test xanh vì chúng dựng alignment bằng
   tay; **adapter MAFFT chưa từng có test**. Bài học: một adapter không có test thì hành vi thật
   của công cụ chưa từng được kiểm — "đã kiểm bằng unit test" không đồng nghĩa với "đã chạy thật".
   Đã sửa: upper-case ở biên, và model **từ chối** chữ thường.
2. **Hai scope, hai con số, không được lẫn**: Ts/Tv **2,6917** (chỉ CDS · 79 gen · bỏ IR — số đo
   lại 04/10, §8.27; trước đó là 2,686 trên 76 gen) so với **1,8416** (cả bộ gen). Cùng một cặp genome. Mọi bảng phải mang nhãn scope — đây là ca chứng minh vì sao
   luật đó tồn tại.

### H. Một con số được ghi vào kế hoạch nhưng không thuộc về dữ liệu của nó (P4-01/D4-04)

1. **Check của pha 4 ghi "21 khối, %IDY 83,7–98,2 (mean 91,5)" cho cặp Iris — không tái lập được
   từ chính đầu vào của probe.** Quét đầy đủ tham số (§8.26): không tổ hợp `nucmer`/`delta-filter`
   nào cho 21 khối, và cặp cùng chi này thẳng hàng ở **98–99 %**, không phải 91,5 %. Dải số ghi
   trong kế hoạch khớp với một so sánh **khác chi** (đối chứng *I. domestica* × *O. diffusa*:
   81,89–98,84, mean 89,54). Con số đến từ dòng khảo sát công cụ §7.1 (đo 26/9, dữ liệu khác) rồi
   bị chép vào §13.4/ROADMAP kèm tên cặp genome — **một số đo đúng nhưng gán sai chỗ**.
   Hệ quả: check của pha sửa thành số **đo được** (7 khối · 89,66–99,11 · mean 97,37), kèm tham số
   và hệ toạ độ; bài học giống hệt luật "số vàng phải đo trên chính dữ liệu của mình" ở CLAUDE.md.
2. **File delta của MUMmer 4 đọc sai thì không báo lỗi** — một header cho cả cặp trình tự (không
   phải mỗi alignment một header), và gap mã hoá bằng `|d| − 1` cột ghép cặp rồi mới tách hướng.
   Đã khoá bằng fixture thật + test kiểm toạ độ và tổng lỗi từng khối (§8.26). Cùng họ với lỗi
   MAFFT chữ thường ở mục G: **định dạng của công cụ phải được đo, không được suy từ tài liệu**.
3. **"Gen không được chia sẻ" có thể chỉ là "gen bị viết tên khác"** (§8.27). `atpI`/`psaI`/`psbB`
   của bản ghi vàng được NC_024936.1 viết là `atpL`/`psaL`/`bspB` — cùng vị trí, cùng `/product`,
   và cặp `atpI` **giống hệt từng base**. Phép so tên làm ba locus biến mất khỏi **mọi** so sánh
   theo gen, và kéo theo **số vàng Ts/Tv §6** (76 → 79 gen, 2,686 → 2,6917). Chỉ mục checklist
   "liệt kê gen không chấm được" làm lộ ra. Bài học: một bảng chỉ in những gì chấm được sẽ giấu
   đúng loại lỗi này.
4. **Hai tên tRNA viết thường làm mất locus trong im lặng** (§8.29). `trnk-UUU`/`trnl-GUA` của
   NC_024936.1 là `trnK-UUU`/`trnI-GAU` — đo bằng trình tự, giống hệt từng base. Bộ chuẩn hoá chỉ
   sửa tiền tố `trn`/`rrn`, không đụng chữ cái axit amin, nên hai cách viết là hai cái tên và locus
   biến mất khỏi mọi bảng khoá theo tên. Đã thêm quy tắc cú pháp (chữ cái axit amin không phân
   biệt hoa/thường, trừ `trnfM-CAU`) và một alias `trnL-GUA → trnI-GAU`. Cùng họ với mục 3: **một
   bảng chỉ in những gì khớp tên sẽ giấu đúng loại lỗi này**.

5. **yn00 đọc ba base một lần và không kiểm khung** (§8.27, `docs/tools/paml.md`). Một gap 1 base
   làm nó báo stop codon giả rồi exit 255 — gặp may; khung lệch mà không đụng stop sẽ được chấm
   như thật. Phép kiểm khung là của mình (`core.align.codon_frame_breaks`) và chạy **trước** khi
   yn00 thấy dữ liệu; `accD` của cặp vàng là ca thật bị chặn (gap sau 1.445 base).
6. **Một cổng kiểm không thể fail thì chưa phải cổng kiểm** (P2B-01, 10/10/2026). `_check_span`
   của `introns` (P2-05) so "exons + gaps" với span của gen — và đẳng thức ấy **đúng với mọi cặp
   parts**, vì khe được tính bằng `min`/`max` nên phần chồng lấn tự triệt tiêu trong tổng. Nó chưa
   từng chạy và không thể chạy; cái nó để lọt đúng là ca nó sinh ra để bắt: parts chồng hoặc lồng
   nhau cho "intron" **âm**, và độ dài âm ấy được xuất ra như một con số. Thay bằng phép kiểm khe
   ≥ 1 bp — đo trên **131 bản ghi / 2.850 gen nhiều phần** trong cây: đúng **1** ca (OX459128.1,
   mito thay thế, parts lồng nhau −138 bp), nay dừng có tên. Cùng họ với mục 1–5: một phép kiểm
   chỉ có giá trị khi đã **đo được nó fail** — và test âm là cách duy nhất chứng minh điều đó.

### I. Chạy cả 34 bản ghi — hai lỗi mà dữ liệu tự chọn không có hình dạng để chạm (P4a-03)

Cả hai lỗi dưới đây **không phải** lỗi tính toán: chúng là những hình dạng dữ liệu thật mà bộ
fixture tự chọn chưa từng có. Genome vàng chạy đúng **vì may**, và đó chính là điều một bộ test
do người dùng tự chọn không kiểm được.

1. **Một cửa sổ 3 bp làm chết cả run `structure`.** *Chasmanthe floribunda* NC_072562.1 có LSC
   **81.003 bp**, nên phép chia cửa sổ 1000/1000 kết thúc bằng một window dài **3 bp** — `TTA`.
   Window đó không có G và không có C, `gc_skew` raise, module raise, **cả run chết** (1/34 bản
   ghi). Genome vàng sống sót chỉ vì LSC của nó dài 83.127 ⇒ window cuối **127 bp**, tình cờ có
   G/C. Sửa: window không có G/C là **không xác định** — ô để trống, cột cumulative cũng trống từ
   đó (một tổng chứa số hạng không xác định thì không xác định), summary thêm
   `n_undefined_windows` và một note nêu lý do. **Mức region vẫn raise**: cả một vùng không có
   G/C là dữ liệu hỏng, không phải chuyện làm tròn của phép chia. Module lên **0.2.0** (H-08).
2. **Một bản ghi chú giải cùng một gen hai lần làm chết cả run `orf`.** *Libertia pulchella*
   NC_072575.1 mang **hai** feature CDS `rps16` — `5741-5780,4684-4895` và `5737-5780,4684-4891`.
   Cả hai quy về cùng span `4684-5780`, nên cùng id `rps16:4684-5780`, và `Table` từ chối key lặp
   ⇒ run chết (1/34). Sửa: id duy nhất theo **feature** (hậu tố `#2`, `#3`…), **cả hai dòng đều
   được kiểm**, summary thêm `n_duplicated_cds` + note nêu tên. Không dòng nào bị bỏ: một chú giải
   trùng lặp là **phát hiện về bản ghi**, không phải rác của bộ kiểm. Module lên **0.2.0**.

Bài học chung, viết ra để lần sau khỏi phải học lại: **hai lỗi này chỉ hiện ra ở bản ghi thứ 5 và
thứ 12 của một bộ 34** — cả hai đều nằm trong dữ liệu công khai đã công bố, và cả hai đều biến
một bản ghi thành không phân tích được. Một bộ fixture do mình tự chọn sẽ không bao giờ có chúng,
vì mình chọn những bản ghi chạy được.

### J. Hai bào quan trong mã nguồn — bốn số đo (05/10/2026)

Đo khi soát lại luật "mỗi lần chạy đúng một bào quan" (D5-06) và D5-07 trong chính mã nguồn. Ba
trong bốn số thuộc **đường ty thể** ⇒ chúng được **ghi tên**: H-09 và H-11 **mở**, H-12 **hoãn**
cùng H-06/H-10 theo D5-07; số về bảng mã **sửa ngay** vì nó đi vào bảng của bản thảo lục lạp.

38. **`gene_names` không chọn từ vựng theo bào quan.** Bảng nhóm ở §13 xếp module này vào "Dùng
    chung" vì "chọn từ vựng theo `/organelle` của bản ghi" — nhưng `modules/gene_names.py` nạp
    `load_gene_aliases()` (bảng plastid) cho **mọi** bản ghi; hàm chọn từ vựng thật
    (`core/gene_names.aliases_for`) chỉ được `cds_by_canonical_name` dùng (substitutions · kaks ·
    popgen). Đo: `plastid run -m gene_names` trên `OX459128.1` + `AP045469.1` → **0** renamed ·
    **0** "two nomenclatures"; cùng lúc `mito_genes` trên đúng `OX459128.1` → **31/51** gen được
    gọi bằng ký hiệu chuẩn. Một module khai "dùng chung" mà trả lời bằng từ vựng plastid là kiểu
    "im lặng" quen thuộc ⇒ **H-12**.
39. **Bảng mã dịch: hai câu trả lời cho "lục lạp dùng bảng nào".** `core/codon.py` viết "Plastids
    use the standard code" và `codon_usage`/`amino_acid` in "NCBI table 1", trong khi
    `modules/annotation.py` khẳng định "plastids translate with NCBI table 11" — và D3-03 đã đo
    cả 1.099 CDS đều khai 11. Đo trên genome vàng: **5/86 CDS** mở đầu bằng mã không-ATG
    (`rps19` GTG ×2 · `rpl2` ACG ×2 · `ndhD` ATC) và `translate()` đọc chúng thành **V · T · I**,
    trong khi protein công bố mở đầu bằng M. Sửa 05/10: nhãn bảng nay ghi **bảng mà bản ghi
    khai** kèm quy ước "mã khởi đầu đọc như chính nó" (**không con số nào đổi**), và dòng D2-03
    vào §3. **Cùng gốc, đo cùng ngày:** `orf` chọn bảng mặc định **1** khi bản ghi không khai
    `/transl_table`, còn đường chú giải lục lạp đóng dấu **11** cho đúng sự im lặng đó — hai đường
    trả lời một câu hỏi theo hai cách (`src/plastid/data/defaults.yaml`; H-09).
40. **Artifact mito không khai `/organelle` (xác nhận lại H-11).** Đo trên
    `.panel/mito-prototype/AP045469.1.gb`: **0** feature `source`, **0** qualifier `/organelle`,
    `Genome.source == {}` — `read_fasta` không có gì để khai và `mito_annotate` chép nguyên
    `source` của input. Hệ quả: `aliases_for` trên chính artifact đó rơi về từ vựng plastid.
41. **Artifact mito không qua nổi chính bộ kiểm của `submission`.** `plastid run -m submission`
    trên artifact đó → **24 first-pass problems**, **23** là "CDS carries no /product" (module
    `mito_annotate` cố ý không viết `/product` — lý do có ghi trong docstring của nó), và các
    file **vẫn được ghi ra** kèm ghi chú "not ready to submit". Tức `submission` **báo, không
    chặn** — chữ "refuses" trong H-11 đã sửa lại cho đúng hành vi. Hai chốt (D3A-05 "chú giải +
    nộp" và "không tự bịa `/product`") va nhau ở đây; theo D5-07 ghi vào H-11 rồi hoãn.

### K. Bước đặt tên (`gene_normalize`) — lỗi đo được (08/10/2026)

1. **Artifact của `gene_normalize` trộn hai từ vựng trong một file.** `GENE_KINDS` là
   `("CDS", "tRNA", "rRNA")` (`core/model.py`) nên quyết định đổi tên chỉ sinh cho ba loại đó;
   feature `gene` — theo mô hình là một *span*, không phải một gen — được ghi lại **nguyên tên
   gốc**. Đo trên artifact của chính run §8.45 (`NC_024936.1.normalized.gb`): **6 feature `gene`
   còn nguyên từ vựng B**, đúng bằng 6 ca đổi tên của §8.44 — `atpL` · `psaL` · `bspB` ·
   `trnk-UUU` · `trnl-GUA` ×2 — nằm **cùng span** với CDS/tRNA đã đổi thành `atpI` · `psaI` ·
   `psbB` · `trnK-UUU` · `trnI-GAU`. Artifact của bản ghi vàng: **0** (nó là tham chiếu, không
   đổi tên gì).
   **Vì sao là lỗi chứ không phải lựa chọn:** cả hai bản ghi công bố đều khớp **132/132** feature
   với feature `gene` cùng span cùng tên (đo cùng ngày: 0 lệch, 0 thiếu) — artifact của toolkit
   là bên lệch chuẩn; CGAS, mô hình mà D5-08 bám theo, **có** đổi feature `gene` (D5-09
   Appendix A: "8 features changed per target — 4 tRNA + 4 `gene`"); và D5-09 §3 cấm hai từ vựng
   trong **một bảng** — ở đây chúng nằm trong **một file**, thứ người dùng mở ra xem.
   **Chưa sửa, và cách sửa phải đo trước:** khớp theo *span* hay theo *locus* là câu hỏi thật —
   feature `gene` có thể ôm trọn intron, và bản ghi vàng có **hai gen khác nhau dưới một tên**
   (`trnG-UCC`, §8.41/H-04) nên khớp theo *tên* sẽ đổi nhầm cả hai. Ghi tên **H-17** (TASKS.md).
   Hệ quả hôm nay là *hiển thị*, không phải *số*: `GENE_KINDS` loại feature `gene` khỏi mọi phép
   đếm và mọi mô-đun so sánh.

### L. Bộ gọi gen dùng chung (P3B-07 · P3B-08 · P3B-09) — sáu lỗi im lặng (08–09/10/2026)

Sáu lỗi dưới đây không do test nào bắt: chúng chỉ hiện ra khi **chạy thật** bộ gọi plastid trên bộ
tham chiếu 34 bản ghi — bốn lỗi đầu khi dựng bộ gọi và bộ trọng tài, hai lỗi sau khi **nối cổng
audit vào lệnh `annotate`** (P3B-09, §8.49). Cả sáu đều thuộc đúng loại "sai mà không báo" mà dự án
này tồn tại để chặn.

1. **`-max_target_seqs 500` cắt cụt bộ tham chiếu trong im lặng.** Hằng số 500 được chọn cho 287
   chuỗi của bộ tham chiếu mito ("nothing may be truncated away" — chú thích cũ nói đúng ý định,
   nhưng con số thì không theo kịp dữ liệu). Bộ plastid có **2.838** chuỗi. Đo 08/10/2026 trên
   genome vàng × 33 bản ghi leave-one-out: cap 500 → **929 HSP / 467 chuỗi**; cap = cỡ bộ →
   **4.382 HSP / 2.603 chuỗi** — cùng một lệnh, cùng một dữ liệu, BLAST **không nói gì**. Hệ quả
   nếu không đo: 63 locus công bố "không tìm thấy alignment nào", tức một bản chú giải thiếu
   gần nửa số gen mà trông vẫn hợp lệ. Sửa: `subject_count()` đếm chuỗi từ chính file bộ tham
   chiếu, `search_argv`/`blastp_argv` **bắt buộc** nhận ngưỡng ấy, và hằng số bị xoá — một ngưỡng
   không thể cắt cụt là ngưỡng duy nhất đáng giữ. Test:
   `test_the_target_cap_is_the_size_of_the_reference_set`.
2. **Một tên gen có khoảng trắng làm header FASTA bị BLAST cắt.** `NC_042827.1` không đặt `/gene`
   cho rRNA; reader rơi về `/product` (đúng như thiết kế), nên tên gen là `23S ribosomal RNA`, và
   header `…|23S ribosomal RNA|rrn23` về tới nơi là `…|23S` — **dưới ngưỡng ba field**. Parser từ
   chối và gọi tên file, thay vì khớp một khoá không bao giờ khớp entry nào. Sửa:
   `ReferenceEntry.key` bỏ khoảng trắng trong từng field và `_header` dùng lại chính `key` đó —
   hai chỗ viết cùng một khoá là hai chỗ có thể lệch nhau. Test:
   `test_a_name_with_a_space_gets_a_header_a_blast_hit_can_come_back_as`.
3. **Tên rRNA dạng câu sản phẩm chưa được chuẩn hoá.** `16S ribosomal RNA` là *cùng locus* với
   `rrn16`, nhưng vào bộ từ vựng như một tên riêng — bốn locus không run nào gọi tên được, và
   bảng `not_called` sẽ báo chúng "no alignment" trong khi alignment vẫn ở đó. Sửa: `gene_names`
   mở rộng luật đã có cho dạng `rrn16S` (cùng một tri thức: lớp kích thước nằm trong chuỗi) cho
   dạng câu; `rrnS`/`rrnL` của ty thể vẫn **không** bị chạm.

4. **Một locus bị chú giải hai lần: bản cụt thắng bản đầy đủ.** `_choose` xếp theo identity
   trước, và một bộ tham chiếu mang `ycf1` **cụt** — bản sao nằm trong IR, ~300 residue — khớp mẫu
   ở **98,83 %**, cao hơn mọi ứng viên đầy đủ (97,2–97,7 %), nên `ycf1` của bản ghi vàng được gọi
   thành **900 bp** ở 108.448–109.347 thay vì 5.421 bp ở 122.984–128.404. Đúng ca mà D5-08 §3 gọi
   tên: **cấu trúc quyết trước, điểm chỉ tie-break**. Sửa: một locus mà **một bản ghi** chú giải
   hai lần, bản ngắn dưới 50 % bản dài, là bản sao bị cụt — không vào chỉ mục tham chiếu, và
   **được nêu tên** (`truncated_annotations`); ngưỡng 0,5 nằm giữa khoảng trống đo được
   0,2992 → 0,9927 trên 243 locus trùng của bộ đã pin. Đo lại trên bản ghi vàng: **đúng một dòng
   đổi**, và nó khớp bản công bố từng base (§8.48). Test:
   `test_a_locus_a_record_annotated_twice_is_searched_once` ·
   `test_ycf1_lands_on_the_published_gene_not_on_a_truncated_reference`.

Hai số dưới đây **cũng** chỉ hiện ra khi chạy thật, lần này là khi nối cổng audit vào lệnh
`annotate` (P3B-09, §8.49). Chúng thuộc cùng họ "sai mà không báo", nhưng khác gốc: cái thứ nhất
là **cấu trúc giữa các ứng viên**, cái thứ hai là **một luật BLOCK quá rộng**.

5. **Vế "cấu trúc trước" giữa các ứng viên chưa được thi hành — và giá của nó đo được.** Khi bản
   ghi bị cổng chặn (`NC_024936.1`) bị bỏ ra khỏi bộ tham chiếu, tập call của bản ghi vàng **đổi
   đúng một dòng**: `rps19` từ **279 bp** ở 83.173–83.451 (khớp bản công bố) thành **225 bp** ở
   83.173–83.397, dù identity vẫn 100,00 %. Nguyên nhân đo trên cả 34 bản ghi: `rps19` là
   **27 × 279 bp/92 aa · 2 × 285 bp/94 aa · 3 × 225 bp/74 aa · 1 × 69 bp/22 aa** — ba bản ghi
   (`NC_081929.1` · `NC_081930.1` · `NC_081932.1`) mang một `rps19` **bị chú giải cụt**, và khi
   bản đầy đủ vắng mặt thì một bản 225 bp thắng rồi kéo call về độ dài của chính nó. H-19 mới sửa
   vế **trong-một-bản-ghi** (bản sao cụt không vào chỉ mục); vế **giữa các bản ghi** — D5-08 §3
   "cấu trúc quyết trước, điểm chỉ tie-break" — thì chưa. **H-20**, việc của bộ trọng tài, và là
   một khác biệt P3B-10 phải phân loại.
6. **Luật BLOCK "toạ độ partial" từ chối một bản ghi có protein đầy đủ.** `NC_024936.1` mang
   **hai** feature `rps19`, mỗi cái **279 bp / 92 aa** với `/translation` trọn vẹn; cờ partial
   (`>83020`, `<153124`) mô tả gene nằm ở ranh giới IR, không phải một protein cụt. Bằng chứng:
   **chính bản ghi ấy** là bản ghi cho ra `rps19` khớp bản công bố từng base (100,00 %, không cờ
   partial nào trên call) — nên lý do của luật ("một reference partial đóng dấu partial lên mẫu",
   §8.15 trap 1) không đúng cho dạng này. **H-21**; sửa một luật BLOCK đã chốt ở P3B-03 cần đo lại
   toàn bộ, không làm lặng lẽ trong một task khác.

**Chưa sửa, đã đặt tên — H-18:** `lhbA` (`NC_042827.1`) và `psbZ` (33 bản ghi còn lại, kể cả bản
ghi vàng) là **một locus**, cùng 62 aa; bộ từ vựng vì thế có 80 tên CDS cho 79 locus, và bộ gọi
gọi cả hai. Gốc là `data/gene_aliases.yaml` thiếu một mục — việc của tầng đặt tên (D5-09), không
phải của bộ gọi; bộ gọi báo đúng như nó thấy.

### M. Kiểm thử chéo (P8-02) — lỗi đo được (10/10/2026)

1. **Cột `Context` của module `ssr` đọc "span ngoài" của gen, nên gen trans-spliced nuốt cả
   genome.** Bản ghi vàng chú giải `rps12` là gen **trans-spliced**: exon 1 ở LSC (68.975–69.089),
   exon 2–3 ở IR, và **cả ba exon được ghi trên mỗi feature CDS** — nên hai feature ấy có span
   ngoài 68.975–98.013 và **68.975–139.634**. Hàm `_gene_context` cũ lấy `min(part.start)` …
   `max(part.end)` làm phạm vi gen rồi trả "Intron" cho bất cứ span nào nằm trong đó mà không nằm
   trọn một exon ⇒ **11/37 dòng SSR của bản vàng bị dán nhãn sai**, trong đó SSR nằm **trong**
   `ycf1` (123.879 · 125.538 · 125.713 · 126.493 · 126.645), `ndhA` (118.930), và ba spacer
   `rpl32-trnL-UAG` (112.532 · 113.049), `ndhF-rpl32` (111.933 · 111.976), `rps3-rpl22` (82.711)
   — tất cả đều thành "Intron". Mười một dòng ấy nằm giữa 68.975 và 139.634, đúng span ngoài của
   feature `rps12` thứ hai. **Phát hiện bởi P8-02**: bảng của CPStools in đúng tên gen cho
   từng dòng, nên 11 dòng lệch lộ ra ngay (§8.61). Sửa theo đúng định nghĩa `core/spans.py` đã có
   từ P2B-01 — khoảng giữa hai part là **intron chỉ khi hai đầu cùng một vùng**, cắt ranh giới
   vùng là **junction** — và kiểm exon **trước**, trên mọi gen, vì intron của gen này có thể chứa
   gen khác (`matK` nằm trong intron `trnK`). Test hồi quy: `tests/unit/test_ssr.py`, 5 ca, viết
   **trước** khi sửa. Sau khi sửa: 37/37 dòng khớp cpstools, trừ đúng 2 dòng của cụm compound
   (quy ước đã ghi ở §8.1/§8.2, không phải lỗi).
   *Vì sao lọt lâu:* cổng IR không bắt được (gen trans-spliced là chuyện của **chú giải**, không
   phải của hai bản IR), và không test nào đọc cột `Context` — bảng vàng tier B chỉ kiểm **số đếm**
   (37 · 31 · 4 · 2). Một cột không ai kiểm là một cột không ai biết nó sai.

## 10. Câu hỏi còn treo

1. ~~Phạm vi con số Ts/Tv~~ — **ĐÃ CHỐT 02/10**: cổng dùng **chỉ vùng CDS · gen dùng chung · căn bằng MAFFT**, **không tính vùng lặp IR**. Bộ gen dùng chung được đo lại 04/10 **bằng locus, không bằng cách viết tên** ⇒ **79 gen, 2,6917** (323 Ts / 120 Tv, 68.448 cột — §6, §8.27). Số Ts/Tv của bộ Eleutherine (1,095–1,777) là phạm vi khác (SNP cả plastome) — chỉ dùng cho nội dung riêng của bộ đó và **phải dán nhãn phạm vi** ở mọi chỗ xuất ra.
2. ~~Ngưỡng `misa.ini`~~ — **ĐÃ ĐÓNG 02/10**: file `misa.ini` đang dùng trên máy anh Điềm **chính là file gốc trong zip của IPK** (md5 `49a1c56048631647c0364532a6964cf9` khớp bản tải lại): `1-10 2-6 3-5 4-5 5-5 6-5`, interruptions 100, GFF true ⇒ **dùng đúng mặc định của MISA**, không phải tuỳ chỉnh.
3. Ai có máy Mac để kiểm.
4. **Dữ liệu còn thiếu:** có lấy được reads Illumina thật của 1 mẫu không (để đo ngưỡng bộ nhớ — hiện mới có read mô phỏng).
5. **Dữ liệu còn thiếu:** có lấy được **genome nhân** của mẫu không — NUPT/MTPT cần genome nhân để chạy (hiện chưa có).
6. ~~Hệ toạ độ sau chuẩn hoá~~ — **ĐÃ CHỐT 02/10**: bảng theo từng genome ghi **toạ độ GenBank**; hình so sánh chéo nhiều genome dùng **hệ đã chuẩn hoá** (các vùng thẳng hàng nhau); mọi bảng ghi rõ hệ đang dùng; `rotation_offset` lưu trong `00_run_info`.
7. ~~Cách sinh report~~ — **ĐÃ CHỐT 02/10**: sinh từ đối tượng `Result`. Mỗi số trong report trỏ về **đúng một ô** bảng gốc; template chỉ trỏ (`{{bảng.cột}}`, `{{bảng[dòng].cột}}`), không gõ số; ref không giải được ⇒ **raise**, không in trống. Yêu cầu "kiểm thử truy vấn số" (§13.5) được thoả bằng cấu trúc, không bằng kỷ luật.
8. ~~Reference Panel: snapshot hay truy vấn sống~~ — **ĐÃ CHỐT 02/10**: engine gợi ý qua **cache RefSeq có ngày** (lệnh làm mới tường minh); người dùng chốt danh sách thì **tải + ghim checksum**; chạy lại offline, tái lập được.
9. Mục **(c)** của Reference Panel còn thiếu — đoạn văn bản nguồn (kế hoạch 22/9) **bị cụt đúng
   ở đây**, cần cô Trang bổ sung — xem §13.1, mục (c)
10. **Tên gen `trnG` in trong bảng bản thảo là tên nào (§8.17).** Bản ghi vàng NC_050833.1
    dùng `trnG-UCC` cho **cả hai** bản tRNA glycine (hai bản khác trình tự), Chloe gọi bản không
    intron là `trnG-GCC`, và bộ tRNA chuẩn của plastome có **cả hai** gen. Đây là chuyện *chữ in
    trong bản thảo / bản nộp*, tức việc của cô Trang (D3-06), không phải việc dev tự quyết — dev
    đã ghi lại đủ bằng chứng để cô Trang quyết trong một lần xem.
11. **Câu do máy sinh được phép khẳng định điều gì (D5-05).** Mặc định code đang chạy: chỉ nhắc
    cột mà bảng nguồn mang, kèm n và tên kiểm định; **không** dùng từ ngữ nhân quả/thích nghi.
    Đây là **phạm vi một tuyên bố khoa học** nên thuộc cô Trang — chờ chốt trước khi P5-03 viết
    câu diễn giải (D5-05).
12. **Chính sách dữ liệu cho đường lắp ráp từ xa (M7.4).** usegalaxy.eu là máy chủ **công khai**;
    việc đưa reads **chưa công bố** của đề tài lên đó **chưa được quyết định** (08/10/2026). Phải
    chốt trước khi dùng đường từ xa — nếu không đồng ý, đường từ xa phải là máy chủ của trường/viện.

## 11. Nguồn tài liệu (ngoài file này)

Nằm trong `D:\business\pts\`:

| File | Nội dung |
|---|---|
| `Ke_hoach_xay_dung_bo_cong_cu_phan_tich_bo_gen_luc_lap_22.9.2026_VuThiHuyenTrang.docx` | Kế hoạch gốc (91.729 ký tự) |
| `Bao_cao_tool_thay_the_NguyenThanhDiem.docx` | 5 + 3 + 12 công cụ bị chặn đều có thay thế đã chạy thật |
| `Danh_sach_tool_khong_cai_duoc_NguyenThanhDiem.docx` | 84 công cụ đã thử · 64 chạy · 17 không chạy |
| `REVIEW_2026-09-26_3_tai_lieu_va_muc_tieu_may_pho_thong.docx` | Bản kỹ thuật đầy đủ: mục 9 (khả thi theo khối) · mục 10–12 (đo rủi ro đường nối, chạy thử liền mạch, đóng vòng bằng công cụ chuẩn) |
| `BAO_CAO_4_CAU_HOI_PlastidToolSuite_2026-09-26.docx` | Bản cho người không chuyên (khả thi / vấn đề / giải quyết được không / giải quyết thế nào) + mục 8 chiến lược tích hợp |
| `KE_HOACH_TRIEN_KHAI_PlastidToolSuite_2026-09-26.docx` | Kế hoạch triển khai (môi trường dev, kiểm thử, đóng gói, lộ trình 8 tuần — **dự kiến**, không phải hạn cam kết; xác nhận 08/10/2026) |
| `HUONG_DAN_SETUP_MOI_TRUONG_CODE_2026-10-01.docx` | Hướng dẫn setup WSL / Claude Code / Docker |
| `thanhdiem_result/` | Trả lời phiếu việc của anh Điềm + kết quả thô (`.misa`, `.statistics`, cpstools, MAFFT per-gene) |
| `probe/` | Script + dữ liệu dùng để đo các lỗi ở §8–9 (`chain_e2e.py`, `chain_tools.py`, `ir_gate.py`, `iridaceae20.gb`) |

## 12. Lệnh cài cho máy mới

**Nguồn duy nhất, có pin phiên bản: [`environment.yml`](environment.yml)** — mọi công cụ ở đó được
ghim đúng bản mà §7.8 ghi, và đó là thứ `make env` cùng image Docker dùng.

*Khối lệnh cũ (nguyên văn báo cáo anh Điềm, **không pin phiên bản** — chính là lỗi §9.B-14 đã đo)
đã được gỡ ngày 07/10/2026 để không ai copy-paste nó. Bản gốc 26/9 giữ ở §8.1 làm lịch sử khảo sát;
nguồn và cách cài từng công cụ: §7.8.*

---

## 13. ROADMAP — thứ tự phát triển

### 13.0. Sơ đồ pipeline & công cụ

*Nguồn: bản này là nguồn chuẩn. `docs/D5-09_reference_naming.md` §1 giữ một bản sao để đọc độc lập — sửa sơ đồ thì sửa ở đây trước.*

```
   NGUỒN THAM CHIẾU
   NCBI RefSeq  ·  bản ghi nhóm tự nộp
                │
                │  plastid panel pin        ← accession · ngày · sha256
                ▼
   ┌──────────────────────────────────────┐
   │  DB THAM CHIẾU      (nhiều file .gb) │
   └──────────────────┬───────────────────┘
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐÁNH GIÁ CHẤT LƯỢNG                 │
   │   · file nào KHÁC so với phần còn lại│
   │   · file nào không dùng được (BLOCK) │
   └──────────────────┬───────────────────┘
                      │  báo cáo
                      ▼
   ┌──────────────────────────────────────┐
   │  NGƯỜI DÙNG ĐỌC VÀ CHỌN              │
   │  1 file .gb tốt nhất                 │
   └──────────────────┬───────────────────┘
                      ▼
                ★ CHUẨN TÊN    →    từ vựng A
                      │
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐỒNG NHẤT TÊN CHO CẢ DB THAM CHIẾU  │   gene_normalize · D5-09
   └──────────────────┬───────────────────┘
                      ▼
        DB THAM CHIẾU ĐÃ CHUẨN HOÁ   (A)
                      │
                      │              reads (Illumina · ONT/HiFi)
                      │                │   fastp · BWA/minimap2 · Bandage
                      │                ▼
                      │       ┌─────────────────────────┐
                      │       │ 1  LẮP RÁP              │  GetOrganelle + SPAdes (short read)
                      │       │                         │  OATK · raven (long read)
                      │       └────────────┬────────────┘
                      │                    │  FASTA
                      │                    ▼
                      └───────────►┌─────────────────────────┐
                                   │ 2  CHÚ GIẢI             │  công cụ tự viết (D5-08)
                                   │                         │  blastx · blastn  [BLAST+]
                                   └────────────┬────────────┘
                                                │  GenBank — tên theo A
                                                ▼
                                   ┌─────────────────────────┐
                                   │ 3  KIỂM TRA TÊN         │  gene_normalize (chế độ kiểm)
                                   └────────────┬────────────┘
                                                ▼
                        ╔═══════════════════════════════════════╗
                        ║   BỘ INPUT CHUẨN HOÁ   (từ vựng A)    ║
                        ╚═══════════════════┬═══════════════════╝
                                            │
                             ┌──────────────┴───────────────┐
                             ▼                              ▼
                    [4]  PHÂN TÍCH                [5]  RNA EDITING (C.3)
                    cấu trúc · gen · codon         ◄── reference RIÊNG (B)
                    SSR (MISA) · repeat (vmatch)        ghép cặp có bằng chứng
                    synteny (MUMmer4)                   → alignment → chưa chọn
                    Ka/Ks (PAML) · Ts/Tv (MAFFT)
                    barcode · popgen · pycirclize
                             │
                             ▼
                    [6]  XUẤT BẢN
                    table2asn (.sbt) · Excel · Supplementary
```

**Đọc sơ đồ:** cột trái là **đường cấp tham chiếu** — không phải bước tuần tự của mẫu, mà là nguồn chảy vào bước 2. Mũi tên dài bên trái chính là *"DB tham chiếu đã chuẩn hoá là thứ chú giải cho các mẫu"*.

Ba điều sơ đồ cốt để không quên: **(1)** đồng nhất tên chạy **trước** chú giải và chạy **trên DB tham chiếu** — mẫu thừa hưởng từ vựng A qua bước chú giải, không bị đổi tên sau đó; **(2)** bước 3 là **cùng module chạy lại ở chế độ kiểm** — nếu annotator lấy tên từ bộ đã chuẩn hoá thì đầu ra phải đã đúng A, nên **mọi thay đổi ở bước 3 là một phát hiện**; **(3)** RNA editing là nhánh **duy nhất có reference riêng**, nằm *sau* bộ input chuẩn hoá, nên reference B phải ghép cặp về A.

| # | Bước | Công cụ | Vào → Ra | Reference cần |
|---|---|---|---|---|
| 0 | Nguồn tham chiếu | `plastid panel refresh` · `propose` · `pin` | NCBI / file nhóm → DB đã ghim (`PINNED.json`, sha256) | — |
| 0b | Đánh giá chất lượng | audit (`docs/D5-08_plastid_annotator.md`) | DB → báo cáo BLOCK · REPORT | — |
| — | **Chuẩn hoá DB tham chiếu** | **`gene_normalize`** (D5-09) | DB + ★ file đã chọn → DB theo từ vựng A | **★ file chuẩn tên** |
| 1 | Lắp ráp | GetOrganelle + SPAdes (short read) · OATK · raven-assembler (long read); fastp · BWA/minimap2 · Bandage (QC) | reads → FASTA | seed DB `embplant_pt` · HMM DB `embryophyta_pltd.fam` |
| 2 | Chú giải | **công cụ tự viết** (D5-08) · BLAST+ 2.17.0 | FASTA → GenBank | **bộ tham chiếu, 1..N bản ghi** |
| 3 | Kiểm tra tên | **`gene_normalize`, chế độ kiểm** | GenBank → GenBank đã kiểm | **★ file chuẩn tên** |
| 4 | Phân tích | MAFFT · MISA · MUMmer4 · PAML · vmatch · pycirclize | A → bảng · hình | loại 1, chọn theo taxonomy (D5-09 phụ lục B) |
| 5 | RNA editing | *chưa chọn* (C.3) | A → vị trí editing | **một reference riêng B** |
| 6 | Xuất bản | table2asn | bảng → `.sqn` | template `.sbt` |

**Sửa 07/10/2026 — kiểm lại sơ đồ ở mức *công cụ*, không phải mức ô** (luật phủ: `ROADMAP.md` → "The pipeline"): bỏ **`Flye`** — đo được là chết vì RAM trên máy 8 GB và đã thay bằng `raven-assembler`/`OATK` (§7.1), hai công cụ này nay thuộc **M7.6** (đường long-read, đang chờ read ONT/HiFi); bỏ **`Sibelia`** — khảo sát 26/9 có chạy thật trên 3 plastome nhưng tuyến synteny của dự án là `nucmer` (D4-02, M4.2) và không mục nào cần khối LCB. `BWA/minimap2` (bước map-back ở ô QC) thuộc **M7.1**; template `.sbt` của nhóm là một dòng "đang chờ" ở `TASKS.md`. Đợt kiểm 07/10 trước đó chỉ đọc ở mức ô nên bỏ sót cả ba.

---

**Thứ tự đã đổi theo quyết định 01/10: phân tích trước, lắp ráp sau.** Làm từ trên xuống; mỗi giai đoạn ra một thứ dùng được ngay.

**Điều kiện chung cho mọi giai đoạn** (đọc trước khi làm bất kỳ giai đoạn nào):

- **Áp dụng cho cả lục lạp và ty thể** — trừ phần ước tính thời gian phân kỳ (chỉ cần cho bài báo công cụ). Riêng **bảy phân tích chuyển tiếp** (bảng nhóm bên dưới) hôm nay **chỉ chạy được trên lục lạp**: chạy trên mitogenome thì dừng ở cổng IR trước khi mô-đun chạy — mở cho ty thể là H-10, **tạm hoãn** theo D5-07.
- **Cổng bắt buộc ở mọi phân tích theo vị trí: hai bản IR phải giống ≥ 99,9 %.** Không đạt ⇒ cờ đỏ, **dừng**, không cho phân tích theo vị trí chạy tiếp. Lý do đo được: cách nhận diện IR thông thường sai âm thầm **2/9 ca (22 %)**, cổng này bắt đúng 2/2 ca sai và báo oan 0/7 ca lành (chi tiết §9 mục A).
- **Cách đếm gen (không ghi trần "132 gen"):** **132 = 86 CDS + 38 tRNA + 8 rRNA**, không dùng feature `gene` (feature `gene` = 133, thừa 1 bản `ycf1` cụt; tên gen duy nhất = 112; gen lặp trong IR = 21). Khớp bài công bố Ai và cs. 2019.
- **Phạm vi con số Ts/Tv dùng làm cổng:** **chỉ vùng CDS · 79 gen dùng chung · căn bằng MAFFT · không tính IR** ⇒ **2,6917** (323 Ts / 120 Tv, 68.448 cột căn). Phạm vi chốt 02/10; **bộ gen dùng chung đo lại 04/10 bằng locus** — ba locus viết tên khác (`atpI`·`psaI`·`psbB` ↔ `atpL`·`psaL`·`bspB`) được tính vào (76 → 79 gen; 2,686 → 2,6917, §8.27). Số của bộ Eleutherine (1,095–1,777) là **khác phạm vi** (SNP của cả plastome) — không được lấy ra đối chiếu trực tiếp.
- Mỗi lần chạy phải xuất **report công cụ + thông số** (tên công cụ · phiên bản · đường dẫn · tham số/ngưỡng · thời gian) và **chia output theo thư mục** để người dùng truy lại — yêu cầu đã chốt 02/10 (§2). Đây là điểm CGAS đang thiếu (§7.7).
- **Chuẩn hoá bộ gen dùng chung**: xác định điểm bắt đầu/chiều và tách 4 vùng LSC/IRa/SSC/IRb — làm ở giai đoạn 1 rồi mọi giai đoạn dùng lại.
- **Chạy ở đâu**: các giai đoạn 1–5 và 8 chạy trên **máy phổ thông**; giai đoạn 6 (thời gian phân kỳ) cần **máy chủ từ xa** (máy 8 GB không đủ — §9.C-28). **Giai đoạn 7 (lắp ráp) tuỳ máy người dùng** (xác nhận 08/10/2026): đủ RAM thì chạy **tại chỗ**, không đủ thì dùng **đường từ xa Galaxy** (§13.7); chính sách đưa dữ liệu chưa công bố lên máy chủ công khai **chưa quyết định** (§10.12).
- **Ca kiểm thử chuẩn**: bộ dữ liệu **Eleutherine bulbosa** + **bộ số vàng** ở §6.
- **Kiểm thử ngay trong từng giai đoạn, không dồn về cuối**: mỗi giai đoạn phải có mục *“Kiểm tra trong giai đoạn”* và chỉ được coi là xong khi số kiểm tra khớp; sai thì **dừng lại sửa ngay trong giai đoạn đó**, không để lỗi trôi sang giai đoạn sau. Giai đoạn 8 chỉ làm phần **kiểm thử chéo**.
- **Chưa chốt (dev tự quyết khi lập kế hoạch với Claude)**: cách tổ chức mô-đun/plugin, định dạng dữ liệu trung gian giữa các bước, và quy trình chạy tự động nhiều bước.

**Hai bào quan — mỗi lần chạy đúng một bào quan (D5-06).** Lục lạp và ty thể dùng chung khung, nhưng một lần chạy chỉ phục vụ **một** bào quan: xác định theo qualifier `/organelle` của chính bản ghi; chuỗi trần không có bản ghi để hỏi thì bộc lộ bằng **đường vào** (`plastid mito annotate` là đường ty thể). **Bộ dữ liệu đầu vào cũng một bào quan: không trộn.** Một lần chạy nhận nhiều bản ghi (`-i` lặp lại) và một panel ghim nhiều bản ghi — cả hai đều phải **cùng một bào quan**. **Ngoại lệ duy nhất, có tên: phân tích chuyển đoạn giữa các ngăn — nhân ↔ lục lạp ↔ ty thể (NUPT/MTPT, M6.1).** Nó đọc cả ba ngăn trong cùng một phân tích một cách có chủ đích (muốn thấy một đoạn lục lạp nằm trong ty thể thì phải đọc cả hai cùng lúc), và là nội dung riêng của Giai đoạn 6 — không phải một bộ dữ liệu trộn. Mô-đun chia nhóm — bảng này là nguồn duy nhất; thêm mô-đun mới thì khai vào đúng nhóm:

| Nhóm | Mô-đun | Vì sao ở nhóm đó (đo được) |
|---|---|---|
| **Dùng chung** | `orf` · `gene_names` · `kaks` · `popgen` · `submission` · `stats` | theo gen/chuỗi, không theo vị trí. `kaks` + `popgen` **đã chạy trên mito thật** (`standin-mito`, §8.30); `orf` chấm chính bản chú giải mito (§8.22); `gene_names` là **tầng tên** dùng chung (hàm `aliases_for` chọn từ vựng theo `/organelle`), nhưng **hôm nay module `gene_names` chỉ nạp bảng plastid** — đo 05/10: 0 renamed / 0 aliased trên bản ghi mito trong khi `mito_genes` gọi tên được 31/51 locus (§9.J) ⇒ **H-12**, tạm hoãn theo D5-07 |
| **Chỉ lục lạp** | `structure` · `genome_stats` · `plastid_annotate` · `gene_table` · `barcodes` · `gene_normalize` | định nghĩa bằng IR/vùng LSC/IRa/SSC/IRb (`structure`, `genome_stats`), bằng công cụ chú giải lục lạp (`plastid_annotate` — engine tự viết, chuyển chú giải từ bản ghi plastome và đọc vùng lặp đảo; ty thể có `mito_annotate`) **[ĐÃ THAY 10/10/2026 — P3B-11: `annotation` ra khỏi bảng cùng chloe · pga; tập 7 mô-đun nhóm plastid nay là 6, xem §8.56]**, hoặc bằng **dữ liệu plastome**: `gene_groups.yaml` → `gene_table` · `barcode_loci.yaml` → `barcodes`. `gene_normalize` ở đây **tạm thời** và có lý do đo được: bản ghi ★ hiện là một plastome, nên chuẩn hoá một bản ghi mito theo nó sẽ in cách viết plastid cho gen ty thể — đúng bài học của H-12 (§8.44, D5-09 §9.8) |
| **Chỉ ty thể** | `mito_annotate` · `mito_genes` | chú giải mitogenome từ bộ tham chiếu + bộ gen ty thể; không cổng IR (§8.21–8.24) |
| **Chuyển tiếp — dùng chung về bản chất, chưa mở cho ty thể** | `ssr` · `introns` · `codon_usage` · `amino_acid` · `substitutions` · `diversity` · `synteny` · `repeats` | SSR · codon usage · amino acid · intron · Ts/Tv · π · mVISTA đều là phân tích chuẩn của bài mitogenome, nhưng mỗi mô-đun báo theo **vùng plastome** hoặc gọi `ctx.require_ir()` ⇒ chạy trên mitogenome thì run **dừng ở cổng, trước khi mô-đun chạy** (đo 03/10: cả hai mitogenome stand-in ném `IRDetectionError`, §8.6). Mở cho ty thể = một chế độ "không vùng" — H-10. `repeats` (M2.12, 10/10/2026) vào hàng này cùng lý do: nó in **vùng plastome** cho từng bản sao match và gọi `ctx.require_ir()`, còn phép tìm repeat thì bào quan nào cũng chạy được |

**[ĐÃ SỬA 09/10/2026 — P3B-09, §8.49]** Nhóm nay được **khai trong registry**: `ModuleSpec.group`
nhận `shared` · `plastid` · `mito`, **bắt buộc** với mọi mô-đun (bỏ sót không phải một lời khai —
thiếu thì import đổ), và cả **22 mô-đun** đã khai theo đúng bảng trên. Luật thi hành nằm ở
`core/organelle.py`: đọc `/organelle` của bản ghi (không đoán khi bản ghi không khai — D5-06: chuỗi
trần tự khai bằng **đường vào**), **hai chiều lệch đều dừng có tên** kèm đường đúng
(`plastid annotate` ↔ `plastid mito annotate`), và một **bộ trộn hai bào quan cũng dừng** — trước
khi mẫu nào được chú giải. Tám mô-đun "chuyển tiếp" khai `shared` (đúng như bảng gọi chúng: *dùng
chung về bản chất*); thứ giữ chúng khỏi mitogenome hôm nay là `requires_ir` — một câu hỏi khác, có
cổng riêng — nên hàng ấy vẫn nguyên hiệu lực và **H-10** vẫn là việc mở chúng ra. Ba mô-đun chưa
có trong bảng (thêm sau 05/10) nay đã vào: `plastid_annotate` → **chỉ lục lạp** (chuyển chú giải
từ bản ghi plastome, đọc vùng lặp đảo), `stats` → **dùng chung** (so sánh nhóm, không theo vị trí),
`reference_audit` → không đăng ký trong registry (nó là bước audit, gọi trực tiếp).

**Hiện trạng trước 09/10/2026, giữ nguyên để đọc lịch sử:** cờ `requires_ir` là tín hiệu duy nhất mô-đun mang, nên hai chiều lệch **không đối xứng**: chạy mô-đun lục lạp trên mitogenome thì **cổng IR đã chặn**, còn chạy `mito_genes` / `mito_annotate` trên plastome thì **không có gì chặn** (không mô-đun nào hỏi bào quan của bản ghi) — và một bộ dữ liệu **trộn hai bào quan** cũng không có gì chặn: không chỗ nào trong run/panel hỏi bào quan của **cả bộ**. Thêm một lỗ của chính luật này, đo 05/10: **artifact do bộ công cụ tự viết chưa khai bào quan** — `plastid mito annotate` xuất GenBank **không có `source` feature và không có `/organelle`**, nên đường chọn từ vựng tên rơi về bảng plastid (H-11). Việc còn lại có tên: **H-09** (khai báo nhóm trong registry + từ chối khi lệch, **kể cả** khi bộ dữ liệu trộn bào quan — **phần khai báo và từ chối đã xong 09/10/2026, §8.49**), **H-10** (mở bảy mô-đun chuyển tiếp), **H-11** (artifact tự viết phải khai `/organelle`) và
**H-12** (`gene_names` chỉ nạp từ vựng plastid — §9.J).

**Ưu tiên (D5-07, cả nhóm chốt 05/10): lục lạp trước.** Dữ liệu ty thể thực vật còn rất ít ⇒ phát triển mô-đun mới tập trung vào **lục lạp**; đường ty thể **giữ nguyên, không mở rộng** (H-10 và H-06 **tạm hoãn**, không phải bỏ — mở lại khi có dữ liệu ty thể). Bảng nhóm ở trên và các điểm mở của D5-06 chính là **cửa để sau này tích hợp phân tích ty thể mà không phải đập lại cấu trúc**.

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
- Cân nhắc dùng engine điều phối có sẵn (Snakemake) thay vì tự viết — **đã chốt 02/10: không dùng ở giai đoạn đầu** (§3); chỉ xem lại khi pipeline vượt ~20 bước hoặc cần chạy lại từng phần

#### Reference Panel — định nghĩa 4 tầng taxonomy (nguyên văn kế hoạch 22/9, mục 0.7 chế độ B)

| Tầng | Là gì | Dùng cho |
|---|---|---|
| Tầng 0 | các accession khác của cùng loài | C.8 di truyền quần thể · phần "trong loài" của C.6 |
| Tầng 1 | cùng chi | nhánh B · cặp Ka/Ks gần (C.2) · lớp mVISTA gần (C.4) |
| Tầng 2 | cùng họ, khác chi | lớp mVISTA xa (C.4) *(ancestral state — C.9 — ra khỏi kế hoạch 10/10/2026, §13.6)* |
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
- Lặp lại phân tán + palindromic (repeat-match · Vmatch) — **[bổ sung 07/10/2026: có ở §1 và §7.1 nhưng thiếu ở đây; milestone M2.12]**
- Đa dạng nucleotide (π) theo gen/intron/vùng liên gen, sliding-window

**Áp dụng cho:** cả **lục lạp** và **ty thể** — nhưng hôm nay phần lớn phân tích dưới đây nằm trong **danh sách chuyển tiếp** (bảng nhóm ở đầu §13): chạy trên mitogenome thì dừng ở cổng IR, việc mở cho ty thể là H-10 (tạm hoãn, D5-07).

**Ra được:** bảng số + hình cho **2 bản thảo Oldenlandia / Coleus** (bản nháp mục Results).

**Kiểm tra trong giai đoạn:** bộ số vàng §6 — độ dài **153.724 bp** · **132 gen = 86 CDS + 38 tRNA + 8 rRNA** (cách đếm ở điều kiện chung) · IR **26.214 bp** (2 bản giống ≥ 99,9 %) · LSC **83.127** / SSC **18.169** · **37 SSR** (mono 31 · di 4 · tri 2) khi có MISA/cpstools (tầng B — §3). Lệch ⇒ dừng, sửa ngay tại đây.

### 13.3. Giai đoạn 3 — Chú giải bộ gen

Chạy được trên bộ gen đã có sẵn.

- Chú giải: PGA (lục lạp), MITOS2/GeSeq (ty thể) — **[ĐÃ THAY — xem §3 (D3-01 · D3A-02 · D5-08): lục lạp = Chloe từ 03/10 (§8.15–8.17), công cụ tự viết theo D5-08 sẽ là giai đoạn mới 3B; ty thể = annotator tự viết từ Giai đoạn 3A]**
- Kiểm tra & chuẩn hoá tên gen — **[ĐÃ THAY 07/10/2026 — D5-09: đồng nhất tên chuyển lên TRƯỚC chú giải, chạy trên DB tham chiếu; sau chú giải chỉ còn `gene_normalize` chế độ kiểm; slot của bước mới còn mở — ROADMAP → Open decisions]**
- Kiểm tra khung đọc (start/stop/internal codon)
- Chuẩn bị hồ sơ nộp NCBI (FASTA/TBL/SQN)

**Kiểm tra trong giai đoạn:** đọc lại file kết quả bằng chính bộ đọc GenBank (không hỏng định dạng) · kiểm tra khung đọc · đếm lại gen theo **cách đếm** đã chốt (CDS + tRNA + rRNA, cả hai bản IR) — **132 là giá trị của genome vàng, không phải hằng số**: P2-10 đo 131 và 133 ở hai loài khác, cả hai đều đúng. *Riêng bộ gen tự lắp ráp* (chỉ có ở giai đoạn 7): phải đối chiếu họ hàng gần hoặc read bắc qua ranh giới LSC–IR / SSC–IR, nếu không thì dán nhãn *chưa xác định* — ở giai đoạn 3 phần này bỏ qua.

### 13.4. Giai đoạn 4 — So sánh nhiều bộ gen

- Áp lực chọn lọc Ka/Ks
- Synteny/collinearity kiểu mVISTA: % tương đồng theo toạ độ
- Ứng viên DNA barcode & vùng biến động cao
- Di truyền quần thể kiểu DnaSP (Hd, θ, Tajima's D, Fu & Li, Fst…)
- Dự đoán RNA editing (cần dữ liệu RNA)

**Áp dụng cho:** cả **lục lạp** và **ty thể** — nhưng hôm nay phần lớn phân tích dưới đây nằm trong **danh sách chuyển tiếp** (bảng nhóm ở đầu §13): chạy trên mitogenome thì dừng ở cổng IR, việc mở cho ty thể là H-10 (tạm hoãn, D5-07).

**Ra được:** nội dung so sánh cho 2 bản thảo.

**Kiểm tra trong giai đoạn:** Ts/Tv = **2,6917** cho cặp *Iris domestica* × *Iris gatesii* (**chỉ CDS · 79 gen dùng chung · MAFFT · không tính IR** — đúng phạm vi ở điều kiện chung; số đo lại 04/10, §8.27) · biểu đồ tương đồng phải ra **7 khối**, %IDY 89,66–99,11 (mean 97,37), coverage 99,09 % — số đo 04/10 trên chính cặp Iris (§8.26). *Con số "21 khối · mean 91,5" trong kế hoạch cũ không tái lập được từ đầu vào của nó và thuộc một so sánh khác chi — ghi thành phát hiện ở §9.H.* Số không khớp ⇒ dừng, kiểm lại dữ liệu vào và tham số.

### 13.5. Giai đoạn 5 — Hình và báo cáo

- Thống nhất theme màu/font; các loại hình: bản đồ plastome vòng · RSCU heatmap · π sliding window (có tô nền LSC/SSC/IR) · sơ đồ IR junction · mVISTA-style %IDY + synteny tuyến tính · phân bố SSR/repeat theo vùng & motif · cây có annotation · NUPT/MTPT theo vị trí gen · mạng haplotype · cấu hình đa hình ty thể · **đánh giá DNA barcode ứng viên**; xuất vector (SVG/PDF/EPS); **quy tắc mở rộng hình theo số bộ gen**
- Thống kê: mô tả, kiểm định phi tham số, tương quan, hiệu chỉnh đa kiểm định (FDR/Bonferroni), diễn giải dạng câu
- Xuất bản: Excel tổng hợp đa sheet · Word (bảng gen/SSR/thống kê) · gói Supplementary Materials · báo cáo PDF 1 trang/bộ gen · báo cáo QC kỹ thuật khách quan
- **Yêu cầu 02/10:** mỗi lần chạy xuất **report công cụ + thông số**, và **output chia theo thư mục** để người dùng truy lại (§2)

**Kiểm tra trong giai đoạn:** sinh lại được bảng/hình từ dữ liệu trung gian (không làm tay) · **kiểm thử truy vấn số**: mỗi số trong report phải trỏ được về **đúng một ô** trong bảng gốc (nếu không, lỗi “số trông hợp lý nhưng sai” vẫn lọt) · report công cụ + thông số có đủ tên, phiên bản, đường dẫn, ngưỡng.

### 13.6. Giai đoạn 6 — Các phân tích cho bài báo công cụ  **[ĐÃ THAY 10/10/2026 — tên gọi, không phải công việc: đây là **năm phân tích** công cụ phải chạy được (NUPT/MTPT · cấu trúc mito · so sánh lục lạp–ty thể · ancestral state · divergence time), tức việc xây công cụ của dev. *Nội dung* bài báo — phân tích gì, khẳng định gì, đưa gì vào bài — là của **cô Trang** (CLAUDE.md → who decides).]**

- NUPT/MTPT — chuyển gen giữa nhân, lục lạp, ty thể (trọng tâm)
- Cấu trúc đa hình ty thể + tái tổ hợp qua trình tự lặp (Chimeric ORF/CMS)
- So sánh song song lục lạp – ty thể trong cùng loài
- **[ĐÃ XOÁ 10/10/2026 — dev quyết: "Ancestral state reconstruction" và "Ước tính thời gian phân
  kỳ" ra khỏi kế hoạch]** — hai cái đó không phải mô-đun của bộ công cụ: chúng cần **cây · hiệu
  chuẩn hoá thạch · một câu hỏi khoa học**, không cần bộ gen, nên là phân tích **cô Trang** chạy cho
  bài báo bằng công cụ của cô. `mcmctree` · `BEAST2` được ghi ở §7 là *đã khảo sát, không đưa vào*.

**Kiểm tra trong giai đoạn:** phần cấu hình ty thể **có mốc công bố để đối chiếu ngay tại đây** — các mitogenome có cấu trúc phức tạp đã công bố: *Damnacanthus indicus* (qua trình tự lặp trực tiếp), *Salvia* (S. officinalis 2023 · S. splendens 2025 · S. miltiorrhiza 2022), *Hemiboea*, *Coffea arabica* (2023) — đều là mốc trong kế hoạch 22/9. Riêng **NUPT/MTPT chưa có số công bố** ⇒ mốc kiểm duy nhất hiện có là đối chiếu với bài Oldenlandia/Coleus, thực hiện ở giai đoạn 8. Giai đoạn này **chạy trên dữ liệu đã công bố / GenBank**; phần cần đồ thị lắp ráp thì chuyển sang giai đoạn 7.

### 13.7. Giai đoạn 7 — Lắp ráp (đã chuyển xuống sau)

Cần dữ liệu đọc thô; nặng về bộ nhớ.

- Lắp ráp + QC hàng loạt (fastp → GetOrganelle)
- Kiểm tra đồ thị assembly (Bandage)
- Cấp dữ liệu cho phần cấu hình ty thể ở giai đoạn 6
- **Đường chạy từ xa** cho máy yếu: usegalaxy.eu (`get_organelle_from_reads`) — cắt bớt read (subsample) trước khi upload, **pin phiên bản công cụ + ghi lại history ID**. *Đưa dữ liệu chưa công bố lên máy chủ công khai: **chưa quyết định** (08/10/2026) — §10.12.*
- **Ca read Illumina thật**: hiện **chưa có** dữ liệu read của đề tài (đang nợ) ⇒ mới chỉ có cổng dựa trên read mô phỏng

**Kiểm tra trong giai đoạn:** so với kết quả đã chạy thật — GetOrganelle + SPAdes cho **127.679 bp** từ 20.000 cặp read mô phỏng · OATK: 126.660 bp, identity **100,00 %** · Bandage ra PNG (103.289 B). Khi có read thật: bổ sung mục đo **đỉnh RAM**.

### 13.8. Giai đoạn 8 — Kiểm thử chéo, hiệu năng và tự động hoá

**Không phải nơi test lần đầu.** Mỗi giai đoạn đã tự kiểm bằng mục *“Kiểm tra trong giai đoạn”*. Giai đoạn này chỉ làm phần không thể làm sớm hơn:

- **Đối chiếu chéo với công cụ khác** trên cùng bộ dữ liệu: CGAS · CPStools · CPGview
  *(CPGview: **khảo sát — không đưa vào**, lý do đo được ở §7.8 · §8.60; nhánh ấy của M8.1 đóng
  bằng phép đo đó, không bằng một lần chạy — P8-01, 10/10/2026)*
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

*Ảnh chụp **26/9/2026** (đợt khảo sát §7.1) đã được gỡ ngày 07/10/2026 — từng dòng của nó nay nằm
ở **§8** (nhật ký sống, từ §8.6 trở đi) hoặc **§7.1** (bản tổng của đợt khảo sát đó), và dòng
"chưa chạy thử lần nào" đã sai từ lâu (Giai đoạn 2–5 chạy hết danh sách ấy). Không còn con số nào
chỉ sống ở đây.*

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
