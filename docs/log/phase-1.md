# Nhật ký đo lường — Giai đoạn 1 — khung & dữ liệu dùng chung

**Phần này là `PROJECT.md` §8** (§8.1–8.5), tách ra ngày 04/10/2026 để `PROJECT.md` đọc được trọn
(§8 khi đó chiếm 1.531/2.379 dòng, 118/216 KB). **Số mục giữ nguyên** như khi còn nằm trong
`PROJECT.md` — mọi tham chiếu kiểu `§8.27` trong code, test, `docs/` và board vẫn trỏ đúng chỗ —
và **không bao giờ đánh số lại**. Mục mới của giai đoạn này ghi tiếp vào **cuối file này**, rồi
thêm một dòng vào bảng mục lục ở `PROJECT.md` §8.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "Nhật ký & tài liệu".

---


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

