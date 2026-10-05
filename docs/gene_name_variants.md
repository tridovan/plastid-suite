# Tên gen — bảng tra cứu (mọi ca đã đo)

**File này để làm gì.** Một chỗ duy nhất để tra khi gặp một tên gen lạ: nó là tên khác của
locus nào, hoặc là ca "một tên nhiều locus", hoặc chưa ai đo. **Mọi dòng đều có nguồn đo** —
file này không chứa phỏng đoán. Cập nhật: 05/10/2026.

**Cách dùng nhanh:**

1. Gặp tên lạ ở **lục lạp** → §1 (đồng nghĩa) → không có thì §2 (cú pháp) → §3 (một tên nhiều locus).
2. Gặp **ty thể** → §4.
3. Muốn biết ca nào **còn thiếu/cần tra thêm** → §7.
4. Nhớ một điều: **"không xếp được" không bao giờ bị đoán** — tên không có trong bảng được in
   nguyên như nó đến, kèm lý do.

**Hướng chuẩn hoá (D3-02, §8.10):** tên **canonical của toolkit = phương ngữ của nhóm** (các bản
ghi của nhóm nộp ra từ GeSeq), tên cổ điển in kèm trong ngoặc — `pafI (ycf3)`. CGAS làm **ngược
lại** (kéo về tên cổ điển); xem §6.

**Tra trong toolkit:** `plastid run -m gene_names` in từng tên kèm **rule** của nó —
`identical` · `syntax` · `semantic` · `sequence` — và cặp canonical/synonym. Từ vựng nằm ở
`src/plastid/data/gene_aliases.yaml` (lục lạp, 9 locus) và `mito_gene_names.yaml` (ty thể, 31 locus).

---

## 1. Lục lạp — một locus, nhiều tên (đồng nghĩa ngữ nghĩa)

| Toolkit dùng | Tên khác đã gặp | Đo được ở đâu (số) | Toolkit xử lý | CGAS có? |
|---|---|---|---|---|
| `pafI` | `ycf3` | §8.9–8.10, 9 bản ghi | alias ✓ | có — **ngược chiều** |
| `pafII` | `ycf4` | §8.9–8.10; corpus: chỉ bản công bố có, 2/34 (§8.33) | alias ✓ | có — ngược chiều |
| `pbf1` | `psbN` | §8.9–8.10 + corpus **21/34** | alias ✓ | có — ngược chiều |
| `clpP1` | `clpP` | §8.9–8.10 + corpus **23/34** | alias ✓ | có — ngược chiều |
| `ndhI` | `ndhL` | §8.17 + corpus **8/34** | alias ✓ | **không có** |
| `atpI` | `atpL` | §8.27 + corpus **7/34** | alias ✓ | **không có** |
| `psaI` | `psaL` | §8.27 + corpus **7/34** | alias ✓ | **không có** |
| `psbB` | `bspB` | §8.27 + corpus **7/34** | alias ✓ | **không có** |
| `trnI-GAU` | `trnL-GUA` | **§8.29** · NC_024936.1 (77 bp trùng khít) | alias ✓ | **không có** |
| `trnI-CAU` | `trnL-CAU` | §8.33 · trong 4 bản ghi corpus có khác biệt từ vựng | **⚠ chưa có alias** | không |
| `psbZ` | `lhbA` | §8.31 · corpus **1/34** | **⚠ chưa có alias** | có |
| `cemA` | `ycf10` | **chưa đo ở mình** | chưa có | có |
| `matK` | `maturaseK` | **chưa đo ở mình** | chưa có | có |
| `accD` | `accD1` | **chưa đo ở mình** | chưa có | có |

Ghi chú: `trnL-GUA`/`trnL-CAU` — không có anticodon leucine `GUA` lẫn `CAU` (bộ ba Leu là
UUA/UUG/CUU/CUC/CUA/CUG), nên hai tên đó là cách viết của tRNA isoleucine; bằng chứng mạnh nhất
là trình tự trùng khít từng base.

## 2. Cú pháp — viết khác, không cần kiến thức ngoài

| Dạng gặp | Chuẩn hoá thành | Đo được | Ghi chú |
|---|---|---|---|
| `trnL(UAA)` | `trnL-UAA` | 30 tên, 1 bản ghi (NC_058328.1) | ngoặc tròn → gạch nối |
| `rrn16S` | `rrn16` | — | **ty thể `rrnS`/`rrnL` cố ý KHÔNG đụng**: ở đó chữ cuối là tên, không phải cỡ |
| `rrn4.5S` | `rrn4.5` | — | họ với `rrn16S` |
| `trnk-UUU` → `trnK-UUU`; `trmA-UGC` → `trnA-UGC` | chuẩn hoá chữ axit amin | §8.29, §8.33 | chữ cái axit amin không phân biệt hoa/thường; **`trnfM-CAU` không bị đụng** (chữ `f` là phần của tên) |
| `RRN16`, khoảng trắng | `rrn16` | — | hoa/thường + khoảng trắng chung |
| `16S ribosomal RNA` (cả câu làm tên) | `rrn16` | §8.31 · **1/34** | câu không phải ký hiệu |
| `trnK`, `trnI`, `trnL` (trần, không anticodon) | giải bằng **trình tự** | §8.31 · **8/34** | giải được khi khớp 100 % một tRNA đã biết; không thì báo *ambiguous* kèm tên các ứng viên. Đọc anticodon bằng offset cố định chỉ khớp **106/293** tRNA (11/38 genome vàng) ⇒ **đã loại** |

## 3. Một tên — nhiều locus (trùng nghĩa hẹp — ca nguy hiểm nhất)

| Ca | Chi tiết | Trạng thái |
|---|---|---|
| `trnG-UCC` | Genome vàng NC_050833.1 dùng **một tên cho hai tRNA glycine khác trình tự**; bộ chuẩn có **cả** `trnG-UCC` và `trnG-GCC`. Trong corpus (§8.33): bản công bố viết `trnG-UCC` nơi **PGA và Chloe cùng** viết `trnG-GCC` | **H-04** — cô Trang quyết tên in trong bản thảo (§8.17, §10.10) |
| `"ribosomal protein S3"` | **Cùng một câu `/product` ở cả lục lạp lẫn ty thể** — hai gen khác nhau | Lý do có hai từ vựng riêng; chỉ đường mito nạp `mito_gene_names.yaml` |
| `trnH` | Bản ghi ghi `trnH-GUG`, ARAGORN đọc `trnH-AUG` (ty thể) | §8.23 — lệch anticodon, 1 ca, đã báo chứ không tự sửa |
| `"5S"` / `"18S"` | *Psychotria* gọi feature **722 bp** là "5S", **1112 bp** là "18S" — sai cỡ thật (5S 107–122 · 18S 1692–1935) | Builder loại + in ra (§8.21) |

## 4. Ty thể — tên là **câu `/product`**, chưa phải ký hiệu

Bản ghi mito thường **không có `/gene`** (đo trên OX459128.1: 0 feature có `/gene`), nên "tên" là
cả câu product ("NADH dehydrogenase subunit 1"). `mito_gene_names.yaml` map câu → ký hiệu chuẩn.

- **31 locus đã map.** File tự khai: phần *câu* là đo được (trích nguyên văn từ OX459128.1), phần
  *ký hiệu* là danh pháp chuẩn **chưa đối chiếu một chú giải được kiểm duyệt** — việc của Phase 8.
  Đây là **bằng chứng yếu nhất trong repo** và file nói rõ điều đó.
- **10 gen chuẩn chưa có tên** (H-06, tạm hoãn theo D5-07): `rpl2` · `rpl10` · `rpl16` · `rpl23` ·
  `rps2` · `rps7` · `rps12` · `rps14` · `rps19` · `sdh3`.
- Hit **chưa có tên** đo trên bản đáp án (§8.22): `rps19` · `rps14` · `rpl16` · `rps2` · `ccmFN2` ·
  `rbcL` · `atpB` · `orfX` · `orfB`.
- **Gen nguồn gốc lục lạp nằm trong ty thể** (đừng nhầm là gen mito — đây chính là nội dung
  MTPT): `rbcL` · `petG` · `petL` · `atpB` · `atpE` (*Salvia*) · `citB` · `epd` · `ltrA` · `rnhA` ·
  `tuf` (*Psychotria*) — §8.21.
- Tên vô danh trong bộ tham chiếu mito (bị loại có chủ đích): **236 ORF** (`orf101a`, `orf297a`…) +
  7 feature tên là cả câu (`hypothetical protein`, `LAGLIDADG endonuclease`) — §8.21.

## 5. Không phải chuyện tên — đừng nhầm khi tra

Các ca này nằm cùng bảng thống kê với các ca tên nhưng **khác kết luận** (gen có/không, số bản sao):

| Ca | Số bản ghi | Nguồn |
|---|---|---|
| `ycf1`: bản công bố có CDS, **cả PGA và Chloe đều không gọi** | **24/34** | §8.31, §8.33 |
| `rps12` (trans-spliced): Chloe gọi **2 bản sao**, hai bên kia gọi 1 | **22/34** | §8.33 |
| `rps12` mất một bản sao (PGA so bản công bố) | 12/34 | §8.31 |
| `pafII` (`ycf4`): chỉ bản công bố có | 2/34 | §8.33 |
| `rps19`: chỉ bản công bố có | 1/34 | §8.33 |
| `rps16`: chỉ bản công bố có (bản ghi *Libertia*, §9.I.2) | 1/34 | §8.33 |

## 6. CGAS — bảng đối chiếu (đọc mã nguồn 05/10/2026)

Repo `github.com/abdullah30/Chloroplast-Genome-Analysis-Suite-CGAS`, commit `52b82e60…`
(26/05/2026), MIT. Bản đọc: `.data/cgas-review/` (gitignored).

- **`GENE_SYNONYMS`** (`cgas/cgas_module5.py:87-121`): 19 mục, **7 mục vô nghĩa**
  (`infA→infA`, `matK→matK`, `accD→accD`, `psbA→psbA`, `rbcL→rbcL`, `atpA→atpA`, `atpB→atpB`);
  12 mục thật: `pafI/paf1→ycf3` · `pafII/paf2→ycf4` · `ycf10→cemA` · `lhbA→psbZ` ·
  `pbf1/psb1→psbN` · `clpP1/clp1→clpP` · `maturaseK→matK` · `accD1→accD`.
  **Hướng: về tên cổ điển — ngược với D3-02 của mình.**
- **Ba engine đổi tên không thống nhất**: module 3 (khớp theo reference + *ghi vào file*;
  tRNA thiếu anticodon → **"return first as best guess"**, `cgas_module3.py:208`), module 5
  (synonym cứng + *chỉ báo cáo*; tRNA thiếu anticodon → **giữ nguyên + cảnh báo**),
  module 13 (dict chuẩn hoá hoa/thường **thứ ba**, viết lại từ đầu).
- **Module 3 khớp tiền tố hai chiều** (`:150-153`): tên lạ có tiền tố tên chuẩn bị gán bừa
  (`ycf1-like` → `ycf1` theo code) — không report.
- **Docstring mâu thuẫn code**: docstring module 5 dạy `trnfM-CAU → trnM-CAU` ("remove 'f'"),
  code thì cấm đúng ca đó (`:121-124`).
- **Module 4** (`:89-104`): đọc `/organelle` của feature source để ghi header nộp NCBI
  `[location=…]`; **thiếu qualifier thì mặc định `chloroplast`** — cùng lỗ hổng H-11 của mình.
- **CGAS thiếu đúng những ca corpus của nhóm cần**: không có `ndhL→ndhI`, không có bộ ba
  `atpL/psaL/bspB`, không có `trnI-GAU/trnL-GUA`.

## 7. Còn thiếu / cần tra cứu thêm

| # | Việc | Bằng chứng đã có? | Ghi chú |
|---|---|---|---|
| 1 | Thêm alias `trnL-CAU` → `trnI-CAU` | **có** (§8.33, không có anticodon Leu CAU) | họ với `trnL-GUA`; cần test canh như các alias khác |
| 2 | Thêm alias `lhbA` → `psbZ` | đo 1/34 (§8.31) | CGAS có; mình chưa |
| 3 | `trnG-UCC` vs `trnG-GCC` in tên nào | có (§8.17, §8.33) | **H-04 — cô Trang quyết** (§10.10) |
| 4 | Xác minh `ycf10`→`cemA`, `maturaseK`→`matK`, `accD1`→`accD` | **chưa đo ở mình** — chỉ CGAS khẳng định | đối chiếu NCBI/literature trước khi thêm |
| 5 | 10 gen mito chưa có tên | có (H-06) | tạm hoãn theo D5-07 (lục lạp trước) |

## 8. Nguồn & một lỗi trích dẫn đã biết

- Từ vựng: `src/plastid/data/gene_aliases.yaml` · `mito_gene_names.yaml`; quy tắc:
  `src/plastid/core/gene_names.py` (rule `syntax`/`semantic`/`sequence`, và số đo 106/293).
- Nhật ký: §8.9 · §8.10 · §8.17 · §8.21 · §8.22 · §8.23 · §8.27 · §8.29 · §8.31 · §8.33;
  câu hỏi mở: §10.10 (trnG), H-04/H-06 trong `TASKS.md`.
- **Lỗi trích dẫn đã biết:** alias `trnL-GUA → trnI-GAU` được đo ở **§8.29** (entry "Barcode +
  vùng siêu biến…"), nhưng `gene_aliases.yaml` (mục `trnI-GAU`) và §8.33 đều ghi **§8.28**
  (entry popgen — không chứa ca này). Sửa ở `gene_aliases.yaml`; §8.33 là nhật ký (append-only)
  nên giữ nguyên, ghi chú tại đây.

*Bảng này KHÔNG thay `gene_aliases.yaml` — file YAML mới là nguồn thi hành (loader từ chối entry
thiếu bằng chứng). Đây là bản để người đọc tra và quyết định thêm gì.*
