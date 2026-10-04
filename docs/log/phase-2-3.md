# Nhật ký đo lường — Giai đoạn 2 & 3 — đọc/gen so sánh, chú giải, nộp NCBI

**Phần này là `PROJECT.md` §8** (§8.6–8.19), tách ra ngày 04/10/2026 để `PROJECT.md` đọc được trọn
(§8 khi đó chiếm 1.531/2.379 dòng, 118/216 KB). **Số mục giữ nguyên** như khi còn nằm trong
`PROJECT.md` — mọi tham chiếu kiểu `§8.27` trong code, test, `docs/` và board vẫn trỏ đúng chỗ —
và **không bao giờ đánh số lại**. Mục mới của giai đoạn này ghi tiếp vào **cuối file này**, rồi
thêm một dòng vào bảng mục lục ở `PROJECT.md` §8.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "Nhật ký & tài liệu".

---

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

