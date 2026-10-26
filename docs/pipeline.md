# Pipeline của PlastidToolSuite

Đây là **định hướng sản phẩm và luồng dữ liệu chuẩn**. Sơ đồ gồm phần đã có và phần dự kiến;
không có nghĩa mọi hộp đã nối thành một command. Bảng đối chiếu bên dưới ghi rõ hiện trạng đọc từ
mã. Lỗi đã tái hiện ở [known issues](known_issues.md). Nhật ký phiên và thứ tự coding không phải
là nguồn định nghĩa pipeline.

## Sơ đồ chuẩn

```text
   NGUỒN THAM CHIẾU
   NCBI RefSeq  ·  bản ghi nhóm tự nộp
                │
                │  plastid panel refresh/propose/pin
                │  accession · nguồn · ngày · sha256
                ▼
   ┌──────────────────────────────────────┐
   │  DB THAM CHIẾU GỐC  (nhiều file .gb) │
   └──────────────────┬───────────────────┘
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐÁNH GIÁ CHẤT LƯỢNG                 │
   │  · file nào KHÁC (REPORT)            │
   │  · file nào không dùng được (BLOCK)  │
   └──────────────────┬───────────────────┘
                      │  báo cáo · sửa hoặc loại file lỗi có ghi nhận
                      ▼
   ┌──────────────────────────────────────┐
   │  NGƯỜI DÙNG ĐỌC VÀ CHỌN              │
   │  1 file .gb phù hợp làm chuẩn tên    │
   └──────────────────┬───────────────────┘
                      ▼
                ★ CHUẨN TÊN    →    từ vựng A
                      │
                      ▼
   ┌──────────────────────────────────────┐
   │  ĐỒNG NHẤT TÊN CHO CẢ DB THAM CHIẾU  │  gene_normalize mode=normalize
   └──────────────────┬───────────────────┘
                      ▼
        DB THAM CHIẾU ĐÃ CHUẨN HOÁ (A)
        artifact mới · checksum mới · truy về DB gốc và ★
                      │
                      │             reads (Illumina · ONT/HiFi)
                      │                         │
                      │                         ▼
                      │       ┌─────────────────────────────────────┐
                      │       │ [1] LẮP RÁP + QC  [HOÃN]             │
                      │       │ short: GetOrganelle + SPAdes         │
                      │       │ long: OATK / raven, cần kiểm chứng    │
                      │       │ fastp: QC short reads                 │
                      │       │ BWA/minimap2: map-back                │
                      │       │ Bandage: đồ thị; kiểm junction/reads  │
                      │       └─────────────────┬───────────────────┘
                      │                         │ FASTA + kết quả QC
                      │                         │
                      │       FASTA đã lắp ráp ─┤  đường vào hiện tại
                      │                         ▼
                      └─────────────►┌───────────────────────────────┐
                                     │ [2] CHÚ GIẢI                   │
                                     │ engine tự viết                │
                                     │ blastx/blastn · ARAGORN · IR  │
                                     └──────────────┬────────────────┘
                                                    │ GenBank — tên theo A
   GenBank đã chú giải bên ngoài ────────────────────┤
   kiểm format/feature/bào quan; giữ nguồn gốc        ▼
                                     ┌───────────────────────────────┐
                                     │ [3] KIỂM TRA TÊN + ANNOTATION │
                                     │ gene_normalize mode=verify    │
                                     │ ORF · vị trí · qualifier      │
                                     └──────────────┬────────────────┘
                                                    ▼
                         ╔══════════════════════════════════════════╗
                         ║ BỘ INPUT ĐÃ KIỂM (từ vựng A)            ║
                         ║ ID mẫu · hash · frame · status từng mẫu  ║
                         ╚══════════════════╤═══════════════════════╝
                                            │
                   ┌────────────────────────┼─────────────────────────┐
                   ▼                        ▼                         ▼
        [4] PHÂN TÍCH             [5] RNA EDITING [CHƯA XÂY] [6b] HỒ SƠ NCBI
        cấu trúc · gen · codon    PREP-Cp gốc · cutoff 0,8   GenBank đã kiểm
        intron · SSR (MISA)       CDS → batch → kết quả gốc  → .fsa + .tbl
        exact repeat             B = reference của PREP-Cp  + template .sbt
          repeat-match           PTS chuẩn bị/import; map   → table2asn
          (Vmatch đối chứng)      CDS → GenBank, tên theo A → .sqn + validator
        synteny (MUMmer4)
        Ka/Ks (PAML yn00)
        Ts/Tv · π · barcode
        popgen · thống kê nhóm
        MAFFT khi cần alignment
                   │                        │ khi M4.5 được triển khai
                   └──────────────┬─────────┘
                                  ▼
                    [6a] BẢNG · HÌNH · BÁO CÁO NHÁP
                    CSV/JSON · Excel · Word · PDF một genome
                    matplotlib · R · cây người dùng cung cấp
                    Supplementary · QC · provenance và checksum
```

**Cách đọc và ranh giới quan trọng:**

- Đường reference cấp DB cho chú giải; không phải bước thay thế assembly. Chuẩn hoá tên trên DB
  chạy **trước** chú giải. Chuẩn ★ là lựa chọn của người dùng, không phải một genome cố định của toolkit.
- Mẫu nhận tên A qua chú giải. Bước verify không phải bước âm thầm sửa tên lần hai: khác biệt là
  phát hiện cần được giải quyết/ghi nhận trước khi pipeline tự nhận input đã chuẩn hoá.
- GenBank ngoài có thể cần một bước normalize có bằng chứng nếu vốn dùng từ vựng khác; không
  đánh đồng khác biệt đó với lỗi của annotator. Giữ bản gốc, mapping tên và hash bản đã chuyển đổi.
- Gate IR chỉ áp dụng cho phần thực sự cần vị trí/vùng IR. Có gate-pass không tự chứng minh biên
  IR đúng; các lỗi hiện hành đang chặn nghiệm thu phần này.
- RNA editing (M4.5) là **việc đang mở, chưa có module**, không phải nhánh bị loại khỏi vòng này.
  **Đã chọn PREP-Cp gốc**, cutoff 0,8: PTS chuẩn bị CDS và import/map kết quả, không viết lại
  predictor. RNA-seq không bắt buộc để chạy dự đoán; dùng bộ công khai để kiểm chứng độc lập.
  Tích hợp đầu tiên là prepare → chạy dịch vụ ngoài → import; chưa xác minh CLI local hoặc API
  để tự động gọi server. Chi tiết và bộ đối chứng đã chọn ở [PREP-Cp](tools/prep_cp.md).
- Hồ sơ NCBI đi từ **sequence và annotation đã kiểm**, không được dựng từ bảng thống kê. `.sbt`
  do người nộp cung cấp. Sinh `.sqn` không đồng nghĩa validator pass hoặc hồ sơ được NCBI chấp nhận.
- Bản đồ vòng đã có renderer matplotlib (`circular_map`); không cần khôi phục `pycirclize` chỉ vì
  nó có trong sơ đồ cũ. Repeat chính dùng `repeat-match`; Vmatch là đối chứng, không phải engine mặc định.
- Đây là luồng mục tiêu. Hiện nhiều bước phải nối thủ công; chưa được quảng bá một lệnh end-to-end.
  `plastid run` chạy nhiều module trên cùng bộ input, không tự truyền GenBank artifact của module
  trước vào module sau. Liệt kê các module trong một lệnh chưa chứng minh đã nối được pipeline.

## Các vai trò reference không được trộn

| Vai trò | Dữ liệu | Lựa chọn/nguồn gốc cần ghi |
|---|---|---|
| DB annotation | Nhiều GenBank đã audit và normalize | Danh sách file, hash gốc/hash chuyển đổi, loại trừ, verdict audit |
| Naming reference ★ | Một GenBank quyết định từ vựng A | Người dùng chọn; file/hash và mapping tên; không tự đổi theo taxonomy mỗi run |
| Reference so sánh | Genome làm hàng đầu/mốc toạ độ cho phân tích | Hiện nhiều module lấy input đầu; nêu rõ record và thứ tự input, không coi mặc nhiên là ★ |
| Assembly seed/HMM | DB seed GetOrganelle hoặc HMM OATK | Phiên bản/checksum DB và chương trình sử dụng; không dùng thay DB annotation |
| Mapping reads | Sequence mẫu đang kiểm | Hash bản sequence và frame tương ứng; không map sang ★ rồi coi là hỗ trợ cho mẫu |
| Reference RNA editing B | Reference/alignment nội bộ của PREP-Cp gốc | PTS không thay DB này; giữ metadata/alignment tool cung cấp, ghi not reported nếu không có revision; pin input/output, map gene ID và vị trí CDS về mẫu A theo [đặc tả](tools/prep_cp.md) |
| Template submission | `.sbt` | Metadata người nộp, không phải genome reference |
| Tree/outgroup | Cây người dùng dựng và taxon được khai | Hash cây, root/drop/order thực hiện; không suy ra outgroup từ tầng taxonomy |

Mode A/B của **panel proposal** là hai cách chọn tập ứng viên, khác với từ vựng A/B ở bước naming.
`panel propose` hiện gợi ý theo taxonomy, không tự chọn ★ hay tự chứng minh reference có chất lượng.

## Hiện trạng và khoảng trống ở từng bước

“Có thành phần” chỉ xác nhận có đường code/kiểm tra, chưa phải đã nghiệm thu toàn pipeline.

| Bước | Thành phần hiện có | Khoảng trống cần nghiệm thu |
|---|---|---|
| 0 — cấp reference | `panel refresh/propose/pin`, `reference verify/loo` | Cache taxonomy/nguồn; dữ liệu bên ngoài; checksum và nguồn gốc qua các bước chuyển đổi |
| 0b — audit/chọn ★ | `reference audit`, `annotate --audit-only`, audit exceptions theo digest | Audit còn chịu lỗi IR; reference mặc định bị chặn; lựa chọn ★ cần đi vào artifact/manifest |
| 0c — normalize DB | `gene_normalize mode=normalize`, blastp cho protein ranking, GenBank artifact và pairing tables | Derived DB phải giữ đủ dữ liệu downstream, có hash mới, không sửa pin gốc; đánh giá lại ngoại lệ audit |
| 1 — assembly/QC | Chỉ có khảo sát tool và khai báo một số dependency | Chưa có module; reads thật, seed/HMM, map-back, graph, junction support, resource budget đều cần tích hợp |
| 2 — annotate | `plastid annotate` gọi `plastid_annotate`, batch mỗi mẫu một run | DB chuẩn hoá chưa nối tự động; feature/locus/IR chênh độ dài và kiểm chứng ngoài họ cần sửa |
| 3 — verify/ORF | `gene_normalize mode=verify`, `orf` | Verify hiện báo would-change, không tự xuất GenBank mới; cần quy tắc dừng/chấp nhận trong workflow và bảo toàn annotation |
| 4 — analyses | Các module liệt kê bên dưới | Cần hợp đồng số liệu từng module và kiểm chứng trước khi tổng hợp; không chỉ nghiệm thu bằng IR/gene count |
| 5 — RNA editing | Chưa có module; đã chọn PREP-Cp gốc, không chọn engine minimap2/samtools mô phỏng | Prepare/import/mapping/provenance; thực chạy engine gốc và đối chiếu Aegilops KJ614412.1 theo [đặc tả](tools/prep_cp.md) |
| 6a — export | CSV/JSON, Word, Excel, `figures`, `package`, `qc`, `genome_page` | Chưa có đầu ra batch thống nhất; `run --figures`/`package`/`qc` hiện chỉ render catalogue matplotlib, không gom R như lệnh `figures` |
| 6b — submission | `submission`: `.fsa/.tbl`, tùy chọn `.sqn` qua table2asn + `.sbt` | Gắn artifact đã kiểm vào đúng run; phân biệt validator chưa chạy với 0 lỗi; kiểm lại loss/partial/qualifier |

Nguồn code để đối chiếu: [CLI](../src/plastid/cli.py),
[panel](../src/plastid/core/panel.py), [audit](../src/plastid/core/reference_audit.py),
[naming](../src/plastid/modules/gene_normalize.py),
[annotator](../src/plastid/modules/plastid_annotate.py),
[Supplementary](../src/plastid/report/package.py), [QC](../src/plastid/report/qc.py),
[submission](../src/plastid/modules/submission.py).

## Phủ các module phân tích và đầu ra

| Module/nhóm | Đầu ra/định nghĩa phải kiểm | Không được suy ra từ tên tính năng |
|---|---|---|
| `genome_stats`, `structure` | Length/GC, số feature và số locus, bốn vùng, frame gốc/chuẩn hoá | Một tổng gene cố định hay một gate-pass không chứng minh mọi vị trí đúng |
| `gene_table`, `gene_names`, `gene_normalize` | Nhóm gene, tên A, synonym/pairing, thông tin mất/đổi khi serialize | Tên giống nhau không tự chứng minh cùng locus; thiếu gene feature không có nghĩa thiếu CDS |
| `orf`, `codon_usage`, `amino_acid` | Transcript strand/order/frame/table, stop/pseudo/partial, RSCU/ENC/CAI/CBI/Fop và reference | Self-reference của bias indices không phải bộ gene biểu hiện cao đã được chọn |
| `introns` | CDS/tRNA/trans-splicing, số lượng và khoảng genomic đúng | Hiện intron tRNA chưa được phủ đầy đủ theo bàn giao |
| `ssr`, `repeats` | Ngưỡng MISA + Perl, motif/context, từng cặp exact-repeat, toạ độ/strand/class | SSR và exact repeats không đồng nghĩa đã có module tandem-repeat tổng quát |
| `substitutions` | Reference/thứ tự, scope CDS/whole, shared loci, xử lý IR/gap/ambiguity, Ts/Tv và loại thay thế | Whole scope mặc định hiện tại khác CDS scope của golden |
| `diversity`, `barcodes` | Alignment columns↔reference↔frame, π/window/step/usable sites, gene/intron/IGS, rule gọi ứng viên | “Ứng viên barcode” không tự chứng minh hiệu năng phân biệt loài |
| `kaks` | Codon alignment, frame/stop, estimator yn00, Ks=0, missing và cách gộp | Không coi tỉ số không tính được là 0, không coi input đầu là naming ★ |
| `synteny` | MUMmer4: nucmer → delta-filter → show-coords, block/identity/coverage và scope tham số | Không đồng nhất mọi phép so sánh genome với một phương pháp LCB khác chỉ vì tool có trong environment |
| `popgen`, `stats` | Nhóm mẫu, n, usable sites, công thức/giả định, đơn vị quan sát, FDR; extraction đúng strand/frame | Cột Fst hiện dùng Nei's Gst; không hứa tương đương mọi estimator. Có số haplotype chưa phải có hình network |
| `submission`, export, figures | Artifact đúng sample/hash/scope, validator, bảng↔hình↔report↔package | Có file không đồng nghĩa đã pass QC; có CSV không chứng minh Excel/Word giữ nguyên số |

Nguồn đầy đủ: [module registry](../src/plastid/modules/__init__.py),
[figure catalogues](../src/plastid/result/figures/__init__.py),
[platform/tool table](platforms.md). `mito_annotate`/`mito_genes` vẫn là prototype ngoài phạm vi
phát triển plastid hiện tại; các kiểm tra hiện có của chúng được giữ.

## Điều kiện nghiệm thu toàn pipeline

1. Có bản input/output của mỗi bước, ID mẫu/hash/frame/phiên bản schema và provenance của chuyển đổi;
   normalized DB là artifact dẫn xuất, không ghi đè raw pin. Biết chính xác file nào được bước sau đọc.
2. Chạy từ cả assembled FASTA và annotated GenBank trên ít nhất hai họ; có ca một/nhiều mẫu và
   mẫu lỗi. Bước thất bại không được tạo trạng thái “đã kiểm” cho downstream; unsupported khác 0.
3. Mỗi chỉ số có đối chứng độc lập cùng định nghĩa, không chỉ cell khớp một bài báo. Đặc biệt kiểm
   sequence extraction, IR copy policy, alignment mapping và missing/unknown.
4. Có đủ nhánh tables/figures/report/QC/Supplementary và submission; R/matplotlib được kê đúng
   phạm vi hỗ trợ. Re-render không chạy lại analysis, không giữ output cũ gây hiểu nhầm.
5. Fresh wheel/checkout có môi trường current-scope dùng được; tool cài trong environment không
   mặc nhiên là chức năng tích hợp. Công cụ/DB của các nhánh hoãn phải được phân loại rõ.
6. M4.5 còn mở phải hiện trong verdict nghiệm thu toàn pipeline. Nếu phát hành một phần trước khi
   RNA editing hoàn thành, phạm vi và phần thiếu phải được quyết định/công bố rõ; không đổi
   “chưa xây” thành “hoãn” chỉ để đóng pipeline. Assembly và các nhánh đã hoãn vẫn có DoD riêng.

Đây là tiêu chí mục tiêu; phiên cập nhật tài liệu chưa triển khai hay chạy nghiệm thu những phần thiếu.
