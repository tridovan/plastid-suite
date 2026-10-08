# Nhật ký đo lường — Giai đoạn 3B — annotator lục lạp tự viết (D5-08) & bước chuẩn tên (D5-09)

**Phần này là `PROJECT.md` §8** (**§8.40–**). Mục mới của giai đoạn này ghi tiếp vào **cuối file
này**, rồi thêm **một dòng** vào bảng mục lục ở `PROJECT.md` §8. **Số mục không bao giờ đánh số
lại** — mọi tham chiếu kiểu `§8.40` trong code, test, `docs/` và board trỏ theo số.

Mục lục: [`PROJECT.md` §8](../../PROJECT.md) · quy tắc ghi: `CLAUDE.md` → "The run log".

---

### 8.40. Bộ tham chiếu lục lạp: 34 bản ghi Iridaceae được pin lại thành một bộ riêng, và luật leave-one-out viết thành dữ liệu (P3B-01, 08/10/2026)

**Việc.** Dựng bộ tham chiếu cho annotator lục lạp (D5-08, D-b) và pin nó: `.panel/plastid-reference/`
— 34 bản ghi, `PINNED.json` (accession · sha256 · md5 · byte · ngày pin) + `REFERENCE_SET.json`
(bộ + luật + thiết kế leave-one-out). Không phải dữ liệu mới: **chính 34 bản ghi của bộ test
anh Điềm** (`.panel/iridaceae-34/`, §8.31), pin lại thành một bộ riêng vì lý do ở dưới.

**Vì sao là họ Iridaceae.** Luật taxonomy của dự án là **cùng họ trước** (bộ `mito-reference` đã
theo đúng luật này). Genome hợp đồng là *Iris domestica* `NC_050833.1` (họ Iridaceae) và **bộ trả
lời** cũng là Iridaceae — 34 plastome có bản chú giải công bố, thứ duy nhất khiến một kết quả
chuyển chú giải kiểm được. Bộ tham chiếu là **dữ liệu theo họ**: mẫu thuộc họ khác sẽ được pin một
bộ theo họ đó, cùng luật, không phải chọn lại mỗi run.

**Số đo.** 34/34 bản ghi qua cổng IR; tổng **5.204.491 bp**; CDS mỗi bản ghi **84–87** — tức
`132 = 86 + 38 + 8` là **số của genome vàng, không phải hằng số**, và một bộ tham chiếu mà các
thành viên bất đồng về sự tồn tại của một gen chính là thứ mà bước kiểm số 3 sinh ra để in ra.
`PX516262.1` (*Eleutherine bulbosa*) là bản nộp GenBank, không phải RefSeq — vẫn nằm trong bộ vì
bộ test có nó. **34/34 sha256 của bộ tham chiếu trùng khít bộ test** (`test_the_reference_set_holds_the_corpus_records_byte_for_byte`):
pin lại không được đổi byte, vì bộ trả lời thuộc về các file đó.

**Luật leave-one-out, viết thành dữ liệu.** Bộ test chứa chính chuỗi của genome vàng, và bản chú
giải PGA của nó đọc **85 CDS** nơi hợp đồng ghi 86 (§8.31) — để bộ trả lời nằm trong bộ tham chiếu
thì annotator chỉ đang tự chấm bài mình. Nên: mỗi bản ghi trong bộ được chuyển chú giải bởi **chính
bộ đó trừ nó ra** (33 bản ghi), và `leave_one_out` là luật ấy ở dạng đọc được bằng máy
(`core/plastid_reference.py`). Bản cam kết `src/plastid/data/plastid_reference.json` để thiết kế đọc
được **mà không cần 28 MB dữ liệu**, kèm `plastid_reference.md` ghi nguồn, luật và cách dựng lại.
Một đích **ngoài** bộ thì được chuyển chú giải bởi cả bộ — không có gì của nó ở trong để phải loại.

**Một điều phải nói rõ vì nó đi vào bài báo nếu bị bỏ quên.** Đây là bộ test của dự án
(D-b: "bộ test là để kiểm, không bao giờ là production"). Nó làm bộ tham chiếu ở đây vì genome hợp
đồng thuộc Iridaceae và không có nguồn Iridaceae nào khác được pin; **không con số nào đo trên bộ
này được vào bảng bản thảo của mẫu thật.**

**Hiện thực.** `core/plastid_reference.py` (luật + kiểm tra pin) · `plastid reference verify`
(băm lại từng bản ghi, in `file — reason`, thoát 1 khi bộ đã đổi) · `plastid reference loo` (sinh
manifest). Một chi tiết của đường dẫn: `pin_files` ghi lại **đúng đường dẫn nó được đưa**, nên một
pin tạo bằng `--out` tương đối mang đường dẫn tương đối — `verify` tìm theo đường dẫn ghi trong pin
trước, rồi tới chính thư mục pin, nên nó chạy đúng từ bất kỳ thư mục nào (có test riêng).
