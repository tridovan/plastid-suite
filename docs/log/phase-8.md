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
