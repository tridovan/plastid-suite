# =============================================================================
# 00_ham_chung.R — hàm dùng chung cho các script vẽ hình
# =============================================================================
# Không chạy trực tiếp; các script R_F*.R sẽ source file này.
#
# Quy ước dòng lệnh của mọi script trong bộ này:
#   Rscript <script>.R --input <file dữ liệu> --outdir <thư mục kết quả> [tham số riêng]
# Mọi script đều tự tạo <outdir>/Figures và ghi PNG + PDF vào đó.
# =============================================================================

# Đọc một tham số dòng lệnh, trả về giá trị mặc định nếu không có.
tham_so <- function(flag, mac_dinh = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (!is.na(i) && length(a) >= i + 1) return(a[i + 1])
  mac_dinh
}

# Cờ boolean: có mặt thì TRUE.
co_co <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)

# Thư mục chứa script đang chạy — để source được 00_ham_chung.R dù gọi từ đâu.
thu_muc_script <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
}

# Khối đầu mỗi script: lấy --input/--outdir, tạo thư mục, đặt wd = outdir.
#   IN_DIR  : nơi đọc dữ liệu
#   OUT_DIR : nơi ghi kết quả (wd, nên mọi đường dẫn "Figures/..." trong script là tương đối)
khoi_tao <- function(mac_dinh_input = NULL) {
  sc <- thu_muc_script()
  IN_DIR  <<- tham_so("--inputdir", mac_dinh_input %||% file.path(sc, "data_mau"))
  OUT_DIR <<- tham_so("--outdir", file.path(sc, "outputs"))
  IN_FILE <<- tham_so("--input", NULL)
  # Chốt đường dẫn dữ liệu vào thành TUYỆT ĐỐI trước setwd() bên dưới: sau khi đổi thư mục làm
  # việc, một đường dẫn tương đối trỏ vào <outdir> chứ không phải chỗ gọi lệnh — đo 17/10/2026,
  # `--inputdir data_mau` làm cả bảy script không tìm thấy dữ liệu của chính nó.
  IN_DIR  <<- normalizePath(IN_DIR, mustWork = FALSE)
  IN_FILE <<- if (is.null(IN_FILE)) NULL else normalizePath(IN_FILE, mustWork = FALSE)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_DIR, "Figures"), showWarnings = FALSE, recursive = TRUE)
  setwd(OUT_DIR)
  invisible(TRUE)
}

# `a %||% b` — giá trị mặc định khi a NULL.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Đường dẫn file dữ liệu: nếu --input là file thì dùng, không thì lấy tên gốc trong data_mau/
duong_dan_vao <- function(ten_goc) {
  if (!is.null(IN_FILE) && file.exists(IN_FILE)) return(IN_FILE)
  p <- file.path(IN_DIR, ten_goc)
  if (!file.exists(p) && !is.null(IN_FILE)) p <- IN_FILE
  p
}

# Gộp nhiều hình ggplot thành một trang theo chiều dọc (thay cowplot::plot_grid /
# gridExtra::arrangeGrob — hai gói này không có sẵn trong R 4.1.2 của máy).
# Ghi thẳng ra PNG + PDF để không phải gọi lại lệnh vẽ.
ghep_doc <- function(ds_hinh, ten_tep, rong, cao, ti_le_cao = NULL, dpi = 300,
                     dev_pdf = TRUE) {
  n <- length(ds_hinh)
  if (is.null(ti_le_cao)) ti_le_cao <- rep(1, n)
  ti_le_cao <- ti_le_cao / sum(ti_le_cao)

  ve <- function() {
    grid::grid.newpage()
    y <- 1
    for (i in seq_len(n)) {
      h <- ti_le_cao[i]
      grid::pushViewport(grid::viewport(x = 0, y = y - h, width = 1, height = h,
                                        just = c("left", "bottom")))
      print(ds_hinh[[i]], newpage = FALSE)
      grid::popViewport()
      y <- y - h
    }
  }

  png_tep <- paste0(ten_tep, ".png")
  png(png_tep, width = rong, height = cao, units = "in", res = dpi, bg = "white")
  ve(); dev.off()
  if (dev_pdf) {
    pdf_tep <- paste0(ten_tep, ".pdf")
    cairo_pdf(pdf_tep, width = rong, height = cao)
    ve(); dev.off()
  }
  cat(sprintf("  -> %s\n", png_tep))
  if (dev_pdf) cat(sprintf("  -> %s\n", paste0(ten_tep, ".pdf")))
  invisible(png_tep)
}

# Ghi một hình ggplot ra PNG + PDF.
luu_hinh <- function(hinh, ten_tep, rong, cao, dpi = 300) {
  ggsave(paste0(ten_tep, ".png"), plot = hinh, width = rong, height = cao,
         dpi = dpi, bg = "white")
  ggsave(paste0(ten_tep, ".pdf"), plot = hinh, width = rong, height = cao,
         device = cairo_pdf)
  cat(sprintf("  -> %s.png\n  -> %s.pdf\n", ten_tep, ten_tep))
  invisible(hinh)
}

# In tiến trình cho dễ theo dõi khi chạy hàng loạt.
buoc <- function(...) cat(sprintf("  %s\n", paste0(...)))
