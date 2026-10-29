# =============================================================================
# 00_common.R — shared helpers for the figure scripts
# =============================================================================
# Not run directly; the R_F*.R scripts source this file.
#
# Command-line convention for every script in this set:
#   Rscript <script>.R --input <data file> --outdir <output directory> [script-specific parameters]
# Every script creates <outdir>/Figures itself and writes PNG + PDF there.
# =============================================================================

# Read one command-line parameter; return the default value when it is absent.
get_param <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (!is.na(i) && length(a) >= i + 1) return(a[i + 1])
  default
}

# Boolean flag: TRUE when present.
has_flag <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)

# A comma-separated list parameter: character(0) when absent or empty.
split_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

# Directory containing the running script — so 00_common.R can be sourced no matter where it is called from.
script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
}

# Opening block of every script: read --input/--outdir, create directories, set wd = outdir.
#   IN_DIR  : where data is read from
#   OUT_DIR : where results are written (wd, so every "Figures/..." path in a script is relative)
init_script <- function(default_input = NULL) {
  sc <- script_dir()
  IN_DIR  <<- get_param("--inputdir", default_input %||% file.path(sc, "sample_data"))
  OUT_DIR <<- get_param("--outdir", file.path(sc, "outputs"))
  IN_FILE <<- get_param("--input", NULL)
  # Freeze the input data path as ABSOLUTE before setwd() below: after changing the working
  # directory, a relative path points into <outdir> rather than the calling location — measured
  # 2026-10-17, `--inputdir sample_data` made all seven scripts fail to find their own data.
  IN_DIR  <<- normalizePath(IN_DIR, mustWork = FALSE)
  IN_FILE <<- if (is.null(IN_FILE)) NULL else normalizePath(IN_FILE, mustWork = FALSE)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_DIR, "Figures"), showWarnings = FALSE, recursive = TRUE)
  setwd(OUT_DIR)
  invisible(TRUE)
}

# `a %||% b` — fallback value when a is NULL.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Data file path: use --input when it is a file, otherwise take base_name from sample_data/
input_path <- function(base_name) {
  if (!is.null(IN_FILE) && file.exists(IN_FILE)) return(IN_FILE)
  p <- file.path(IN_DIR, base_name)
  if (!file.exists(p) && !is.null(IN_FILE)) p <- IN_FILE
  p
}

# Stack several ggplot figures vertically on one page (replaces cowplot::plot_grid /
# gridExtra::arrangeGrob — neither package is available in this machine's R 4.1.2).
# Write PNG + PDF directly so the drawing command does not have to be called again.
combine_figures <- function(plot_list, file_stem, width, height, height_ratios = NULL, dpi = 300,
                            dev_pdf = TRUE) {
  n <- length(plot_list)
  if (is.null(height_ratios)) height_ratios <- rep(1, n)
  height_ratios <- height_ratios / sum(height_ratios)

  draw <- function() {
    grid::grid.newpage()
    y <- 1
    for (i in seq_len(n)) {
      h <- height_ratios[i]
      grid::pushViewport(grid::viewport(x = 0, y = y - h, width = 1, height = h,
                                        just = c("left", "bottom")))
      print(plot_list[[i]], newpage = FALSE)
      grid::popViewport()
      y <- y - h
    }
  }

  png_file <- paste0(file_stem, ".png")
  png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white")
  draw(); dev.off()
  if (dev_pdf) {
    pdf_file <- paste0(file_stem, ".pdf")
    cairo_pdf(pdf_file, width = width, height = height)
    draw(); dev.off()
  }
  cat(sprintf("  -> %s\n", png_file))
  if (dev_pdf) cat(sprintf("  -> %s\n", paste0(file_stem, ".pdf")))
  invisible(png_file)
}

# Write one ggplot figure to PNG + PDF.
save_figure <- function(plot, file_stem, width, height, dpi = 300) {
  ggsave(paste0(file_stem, ".png"), plot = plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(file_stem, ".pdf"), plot = plot, width = width, height = height,
         device = cairo_pdf)
  cat(sprintf("  -> %s.png\n  -> %s.pdf\n", file_stem, file_stem))
  invisible(plot)
}

# Print progress so batch runs are easy to follow.
log_step <- function(...) cat(sprintf("  %s\n", paste0(...)))
