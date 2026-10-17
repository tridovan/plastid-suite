#!/usr/bin/env Rscript
# =============================================================================
# R_F6_nucleotide_diversity_pi.R — Nucleotide diversity (pi): (A) gene mã hoá + RNA, (B) vùng gian gene + intron
# Hình trong bài: Figure 6
# Gọi:
#   Rscript R_F6_nucleotide_diversity_pi.R --input <file dữ liệu> --outdir <thư mục kết quả>
# Mặc định: đọc data_mau/genes_for_R_plot.txt / noncoding_for_R_plot.txt, ghi <thư mục>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_ham_chung.R"))
khoi_tao()
cat("\n=== R_F6_nucleotide_diversity_pi.R ===\n")

# CGAS Module 13: Nucleotide Diversity Visualization
# Generated automatically

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

cat("\n========================================\n")
cat("CGAS MODULE 13: DIVERSITY VISUALIZATION\n")
cat("========================================\n\n")

# Read data
cat("Reading nucleotide diversity data...\n")
df_A <- read.delim(duong_dan_vao("genes_for_R_plot.txt"), sep = "\t", check.names = FALSE)
df_B <- read.delim(duong_dan_vao("noncoding_for_R_plot.txt"), sep = "\t", check.names = FALSE)

cat(paste("  ✓ Loaded genes:", nrow(df_A), "regions\n"))
cat(paste("  ✓ Loaded non-coding:", nrow(df_B), "regions\n\n"))

# Smart y-axis interval function
# Determines appropriate interval based on max value
# Each panel calculates independently since genes and non-coding have different ranges
get_y_axis_breaks <- function(max_val) {
  if (max_val <= 0.04) {
    interval <- 0.005
  } else if (max_val <= 0.08) {
    interval <- 0.01
  } else if (max_val <= 0.12) {
    interval <- 0.02
  } else if (max_val <= 0.3) {
    interval <- 0.05
  } else if (max_val <= 0.6) {
    interval <- 0.1
  } else if (max_val <= 1.2) {
    interval <- 0.2
  } else if (max_val <= 2.0) {
    interval <- 0.4
  } else {
    interval <- 0.5
  }
  max_break <- ceiling(max_val / interval) * interval
  return(list(breaks = seq(0, max_break, by = interval), interval = interval))
}

# Calculate optimal y-axis breaks for both panels
y_axis_A <- get_y_axis_breaks(max(df_A$Value, na.rm = TRUE))
y_axis_B <- get_y_axis_breaks(max(df_B$Value, na.rm = TRUE))

cat(sprintf("Y-axis intervals:\n"))
cat(sprintf("  Panel A: max=%.4f, using %.3f intervals\n", 
            max(df_A$Value, na.rm = TRUE), y_axis_A$interval))
cat(sprintf("  Panel B: max=%.4f, using %.3f intervals\n\n", 
            max(df_B$Value, na.rm = TRUE), y_axis_B$interval))

# Create figures directory
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# ============================================================================
# PANEL A: GENES  (blue line + points, your original fine scaling)
# ============================================================================

cat("Generating Panel A: Genes...\n")

plot_A <- ggplot(df_A, aes(x = factor(Region, levels = Region), y = Value, group = 1)) +
  geom_line(color = "#2c7fb8", size = 0.5) +
  geom_point(size = 1.5, color = "#2c7fb8", fill = "white",
             shape = 21, stroke = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "A: Genes", x = NULL,
       y = expression(paste("Nucleotide diversity (", pi, ")"))) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 9, face = "italic",
      margin = margin(t = 2, b = 20)
    ),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),
    axis.line    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(
      size = 11, face = "bold", hjust = 0,
      margin = margin(b = 10)
    ),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  ) +
  scale_y_continuous(
    breaks = y_axis_A$breaks,
    limits = c(0, max(df_A$Value, na.rm = TRUE) * 1.1),
    labels = scales::label_number(accuracy = y_axis_A$interval),
    expand = expansion(mult = c(0.05, 0.05))
  )

# ============================================================================
# PANEL B: INTERGENIC SPACERS AND INTRONS
# ============================================================================

cat("Generating Panel B: Non-coding regions...\n")

plot_B <- ggplot(df_B, aes(x = factor(Region, levels = Region), y = Value, group = 1)) +
  geom_line(color = "#2c7fb8", size = 0.5) +
  geom_point(size = 1.5, color = "#2c7fb8", fill = "white", shape = 21, stroke = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "B: Intergenic spacers and introns", x = NULL,
       y = expression(paste("Nucleotide diversity (", pi, ")"))) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               size = 9, face = "italic",
                               margin = margin(t = 2, b = 20)),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),
    axis.line    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, 
                              margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 10, b = 50, l = 10)
  ) +
  scale_y_continuous(breaks = y_axis_B$breaks,
                     limits = c(0, max(df_B$Value, na.rm = TRUE) * 1.1),
                     labels = scales::number_format(accuracy = y_axis_B$interval, decimal.mark = "."),
                     expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(clip = "off")

# ============================================================================
# COMBINE PANELS AND SAVE
# ============================================================================

cat("Combining panels and saving...\n")

ghep_doc(list(plot_A, plot_B),
         file.path("Figures", "nucleotide_diversity_plot"),
         rong = as.numeric(tham_so("--rong", "20")),
         cao  = as.numeric(tham_so("--cao",  "10")),
         ti_le_cao = c(1, 1),
         dpi  = as.numeric(tham_so("--dpi", "600")))

cat("\nXong: Figures/nucleotide_diversity_plot.png + .pdf\n")
