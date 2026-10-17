#!/usr/bin/env Rscript
# =============================================================================
# R_F7_pairwise_substitution.R — Thay thế nucleotide theo cặp (E. bulbosa vs 33 loài) + tỉ lệ Ts/Tv
# Hình trong bài: Figure 7
# Gọi:
#   Rscript R_F7_pairwise_substitution.R --input <file dữ liệu> --outdir <thư mục kết quả>
# Mặc định: đọc data_mau/Complete_Substitution_Analysis_filtered.csv, ghi <thư mục>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_ham_chung.R"))
khoi_tao()
cat("\n=== R_F7_pairwise_substitution.R ===\n")

# CGAS Module 10: Enhanced SNP/Substitution Visualization Script
# Generated automatically with dynamic plotting capabilities

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

cat("\n========================================\n")
cat("CGAS MODULE 10: ENHANCED VISUALIZATION\n")
cat("========================================\n\n")

# Read data
cat("Reading SNP data from CSV...\n")
data <- read.csv(duong_dan_vao("Complete_Substitution_Analysis_filtered.csv"), check.names=FALSE)

# Safety check: Empty CSV
if (nrow(data) == 0) {
  stop("ERROR: CSV file is empty - no data to plot")
}

cat(paste("  ✓ Loaded:", nrow(data), "rows\n\n"))

# Create figures directory
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# ============================================================================
# DYNAMIC PLOT SIZING
# ============================================================================

# Get sample columns
sample_cols <- names(data)[!(names(data) %in% c("Substitution_Type", "Category"))]
num_samples <- length(sample_cols)

cat(paste("  → Number of samples:", num_samples, "\n"))
cat(paste("  → Samples are alphabetically ordered\n\n"))

# Calculate dynamic dimensions
if (num_samples == 1) {
  plot_width_vert <- 4
  plot_width_horiz <- 8
  plot_height_vert <- 6
  plot_height_horiz <- 3
  base_font <- 14
  axis_font <- 12
  axis_title_font <- 13
  title_font <- 16
  legend_font <- 11
  bar_width <- 0.5
} else if (num_samples == 2) {
  plot_width_vert <- 5
  plot_width_horiz <- 8
  plot_height_vert <- 6
  plot_height_horiz <- 4
  base_font <- 13
  axis_font <- 11
  axis_title_font <- 12
  title_font <- 15
  legend_font <- 10
  bar_width <- 0.6
} else if (num_samples <= 5) {
  plot_width_vert <- 6 + (num_samples - 2) * 0.8
  plot_width_horiz <- 9
  plot_height_vert <- 6
  plot_height_horiz <- 4 + (num_samples - 2) * 0.3
  base_font <- 12
  axis_font <- 10
  axis_title_font <- 11
  title_font <- 14
  legend_font <- 10
  bar_width <- 0.7
} else if (num_samples <= 10) {
  plot_width_vert <- 8.4 + (num_samples - 5) * 0.7
  plot_width_horiz <- 10
  plot_height_vert <- 6.5
  plot_height_horiz <- 6 + (num_samples - 5) * 0.35
  base_font <- 11
  axis_font <- 9
  axis_title_font <- 10
  title_font <- 13
  legend_font <- 9
  bar_width <- 0.7
} else if (num_samples <= 15) {
  plot_width_vert <- 11.9 + (num_samples - 10) * 0.5
  plot_width_horiz <- 11
  plot_height_vert <- 7
  plot_height_horiz <- 7.75 + (num_samples - 10) * 0.3
  base_font <- 10
  axis_font <- 8.5
  axis_title_font <- 9.5
  title_font <- 12
  legend_font <- 8.5
  bar_width <- 0.7
} else if (num_samples <= 20) {
  plot_width_vert <- 14.4 + (num_samples - 15) * 0.4
  plot_width_horiz <- 12
  plot_height_vert <- 7.5
  plot_height_horiz <- 9.25 + (num_samples - 15) * 0.25
  base_font <- 9.5
  axis_font <- 8
  axis_title_font <- 9
  title_font <- 11.5
  legend_font <- 8
  bar_width <- 0.7
} else {
  # Fallback for >20 (should not occur with filtering)
  plot_width_vert <- 16
  plot_width_horiz <- 13
  plot_height_vert <- 8
  plot_height_horiz <- 12
  base_font <- 9
  axis_font <- 7.5
  axis_title_font <- 8.5
  title_font <- 11
  legend_font <- 7.5
  bar_width <- 0.65
}

cat(paste("  → Vertical plot dimensions:", round(plot_width_vert, 1), "×", 
          plot_height_vert, "inches\n"))
cat(paste("  → Horizontal plot dimensions:", round(plot_width_horiz, 1), "×", 
          round(plot_height_horiz, 1), "inches\n\n"))

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Filter substitution data
sub_data <- data[data$Category == "Substitution", ]

if (nrow(sub_data) == 0) {
  stop("ERROR: No substitution data found in CSV")
}

# Define complementary pair mapping
pair_map <- data.frame(
  Original = c("A_to_G", "G_to_A", "C_to_T", "T_to_C", "C_to_G", "G_to_C",
               "A_to_T", "T_to_A", "A_to_C", "C_to_A", "T_to_G", "G_to_T"),
  Combined = c("A↔G", "A↔G", "C↔T", "C↔T", "C↔G", "C↔G",
               "A↔T", "A↔T", "A↔C", "A↔C", "T↔G", "T↔G"),
  Type = c("Transition", "Transition", "Transition", "Transition",
           "Transversion", "Transversion", "Transversion", "Transversion",
           "Transversion", "Transversion", "Transversion", "Transversion"),
  stringsAsFactors = FALSE
)

# Merge and combine pairs
sub_data <- merge(sub_data, pair_map,
                  by.x = "Substitution_Type", by.y = "Original", all.x = TRUE)

# Sum complementary pairs
combined_data <- sub_data %>%
  group_by(Combined, Type) %>%
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

# Convert to long format
plot_data <- combined_data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Count")

# Order substitution types (transitions first)
sub_order <- c("A↔G", "C↔T", "A↔C", "T↔G", "C↔G", "A↔T")
plot_data$Combined <- factor(plot_data$Combined, levels = sub_order)

# Samples should already be alphabetically ordered from CSV
# Preserve that order
plot_data$Sample <- factor(plot_data$Sample, levels = sample_cols)

if (nrow(plot_data) == 0) {
  stop("ERROR: No data to plot after processing")
}

# ============================================================================
# COLOR SCHEME
# ============================================================================

# Publication-quality color scheme for 6 substitution types
sub_colors <- c(
  "A↔G" = "#3498DB",    # Blue (Transition)
  "C↔T" = "#E74C3C",    # Red (Transition)
  "A↔C" = "#95A5A6",    # Gray (Transversion)
  "T↔G" = "#2ECC71",    # Green (Transversion)
  "C↔G" = "#F39C12",    # Orange (Transversion)
  "A↔T" = "#9B59B6"     # Purple (Transversion)
)

# ============================================================================
# PUBLICATION-QUALITY THEME
# ============================================================================

pub_theme <- theme_bw(base_size = base_font) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legend_font),
    legend.text = element_text(face = "italic", size = legend_font),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = title_font, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ============================================================================
# FIGURE 1: VERTICAL BAR PLOT (GROUPED BARS)
# ============================================================================

cat("Generating vertical grouped bar plot...\n")

p_vertical <- ggplot(plot_data, aes(x = Combined, y = Count, fill = Sample)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.85),
    width = bar_width,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = colorRampPalette(c("#3498DB", "#E74C3C", "#2ECC71", "#F39C12", 
                                 "#9B59B6", "#1ABC9C", "#34495E", "#E67E22"))(num_samples),
    name = "Species"
  ) +
  labs(
    x = "Substitution type",
    y = "Number of substitutions",
    title = "Nucleotide Substitution Patterns"
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_font),
    axis.text.y = element_text(size = axis_font),
    axis.title = element_text(size = axis_title_font)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save vertical plot
pdf("Figures/Substitution_Types_vertical.pdf", width = plot_width_vert, height = plot_height_vert)
print(p_vertical)
dev.off()

ggsave("Figures/Substitution_Types_vertical.png", plot = p_vertical, 
       width = plot_width_vert, height = plot_height_vert, dpi = 600)

cat("  ✓ Vertical plot saved\n\n")

# ============================================================================
# FIGURE 2: HORIZONTAL BAR PLOT (GROUPED BARS)
# ============================================================================

cat("Generating horizontal grouped bar plot...\n")

p_horizontal <- ggplot(plot_data, aes(x = Combined, y = Count, fill = Sample)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.85),
    width = bar_width,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = colorRampPalette(c("#3498DB", "#E74C3C", "#2ECC71", "#F39C12", 
                                 "#9B59B6", "#1ABC9C", "#34495E", "#E67E22"))(num_samples),
    name = "Species"
  ) +
  coord_flip() +
  labs(
    x = "Substitution type",
    y = "Number of substitutions",
    title = "Nucleotide Substitution Patterns"
  ) +
  pub_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_font),
    axis.text.y = element_text(size = axis_font),
    axis.title = element_text(size = axis_title_font)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save horizontal plot
pdf("Figures/Substitution_Types_horizontal.pdf", width = plot_width_horiz, height = plot_height_horiz)
print(p_horizontal)
dev.off()

ggsave("Figures/Substitution_Types_horizontal.png", plot = p_horizontal, 
       width = plot_width_horiz, height = plot_height_horiz, dpi = 600)

cat("  ✓ Horizontal plot saved\n\n")

# ============================================================================
# FIGURE 3: TS/TV RATIO
# ============================================================================

cat("Generating Ts/Tv ratio graph...\n")

# Extract Ts/Tv data
tstv_data <- data[data$Substitution_Type == "Ts/Tv Ratio", ]

if (nrow(tstv_data) == 0) {
  stop("ERROR: No Ts/Tv Ratio found in CSV")
}

# Convert to long format
tstv_long <- tstv_data %>%
  select(-Category, -Substitution_Type) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Ratio")

# Convert Inf to numeric (will be handled in plotting)
tstv_long$Ratio <- as.numeric(tstv_long$Ratio)

# Remove NA values
tstv_long <- tstv_long[!is.na(tstv_long$Ratio) & is.finite(tstv_long$Ratio), ]

if (nrow(tstv_long) == 0) {
  cat("  ⚠ Warning: All Ts/Tv ratios are NA or Inf - skipping Ts/Tv plot\n")
} else {
  # Preserve alphabetical order
  tstv_long$Sample <- factor(tstv_long$Sample, levels = sample_cols)
  
  # Create bar plot
  p_tstv <- ggplot(tstv_long, aes(x = Sample, y = Ratio)) +
    geom_bar(stat = "identity", fill = "#2ECC71", color = "black", 
             linewidth = 0.5, width = bar_width) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(
      x = "",
      y = "Ts/Tv Ratio",
      title = "Transition/Transversion Ratio"
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                 size = axis_font, face = "italic"),
      axis.text.y = element_text(size = axis_font),
      axis.title.y = element_text(size = axis_title_font),
      plot.title = element_text(size = title_font),
      legend.position = "none"
    )
  
  # Safe Y-axis scaling
  max_y <- max(tstv_long$Ratio, na.rm = TRUE)
  if (!is.finite(max_y) || max_y <= 0) {
    max_y <- 2.0
  }
  
  p_tstv <- p_tstv + scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = seq(0, ceiling(max_y * 5) / 5, by = 0.2)
  )
  
  # Save Ts/Tv plot
  pdf("Figures/Ts_Tv_Ratio.pdf", width = plot_width_vert, height = plot_height_vert)
  print(p_tstv)
  dev.off()
  
  ggsave("Figures/Ts_Tv_Ratio.png", plot = p_tstv, 
         width = plot_width_vert, height = plot_height_vert, dpi = 600)
  
  cat("  ✓ Ts/Tv ratio plot saved\n\n")
}

cat("========================================\n")
cat("All visualizations completed!\n")
cat("========================================\n")
cat("\nSummary:\n")
cat("  - Complementary pairs combined (e.g., C→G + G→C = C↔G)\n")
cat("  - Total substitution pairs: 6\n")
cat("  - Transitions: A↔G, C↔T\n")
cat("  - Transversions: A↔T, C↔G, A↔C, T↔G\n")
cat("  - Samples ordered alphabetically\n")
cat(paste("  - Plot dimensions adapted for", num_samples, "samples\n"))
