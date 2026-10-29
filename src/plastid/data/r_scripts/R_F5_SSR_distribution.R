#!/usr/bin/env Rscript
# =============================================================================
# R_F5_SSR_distribution.R — SSR distribution across 34 plastomes: (A) gene location, (B) LSC/SSC/IR region, (C) motif type
# Figure in the paper: Figure 5
# Usage:
#   Rscript R_F5_SSR_distribution.R --input <data file> --outdir <output directory>
# Default: reads sample_data/ssr_summary_for_r.csv, writes <directory>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_common.R"))
init_script()
cat("\n=== R_F5_SSR_distribution.R ===\n")

# CGAS Module 12: SSR Visualization Script with Formula-Based Dynamic Sizing
# Generated automatically

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

cat("\n========================================\n")
cat("CGAS MODULE 12: SSR VISUALIZATION\n")
cat("========================================\n\n")

# Read data
cat("Reading SSR summary data from CSV...\n")
data <- read.csv(input_path("ssr_summary_for_r.csv"), check.names=FALSE)

if (nrow(data) == 0) {
  stop("ERROR: CSV file is empty - no data to plot")
}

num_species <- nrow(data)
cat(paste("  ✓ Loaded:", nrow(data), "rows\n"))
cat(paste("  ✓ Number of species:", num_species, "\n\n"))

# Create figures directory
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# Dynamic sizing formula - publication quality
# Optimized for 2-30 species range
if (num_species <= 3) {
  plot_width <- 4 + num_species * 1.2
} else if (num_species <= 5) {
  plot_width <- 7.6 + (num_species - 3) * 0.9
} else if (num_species <= 10) {
  plot_width <- 9.4 + (num_species - 5) * 0.6
} else if (num_species <= 20) {
  plot_width <- 12.4 + (num_species - 10) * 0.36
} else {
  plot_width <- 16
}

plot_width <- min(18, max(5, plot_width))
plot_height <- 6

# Bar width - tighter for fewer species
if (num_species <= 5) {
  bar_width_proportion <- 0.7
} else if (num_species <= 10) {
  bar_width_proportion <- 0.75
} else {
  bar_width_proportion <- 0.8
}

# Font sizes
if (num_species <= 5) {
  base_font_size <- 12
  axis_text_size <- 10
  title_size <- 14
  legend_text_size <- 10
} else if (num_species <= 10) {
  base_font_size <- 11
  axis_text_size <- 9
  title_size <- 13
  legend_text_size <- 9
} else if (num_species <= 20) {
  base_font_size <- 9
  axis_text_size <- 7
  title_size <- 11
  legend_text_size <- 8
} else {
  base_font_size <- 8
  axis_text_size <- 6
  title_size <- 10
  legend_text_size <- 7
}

cat(paste("  ✓ Plot dimensions:", round(plot_width, 1), "x", plot_height, "inches\n"))
cat(paste("  ✓ Bar width:", bar_width_proportion, "\n"))
cat(paste("  ✓ Fonts: base=", base_font_size, ", axis=", axis_text_size, "\n\n"))

# ============================================================================
# FIGURE 1: SSR DISTRIBUTION BY GENOMIC REGIONS (LSC, SSC, IR)
# ============================================================================

cat("Generating Figure 1: SSR Distribution by Genomic Regions...\n")

# Prepare data for stacked bar plot
region_data <- data %>%
  select(Species, LSC, SSC, IR) %>%
  pivot_longer(cols = c(LSC, SSC, IR), names_to = "Region", values_to = "Count")

# Order regions
region_data$Region <- factor(region_data$Region, levels = c("LSC", "SSC", "IR"))

# Define colors for regions
region_colors <- c(
  "LSC" = "#4472C4",    # Blue
  "SSC" = "#ED7D31",    # Orange
  "IR" = "#A5A5A5"      # Gray
)

# Create stacked bar plot with formula-based dynamic sizing
p1 <- ggplot(region_data, aes(x = Species, y = Count, fill = Region)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3, width = bar_width_proportion) +
  scale_fill_manual(values = region_colors, name = "Genomic Region") +
  labs(x = "", y = "SSR Count", title = "SSR Distribution by Genomic Regions") +
  theme_minimal(base_size = base_font_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_text_size, face = "italic"),
    axis.text.y = element_text(size = axis_text_size),
    axis.title.y = element_text(size = base_font_size, face = "plain", margin = margin(r = 10)),
    plot.title = element_text(size = title_size, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(size = base_font_size, face = "plain"),
    legend.text = element_text(size = legend_text_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05)))

# Save with formula-based dynamic dimensions
pdf("Figures/SSR_Genomic_Regions.pdf", width = plot_width, height = plot_height)
print(p1)
dev.off()

ggsave("Figures/SSR_Genomic_Regions.png", plot = p1, width = plot_width, height = plot_height, dpi = 600)

cat("  ✓ Figure 1 saved\n\n")

# ============================================================================
# FIGURE 2: SSR DISTRIBUTION BY MOTIF TYPE
# ============================================================================

cat("Generating Figure 2: SSR Distribution by Motif Type...\n")

# Prepare data for stacked bar plot
motif_data <- data %>%
  select(Species, Mono, Di, Tri, Tetra, Penta, Hexa) %>%
  pivot_longer(cols = c(Mono, Di, Tri, Tetra, Penta, Hexa), names_to = "Motif", values_to = "Count")

# Order motifs
motif_data$Motif <- factor(motif_data$Motif, levels = c("Mono", "Di", "Tri", "Tetra", "Penta", "Hexa"))

# Define colors for motifs
motif_colors <- c(
  "Mono" = "#4472C4",    # Blue
  "Di" = "#ED7D31",      # Orange
  "Tri" = "#A5A5A5",     # Gray
  "Tetra" = "#FFC000",   # Yellow
  "Penta" = "#70AD47",   # Green
  "Hexa" = "#5B9BD5"     # Light Blue
)

# Create stacked bar plot with formula-based dynamic sizing
p2 <- ggplot(motif_data, aes(x = Species, y = Count, fill = Motif)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3, width = bar_width_proportion) +
  scale_fill_manual(values = motif_colors, name = "Motif Type") +
  labs(x = "", y = "SSR Count", title = "SSR Distribution by Motif Type") +
  theme_minimal(base_size = base_font_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_text_size, face = "italic"),
    axis.text.y = element_text(size = axis_text_size),
    axis.title.y = element_text(size = base_font_size, face = "plain", margin = margin(r = 10)),
    plot.title = element_text(size = title_size, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(size = base_font_size, face = "plain"),
    legend.text = element_text(size = legend_text_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05)))

# Save with formula-based dynamic dimensions
pdf("Figures/SSR_Motif_Types.pdf", width = plot_width, height = plot_height)
print(p2)
dev.off()

ggsave("Figures/SSR_Motif_Types.png", plot = p2, width = plot_width, height = plot_height, dpi = 600)

cat("  ✓ Figure 2 saved\n\n")

# ============================================================================
# FIGURE 3: SSR DISTRIBUTION BY GENOMIC LOCATION
# ============================================================================

cat("Generating Figure 3: SSR Distribution by Genomic Location...\n")


# Prepare data for stacked bar plot (excluding zero columns)
location_data <- data %>%
  select(Species, Gene, tRNA, Intron, IGS) %>%
  pivot_longer(cols = c("Gene", "tRNA", "Intron", "IGS"), names_to = "Location", values_to = "Count")

# Order locations
location_data$Location <- factor(location_data$Location, levels = c("Gene", "tRNA", "Intron", "IGS"))

# Define colors for locations
location_colors <- c(
  "Gene" = "#4472C4",    # Blue
  "tRNA" = "#ED7D31",    # Orange
  "rRNA" = "#A5A5A5",    # Gray
  "Intron" = "#FFC000",  # Yellow
  "IGS" = "#70AD47"      # Green
)

# Create stacked bar plot with formula-based dynamic sizing
p3 <- ggplot(location_data, aes(x = Species, y = Count, fill = Location)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3, width = bar_width_proportion) +
  scale_fill_manual(values = location_colors, name = "Genomic Location") +
  labs(x = "", y = "SSR Count", title = "SSR Distribution by Genomic Location") +
  theme_minimal(base_size = base_font_size) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_text_size, face = "italic"),
    axis.text.y = element_text(size = axis_text_size),
    axis.title.y = element_text(size = base_font_size, face = "plain", margin = margin(r = 10)),
    plot.title = element_text(size = title_size, face = "plain", hjust = 0.5, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.title = element_text(size = base_font_size, face = "plain"),
    legend.text = element_text(size = legend_text_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05)))

# Save with formula-based dynamic dimensions
pdf("Figures/SSR_Genomic_Locations.pdf", width = plot_width, height = plot_height)
print(p3)
dev.off()

ggsave("Figures/SSR_Genomic_Locations.png", plot = p3, width = plot_width, height = plot_height, dpi = 600)

cat("  ✓ Figure 3 saved\n\n")

cat("\n========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("========================================\n")
cat(paste("  ✓ Generated figures for", num_species, "species\n"))
cat(paste("  ✓ Plot dimensions:", round(plot_width, 1), "x", round(plot_height, 1), "inches\n"))
cat(paste("  ✓ Bar width proportion:", bar_width_proportion, "\n"))
cat(paste("  ✓ Font size:", round(base_font_size, 1), "pt\n"))
cat("\n========================================\n")
