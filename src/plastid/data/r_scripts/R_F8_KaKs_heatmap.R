#!/usr/bin/env Rscript
# =============================================================================
# R_F8_KaKs_heatmap.R — Pairwise Ka/Ks heatmap (gene × species), log2 scale
# Figure in the paper: Figure 8
# Usage:
#   Rscript R_F8_KaKs_heatmap.R --input <data file> --outdir <output directory>
# Default: reads sample_data/KaKs_matrix_79genes_BO_D.tsv, writes <directory>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_common.R"))
init_script()
cat("\n=== R_F8_KaKs_heatmap.R ===\n")
# Heatmap Ka/Ks: 33 Iridaceae species vs Eleutherine reference
# Bright vibrant version

library(ggplot2)
library(reshape2)

tsv_path <- input_path("KaKs_matrix_79genes_BO_D.tsv")
raw <- read.delim(tsv_path, header = TRUE, row.names = 1, check.names = FALSE)

mat <- as.matrix(raw)
mode(mat) <- "numeric"

# Filter ≤20% NA
na_per_gene <- rowSums(is.na(mat))
na_per_spp  <- colSums(is.na(mat))
max_na_threshold <- as.numeric(get_param("--max-na", "1"))   # 1 = no filtering (matches the published figure)
keep_genes <- na_per_gene / ncol(mat) <= max_na_threshold
keep_spp   <- na_per_spp / nrow(mat) <= max_na_threshold
mat <- mat[keep_genes, keep_spp]

clean_name <- function(x) gsub("_", " ", x)

df <- melt(mat, varnames = c("Gene", "Species"), value.name = "KaKs")
df$Gene    <- factor(clean_name(df$Gene),    levels = rev(clean_name(rownames(mat))))
df$Species <- factor(clean_name(df$Species), levels = clean_name(colnames(mat)))

# Replace 0 for log scale
df$KaKs[!is.na(df$KaKs) & df$KaKs == 0] <- 0.001

# ── Bright vibrant palette ──
pal_colors <- c("#00B0F0", "#40C4FF", "#80D8FF", "#B3E5FC",
                "#FFFFFF",
                "#FFCDD2", "#EF9A9A", "#F44336", "#D50000")

p <- ggplot(df, aes(x = Species, y = Gene, fill = KaKs)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(
        colors = pal_colors,
        name = "Ka/Ks",
        na.value = "grey90",
        trans = "log2",
        breaks = c(0.01, 0.05, 0.2, 0.5, 1, 2, 5),
        labels = c("0.01", "0.05", "0.2", "0.5", "1", "2", "5+"),
        limits = c(0.005, NA),
        oob = scales::squish
    ) +
    theme_minimal(base_family = "sans") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                    size = 9, face = "bold.italic", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text  = element_text(size = 8),
        legend.position = "right",
        legend.key.height = unit(1.2, "cm"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey95", color = NA),
        plot.margin = margin(10, 10, 10, 10),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 12))
    ) +
    labs(title = get_param("--title", NULL), fill = "Ka/Ks")

save_figure(p, file.path("Figures", get_param("--tag", "Figure8_KaKs_heatmap")),
         width = as.numeric(get_param("--width", "18")),
         height = as.numeric(get_param("--height",  "16")),
         dpi  = as.numeric(get_param("--dpi",  "300")))

cat(sprintf("Saved: Figure8 KaKs heatmap (%d genes x %d species)\n",
    nrow(mat), ncol(mat)))
cat(sprintf("Mean Ka/Ks: %.3f | Median: %.3f\n",
    mean(mat, na.rm = TRUE), median(mat, na.rm = TRUE)))
