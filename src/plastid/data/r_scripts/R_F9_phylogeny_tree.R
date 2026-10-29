#!/usr/bin/env Rscript
# =============================================================================
# R_F9_phylogeny_tree.R — draw the paper's ML tree (ape): circular fan + rectangular
# Figure in the paper: Figure 9 (34 taxa, whole-plastome supermatrix, rooted with Orchidaceae)
#                      and Figures S1–S3 (32/36 taxa, coding-gene supermatrix, midpoint-rooted)
# Usage:
#   Rscript R_F9_phylogeny_tree.R --tree sample_data/tree_34taxa_fullgenome.treefile --outdir outputs
#   Rscript R_F9_phylogeny_tree.R --tree sample_data/tree_32taxa_cds.treefile \
#       --tag FigureS3 --drop-tips "" --root midpoint --width 4200 --height 4200
# Parameters:
#   --tag <prefix>             output file name (default "phylo")
#   --drop-tips "a,b"          taxa to remove from the tree (default: drop no taxa)
#   --outgroup "a,b"           outgroup used to root (default: no outgroup selected)
#   --outgroup-label "(...) "  family annotation appended after outgroup labels ("" = none)
#   --root outgroup|midpoint|none   rooting mode (default: midpoint)
#   --plot-type both|fan|rect       default both
#   --legend bottom|topleft|none     default bottom (separate legend strip below the figure)
#   --width/--height/--dpi         dimensions (default 4200/4200/300)
#   --title ""                 default: no title
#
# PTS does not default to dropping taxa or selecting a paper's outgroup: both choices must be
# passed explicitly (--drop-tips/--outgroup/--root). The midpoint default is a tree that is not
# assigned any paper's outgroup.
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_common.R"))
init_script()
cat("\n=== R_F9_phylogeny_tree.R ===\n")

suppressPackageStartupMessages(library(ape))

TREE_FILE <- input_path("tree_34taxa_fullgenome.treefile")
TAG       <- get_param("--tag", "phylo")
PLOT_TYPE <- get_param("--plot-type", "both")
LEGEND    <- get_param("--legend", "bottom")
WIDTH     <- as.numeric(get_param("--width", "4200"))
HEIGHT    <- as.numeric(get_param("--height", "4200"))
DPI       <- as.numeric(get_param("--dpi", "300"))
TITLE     <- get_param("--title", NULL)
ROOT_MODE <- get_param("--root", "midpoint")

split_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}
DROP_TIPS <- split_arg(get_param("--drop-tips", ""))
OUTGROUP  <- split_arg(get_param("--outgroup", ""))
OG_LABEL  <- get_param("--outgroup-label", "")

# ── 1. Read the tree ─────────────────────────────────────────────────────────
cat(sprintf("  reading tree: %s\n", TREE_FILE))
tree <- read.tree(TREE_FILE)
tree$tip.label <- gsub("_NC_[0-9]+.*", "", tree$tip.label)   # strip accession
tree$tip.label <- gsub("_", " ", tree$tip.label)
cat(sprintf("  %d tips, %d node\n", Ntip(tree), Nnode(tree)))

# ── 2. Drop unused taxa ──────────────────────────────────────────────────────
if (length(DROP_TIPS)) {
  matched <- intersect(DROP_TIPS, tree$tip.label)
  if (length(matched)) { tree <- drop.tip(tree, matched); cat(sprintf("  dropped: %s\n", paste(matched, collapse = ", "))) }
}

# ── 3. Annotate the family for the outgroup ──────────────────────────────────
og_labeled <- character(0)
if (length(OUTGROUP) && nzchar(OG_LABEL)) {
  matched <- intersect(OUTGROUP, tree$tip.label)
  if (length(matched)) {
    tree$tip.label[tree$tip.label %in% matched] <- paste0(tree$tip.label[tree$tip.label %in% matched], "  ", OG_LABEL)
    og_labeled <- paste0(matched, "  ", OG_LABEL)
  }
} else if (length(OUTGROUP)) {
  og_labeled <- intersect(OUTGROUP, tree$tip.label)
}

# ── 4. Root ──────────────────────────────────────────────────────────────────
root_midpoint <- function(tree) {          # root at the midpoint of the longest path
  n <- Ntip(tree); dm <- dist.nodes(tree); d <- dm[1:n, 1:n]; mx <- max(d)
  idx <- which(d == mx, arr.ind = TRUE)[1, ]; path <- nodepath(tree, idx[1], idx[2])
  half <- mx / 2; cum <- 0
  for (k in seq_len(length(path) - 1)) {
    u <- path[k]; v <- path[k + 1]
    e <- which((tree$edge[, 1] == u & tree$edge[, 2] == v) |
               (tree$edge[, 1] == v & tree$edge[, 2] == u))
    el <- tree$edge.length[e]
    if (cum + el >= half) {
      pos <- half - cum
      if (tree$edge[e, 1] == u) { pu <- u; pv <- v; off <- pos } else { pu <- v; pv <- u; off <- el - pos }
      newnode <- Ntip(tree) + Nnode(tree) + 1
      tree$edge[e, 2] <- newnode; tree$edge.length[e] <- off
      tree$edge <- rbind(tree$edge, c(newnode, pv)); tree$edge.length <- c(tree$edge.length, el - off)
      tree$Nnode <- tree$Nnode + 1; tree$node.label <- c(tree$node.label, "")
      return(root(tree, node = newnode, resolve.root = TRUE))
    }
    cum <- cum + el
  }
  tree
}

if (ROOT_MODE == "outgroup" && length(og_labeled) >= 2) {
  tree <- root(tree, og_labeled, resolve.root = TRUE); cat("  rooted by outgroup\n")
} else if (ROOT_MODE == "none") {
  cat("  kept the tree's original root\n")
} else {
  tree <- root_midpoint(tree); cat("  rooted at the midpoint\n")
}
tree <- ladderize(tree)
# Write the tree that was actually drawn: the transformation (tip drops, root, order) is part of
# the result, and a reader must be able to check it instead of inferring it from the image.
USED_TREE <- paste0(TAG, "_used.treefile")
write.tree(tree, USED_TREE)
cat(sprintf("  -> %s (the tree as drawn)\n", USED_TREE))

# ── 5. Parse support labels "SH-aLRT/UFBoot" ─────────────────────────────────
node_labels <- as.character(tree$node.label)
parse_number <- function(x, part) {
  if (is.na(x) || x == "" || x == "Root") return(NA)
  p <- strsplit(x, "/")[[1]]
  if (length(p) >= part && nzchar(p[part])) suppressWarnings(as.numeric(p[part])) else NA
}
shalrt <- vapply(node_labels, parse_number, numeric(1), part = 1)
ufboot <- vapply(node_labels, parse_number, numeric(1), part = 2)

n <- length(tree$tip.label)
tip_cols <- rep("black", n)
tip_cols[tree$tip.label %in% og_labeled] <- "gray35"
fmt_shalrt <- function(v) if (is.na(v)) NA else format(v, nsmall = 1)

# ── 6. Drawing function ──────────────────────────────────────────────────────
draw_tree <- function(plot_type, file_stem, w, h) {
  png(paste0(file_stem, ".png"), width = w, height = h, res = DPI)
  if (LEGEND == "bottom") {
    layout(matrix(1:2, nrow = 2), heights = c(1, 0.12))
    par(mar = if (plot_type == "fan") c(0, 0, 0, 0) else c(1, 1, 1, 10))
  } else {
    par(mar = if (plot_type == "fan") c(0, 0, 0, 0) else c(1, 1, 1, 12))
  }
  if (plot_type == "fan") {
    plot(tree, type = "fan", font = 2, cex = 0.85, edge.width = 1.5, tip.color = tip_cols,
         label.offset = 0.005, no.margin = TRUE, main = TITLE)
    adj_sh <- c(0.5, -0.2); adj_uf <- c(0.5, 1.6)
  } else {
    plot(tree, type = "phylogram", font = 2, cex = 0.9, edge.width = 1.5, tip.color = tip_cols,
         label.offset = 0.001, no.margin = FALSE, align.tip.label = FALSE, main = TITLE)
    adj_sh <- c(0.5, 1.6); adj_uf <- c(0.5, -0.4)
  }
  for (i in seq_along(shalrt)) if (!is.na(shalrt[i]))
    nodelabels(text = fmt_shalrt(shalrt[i]), node = i + n, frame = "none", cex = 0.85,
               font = 2, adj = adj_sh, col = "red")
  for (i in seq_along(ufboot)) if (!is.na(ufboot[i]))
    nodelabels(text = ufboot[i], node = i + n, frame = "none", cex = 0.85,
               font = 2, adj = adj_uf, col = "black")
  add.scale.bar(cex = 0.8, font = 2)
  if (LEGEND == "topleft") {
    legend("topleft", legend = c("SH-aLRT", "UFBoot"), text.col = c("red", "black"),
           text.font = 2, bty = "n", cex = 1.6, title = "Support")
  } else if (LEGEND == "bottom") {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("bottom", legend = c("SH-aLRT", "UFBoot"), text.col = c("red", "black"),
           text.font = 2, bty = "n", horiz = TRUE, cex = 1.4)
    text(0.5, 0.9, "Support", font = 2, cex = 1.4)
  }
  dev.off()
  cat(sprintf("  -> %s.png\n", file_stem))
}

if (PLOT_TYPE %in% c("both", "fan"))
  draw_tree("fan", file.path("Figures", paste0(TAG, "_circular")), WIDTH, WIDTH)
if (PLOT_TYPE %in% c("both", "rect"))
  draw_tree("rect", file.path("Figures", paste0(TAG, "_rectangular")), WIDTH, round(HEIGHT * 0.88))

cat(sprintf("  tips = %d | nodes with support = %d | SH-aLRT red, UFBoot black\n", n, sum(!is.na(shalrt))))
