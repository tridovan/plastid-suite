#!/usr/bin/env Rscript
# =============================================================================
# R_F9_phylogeny_tree.R — vẽ cây ML của bài (ape): dạng toả tròn (fan) + hình chữ nhật
# Hình trong bài: Figure 9 (34 taxa, supermatrix toàn plastome, root bằng Orchidaceae)
#                 và Figure S1–S3 (32/36 taxa, supermatrix gen mã hoá, root giữa)
# Gọi:
#   Rscript R_F9_phylogeny_tree.R --tree data_mau/cay_34taxa_fullgenome.treefile --outdir outputs
#   Rscript R_F9_phylogeny_tree.R --tree data_mau/cay_32taxa_cds.treefile \
#       --tag FigureS3 --drop-tips "" --root midpoint --rong 4200 --cao 4200
# Tham số:
#   --tag <tiền tố>            tên file kết quả (mặc định "phylo")
#   --drop-tips "a,b"          loài bỏ khỏi cây (mặc định 2 loài Asparagus)
#   --outgroup "a,b"           nhóm ngoài dùng để root (mặc định 2 loài Orchidaceae)
#   --og-nhan "(Orchidaceae)"  chú thích họ gắn sau nhãn nhóm ngoài ("" = không gắn)
#   --root outgroup|midpoint|none   cách root (mặc định: outgroup)
#   --kieu both|fan|rect       mặc định both
#   --legend duoi|topleft|khong     mặc định duoi (dải chú thích riêng dưới hình)
#   --rong/--cao/--dpi         kích thước (mặc định 4200/4200/300)
#   --tieude ""                mặc định không có tiêu đề
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_ham_chung.R"))
khoi_tao()
cat("\n=== R_F9_phylogeny_tree.R ===\n")

suppressPackageStartupMessages(library(ape))

TEP_CAY  <- duong_dan_vao("cay_34taxa_fullgenome.treefile")
TAG      <- tham_so("--tag", "phylo")
KIEU     <- tham_so("--kieu", "both")
LEGEND   <- tham_so("--legend", "duoi")
RONG     <- as.numeric(tham_so("--rong", "4200"))
CAO      <- as.numeric(tham_so("--cao", "4200"))
DPI      <- as.numeric(tham_so("--dpi", "300"))
TIEU_DE  <- tham_so("--tieude", NULL)
CACH_ROOT <- tham_so("--root", "outgroup")

tach <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}
BO_LOAI  <- tach(tham_so("--drop-tips", "Asparagus officinalis,Asparagus schoberioides"))
NHOM_NGOAI <- tach(tham_so("--outgroup", "Dendrobium officinale,Phalaenopsis equestris"))
OG_NHAN  <- tham_so("--og-nhan", "(Orchidaceae)")

# ── 1. Đọc cây ───────────────────────────────────────────────────────────────
cat(sprintf("  đọc cây: %s\n", TEP_CAY))
tree <- read.tree(TEP_CAY)
tree$tip.label <- gsub("_NC_[0-9]+.*", "", tree$tip.label)   # bỏ accession
tree$tip.label <- gsub("_", " ", tree$tip.label)
cat(sprintf("  %d tips, %d node\n", Ntip(tree), Nnode(tree)))

# ── 2. Bỏ loài không dùng ────────────────────────────────────────────────────
if (length(BO_LOAI)) {
  co <- intersect(BO_LOAI, tree$tip.label)
  if (length(co)) { tree <- drop.tip(tree, co); cat(sprintf("  bỏ: %s\n", paste(co, collapse = ", "))) }
}

# ── 3. Chú thích họ cho nhóm ngoài ───────────────────────────────────────────
og_labeled <- character(0)
if (length(NHOM_NGOAI) && nzchar(OG_NHAN)) {
  co <- intersect(NHOM_NGOAI, tree$tip.label)
  if (length(co)) {
    tree$tip.label[tree$tip.label %in% co] <- paste0(tree$tip.label[tree$tip.label %in% co], "  ", OG_NHAN)
    og_labeled <- paste0(co, "  ", OG_NHAN)
  }
} else if (length(NHOM_NGOAI)) {
  og_labeled <- intersect(NHOM_NGOAI, tree$tip.label)
}

# ── 4. Root ──────────────────────────────────────────────────────────────────
root_giua <- function(tree) {          # root tại trung điểm đường đi dài nhất
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

if (CACH_ROOT == "outgroup" && length(og_labeled) >= 2) {
  tree <- root(tree, og_labeled, resolve.root = TRUE); cat("  root theo nhóm ngoài\n")
} else if (CACH_ROOT == "none") {
  cat("  giữ nguyên gốc cây\n")
} else {
  tree <- root_giua(tree); cat("  root tại trung điểm\n")
}
tree <- ladderize(tree)

# ── 5. Tách nhãn hỗ trợ "SH-aLRT/UFBoot" ─────────────────────────────────────
node_labels <- as.character(tree$node.label)
lay_so <- function(x, phan) {
  if (is.na(x) || x == "" || x == "Root") return(NA)
  p <- strsplit(x, "/")[[1]]
  if (length(p) >= phan && nzchar(p[phan])) suppressWarnings(as.numeric(p[phan])) else NA
}
shalrt <- vapply(node_labels, lay_so, numeric(1), phan = 1)
ufboot <- vapply(node_labels, lay_so, numeric(1), phan = 2)

n <- length(tree$tip.label)
tip_cols <- rep("black", n)
tip_cols[tree$tip.label %in% og_labeled] <- "gray35"
fmt_shalrt <- function(v) if (is.na(v)) NA else format(v, nsmall = 1)

# ── 6. Hàm vẽ ────────────────────────────────────────────────────────────────
ve <- function(loai_hinh, ten_tep, w, h) {
  png(paste0(ten_tep, ".png"), width = w, height = h, res = DPI)
  if (LEGEND == "duoi") {
    layout(matrix(1:2, nrow = 2), heights = c(1, 0.12))
    par(mar = if (loai_hinh == "fan") c(0, 0, 0, 0) else c(1, 1, 1, 10))
  } else {
    par(mar = if (loai_hinh == "fan") c(0, 0, 0, 0) else c(1, 1, 1, 12))
  }
  if (loai_hinh == "fan") {
    plot(tree, type = "fan", font = 2, cex = 0.85, edge.width = 1.5, tip.color = tip_cols,
         label.offset = 0.005, no.margin = TRUE, main = TIEU_DE)
    adj_sh <- c(0.5, -0.2); adj_uf <- c(0.5, 1.6)
  } else {
    plot(tree, type = "phylogram", font = 2, cex = 0.9, edge.width = 1.5, tip.color = tip_cols,
         label.offset = 0.001, no.margin = FALSE, align.tip.label = FALSE, main = TIEU_DE)
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
  } else if (LEGEND == "duoi") {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("bottom", legend = c("SH-aLRT", "UFBoot"), text.col = c("red", "black"),
           text.font = 2, bty = "n", horiz = TRUE, cex = 1.4)
    text(0.5, 0.9, "Support", font = 2, cex = 1.4)
  }
  dev.off()
  cat(sprintf("  -> %s.png\n", ten_tep))
}

if (KIEU %in% c("both", "fan"))
  ve("fan", file.path("Figures", paste0(TAG, "_circular")), RONG, RONG)
if (KIEU %in% c("both", "rect"))
  ve("rect", file.path("Figures", paste0(TAG, "_rectangular")), RONG, round(CAO * 0.88))

cat(sprintf("  tips = %d | node có hỗ trợ = %d | SH-aLRT đỏ, UFBoot đen\n", n, sum(!is.na(shalrt))))
