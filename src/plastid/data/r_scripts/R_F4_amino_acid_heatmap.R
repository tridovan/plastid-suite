#!/usr/bin/env Rscript
# =============================================================================
# R_F4_amino_acid_heatmap.R — Amino acid composition (%), one column per species in the input
# Figure in the paper: Figure 4 (34 Iridaceae species; any input draws its own species)
# Usage:
#   Rscript R_F4_amino_acid_heatmap.R --input <data file> --outdir <output directory>
# Default: reads sample_data/Complete_Amino_Acid_Analysis.csv, writes <directory>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_common.R"))
init_script()
cat("\n=== R_F4_amino_acid_heatmap.R ===\n")
# CGAS Module 9: Amino Acid Comparison — ALL species in EXACT CGAS style
# Dodged bars, x=amino acid (3-letter, sorted by avg frequency), fill=species
# NO percentage labels (unreadable with many species)
# Draws every species column in the input, in the input order unless --species-order is given.
#
# USAGE:
#   cd Module9_Amino_Acid_Analysis/
#   Rscript /path/to/aa_34sp_cgas_style.R
#
# PREREQUISITES: Module 9 must be run first to generate Complete_Amino_Acid_Analysis.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(RColorBrewer)
})

cat("\n========================================\n")
cat("CGAS MODULE 9: AMINO ACID COMPARISON — ALL SPECIES\n")
cat("========================================\n\n")

# === SPECIES COME FROM THE INPUT ===
# One row per residue, then one column per sample; every column after the two metadata columns is
# a species, in the input's own order. `--species-order` fixes another order and
# `--exclude-species` drops named columns. No species list is built in.
META_COLS <- c("AA_Code", "Amino_Acid")
SPECIES_ORDER  <- split_list(get_param("--species-order", ""))
EXCLUDE_SPECIES <- split_list(get_param("--exclude-species", ""))

data <- read.csv(input_path("Complete_Amino_Acid_Analysis.csv"), check.names=FALSE, stringsAsFactors=FALSE)
available_sp <- setdiff(colnames(data), META_COLS)
if (length(EXCLUDE_SPECIES)) available_sp <- setdiff(available_sp, EXCLUDE_SPECIES)
if (length(SPECIES_ORDER)) {
  missing_sp <- setdiff(SPECIES_ORDER, available_sp)
  if (length(missing_sp)) {
    stop(sprintf("--species-order names %d species not in the input: %s",
                 length(missing_sp), paste(missing_sp, collapse = ", ")))
  }
  available_sp <- SPECIES_ORDER
}
if (!length(available_sp)) stop("the input has no species column after --exclude-species")
cat(paste("Species:", length(available_sp), "\n"))

# Clean names
colnames(data)[3:ncol(data)] <- gsub("_"," ",colnames(data)[3:ncol(data)])
sp_display <- gsub("_"," ",available_sp)

plot_data <- data[,c("AA_Code","Amino_Acid",sp_display)]

aa_data <- plot_data %>%
  pivot_longer(cols=-c(AA_Code,Amino_Acid),names_to="Species",values_to="Percentage") %>%
  mutate(Species=factor(Species,levels=sp_display),AA_Code=factor(AA_Code))
colnames(aa_data)[colnames(aa_data)=="AA_Code"] <- "AminoAcid"

# 1-letter to 3-letter (exactly like CGAS)
aa_1to3 <- c("A"="Ala","R"="Arg","N"="Asn","D"="Asp","C"="Cys","Q"="Gln","E"="Glu",
  "G"="Gly","H"="His","I"="Ile","L"="Leu","K"="Lys","M"="Met","F"="Phe","P"="Pro",
  "S"="Ser","T"="Thr","W"="Trp","Y"="Tyr","V"="Val")
aa_data$AminoAcid <- aa_1to3[as.character(aa_data$AminoAcid)]

# Reorder by average percentage (exactly like CGAS)
avg_order <- aa_data %>% group_by(AminoAcid) %>% summarise(avg=mean(Percentage)) %>%
  arrange(desc(avg)) %>% pull(AminoAcid)
aa_data$AminoAcid <- factor(aa_data$AminoAcid,levels=avg_order)

n_sp <- length(available_sp)
# Scale to the input: the fixed brewer stacks ran out past ~41 colours and returned NA.
species_colors <- grDevices::hcl.colors(n_sp, "Dark 3")

# CGAS-EXACT plot: dodged bars, NO percentage labels (unreadable with 34 species)
p <- ggplot(aa_data,aes(x=AminoAcid,y=Percentage,fill=Species))+
  geom_bar(stat="identity",position=position_dodge(width=0.9),width=0.8,
           color="black",linewidth=0.2)+
  scale_fill_manual(values=species_colors,name="")+
  scale_y_continuous(expand=expansion(mult=c(0,0.12)),limits=c(0,14),breaks=seq(0,14,2))+
  labs(x="Amino acid",y="Amino acid frequency")+
  theme_classic(base_size=12)+
  theme(axis.line=element_line(color="black",size=0.8),
        axis.ticks=element_line(color="black",size=0.8),
        axis.title.x=element_text(size=14,face="bold",margin=margin(t=15)),
        axis.title.y=element_text(size=14,face="bold",margin=margin(r=10)),
        axis.text.y=element_text(size=10,color="black"),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,size=18,color="black"),
        legend.position="bottom",legend.direction="horizontal",
        legend.text=element_text(size=6,face="italic"),legend.title=element_blank(),
        legend.key.size=unit(0.4,"cm"),legend.background=element_rect(fill="white",color=NA),
        panel.grid=element_blank(),
        plot.margin=margin(t=0.5,r=0.5,b=2,l=0.5,unit="cm"),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.8))+
  guides(fill=guide_legend(nrow=6))

dir.create("Figures",showWarnings=FALSE)
# Width follows the number of species columns for the same reason as Figure 3.
FIG_WIDTH  <- as.numeric(get_param("--width", max(22, 11 + 0.35 * n_sp)))
FIG_HEIGHT <- as.numeric(get_param("--height", 12))
# A large panel is a legitimate wide figure, but a raster device has a pixel ceiling: cap the
# resolution so the long side stays at 20,000 px.
FIG_DPI <- min(600, floor(20000 / max(FIG_WIDTH, FIG_HEIGHT)))
ggsave("Figures/Amino_Acid_Comparison.png",p,width=FIG_WIDTH,height=FIG_HEIGHT,dpi=FIG_DPI,
       limitsize=FALSE)
pdf("Figures/Amino_Acid_Comparison.pdf",width=FIG_WIDTH,height=FIG_HEIGHT); print(p); dev.off()
cat(sprintf("\nSaved: Figures/Amino_Acid_Comparison.pdf and .png (%d species, %.0f x %.0f in)\n",
            n_sp, FIG_WIDTH, FIG_HEIGHT))
cat("========================================\nDONE\n========================================\n")
