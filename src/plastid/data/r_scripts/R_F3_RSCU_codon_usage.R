#!/usr/bin/env Rscript
# =============================================================================
# R_F3_RSCU_codon_usage.R — RSCU of the 64 codons, one column per species in the input
# Figure in the paper: Figure 3 (34 Iridaceae species; any input draws its own species)
# Usage:
#   Rscript R_F3_RSCU_codon_usage.R --input <data file> --outdir <output directory>
# Default: reads sample_data/Complete_Codon_Usage_Analysis.csv, writes <directory>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_common.R"))
init_script()
cat("\n=== R_F3_RSCU_codon_usage.R ===\n")
# CGAS Module 8: Codon Usage Comparison — ALL species, EXACT CGAS style
# Dodged bars, facet by amino acid, codon labels inside bars
# Draws every species column in the input, in the input order unless --species-order is given.
# Species names appear under EACH facet (not just bottom strip)
#
# USAGE:
#   cd Module8_Codon_Usage_Analysis/
#   Rscript /path/to/codon_34sp_cgas_style.R
#
# PREREQUISITES: Module 8 must be run first. Requires R package 'zoo'.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(RColorBrewer); library(zoo)
})

cat("\n========================================\n")
cat("CGAS MODULE 8: CODON USAGE COMPARISON — ALL SPECIES\n")
cat("========================================\n\n")

# === SPECIES COME FROM THE INPUT ===
# One row per codon, then one column per sample. Every column after the two metadata columns is a
# species, drawn in the input's own order. `--species-order "a,b,c"` fixes another order and
# `--exclude-species "x,y"` drops named columns; no species list is built in, so a panel outside
# the manuscript's 34 Iridaceae names draws its own species instead of an empty panel.
META_COLS <- c("Amino_Acid", "Codon")
SPECIES_ORDER  <- split_list(get_param("--species-order", ""))
EXCLUDE_SPECIES <- split_list(get_param("--exclude-species", ""))

# Read data
data <- read.csv(input_path("Complete_Codon_Usage_Analysis.csv"), check.names=FALSE, stringsAsFactors=FALSE)
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

# AA 3-letter codes (exactly like CGAS)
aa_name_to_code <- c("Alanine"="Ala","Arginine"="Arg","Asparagine"="Asn","Aspartate"="Asp",
  "Cysteine"="Cys","Glutamine"="Gln","Glutamate"="Glu","Glycine"="Gly","Histidine"="His",
  "Isoleucine"="Ile","Leucine"="Leu","Lysine"="Lys","Methionine"="Met","Phenylalanine"="Phe",
  "Proline"="Pro","Serine"="Ser","Threonine"="Thr","Tryptophan"="Trp","Tyrosine"="Tyr",
  "Valine"="Val","Stop"="Stop")

amino_acids <- data$Amino_Acid; amino_acids[amino_acids==""] <- NA
amino_acids <- na.locf(amino_acids)
aa3 <- aa_name_to_code[amino_acids]; aa3[is.na(aa3)] <- amino_acids[is.na(aa3)]

data_subset <- data[,c("Amino_Acid","Codon",available_sp)]
colnames(data_subset)[3:ncol(data_subset)] <- gsub("_"," ",colnames(data_subset)[3:ncol(data_subset)])
sp_display <- gsub("_"," ",available_sp)

data_long <- data_subset %>%
  pivot_longer(cols=-c(Amino_Acid,Codon),names_to="Species",values_to="RSCU") %>%
  mutate(Species=factor(Species,levels=sp_display),
         AA=aa3[match(Codon,data$Codon)],AA=factor(AA,levels=unique(aa3)),
         Codon=factor(Codon,levels=unique(Codon)))
data_long <- data_long[!is.na(data_long$AA),]
data_long <- data_long %>% group_by(AA,Species) %>% arrange(Codon) %>%
  mutate(cum_RSCU=cumsum(RSCU),ypos=cum_RSCU-RSCU/2) %>% ungroup()

n_sp <- length(available_sp)
# Scale to the input, not to 34: the fixed brewer stacks ran out past ~41 colours and returned
# NA for the rest, which draws an empty fill without a word.
species_colors <- grDevices::hcl.colors(n_sp, "Dark 3")

# CGAS-EXACT: dodged bars, facet by AA, codon labels, species names under EVERY facet
p <- ggplot(data_long,aes(x=Species,y=RSCU,fill=Species))+
  geom_bar(stat="identity",width=0.7,color="black",linewidth=0.3)+
  geom_text(aes(label=Codon,y=ypos),size=2.2,angle=90,vjust=0.5,fontface="bold")+
  facet_wrap(~AA,nrow=4,scales="free_x")+  # free_x = names under each facet
  scale_fill_manual(values=species_colors,guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.05)))+
  labs(y="Relative Synonymous Codon Usage",x="")+
  theme_minimal(base_size=11)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5.5,face="italic"),
        strip.text=element_text(size=10,face="bold"),
        panel.grid.major.y=element_line(color="gray90",linewidth=0.5),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.8),
        axis.line=element_line(color="black",linewidth=0.8),
        axis.ticks.y=element_line(color="black",linewidth=0.8),
        legend.position="none",plot.background=element_rect(fill="white",color=NA),
        plot.margin=margin(10,10,10,10))+
  geom_hline(yintercept=0,color="black",linewidth=0.8)

dir.create("Figures",showWarnings=FALSE)
# The drawing width follows the number of species columns, because a dodged bar per species at a
# fixed width stops being readable as the panel grows; --width/--height override it.
FIG_WIDTH  <- as.numeric(get_param("--width", max(34, 17 + 0.5 * n_sp)))
FIG_HEIGHT <- as.numeric(get_param("--height", 20))
# A large panel is a legitimate wide figure, but a raster device has a pixel ceiling: cap the
# resolution so the long side stays at 20,000 px (the paper's 34x20 in at 600 dpi was 20,400).
FIG_DPI <- min(600, floor(20000 / max(FIG_WIDTH, FIG_HEIGHT)))
ggsave("Figures/Codon_Usage_Comparison.png",p,width=FIG_WIDTH,height=FIG_HEIGHT,dpi=FIG_DPI,
       limitsize=FALSE)
pdf("Figures/Codon_Usage_Comparison.pdf",width=FIG_WIDTH,height=FIG_HEIGHT); print(p); dev.off()
cat(sprintf("\nSaved: Figures/Codon_Usage_Comparison.pdf and .png (%d species, %.0f x %.0f in)\n",
            n_sp, FIG_WIDTH, FIG_HEIGHT))
cat("========================================\nDONE\n========================================\n")
