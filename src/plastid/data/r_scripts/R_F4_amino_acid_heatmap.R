#!/usr/bin/env Rscript
# =============================================================================
# R_F4_amino_acid_heatmap.R — Heatmap thành phần acid amin (%) của 34 loài Iridaceae
# Hình trong bài: Figure 4
# Gọi:
#   Rscript R_F4_amino_acid_heatmap.R --input <file dữ liệu> --outdir <thư mục kết quả>
# Mặc định: đọc data_mau/Complete_Amino_Acid_Analysis.csv, ghi <thư mục>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_ham_chung.R"))
khoi_tao()
cat("\n=== R_F4_amino_acid_heatmap.R ===\n")
# CGAS Module 9: Amino Acid Comparison — ALL species in EXACT CGAS style
# Dodged bars, x=amino acid (3-letter, sorted by avg frequency), fill=species
# NO percentage labels (unreadable with many species)
# Species ordered by genetic distance from focal species
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

# === CONFIGURATION (edit for your dataset) ===
OUTGROUP_COLS <- c("Curculigo_orchioides", "Molineria_capitulata")
SP_ORDER <- c(  # order by genetic distance from focal (closest first)
  "Eleutherine_bulbosa", "Phalocallis_coelestis", "Tigridia_pavonia",
  "Trimezia_steyermarkii", "Neomarica_candida", "Diplarrena_moraea",
  "Iris_domestica", "Iris_afghanica", "Iris_lortetii",
  "Iris_bismarckiana", "Iris_gatesii", "Iris_hermona",
  "Libertia_pulchella", "Iris_orchioides", "Isophysis_tasmanica",
  "Moraea_spathulata", "Olsynium_douglasii", "Moraea_polystachya",
  "Sisyrinchium_idahoense", "Gladiolus_communis", "Watsonia_pillansii",
  "Nivenia_stokoei", "Chasmanthe_floribunda", "Aristea_ecklonii",
  "Freesia_laxa", "Crocus_cartwrightianus", "Crocus_sativus",
  "Tritonia_disticha", "Patersonia_fragilis", "Witsenia_maura",
  "Dierama_igneum", "Hesperantha_coccinea", "Gladiolus_x_gandavensis",
  "Iris_missouriensis"
)

data <- read.csv(duong_dan_vao("Complete_Amino_Acid_Analysis.csv"), check.names=FALSE, stringsAsFactors=FALSE)
for (col in OUTGROUP_COLS) { if (col %in% colnames(data)) data <- data %>% select(-all_of(col)) }

available_sp <- intersect(SP_ORDER, colnames(data))
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
species_colors <- c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),
                     brewer.pal(8,"Set2"),brewer.pal(9,"Set1"))[1:n_sp]

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
ggsave("Figures/Amino_Acid_Comparison_34sp.png",p,width=22,height=12,dpi=600)
pdf("Figures/Amino_Acid_Comparison_34sp.pdf",width=22,height=12); print(p); dev.off()
cat("\nSaved: Figures/Amino_Acid_Comparison_34sp.pdf and .png\n")
cat("========================================\nDONE\n========================================\n")
