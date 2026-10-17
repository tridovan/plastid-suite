#!/usr/bin/env Rscript
# =============================================================================
# R_F3_RSCU_codon_usage.R — RSCU 64 codon của 34 loài Iridaceae, nhóm theo acid amin
# Hình trong bài: Figure 3
# Gọi:
#   Rscript R_F3_RSCU_codon_usage.R --input <file dữ liệu> --outdir <thư mục kết quả>
# Mặc định: đọc data_mau/Complete_Codon_Usage_Analysis.csv, ghi <thư mục>/Figures/*.png + *.pdf
# =============================================================================
.local_dir <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.local_dir, "00_ham_chung.R"))
khoi_tao()
cat("\n=== R_F3_RSCU_codon_usage.R ===\n")
# CGAS Module 8: Codon Usage Comparison — ALL species, EXACT CGAS style
# Dodged bars, facet by amino acid, codon labels inside bars
# Species ordered by genetic distance from focal species (Eleutherine for Iridaceae)
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

# Read data
data <- read.csv(duong_dan_vao("Complete_Codon_Usage_Analysis.csv"), check.names=FALSE, stringsAsFactors=FALSE)
for (col in OUTGROUP_COLS) { if (col %in% colnames(data)) data <- data %>% select(-all_of(col)) }
available_sp <- intersect(SP_ORDER, colnames(data))
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
species_colors <- c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),
                     brewer.pal(8,"Set2"),brewer.pal(9,"Set1"))[1:n_sp]

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
ggsave("Figures/Codon_Usage_Comparison_34sp.png",p,width=34,height=20,dpi=600)
pdf("Figures/Codon_Usage_Comparison_34sp.pdf",width=34,height=20); print(p); dev.off()
cat("\nSaved: Figures/Codon_Usage_Comparison_34sp.pdf and .png\n")
cat("========================================\nDONE\n========================================\n")
