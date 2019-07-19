rm(list = ls())
graphics.off()

library(data.table)
library(reshape2)

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### ClinVar
CV <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf")
CV <- subset(CV, CV$CAT != "")
CV <- CV[,c("CHR", "POS", "REF", "CAT")]
CV <- CV[!duplicated(CV),] # removes SNVs with multiple ALT

CV_align <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_human_to_mouse_alignment.txt")
CV_align <- CV[CV_align, on = c("CHR", "POS", "REF")]
CV_align <- CV_align[,c("CHR", "POS", "REF", "CAT")]
CV_align <- CV_align[!duplicated(CV_align),]

CV <- as.data.frame(table(CV$CAT))
colnames(CV) <- c("annotation", "n_SNV")
CV_align <- as.data.frame(table(CV_align$CAT))
colnames(CV_align) <- c("annotation", "n_SNV_align")

h_dt_total_out <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")
CV_out <- merge(CV, CV_align)
CV_out$n_SNV_align_frac <- CV_out$n_SNV_align/CV_out$n_SNV
CV_out <- merge(CV_out, h_dt_total_out[,c("annotation", "ann_aligned_frac")])

### OUTPUT
fwrite(CV_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv")


### GWAS
GWAS <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_formatted.txt")
GWAS <- subset(GWAS, GWAS$CAT != "")
GWAS <- GWAS[,c("CHR", "POS", "CAT")]
GWAS <- GWAS[!duplicated(GWAS),] # removes SNVs with multiple ALT

GWAS_align <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_human_to_mouse_alignment.txt")
GWAS_align <- GWAS[GWAS_align, on = c("CHR", "POS")]
GWAS_align <- GWAS_align[,c("CHR", "POS", "REF", "CAT")]
GWAS_align <- GWAS_align[!duplicated(GWAS_align),]

GWAS <- as.data.frame(table(GWAS$CAT))
colnames(GWAS) <- c("annotation", "n_SNV")
GWAS_align <- as.data.frame(table(GWAS_align$CAT))
colnames(GWAS_align) <- c("annotation", "n_SNV_align")

h_dt_total_out <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")
GWAS_out <- merge(GWAS, GWAS_align)
GWAS_out$n_SNV_align_frac <- GWAS_out$n_SNV_align/GWAS_out$n_SNV
GWAS_out <- merge(GWAS_out, h_dt_total_out[,c("annotation", "ann_aligned_frac")])

### OUTPUT
fwrite(GWAS_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv")



#########

GWAS_out$CAT <- "GWAS"
CV_out$CAT <- "ClinVar"

tmp <- GWAS_out[,c("annotation", "ann_aligned_frac")]
tmp1 <- GWAS_out[,c("n_SNV", "n_SNV_align", "n_SNV_align_frac")]
colnames(tmp1) <- paste0("GWAS_", colnames(tmp1))
tmp2 <- CV_out[,c("n_SNV", "n_SNV_align", "n_SNV_align_frac")]
colnames(tmp2) <- paste0("ClinVar_", colnames(tmp2))
tmp_df <- cbind(tmp, tmp1, tmp2)
tmp_df <- melt(tmp_df, id.vars = c("annotation", "GWAS_n_SNV", "GWAS_n_SNV_align", "ClinVar_n_SNV", "ClinVar_n_SNV_align"))
tmp_df$variable <- as.character(tmp_df$variable)
tmp_df$variable[tmp_df$variable == "ann_aligned_frac"] <- "Total alignment"
tmp_df$variable[tmp_df$variable == "GWAS_n_SNV_align_frac"] <- "GWAS"
tmp_df$variable[tmp_df$variable == "ClinVar_n_SNV_align_frac"] <- "ClinVar"
tmp_df <- subset(tmp_df, tmp_df$annotation != "Exon - CDS")

p1 <- ggplot(data=tmp_df, aes(x=annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Huaman annotation") +
  ylab("Fraction") +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "black", size=1.5) +
  # geom_text(aes(label=ClinVar_n_SNV_align), vjust=1.6, color="white",
  #           position = position_dodge(0.9), size=3.5) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p1
# add n values (melt value2)
ggsave("~/Dropbox/Figure_SNVs_by_annotation.jpg", plot = p1, height = 6, width = 8)



tmp_x <- data.frame(
                    SNV_cat = c("GWAS", "GWAS", "ClinVar", "ClinVar"),
                    annotation_cat = c("protein-coding", "non-coding", "protein-coding", "non-coding"),
                    n_SNV = c(2747, 83313, 42033, 8651)
                    )
# sum(GWAS_out$n_SNV) - 2747
# sum(CV_out$n_SNV) - 42033

p2 <- ggplot(data=tmp_x, aes(x=SNV_cat, y=n_SNV, fill=annotation_cat)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("") +
  ylab("N SNVs") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p2
# add n values
ggsave("~/Dropbox/Figure_SNVs_by_coding.jpg", plot = p2, height = 6, width = 8)



### PLOT 

GWAS_plot <- GWAS_out
GWAS_plot$n_SNV_align_exp <- GWAS_plot$n_SNV * GWAS_plot$ann_aligned_frac
GWAS_plot$n_SNV_frac <- 1

# GWAS_plot <- GWAS_plot[,c("annotation", "n_SNV", "n_SNV_align", "n_SNV_align_exp")]
# colnames(GWAS_plot) <- c("annotation", "Human", "Mouse", "Mouse expected")
# GWAS_plot <- melt(GWAS_plot, by = "annotation")

GWAS_plot <- GWAS_plot[,c("annotation", "n_SNV_frac", "n_SNV_align_frac", "ann_aligned_frac")]
colnames(GWAS_plot) <- c("annotation", "Human", "Mouse", "Syntenic SNVs")
GWAS_plot <- melt(GWAS_plot, by = "annotation")

p_gwas <- ggplot(data=GWAS_plot, aes(x=annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("GWAS") +
  xlab("Huaman annotation") +
  coord_flip()
p_gwas

CV_plot <- CV_out
CV_plot$n_SNV_align_exp <- CV_plot$n_SNV * CV_plot$ann_aligned_frac
CV_plot$n_SNV_frac <- 1

# CV_plot <- CV_plot[,c("annotation", "n_SNV", "n_SNV_align", "n_SNV_align_exp")]
# colnames(CV_plot) <- c("annotation", "Human", "Mouse", "Mouse expected")
# CV_plot <- melt(CV_plot, by = "annotation")

CV_plot <- CV_plot[,c("annotation", "n_SNV_frac", "n_SNV_align_frac", "ann_aligned_frac")]
colnames(CV_plot) <- c("annotation", "Human", "Mouse", "Syntenic SNVs")
CV_plot <- melt(CV_plot, by = "annotation")

p_cv <- ggplot(data=CV_plot, aes(x=annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("ClinVar") +
  xlab("Huaman annotation") +
  coord_flip()
p_cv

pout <- grid.arrange(p_cv, p_gwas, nrow = 1, widths = c(5, 5))


GWAS_plot <- GWAS_out
GWAS_plot$n_SNV_align_exp <- GWAS_plot$n_SNV * GWAS_plot$ann_aligned_frac
GWAS_plot$n_SNV_frac <- 1
CV_plot <- CV_out
CV_plot$n_SNV_align_exp <- CV_plot$n_SNV * CV_plot$ann_aligned_frac
CV_plot$n_SNV_frac <- 1