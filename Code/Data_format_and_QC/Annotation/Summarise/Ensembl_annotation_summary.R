# Table -- Summary details for unranked annotation categories. 
# (Number of sequences, sequence length (bp) IQR, total bp, genome coverage (%)). 

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### IMPORT

mann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv")
hann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv")

h_dt <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_annotation_overlap.csv")
m_dt <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_annotation_overlap.csv")

### FORMAT

colnames(mann) <- c("category", "chromosome", "start", "end")
colnames(hann) <- c("category", "chromosome", "start", "end")
colnames(h_dt) <- c("category", "N")
colnames(m_dt) <- c("category", "N")

cats <- c("Exon - CDS",
          "Exon - UTR",
          "Exon - other",
          "Promoter",
          "Enhancer - proximal",
          "Enhancer - distal",
          "CTCF binding",
          "Open chromatin",
          "TF binding",
          "Intron - proximal",
          "Intron - distal",
          "Unannotated")
letter <- c("A","B","C","D","E","F","G","H","I","J", "K", "L")

mann$length <- (mann$end - mann$start) + 1
hann$length <- (hann$end - hann$start) + 1

# GRCh38 chromosome lengths
H_chr_length <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                  145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                  101991189,  90338345,  83257441,  80373285,  58617616,  64444167,  46709983,
                  50818468)
M_chr_length <- c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 
                  129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 
                  104043685, 98207768, 94987271, 90702639, 61431566)


### Sequence IQR

n_sequences <- rep(NA, length(cats))
Q1 <- rep(NA, length(cats))
Q2 <- rep(NA, length(cats))
Q3 <- rep(NA, length(cats))
no_overlap <- rep(NA, length(letter))
overlap <- rep(NA, length(letter))
CDS_overlap <- rep(NA, length(letter))
n_bp <- rep(NA, length(cats))
genomic_percentage <- rep(NA, length(cats))

for (cat in 1:length(cats)){
  sub <- subset(hann, hann$category == cats[cat])
  n_sequences[cat] <- as.integer(nrow(sub))
  Q1[cat] <- as.integer(quantile(sub$length, 0.25))
  Q2[cat] <- as.integer(quantile(sub$length, 0.5))
  Q3[cat] <- as.integer(quantile(sub$length, 0.75))
  # n_bp[cat] <- as.integer(sum(sub$length))
  # genomic_percentage[cat] <- (sum(sub$length)/sum(H_chr_length))*100
}
for (i in 1:length(letter)){
  tmp1 <- h_dt[h_dt$category %like% letter[i],]
  n_bp[i] <- sum(tmp1$N)
  genomic_percentage[i] <- (n_bp[i]/sum(h_dt$N))*100
  no_overlap[i] <- ((sum(h_dt$N[h_dt$category == letter[i]]))/n_bp[i])*100
  tmp2 <- tmp1[tmp1$category != letter[i],]
  overlap[i] <- (sum(tmp2$N)/n_bp[i])*100
  tmp3 <- tmp1[tmp1$category %like% "A",]
  CDS_overlap[i] <- (sum(tmp3$N)/n_bp[i])*100
  print(i)
}
hann_summary <- data.table(species = rep("Human", length(cats)),
                           category = cats,
                           n_sequences = n_sequences,
                           length_Q1 = Q1,
                           length_Q2 = Q2,
                           length_Q3 = Q3,
                           n_bp = n_bp,
                           no_overlap_percentage = round(no_overlap, digits = 2),
                           overlap_percentage = round(overlap, digits = 2),
                           CDS_overlap_percentage = round(CDS_overlap, digits = 2),
                           genomic_percentage = round(genomic_percentage, digits = 2))

for (cat in 1:length(cats)){
  sub <- subset(mann, mann$category == cats[cat])
  n_sequences[cat] <- as.integer(nrow(sub))
  Q1[cat] <- as.integer(quantile(sub$length, 0.25))
  Q2[cat] <- as.integer(quantile(sub$length, 0.5))
  Q3[cat] <- as.integer(quantile(sub$length, 0.75))
  # n_bp[cat] <- as.integer(sum(sub$length))
  # genomic_percentage[cat] <- (sum(sub$length)/sum(M_chr_length))*100 
}
for (i in 1:length(letter)){
  tmp1 <- m_dt[m_dt$category %like% letter[i],]
  n_bp[i] <- sum(tmp1$N)
  genomic_percentage[i] <- (n_bp[i]/sum(m_dt$N))*100
  no_overlap[i] <- ((sum(m_dt$N[m_dt$category == letter[i]]))/n_bp[i])*100
  tmp2 <- tmp1[tmp1$category != letter[i],]
  overlap[i] <- (sum(tmp2$N)/n_bp[i])*100
  tmp3 <- tmp1[tmp1$category %like% "A",]
  CDS_overlap[i] <- (sum(tmp3$N)/n_bp[i])*100
  print(i)
}
mann_summary <- data.table(species = rep("Mouse", length(cats)),
                           category = cats,
                           n_sequences = n_sequences,
                           length_Q1 = Q1,
                           length_Q2 = Q2,
                           length_Q3 = Q3,
                           n_bp = n_bp,
                           no_overlap_percentage = round(no_overlap, digits = 2),
                           overlap_percentage = round(overlap, digits = 2),
                           CDS_overlap_percentage = round(CDS_overlap, digits = 2),
                           genomic_percentage = round(genomic_percentage, digits = 2))

### Fraction genome coverage
sum(hann_summary$genomic_percentage*(hann_summary$no_overlap)/100)-sum(hann_summary$genomic_percentage[10:11])
sum(mann_summary$genomic_percentage*(mann_summary$no_overlap)/100)-sum(mann_summary$genomic_percentage[10:11])

output <- rbind(hann_summary, mann_summary)
output$category <- factor(output$category, levels = cats)
output <- output[order(category),]

### OUTPUT

fwrite(output, "~/Dropbox/PhD/Data/Ensembl/Annotation/HM_annotation_coverage_summary.csv")

#####







