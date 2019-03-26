rm(list=ls())
graphics.off()

library(ggplot2)

list.files("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/", pattern = "FRACTION")

PC_exon <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_PC_exon_POS_short.csv", header = F)
enhancer <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_enhancer_POS_short.csv", header = F)
TF <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_TF_binding_POS_short.csv", header = F)
NC_exon <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_NC_exon_POS_short.csv", header = F)
promoter <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_promoter_POS_short.csv", header = F)
UTR <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_UTR_POS_short.csv", header = F)
OC <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_open_chromatin_POS_short.csv", header = F)


colnames(PC_exon) <- c("percentile", "fraction")
colnames(enhancer) <- c("percentile", "fraction")
colnames(TF) <- c("percentile", "fraction")
colnames(NC_exon) <- c("percentile", "fraction")
colnames(promoter) <- c("percentile", "fraction")
colnames(UTR) <- c("percentile", "fraction")
colnames(OC) <- c("percentile", "fraction")

PC_exon$annotation <- "Exon, CDS"
enhancer$annotation <- "Enhancer"
TF$annotation <- "TF binding"
NC_exon$annotation <- "Exon, ncRNA"
promoter$annotation <- "Promoter"
UTR$annotation <- "Exon, UTR"
OC$annotation <- "Open chromatin"


df <- rbind(NC_exon, promoter, UTR, OC, enhancer, PC_exon, TF)

# df$Outlier <- 0
# df$Outlier[df$Percentile == 1] <- 1
# df$Outlier <- as.factor(df$Outlier)

ggplot(df, aes(x=percentile, y=fraction, color=annotation)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 


ggplot(NC_exon, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(OC, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(promoter, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(UTR, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(PC_exon, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(enhancer, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 

ggplot(TF, aes(x=percentile, y=fraction)) + 
  geom_point() +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Genomic territory (fraction)") 


#######




