rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

### SCRIPT 1:
### BY CHR
  ## Import alignment annotations:
    # calculate kmers (vectorise)
  ## Import human annotations
    # calculate annotations by base (data.table)
  ## Import ClinVar
    # calculate annotations by base (data.table)
    # calculate alignment, conservation, and kmer conservation
    # calculate total alignment, conservation, and kmer conservation for each annotation

# CLinvar columns:
#   H_CHR
#   H_POS
#   H_REF
#   H_KMER
#   H_ANN
#   M_CHR
#   M_POS
#   M_REF
#   M_KMER
#   M_ANN

### SET VARS
CV_file <- args[1]
HM_align_file <- args[2]
H_ann_file <- args[3]
out_file <- args[4]
chr <- as.integer(args[5])

# CV_file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf"
# HM_align_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr1.csv"
# H_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_mapping_chr1.csv"
# chr <- as.integer("1")

### IMPORT 

CV <- fread(CV_file)
HM_align <- fread(HM_align_file)
H_ann <- fread(H_ann_file)


### FORMAT ClinVar

CV <- subset(CV, CV$CHR == chr)
CV <- CV[,c("CHR", "POS", "REF")]
colnames(CV) <- c("H_CHR", "H_POS", "H_REF")
CV <- CV[!duplicated(CV),] # removes SNVs with multiple ALT


### HUMAN_MOUSE ALIGNMENT VARIABLES

colnames(HM_align) <- c("H_CHR", "H_POS", "H_REF", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_ANN")
HM_align <- subset(HM_align, HM_align$M_CHR %in% c(1:19, "X"))

## calculate alignment kmers (vectorise)
# add "-" to each non sequential POS
# get bases in kmers (replace non-sequential positions with "-")
H_base1 <- c("-", HM_align$H_REF[1:(length(HM_align$H_REF)-1)])
H_base1[c(1, (which(abs(diff(HM_align$H_POS))!=1))+1)] <- "-"
H_base2 <- HM_align$H_REF
H_base3 <- c(HM_align$H_REF[2:(length(HM_align$H_REF))], "-")
H_base3[c(which(abs(diff(HM_align$H_POS))!=1), length(HM_align$H_POS))] <- "-"
HM_align$H_Kmer <- paste0(H_base1, H_base2, H_base3)
rm(H_base1, H_base2, H_base3)

M_base1 <- c("-", HM_align$M_REF[1:(length(HM_align$M_REF)-1)])
M_base1[c(1, (which(abs(diff(HM_align$M_POS))!=1))+1)] <- "-"
M_base2 <- HM_align$M_REF
M_base3 <- c(HM_align$M_REF[2:(length(HM_align$M_REF))], "-")
M_base3[c(which(abs(diff(HM_align$M_POS))!=1), length(HM_align$M_POS))] <- "-"
HM_align$M_Kmer <- paste0(M_base1, M_base2, M_base3)
rm(M_base1, M_base2, M_base3)


### HUMAN ANNOTATIONS

# number code annotations (plyr)
H_ann$category <- mapvalues(H_ann$category, from=c("Exon - CDS",
                                                   "Exon - UTR",
                                                   "Exon - non-coding",
                                                   "Promoter",
                                                   "Enhancer",
                                                   "Open chromatin",
                                                   "TF binding",
                                                   "Promoter flanking",
                                                   "Intron",
                                                   "Unannotated"), 
                            to=c("A","B","C","D","E","F","G","H","I","J"))
# ensure start <= end 
tmp1 <- subset(H_ann, H_ann$end < H_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(H_ann, H_ann$end >= H_ann$start)
H_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)
# data.table to get annotations by POS
h_dt1 <- data.table(
  CHR = as.character(CV$H_CHR),
  INT = CV$H_POS)
h_dt2 <- data.table(
  CHR = as.character(H_ann$chromosome),
  INT_START = H_ann$start,
  INT_END = H_ann$end,
  CAT = H_ann$category)
# str(h_dt1)
# str(h_dt2)
# which(h_dt2$INT_START > h_dt2$INT_END)
h_dt1 <- h_dt2[h_dt1, on=.(CHR=CHR, INT_START<=INT, INT_END>=INT), allow.cartesian = T]
if (identical(h_dt1$INT_END, h_dt1$INT_START)){
  h_dt1 <- h_dt1[,c("CHR", "INT_END", "CAT")]
  colnames(h_dt1) <- c("CHR", "INT", "CAT")
  h_dt1 <- h_dt1[!duplicated(h_dt1),]
}
h_dt1 <- dcast(h_dt1, CHR + INT ~ CAT, value.var = "CAT")
h_dt1[is.na(h_dt1)] <- ""
if (ncol(h_dt1) > 3){
  h_dt1 <- data.table(
    H_POS = h_dt1$INT,
    H_ANN = apply(h_dt1[,3:ncol(h_dt1)], 1, paste0, collapse = ""))
} else {
  h_dt1 <- as.data.table(h_dt1)
  h_d1 <- h_dt1[,2:3]
  colnames(h_dt1) <- c("H_POS", "H_ANN")
}
# table(h_dt1$H_ANN)
CV <- CV[,c("H_CHR", "H_POS", "H_REF")]
CV <- h_dt1[CV, on=c("H_POS")]
rm(H_ann, h_dt1, h_dt2)


### FORMAT OUTPUT

CV_nonHM <- subset(CV, !CV$H_POS %in% HM_align$H_POS)
HM_CV <- subset(HM_align,  HM_align$H_POS %in% CV$H_POS)
CV <- rbind.fill(HM_CV, CV_nonHM)
CV <- CV[,c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")]
rm(HM_CV, CV_nonHM, HM_align)

### OUTPUT
fwrite(CV, out_file, col.names = F)



##### 

# # identify ind for bases at the begining of human sequences
# H_begin <- c(1, ifelse(abs(diff(test$H_POS))==1, 0, 1))
# # identify ind for bases at the end of human sequences
# H_end <- c(ifelse(abs(diff(test$H_POS))==1, 0, 1), 1)
# 
# H_REF_begin <- test$H_REF
# H_REF_begin[which(H_begin == 1)] <- "-"
# 
# H_REF_end <- test$H_REF
# H_REF_end[which(H_end == 1)] <- "-"
# 
# base1 <- H_REF_begin
# base2 <- test$H_REF
# base3 <- H_REF_end
# 
# test$H_Kmer <- paste0(base1, base2, base3)
