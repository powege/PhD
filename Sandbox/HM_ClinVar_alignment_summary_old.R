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

### GENERATE NULL DATASET! 
### ENSURE REF AND KMER MATCH!
### MOUSE KMERS MAY BE REVERSED -- rerun subsetting on H POS

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

# CV_file <- "/well/lindgren/George/Data/ClinVar/formatted/ClinVar_pathogenic_snps_QCed.vcf"
# HM_align_file <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr1.csv"
# H_ann_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "/well/lindgren/George/Data/ClinVar/formatted/ClinVar_mouse_mapping.csv"
# chr <- as.integer("1")

# CV_file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf"
# HM_align_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr1.csv"
# H_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# # out_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv"
# chr <- as.integer("1")

### FORMAT ClinVar
CV <- fread(CV_file)
CV <- subset(CV, CV$CHR == chr)
CV <- CV[,c("CHR", "POS", "REF")]
colnames(CV) <- c("H_CHR", "H_POS", "H_REF")
CV <- CV[!duplicated(CV),] # removes SNVs with multiple ALT


### HUMAN_MOUSE ALIGNMENT VARIABLES
HM_align <- fread(HM_align_file)
colnames(HM_align) <- c("H_CHR", "H_POS", "H_REF", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_ANN")
HM_align <- subset(HM_align, HM_align$M_CHR %in% c(1:19, "X"))
CV_nonHM <- subset(CV, !CV$H_POS %in% HM_align$H_POS)
HM_CV <- subset(HM_align,  HM_align$H_POS %in% CV$H_POS)


### KMERS (aligned)
# for loop by chromosome
dt_out <- data.table()
for (m_chr in unique(HM_CV$M_CHR)){
  
  HM_align_sub <- subset(HM_align, HM_align$M_CHR == m_chr)
  HM_CV_sub <- subset(HM_CV, HM_CV$M_CHR == m_chr)

  HM_CV_sub$M_Kmer <- NA
  HM_CV_sub$H_Kmer <- NA
  
  for (i in 1:nrow(HM_CV_sub)){
    HK1 <- HM_align_sub$H_REF[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]-1))][1]
    HK2 <- HM_align_sub$H_REF[which(HM_align_sub$H_POS == HM_CV_sub$H_POS[i])][1]
    HK3 <- HM_align_sub$H_REF[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]+1))][1]
    
    cond_neg <- abs(HM_align_sub$M_POS[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]-1))][1] - 
                      HM_align_sub$M_POS[which(HM_align_sub$H_POS == HM_CV_sub$H_POS[i])][1])
    cond_pos <- abs(HM_align_sub$M_POS[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]+1))][1] - 
                      HM_align_sub$M_POS[which(HM_align_sub$H_POS == HM_CV_sub$H_POS[i])][1])
    
    if ( cond_neg == 1 & !is.na(cond_neg) ){
      MK1 <- HM_align_sub$M_REF[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]-1))][1]
    }
    MK2 <- HM_align_sub$M_REF[which(HM_align_sub$M_POS == HM_CV_sub$M_POS[i])][1]
    if ( cond_neg == 1 & !is.na(cond_neg) ){
      MK3 <- HM_align_sub$M_REF[which(HM_align_sub$H_POS == (HM_CV_sub$H_POS[i]+1))][1]
    }
  
    if (is.na(HK1)){ HK1 = "-"}
    if (is.na(HK3)){ HK3 = "-"}
    if (is.na(MK1)){ MK1 = "-"}
    if (is.na(MK3)){ MK3 = "-"}
    if (length(HK1) == 0){ HK1 = "-"}
    if (length(HK3) == 0){ HK3 = "-"}
    if (length(MK1) == 0){ MK1 = "-"}
    if (length(MK3) == 0){ MK3 = "-"}
    
    HM_CV_sub$H_Kmer[i] <- paste0(HK1, HK2, HK3)
    HM_CV_sub$M_Kmer[i] <- paste0(MK1, MK2, MK3)
    
    print(i)
  }
  
  dt_out <- rbind(dt_out, HM_CV_sub)
}


# HM_CV <- cbind(H_dt_out, M_dt_out) ### CANNOT CBIND AS REORDERED
CV <- rbind.fill(dt_out, CV_nonHM)
rm(HM_CV, CV_nonHM, dt_out, HM_align)

### HUMAN ANNOTATIONS
H_ann <- fread(H_ann_file)
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
CV <- CV[,c("H_CHR", "H_POS", "H_REF", "H_Kmer", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")]
CV <- h_dt1[CV, on=c("H_POS")]
rm(H_ann, h_dt1, h_dt2)


### OUTPUT
CV <- CV[,c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")]
fwrite(CV, out_file, col.names = F)


#####

# ### KMERS (aligned)
# 
# # CV <- CV[which(CV$H_POS %in% HM_align$H_POS),] # to test
# H_CV <- CV[,c("H_POS", "H_ANN", "H_CHR", "H_REF")]
# H_CV <- H_CV[!duplicated(H_CV),]
# H_CV <- H_CV[complete.cases(H_CV),]
# H_CV <- H_CV[order(H_CV$H_POS),]
# 
# # remove human aligned duplicates
# HM_align_temp <- HM_align[!duplicated(HM_align$H_POS),]
# # identify which POS do not have an aligned k3mer
# noK_b1_int <- which( (HM_align_temp$H_POS[which(HM_align_temp$H_POS %in% H_CV$H_POS)] -1) 
#                      != 
#                       HM_align_temp$H_POS[which(HM_align_temp$H_POS %in% H_CV$H_POS) -1] )
# base1 <- HM_align_temp$H_REF[which(HM_align_temp$H_POS %in% H_CV$H_POS) - 1]
# 
# base2 <- HM_align_temp$H_REF[which(HM_align_temp$H_POS %in% H_CV$H_POS)]
# 
# noK_b3_int <- which( (HM_align_temp$H_POS[which(HM_align_temp$H_POS %in% H_CV$H_POS)] +1) 
#                      != 
#                       HM_align_temp$H_POS[which(HM_align_temp$H_POS %in% H_CV$H_POS) +1] ) 
# base3 <- HM_align_temp$H_REF[which(HM_align_temp$H_POS %in% H_CV$H_POS) + 1]
# base3[noK_b3_int] <- "-"
# 
# H_CV$H_k3mer <- paste0(base1, base2, base3)
# 
# M_CV <- CV[,c("M_POS", "M_ANN", "M_CHR", "M_REF")]
# M_CV <- M_CV[!duplicated(M_CV),]
# M_CV <- M_CV[complete.cases(M_CV),]
# M_CV <- M_CV[order(M_CV$M_POS),]
# 
# # remove human aligned duplicates
# HM_align_temp <- HM_align[!duplicated(HM_align$M_POS),]
# # identify which POS do not have an aligned k3mer 
# # (mouse seq can be in reverse)
# # (ensure chr match)
# noK_b1_int <-  which( abs( (HM_align_temp$M_POS[which(HM_align_temp$M_POS %in% M_CV$M_POS) -1]) -
#                            HM_align_temp$M_POS[which(HM_align_temp$M_POS %in% M_CV$M_POS)] ) != 1 )
# base1 <- HM_align_temp$M_REF[which(HM_align_temp$M_POS %in% M_CV$M_POS) - 1]
# base1[noK_b1_int] <- "-"
# 
# base2 <- HM_align_temp$M_REF[which(HM_align_temp$M_POS %in% M_CV$M_POS)]
# 
# noK_b3_int <-  which( abs( (HM_align_temp$M_POS[which(HM_align_temp$M_POS %in% M_CV$M_POS) +1]) -
#                              HM_align_temp$M_POS[which(HM_align_temp$M_POS %in% M_CV$M_POS)] ) != 1 )
# base3 <- HM_align_temp$M_REF[which(HM_align_temp$M_POS %in% M_CV$M_POS) + 1]
# base3[noK_b3_int] <- "-"
# 
# M_CV$M_k3mer <- paste0(base1, base2, base3)
# 
# H_CV <- as.data.table(H_CV)
# M_CV <- as.data.table(M_CV)
# CV <- H_CV[CV, on = c("H_POS", "H_ANN" , "H_CHR",  "H_REF")]
# CV <- M_CV[CV, on = c("M_POS", "M_ANN" , "M_CHR",  "M_REF")]
# 
# rm(HM_align, HM_align_temp, CV_sub, HM_sub, H_CV, M_CV, base1, base2, base3, noK_b1_int, noK_b3_int)


# ### KMERS # works if no duplicated POS
# H_ref <- fread(H_ref_file)
# H_ref <- H_ref[,c("POS", "REF")]
# # CV <- CV[which(CV$H_POS %in% H_ref$POS),] # to test
# H_CV <- CV[,c("H_POS", "H_ANN", "H_CHR", "H_REF")]
# H_CV <- H_CV[!duplicated(H_CV),]
# H_CV <- H_CV[complete.cases(H_CV),]
# H_CV <- H_CV[order(H_CV$H_POS),]
# base1 <- H_ref$REF[which(H_ref$POS %in% H_CV$H_POS) - 1]
# base2 <- H_ref$REF[which(H_ref$POS %in% H_CV$H_POS)]
# base3 <- H_ref$REF[which(H_ref$POS %in% H_CV$H_POS) + 1]
# H_CV$H_k3mer <- paste0(base1, base2, base3)
# rm(H_ref)
# 
# ### cannot subset by chr for mouse due to alignment
# M_ref <- fread(M_ref_file)
# M_ref <- M_ref[,c("POS", "REF")]
# M_CV <- CV[,c("M_POS", "M_ANN", "M_CHR", "M_REF")]
# M_CV <- M_CV[!duplicated(M_CV),]
# M_CV <- M_CV[complete.cases(M_CV),]
# M_CV <- M_CV[order(M_CV$M_POS),]
# base1 <- M_ref$REF[which(M_ref$POS %in% M_CV$M_POS) - 1]
# base2 <- M_ref$REF[which(M_ref$POS %in% M_CV$M_POS)]
# base3 <- M_ref$REF[which(M_ref$POS %in% M_CV$M_POS) + 1]
# M_CV$M_k3mer <- paste0(base1, base2, base3)
# rm(M_ref)
# 
# CV <- H_CV[CV, on = c("H_POS", "H_ANN" , "H_CHR",  "H_REF")]
# CV <- M_CV[CV, on = c("M_POS", "M_ANN" , "M_CHR",  "M_REF")]

### STACK OVERFLOW

# set.seed(1)
# dt1 <- data.table(INT = c(9, 10, 45, 50, 3, 37, 9, 15, 48, 82),
#                   CAT = c(rep("X", 4), rep("Y", 6)))
# dt2 <- data.table(INT = c(1:10, 15:1, 31:55, 41:60, 90:81),
#                   CAT = c(rep("X", 10), rep("Y", 15), rep("Y", 25), rep("X", 20), rep("Y", 10)),
#                   BASE = sample (c("A", "G", "C", "T"), 80, replace = T))
# 
# dt1$Kmer <- NA
# for (i in 1:nrow(dt1)){
#   dt2_sub <- subset(dt2, dt2$CAT == dt1$CAT[i])
#   K1 <- dt2_sub$BASE[which(dt2_sub$INT == (dt1$INT[i]-1))]
#   K2 <- dt2_sub$BASE[which(dt2_sub$INT == dt1$INT[i])]
#   K3 <- dt2_sub$BASE[which(dt2_sub$INT == (dt1$INT[i]+1))]
#   if (length(K1) == 0){ K1 = "-"}
#   if (length(K3) == 0){ K3 = "-"}
#   dt1$Kmer[i] <- paste0(K1, K2, K3)
#   print(i)
# }
# 
# 
# # for loop by chromosome
# dt_out <- data.table()
# for (m_chr in unique(HM_CV$M_CHR)){
# 
#   HM_align_sub <- subset(HM_align, HM_align$M_CHR == m_chr)
#   HM_CV_sub <- subset(HM_CV, HM_CV$M_CHR == m_chr)
# 
#   HM_CV_sub$Kmer <- NA
#   for (i in 1:nrow(HM_CV_sub)){
#     K1 <- HM_align_sub$M_REF[which(HM_align_sub$M_POS == (HM_CV_sub$M_POS[i]-1))][1]
#     K2 <- HM_align_sub$M_REF[which(HM_align_sub$M_POS == HM_CV_sub$M_POS[i])][1]
#     K3 <- HM_align_sub$M_REF[which(HM_align_sub$M_POS == (HM_CV_sub$M_POS[i]+1))][1]
#     if (length(K1) == 0){ K1 = "-"}
#     if (length(K3) == 0){ K3 = "-"}
#     HM_CV_sub$Kmer[i] <- paste0(K1, K2, K3)
#     print(i)
#   }
# 
#   dt_out <- rbind(dt_out, HM_CV_sub)
# }
