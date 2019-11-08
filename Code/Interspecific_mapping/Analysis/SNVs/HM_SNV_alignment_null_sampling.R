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

### SCRIPT:
### BY CHR
### NULL SAMPLING
## Import alignment annotations:
# calculate kmers (vectorise)
## Import human reference and annotations
# calculate annotations by base (data.table)
## Import ClinVar formatted
## Calculate number of ClinVar SNVs per annotation grouping
## stratified sampling of human POS by chromosome and annotation grouping
## calculate total alignment, conservation, and kmer conservation by annotation for each sample

### SET VARS
CV_file <- args[1]
HM_align_file <- args[2]
H_ref_file <- args[3]
H_ann_file <- args[4]
out_file <- args[5]
chr <- as.integer(args[6])
n_repeats <- as.integer(args[7])

# CV_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_DDC_mouse_mapping_Hchr1.csv"
# HM_align_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr1.csv"
# H_ref_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr1.txt"
# H_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_DDC_GWAS_null_sample_chr1.csv"
# chr <- as.integer("1")

# CV_file <- "/well/lindgren/George/Data/GWAS/formatted/GWAS_DDC_mouse_mapping_Hchr22.csv"
# HM_align_file <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr22.csv"
# H_ref_file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr22.txt"
# H_ann_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "/well/lindgren/George/Data/GWAS/formatted/GWAS_DDC_mouse_mapping_null_sampling_Hchr22.csv"
# chr <- as.integer("22")

### IMPORT 

CV <- fread(CV_file)
HM_align <- fread(HM_align_file)
H_ref <- fread(H_ref_file)
H_ann <- fread(H_ann_file)

###########
### FORMAT
#########

colnames(CV) <- c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")
colnames(HM_align) <- c("H_CHR", "H_POS", "H_REF", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_ANN")
HM_align <- subset(HM_align, HM_align$M_CHR %in% c(1:19, "X"))
H_ref <- H_ref[,c("POS", "REF")]
colnames(H_ref) <- c("H_POS", "H_REF")
H_ref <- subset(H_ref, H_ref$H_REF %in% c("A", "T", "C", "G"))

### ALIGNMENT KMERS

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


### REFERENCE ANNOTATIONS

# colnames
colnames(H_ann) <- c("category", "chromosome", "start", "end")

# categories
H_ann$category <- as.character(H_ann$category)
H_ann$category[H_ann$category == "TF binding" | H_ann$category == "Open chromatin"] <- "Miscellaneous"

# number code annotations (plyr)
H_ann$category <- mapvalues(H_ann$category, from=c("Exon - CDS",
                                                   "Exon - UTR",
                                                   "Exon - other",
                                                   "Promoter",
                                                   "Enhancer - proximal",
                                                   "Enhancer - distal",
                                                   "CTCF binding",
                                                   "Miscellaneous",
                                                   "Intron - proximal",
                                                   "Intron - distal",
                                                   "Unannotated"), 
                            to=c("A","B","C","D","E","F","G","H","I","J","K"))
# ensure start <= end 
tmp1 <- subset(H_ann, H_ann$end < H_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(H_ann, H_ann$end >= H_ann$start)
H_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)


### Split into subsets as chr 1 and 2 too large to format as one
bases <- c("A", "T", "C", "G")
H_ref_ann <- data.table()
for(i in 1:4){
  
  # subset by base
  H_ref_split <- subset(H_ref, H_ref$H_REF == bases[i])
  
  # data.table to get annotations by POS
  h_dt1 <- data.table(
    CHR = rep(as.character(chr), nrow(H_ref_split)), 
    INT = H_ref_split$H_POS)
  h_dt2 <- data.table(
    CHR = as.character(H_ann$chromosome),
    INT_START = H_ann$start,
    INT_END = H_ann$end,
    CAT = as.character(H_ann$category))
  
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
    h_dt1$H_REF <- bases[i]
    H_ref_ann <- rbind(H_ref_ann, h_dt1)
}
H_ref_ann <- H_ref_ann[order(H_POS),]
rm(H_ann, H_ref_split, h_dt1, h_dt2, H_ref)


### STRATIFIED SAMPLING

# Calculate number of ClinVar SNVs per annotation grouping
CV_cats <- as.data.table(table(CV$H_ANN))
colnames(CV_cats) <- c("cat", "N")
# ref_cats <- as.data.table(table(H_ref_ann$H_ANN))
# colnames(ref_cats) <- c("cat", "N")

### FUNCTION that samples H_POS by category and returns variable counts
sample_sample <- function(i, H_ref_cat, CV_cats){
  
  # sample H_POS from H_ref
  H_sample <- subset(H_ref_cat, H_ref_cat$H_POS %in% sample(H_ref_cat$H_POS, size = CV_cats$N[i], replace = FALSE))
  
  # merge with HM_align
  H_sample <- H_sample[,c("H_POS")]
  H_sample <- HM_align[H_sample, on = "H_POS"]
  
  category <- H_ref_cat$H_ANN[1]
  # count the total SNVs in cat
  total_H <- nrow(H_sample)
  # count the number of SNVs that align to mouse
  total_M_align <- nrow(H_sample[!is.na(H_sample$M_POS),])
  # count the number of SNVs that are conserved in mouse
  total_M_conserv <- nrow(H_sample[H_sample$H_REF == H_sample$M_REF,])
  # count the number of SNVs with the 3-mer conserved in mouse
  total_M_conserv_K <- nrow(H_sample[H_sample$H_Kmer == H_sample$M_Kmer,])
  
  # get total for each annotation in mouse
  total_M_align_A <- nrow(H_sample[H_sample$M_ANN %like% "A",])
  total_M_align_B <- nrow(H_sample[H_sample$M_ANN %like% "B",])
  total_M_align_C <- nrow(H_sample[H_sample$M_ANN %like% "C",])
  total_M_align_D <- nrow(H_sample[H_sample$M_ANN %like% "D",])
  total_M_align_E <- nrow(H_sample[H_sample$M_ANN %like% "E",])
  total_M_align_F <- nrow(H_sample[H_sample$M_ANN %like% "F",])
  total_M_align_G <- nrow(H_sample[H_sample$M_ANN %like% "G",])
  total_M_align_H <- nrow(H_sample[H_sample$M_ANN %like% "H",])
  total_M_align_I <- nrow(H_sample[H_sample$M_ANN %like% "I",])
  total_M_align_J <- nrow(H_sample[H_sample$M_ANN %like% "J",])
  total_M_align_K <- nrow(H_sample[H_sample$M_ANN %like% "K",])
  
  
  return(c(category, total_H, total_M_align, total_M_conserv, total_M_conserv_K,
           total_M_align_A, total_M_align_B, total_M_align_C, total_M_align_D, total_M_align_E,
           total_M_align_F, total_M_align_G, total_M_align_H, total_M_align_I, total_M_align_J, total_M_align_K))
}

# for each annotation category sample 
out_list <- list()
for (i in 1:nrow(CV_cats)){
  
  # subset H_ref by CV_cats[i]
  H_ref_cat <- subset(H_ref_ann, H_ref_ann$H_ANN == CV_cats$cat[i])
  # run function X times
  AA <- replicate(n_repeats, sample_sample(i = i, H_ref_cat = H_ref_cat, CV_cats = CV_cats), simplify = FALSE)
  # index output
  out_matrix <- do.call("rbind", AA)
  out_matrix <- cbind(out_matrix, c(1:nrow(out_matrix)))
  out_list[[i]] <- out_matrix
  print(i)
}

out_dt <- as.data.table(do.call("rbind", out_list))
colnames(out_dt) <- c("sample_ind", "H_ANN", "total_H", "total_M_align", "total_M_conserv", "total_M_conserv_K",
                      "total_M_align_A", "total_M_align_B", "total_M_align_C", "total_M_align_D", "total_M_align_E",
                      "total_M_align_F", "total_M_align_G", "total_M_align_H", "total_M_align_I", "total_M_align_J", "total_M_align_K")
out_dt$H_CHR <- chr


### EXPORT
fwrite(out_dt, out_file, col.names = F)

#####

# Error in sample.int(length(x), size, replace, prob) : 
#   invalid first argument
# Calls: replicate ... subset.data.table -> eval -> eval -> %in% -> sample -> sample.int
# 13

