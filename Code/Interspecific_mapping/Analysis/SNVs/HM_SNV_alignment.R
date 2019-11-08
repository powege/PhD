rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(stringi)

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
## Import human SNVs
# calculate annotations by base (data.table)
# calculate alignment, conservation, and kmer conservation
# calculate total alignment, conservation, and kmer conservation for each annotation

# OUTPUT columns:
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

# CV_file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_nonDDC_QCed.txt"
# HM_align_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr11.csv"
# H_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_mouse_mapping_chr11.csv"
# chr <- as.integer("11")

# CV_file <- "/well/lindgren/George/Data/GWAS/formatted/GWAS_catalog_nonDDC_QCed.txt"
# HM_align_file <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr11.csv"
# H_ann_file <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "/well/lindgren/George/Data/GWAS/formatted/GWAS_nonDDC_mouse_mapping_Hchr11.csv"
# chr <- as.integer("11")

### FUNCTIONS

# FUNCTION that sorts the caracters in a string, or vercor of strings
sort_cat3 <- function(strings){
  apply(stri_extract_all_regex(strings, "\\p{L}", simplify = TRUE), 1, function(i){
    stri_join(stri_sort(i), collapse = "")
  })
}

### IMPORT 

CV <- fread(CV_file)
HM_align <- fread(HM_align_file)
H_ann <- fread(H_ann_file)


### FORMAT 

CV <- subset(CV, CV$CHR == chr)
CV <- CV[,c("CHR", "POS")]
colnames(CV) <- c("H_CHR", "H_POS")
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

# data.table to get annotations by POS
h_dt1 <- data.table(
  CHR = as.character(CV$H_CHR),
  INT = CV$H_POS)
h_dt2 <- data.table(
  CHR = as.character(H_ann$chromosome),
  INT_START = H_ann$start,
  INT_END = H_ann$end,
  CAT = as.character(H_ann$category))
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
# h_dt1$H_ANN <- sort_cat3(h_dt1$H_ANN)
# table(h_dt1$H_ANN)
CV <- CV[,c("H_CHR", "H_POS")]
CV <- h_dt1[CV, on=c("H_POS")]
# table(CV$H_ANN)
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

# test <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_nonDDC_mouse_mapping_Hchr11.csv")
# colnames(test) <- c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")
# table(CV$H_ANN)
# table(test$H_ANN)
# 
# test$H_ANN <- sort_cat(test$H_ANN)
# CV$H_ANN <- sort_cat(CV$H_ANN)
# 
# # x <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_nonDDC_mouse_mapping_Hchr11.csv")
# # table(x$V5)

###################
### STACK OVERFLOW
##################


# set.seed(1)
# strings <- c(do.call(paste0, replicate(4, sample(LETTERS, 10000, TRUE), FALSE)),
#              do.call(paste0, replicate(3, sample(LETTERS, 10000, TRUE), FALSE)),
#              do.call(paste0, replicate(2, sample(LETTERS, 10000, TRUE), FALSE)))
# tableX <- as.data.table(table(strings))
# 
# 
# sort_cat <- function(strings){
#   tmp <- strsplit(strings, split="")
#   tmp <- lapply(tmp, sort)
#   tmp <- lapply(tmp, paste0, collapse = "")
#   tmp <- unlist(tmp)
#   return(tmp)
# }
# test1 <- sort_cat(strings)
# table1 <- as.data.table(table(sort_cat(strings)))
# 
# library(stringi)
# sort_stringi = function(s) {
#   s = stri_split_boundaries(s, type = "character")
#   s = lapply(s, stri_sort)
#   s = lapply(s, stri_join, collapse = "")
#   unlist(s)
# }
# test2 <- sort_stringi(strings)
# table2 <- as.data.table(table(sort_stringi(strings)))
# 
# sort_cat3 <- function(strings){
#   apply(stri_extract_all_regex(strings, "\\p{L}", simplify = TRUE), 1, function(i){
#     stri_join(stri_sort(i), collapse = "")
#   })
# }
# test3 <- sort_cat3(strings)
# table3 <- as.data.table(table(sort_cat3(strings)))
# 
# 
# library(microbenchmark)
# microbenchmark::microbenchmark(
#        one = sort_cat(strings),
#        two = sort_stringi(strings),
#        three = sort_cat3(strings),
#        times = 10)









