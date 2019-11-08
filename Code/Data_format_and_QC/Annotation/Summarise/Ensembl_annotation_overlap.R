# SCRIPT that totals ovelapping annotation categories

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
## Import reference and annotations
# calculate annotations by base (data.table)
# output table annotation overlap


### SET VARS
ref_file <- args[1]
ann_file <- args[2]
out_file <- args[3]
chr <- as.integer(args[4])

# ref_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr1.txt"
# ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_GWAS_null_sample_chr1.csv"
# chr <- as.integer("1")

### IMPORT 

H_ref <- fread(ref_file)
H_ann <- fread(ann_file)

###########
### FORMAT
#########

H_ref <- H_ref[,c("POS", "REF")]
colnames(H_ref) <- c("H_POS", "H_REF")
H_ref <- subset(H_ref, H_ref$H_REF %in% c("A", "T", "C", "G"))

colnames(H_ann) <- c("category", "chromosome", "start", "end")
# number code annotations (plyr)
H_ann$category <- mapvalues(H_ann$category, from=c("Exon - CDS",
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
                                                   "Unannotated"), 
                            to=c("A","B","C","D","E","F","G","H","I","J", "K", "L"))
# ensure start <= end 
tmp1 <- subset(H_ann, H_ann$end < H_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(H_ann, H_ann$end >= H_ann$start)
H_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)


### Split into subsets as chr 1 and 2 too large to format as one
bases <- c("A", "T", "C", "G")
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
    CAT = H_ann$category)
  
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
  H_ref <- h_dt1[H_ref, on=c("H_POS")]
}
rm(H_ann, H_ref_split, h_dt1, h_dt2)


### ANNOTATION TABLE
output <- as.data.table(table(H_ref$H_ANN))


### EXPORT
fwrite(output, out_file, col.names = F)



