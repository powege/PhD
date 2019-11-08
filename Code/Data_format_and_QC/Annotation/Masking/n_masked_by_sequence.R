### SCRIPT that calculates the total of each sequence in a bed file that is masked
### (ie < 90% of individuals/strains with > 10X coverage in gnomAD and MGP)

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 


### SET ARGS
bed_file <- args[1]
mask_file <- args[2]
out_file <- args[3]
CHR <- args[4]

# bed_file <- "/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv" 
# mask_file <- "/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_StrictMask_chr1.txt"
# out_file <- "/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_1KGP_mask_chr1.csv"
# CHR <- "1"

# bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv"
# mask_file <- "~/Dropbox/PhD/Data/MGP/bam/MGP_fraction90_depth_less10X_chr19.txt"
# out_file <- ""
# CHR <- 19


### FUNCTIONS

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT
dt <- fread(bed_file)
mask_pos <- fread(mask_file)


### FORMAT

# remove exons with NA
dt <- subset(dt, !is.na(dt$V2) & !is.na(dt$V3))

# vector of mask POS
colnames(mask_pos) <- "V1"
mask_pos <- mask_pos$V1

# subset exons by CHR
dt_chr <- dt[V1 == CHR]
# dt_chr <- dt_chr[1:98,] # to test

# get vector of bed POS
seq_pos <- seq2(from = dt_chr$V2, to = dt_chr$V3)

# get total of CDS POS that are masked
DT <- data.table(repid=rep(1:length(seq_pos), sapply(seq_pos, length)),
                 val = unlist(seq_pos))
DT <- setDT(stack(setNames(seq_pos, 1:length(seq_pos))))
seq_total_masked <- DT[, x := +(values %in% mask_pos)][, sum(x), keyby=.(ind)]$V1
# seq_total_masked2 <- sapply(seq_pos, function(x) sum(x %in% mask_pos))
# seq_total_masked3 <- rep(NA, length(seq_pos))
# for (i in 1:length(seq_pos)){
#   seq_total_masked2[i] <- length(which(seq_pos[[i]] %in% mask_pos))
#   print(i)
# }
dt_chr$seq_total_masked <- seq_total_masked

# merge data tables
out <- as.data.table(out)
out <- out[dt_chr, on = c("ensembl_transcript_id")]

# subset output columns
out <- out[,c("chromosome_name", 
              "external_gene_name", 
              "ensembl_gene_id", 
              "ensembl_transcript_id", 
              "seq_length", 
              "seq_fraction_masked")]

# remove duplicated rows
out <- out[!duplicated(out),]


### EXPORT
fwrite(out, out_file, col.names = F)


#####

### STACK OVERFLOW

# set.seed(1)
# 
# vec_list <- replicate(100, sample(1:10000000, size=sample(1:10000, 100)), simplify=FALSE)
# vec <- sample(1:10000000, size=10000)
# 
# system.time({
# total_match <- rep(NA, length(vec_list))
# for (i in 1:length(vec_list)){
#   total_match[i] <- length(which(vec_list[[i]] %in% vec))
#   print(i)
# }
# })
# 
# system.time({
# total_match2 <- sapply(vec_list, function(x) sum(x %in% vec))
# })
# 
# system.time({
# DT <- data.table(repid=rep(1:length(vec_list), sapply(vec_list, length)), val=unlist(vec_list))
# DT <- setDT(stack(setNames(vec_list, 1:length(vec_list))))
# DT[, x := +(values %in% vec)][, sum(x), keyby=.(ind)]$V1
# })


