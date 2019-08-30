### SCRIPT that calculates the fraction of each Ensembl canonical transcript that is masked
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
in_file <- args[1]
mask_file <- args[2]
out_file <- args[3]
CHR <- args[4]

# in_file <- "/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv" 
# mask_file <- "/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_StrictMask_chr1.txt"
# out_file <- "/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_1KGP_mask_chr1.csv"
# CHR <- "1"

# in_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv"
# mask_file <- "~/Dropbox/PhD/Data/MGP/bam/MGP_fraction90_depth_less10X_chr19.txt"
# out_file <- ""
# CHR <- 19


### FUNCTIONS

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT
dt <- fread(in_file)
mask_pos <- fread(mask_file)


### FORMAT 

# remove exons with no CDS POS 
dt <- subset(dt, !is.na(dt$genomic_coding_start) & !is.na(dt$genomic_coding_end))

# vector of mask POS
colnames(mask_pos) <- "V1"
mask_pos <- mask_pos$V1

# subset exons by CHR
dt_chr <- subset(dt, dt$chromosome_name == CHR)
# dt_chr <- dt_chr[1:98,] # to test

#get CDS length
# cds_length <- (dt_chr$genomic_coding_end - dt_chr$genomic_coding_start) + 1

# get vector of CDS POS
cds_pos <- seq2(from = dt_chr$genomic_coding_start, to = dt_chr$genomic_coding_end)

# get total of CDS POS that are masked
DT <- data.table(repid=rep(1:length(cds_pos), sapply(cds_pos, length)),
                 val = unlist(cds_pos))
DT <- setDT(stack(setNames(cds_pos, 1:length(cds_pos))))
cds_total_masked <- DT[, x := +(values %in% mask_pos)][, sum(x), keyby=.(ind)]$V1
# cds_total_masked2 <- sapply(cds_pos, function(x) sum(x %in% mask_pos))
# cds_total_masked3 <- rep(NA, length(cds_pos))
# for (i in 1:length(cds_pos)){
#   cds_total_masked2[i] <- length(which(cds_pos[[i]] %in% mask_pos))
#   print(i)
# }
dt_chr$cds_total_masked <- cds_total_masked

# calculate fraction masked by transcript
frac_masked <- function(sub){ sum(sub$cds_total_masked)/sub$cds_length[1] }
out <- ddply(dt_chr, "ensembl_transcript_id", frac_masked)
colnames(out) <- c("ensembl_transcript_id", "cds_fraction_masked")

# merge data tables
out <- as.data.table(out)
out <- out[dt_chr, on = c("ensembl_transcript_id")]

# subset output columns
out <- out[,c("chromosome_name", 
              "external_gene_name", 
              "ensembl_gene_id", 
              "ensembl_transcript_id", 
              "cds_length", 
              "cds_fraction_masked")]

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


