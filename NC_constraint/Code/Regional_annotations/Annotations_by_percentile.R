rm(list=ls())
graphics.off()

# library(dplyr)
# library(tidyr)
library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args <- 2

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

X <- as.integer(args[1])

# list data files to index
files <- c("PC_exon_POS_short.csv", "NC_exon_POS_short.csv", "UTR_POS_short.csv", "intron_POS_short.csv",
           "promoter_POS_short.csv", "promoter_flanking_POS_short.csv", "enhancer_POS_short.csv",
           "TF_binding_POS_short.csv", "open_chromatin_POS_short.csv")

# import data
dt <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/", files[X]))
# dt <- fread(paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/", files[X]))

# get vector of all annotation POS ID 
out <- list()
chr <- c(1:19, "X")
for(i in 1:length(chr)){
  sub <- subset(dt, dt$chromosome == chr[i])
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(chr[i], length(all.POS))
  out[[i]] <- paste(all.CHR, all.POS, sep = "_")
  print(i)
}
out <- unlist(out)
out <- unique(out)


# get proportion of annotation POS in each percentile
prop <- rep(NA, 100)
for (i in c(1:100)){
  # p <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  prop[i] <- table(p$V1 %in% out)[2]/length(p$V1)
  print(i)
}

percentile <- c(1:100)
df <- data.frame(Percentile = percentile, fraction = prop)
fwrite(df, paste0("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_", files[X]), col.names = F)
# fwrite(out, paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/LONG_", files[X]))


