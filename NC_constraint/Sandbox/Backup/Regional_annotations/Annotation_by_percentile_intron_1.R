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

# import data
dt <- fread("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/intron_POS_short.csv")
# dt <- fread(paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/", files[X]))

# get vector of all annotation POS ID 
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  sub <- subset(dt, dt$chromosome == X)
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(X, length(all.POS))
  out <- paste0(all.CHR, "_", all.POS)
  out <- as.data.table(unique(out))
  fwrite(
         out, 
         "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/intron_POS_ID_long.csv",
         append = T,
         col.names = F
         )


# # get proportion of annotation POS in each percentile
# prop <- rep(NA, 100)
# for (i in c(1:100)){
#   # p <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
#   p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
#   prop[i] <- table(p$V1 %in% out)[2]/length(p$V1)
#   print(i)
# }
# 
# percentile <- c(1:100)
# df <- data.frame(Percentile = percentile, fraction = prop)
# fwrite(df, paste0("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_", files[X]), col.names = F)
# # fwrite(out, paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/LONG_", files[X]))


