rm(list=ls())
graphics.off()

# library(dplyr)
# library(tidyr)
library(data.table)

# import
out <- fread("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/unannotated_POS_ID_long.csv", header = F)

# get proportion of annotation POS in each percentile
prop <- rep(NA, 100)
for (i in 1:100){
  # p <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  prop[i] <- table(p$V1 %in% out$V1)[2]/length(p$V1)
  print(i)
}

percentile <- c(1:100)
df <- data.frame(Percentile = percentile, fraction = prop)
fwrite(df, paste0("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_unannotated_POS_short.csv"), col.names = F)
# fwrite(out, paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/LONG_", files[X]))


