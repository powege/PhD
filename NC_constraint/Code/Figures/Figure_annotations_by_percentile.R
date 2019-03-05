rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

gene.pos <- fread("/well/lindgren/George/Data/Ensembl/BioMart/M_exon_POS.csv")
# gene.pos <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/M_exon_POS.csv")


CDS_out <- list()
chr <- c(1:19, "X")
for(i in 1:length(chr)){
  sub <- subset(gene.pos, gene.pos$chromosome_name == chr[i])
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  all.POS <- unlist(seq2(from = sub$exon_chrom_start, to = sub$exon_chrom_end))
  all.CHR <- rep(chr[i], length(all.POS))
  CDS_out[[i]] <- paste(all.CHR, all.POS, sep = "_")
  print(i)
}
CDS_out <- unlist(CDS_out)

prop <- rep(NA, 100)
for (i in c(1:100)){
  # p <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", i, ".csv"), header = F)
  prop[i] <- table(p$V1 %in% CDS_out)[2]/length(p$V1)
  print(i)
}

percentile <- c(1:100)
df <- data.frame(Percentile = percentile, CDS_fraction = prop)
fwrite(df, "/well/lindgren/George/Data/NC_constraint/tmp_CDS_by_prcentile.csv")
