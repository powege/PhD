rm(list = ls())
graphics.off()

library(data.table)

# rmsk <- fread("~/Dropbox/PhD/Data/UCSC/M_rmsk.txt")
rmsk <- fread("/well/lindgren/George/Data/UCSC/rmsk.txt")

colnames(rmsk) <- c("bin",	"swScore",	"milliDiv",	"milliDel",	"milliIns",	"genoName",
                    "genoStart",	"genoEnd", "genoLeft",	"strand",	"repName",	"repClass",	"repFamily",	"repStart",
                    "repEnd",	"repLeft",	"id")
# table(rmsk$genoName)
rmsk <- rmsk[,c("genoName", "genoStart", "genoEnd", "repClass")]
chr <- paste0(rep("chr", 20), c(1:19, "X"))
rmsk <- rmsk[rmsk$genoName %in% chr,]
rmsk$genoName <- gsub("chr", "", rmsk$genoName)
# table(rmsk$genoName)
colnames(rmsk) <- c("chromosome", "start", "end", "repeat_class")
rmsk <- rmsk[complete.cases(rmsk),]

# get vector of all annotation POS ID 
chr <- c(1:19, "X")
for(i in 1:length(chr)){
  sub <- subset(rmsk, rmsk$chromosome == chr[i])
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(chr[i], length(all.POS))
  out <- paste0(all.CHR, "_", all.POS)
  out <- unique(out)
  out <- data.table(ID = out)
  fwrite(out, paste0("/well/lindgren/George/Data/Sandbox/tmp_rmsk_POS_", i, ".csv"), col.names = F)
}


# # get vector of all annotation POS ID 
# out <- list()
# chr <- c(1:19, "X")
# for(i in 1:length(chr)){
#   sub <- subset(rmsk, rmsk$chromosome == chr[i])
#   seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
#   all.POS <- unlist(seq2(from = sub$start, to = sub$end))
#   all.CHR <- rep(chr[i], length(all.POS))
#   out[[i]] <- paste0(all.CHR, "_", all.POS)
#   print(i)
# }
# out <- unlist(out)
# out <- unique(out)


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
# fwrite(df, "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_repeats.csv", col.names = F)
# # fwrite(out, paste0("~/Dropbox/PhD//Data/NC_constraint/Ensembl_annotation_POS/LONG_", files[X]))


