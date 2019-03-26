rm(list=ls())
graphics.off()

library(data.table)

rmsk <- fread("~/Dropbox/PhD/Data/UCSC/M_rmsk.txt")
# rmsk <- fread("/well/lindgren/George/Data/UCSC/rmsk.txt")

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

table(rmsk$repeat_class)
rmsk$length <- rmsk$end - rmsk$start
rmsk_long <- subset(rmsk, rmsk$length > 150)
table(rmsk_long$repeat_class)
rmsk_long$length <- NULL

fwrite(rmsk_long, "~/Dropbox/PhD/Data/UCSC/M_rmsk_200bp_plus.txt")

# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# all.POS <- unlist(seq2(from = rmsk_long$start, to = rmsk_long$end))


