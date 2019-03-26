rm(list = ls())
graphics.off()

library(data.table)

dt1 <- fread("~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_by_window.csv")
dt2 <- fread("~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_variables_by_window.csv")

dt <- dt2[dt1, on = c("CHR", "POS_from", "POS_to")]

per1 <- subset(dt, dt$Constraint_percentile == 1)
per2 <- subset(dt, dt$Constraint_percentile == 2)
per3 <- subset(dt, dt$Constraint_percentile == 3)

table(per1$CHR)
table(per2$CHR)
table(per3$CHR)

hist(dt$Constraint_score)
hist(per1$Constraint_score)
hist(per2$Constraint_score)
tmp <- rbind(per1, per2, per3)
hist(tmp$Constraint_score)

table(per1$n_SNV)
table(per2$n_SNV)
hist(dt$n_SNV, breaks = 100)


hist(dt$p_SNV_given_kmers)
mean(dt$p_SNV_given_kmers)
hist(per1$p_SNV_given_kmers)
mean(per1$p_SNV_given_kmers)
hist(per2$p_SNV_given_kmers)
mean(per2$p_SNV_given_kmers)

hist(dt$Read_depth)
hist(per1$Read_depth)
hist(per2$Read_depth)
hist(per3$Read_depth)

tmp <-  subset(per1, per1$CHR == 15)
hist(tmp$POS_from)
hist(dt$Repeats)
hist(per1$Repeats)
hist(per2$Repeats)
hist(per3$Repeats)
summary(dt$Repeats)
summary(per1$Repeats)
summary(per2$Repeats)
summary(per3$Repeats)



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

head(rmsk)
table(rmsk$repeat_class)

rmsk$length <- rmsk$end - rmsk$start
hist(rmsk$length)
rmsk.sub <- subset(rmsk, rmsk$length > 150)

rmsk.1 <- subset(rmsk.sub, rmsk.sub$chromosome == 1)
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
rmsk.POS <- unlist(seq2(from = rmsk.1$start, to = rmsk.1$end))
out <- rep(NA, 10)
for(i in c(1:10)){
  per <- subset(dt, dt$Constraint_percentile == i)
  per <- subset(per, per$CHR == 1)
  per.POS <- as.vector(seq2(from = per$POS_from, to = per$POS_to))
  x <- summary(per.POS %in% rmsk.POS)
  out[i] <- as.integer(x[3][[1]]) / ( as.integer(x[2][[1]]) + as.integer(x[3][[1]]))
  print(i)
}
out



