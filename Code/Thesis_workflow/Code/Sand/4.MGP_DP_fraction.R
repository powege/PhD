### Script that identifies fraction of strains with converage > 10X for each bp
### INPUT


rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

# set args variables
in.file <- args[1]
out.file1 <- args[2]
out.file2 <- args[3]
CHR <- args[4]

# fread file
dt <- fread(in.file)

# colnames
# colnames(dt) <- c("CHR", "POS")

# tabulate
dt <- as.data.table(table(dt$V2))
colnames(dt) <- c("POS", "N")
dt <- subset(dt, dt$N != 0)

# calculate fraction
n_strain = 35
dt$fraction_10X <- (n_strain - dt$N)/n_strain

# add chr
# dt$CHR <- CHR

# subset columns
dt <- dt[, c("POS", "fraction_10X")]

# fwrite 
fwrite(dt, out.file1, col.names = F)

# subset CHR POS for fraction_10X < 0.9
dt <- subset(dt, dt$fraction_10X < 0.9)

# subset columns
dt <- dt[, c("POS")]

# fwrite 
fwrite(dt, out.file2, col.names = F)







######

# dt <- data.table(CHR = c(1,1,1,2,2,2,3,3,3),
#                  POS = c(1,2,3,1,2,3,1,2,3),
#                  S1 = sample(c(1:100), 9, replace=T),
#                  S2 = sample(c(1:100), 9, replace=T),
#                  S3 = sample(c(1:100), 9, replace=T),
#                  S4 = sample(c(1:100), 9, replace=T),
#                  S5 = sample(c(1:100), 9, replace=T),
#                  S6 = sample(c(1:100), 9, replace=T),
#                  S7 = sample(c(1:100), 9, replace=T)
#                  )
# dt1 <- dt[,1:2]
# dt2 <- as.matrix(dt[,3:9])
# frac_above_10 <- function(vec){ length(vec[vec>=10]) / length(vec) }
# dt1$frac_10X <- apply(dt2, 1, frac_above_10)

#######

# dt <- data.table(CHR = c(1,1,1,2,2,2,3,3,3),
#                  POS = c(1,1,6,1,1,3,1,1,6))
# x <- table(dt)
# x <- as.data.table(x)
# x <- subset(x, x$N != 0)
# 
# n_strain = 35
# x$fraction_10X <- (n_strain - x$N)/n_strain
# 
# x <- x[, c(1,2,4)]



