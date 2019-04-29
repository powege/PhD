### SCRIPT that formats mask from positions defined by RepeatMasker and with MGP read depth < 10

rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Exactly three arguments must be supplied", call.=FALSE)
} 

# set args variables
REF.sm.file <- args[1]
RD.file <- args[2]
out.file <- args[3]

# fread sm REF
dt <- fread(REF.sm.file)

# fread RD
rd <- fread(RD.file)

# subset POS and REP_MASK
dt <- dt[,c("POS", "REP_MASK")]

# set colnames (POS, StrictMask)
colnames(dt) <- c("POS", "Mask")

# Mask POS with low coverage
rd <- rd$V2
dt$Mask[dt$POS %in% rd] <- 1

# fwrite 
fwrite(dt, out.file)




