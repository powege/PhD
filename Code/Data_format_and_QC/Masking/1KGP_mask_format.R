rm(list=ls())
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
StrictMask.file <- args[2]
out.file <- args[3]

### IMPORT
# tmp <- fread(in.file)
fasta <- fread(StrictMask.file, skip=1, header=F)
# REF <- fread(REF.sm.file)

### FORMAT 

# subset POS and REP_MASK
# REF <- REF[,c("POS", "REP_MASK")]

# Format StrictMask fasta 
StrictMask <- paste(fasta$V1, collapse = '')
StrictMask <- strsplit(StrictMask, "")[[1]]
POS <- seq(from=1, to=length(StrictMask))
out <- data.table(POS = POS, 
                  StrictMask = StrictMask)

# merge 
# out <- out[REF, on="POS"]

# MASK bases StrictMask != "P" | REP_MASK == 1
out$Mask <- 1
out$Mask[out$StrictMask == "P"] <- 0
# out$Mask[out$REP_MASK == 1] <- 1

out$StrictMask <- NULL
# out$REP_MASK <- NULL

### Output
fwrite(out, out.file)