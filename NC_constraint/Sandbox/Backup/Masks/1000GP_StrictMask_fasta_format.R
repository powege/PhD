rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Exactly two arguments must be supplied", call.=FALSE)
} 

# set args variables
in.file <- args[1]
out.file <- args[2]

### Import fasta 
# tmp <- fread(in.file)
fasta <- fread(in.file, skip=1, header=F)

### Format
StrictMask <- paste(fasta$V1, collapse = '')
StrictMask <- strsplit(StrictMask, "")[[1]]
POS <- seq(from=1, to=length(StrictMask))
out <- data.table(POS = POS, 
                  StrictMask = StrictMask)

### Output
fwrite(out, out.file)