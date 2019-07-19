rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

gwas.file <- args[1]
alignment.file <- args[2]
out.file <- args[3]
# gwas.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_formatted.txt"

cv <- fread(gwas.file)
out <- list()

for (i in 1:22){

### IMPORT
chr <- fread(paste0(alignment.file, i, ".txt"))
### FORMAT
colnames(chr) <- c("CHR", "POS", "REF", "M_CHR", "M_POS", "M_REF")
# subset pathogenic SNVs
sub_cv <- subset(cv, cv$CHR == chr$CHR[1])
sub_cv$CHR <- as.integer(sub_cv$CHR)
chr$CHR <- as.integer(chr$CHR)
out_sub <- chr[sub_cv, on = c("CHR", "POS")]
out_sub <- out_sub[complete.cases(out_sub),]
out[[i]] <- out_sub
}

out <- do.call("rbind", out)
fwrite(out, file = out.file)


