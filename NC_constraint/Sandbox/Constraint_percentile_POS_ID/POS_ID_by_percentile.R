rm(list=ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

percentile <- args[1]
species <- args[2]

### FUNCTIONS 
get_ID_all <- function(sub){
  all.POS <- as.character(rep(sub$POS_from, each = 100) + 0:99)
  all.CHR <- rep(sub$CHR, each = 100)
  all.ID <- paste0(all.CHR, "_", all.POS)
  return(all.ID)
}

### IMPORT CONSTRAINT SCORES
cs <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Constraint/Constraint_by_window_", species, "_650_50.csv"))
### Subset percentile
sub <- subset(cs, cs$Constraint_percentile_CpG == percentile)
### Get POS_ID
ID <- as.data.table(get_ID_all(sub))
### EXPORT
fwrite(ID, paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", percentile, "_650_50.csv"), col.names = F)
