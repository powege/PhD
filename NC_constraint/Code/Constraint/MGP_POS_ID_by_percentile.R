rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

### FUNCTIONS 
get_ID_all <- function(sub){
  all.POS <- as.character(rep(sub$POS_from, each = 100) + 0:99)
  all.CHR <- rep(sub$CHR, each = 100)
  all.ID <- paste(all.CHR, all.POS, sep = "_")
  return(all.ID)
}

### IMPORT CONSTRAINT SCORES
# cs <- fread("~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_by_window.csv")
cs <- fread("/well/lindgren/George/Data/NC_constraint/MGP_constraint_by_window.csv")

sub <- subset(cs, cs$Constraint_percentile == args[1])
ID <- as.data.table(get_ID_all(sub))
# fwrite(ID, paste0("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", args[1], ".csv"))
fwrite(ID, paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/MGP_POS_ID_percentile_", args[1], ".csv"), col.names = F)


#####

