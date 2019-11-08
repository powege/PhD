### SCRIPT that defines the bases under consideration for each chr
# returns for each chromosome: 
# all bases under consideration (short format)
# the number of each base in the reference

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
in_file <- args[1] 
counts_out_file <- args[2] 
short_out_file <- args[3] 
species <- args[4]

# in_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"
# out_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_length_summary.csv"
# species <- "mouse"

# in_file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr" 
# counts_out_file <- "/well/lindgren/George/Data/Ensembl/Reference/Human_REF_n_bp.csv" 
# short_out_file <- "/well/lindgren/George/Data/Ensembl/Reference/Human_REF_ATGC_POS.csv" 
# species <- "human"

if (species == "mouse"){chromosomes <- c(1:19)}
if (species == "human"){chromosomes <- c(1:22)}

### FUNCTIONS

# Function that converts long to short by ind
# sub <- subset(dt_long, dt_long$ID == "1427119318314353")
long_to_short <- function(sub){
  data.table(
    POS_START = sub$POS[1],
    POS_END = sub$POS[nrow(sub)]
  )
}

counts_list <- list()
short_list <- list()

for (chr in chromosomes){

### IMPORT

dt <- fread(paste0(in_file, chr, ".txt"))


### FORMAT

# count number of each base
counts <- as.data.table(table(dt$REF))
colnames(counts) <- c("Base", "N")
counts$CHR <- chr

# subset A, T, G, C bases
dt_long <- subset(dt, dt$REF %in% c("A", "T", "G", "C"))

# convert long to short
dt_long$ID <- cumsum(c(TRUE, abs(diff(dt_long$POS))!=1))
dt_short <- ddply(dt_long, "ID", long_to_short)
dt_short$ID <- NULL
dt_short$CHR <- chr

counts_list[[chr]] <- counts
short_list[[chr]] <- dt_short

print(chr)

}

counts_out <- do.call("rbind", counts_list)
short_out <- do.call("rbind", short_list)


### OUTPUT
fwrite(counts_out, counts_out_file)
fwrite(short_out, short_out_file)



