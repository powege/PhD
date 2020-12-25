rm(list = ls())
graphics.off()

library(data.table)

### SCRIPT:

## INPUT:
# strain prefix file

## FORMAT
# for each strain: fread samtools_DP.txt.gz; subset DP < 10; add to dt
# count DP < 10 and merge with template (all chr and pos); NA to 0
# calculate % DP < 10

## OUTPUT
# chromosome; pos; n_less10X; p_less10X
# chromosome; start; end (for p_less10X > 0.1)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

raw.path <- args[1]
prefix.file <- args[2]
functions_file <- args[3]
chr <- as.character(args[4])
out.long.file <- args[5]
out.bed.file <- args[6]

# prefix.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Raw/wildSTRAINS.txt"
# functions_file <- "/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
# chr <- as.character(19)
# out.long.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Raw/MGP_wild_derived_n_DP_less10X_chr19.csv.gz"
# out.bed.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/MGP_wild_derived_DP_less10X_more0.1_chr19.bed"

### FUNCTIONS
source(functions_file)

### IMPORT
dt_prefix <- fread(prefix.file, header = F)

### FORMAT
dt_prefix <- dt_prefix$V1 # as vector
list_less10X <- list()
for (i in 1:length(dt_prefix)){ # for each prefix
  tmp <- fread(paste0("gunzip -cq ", raw.path,  dt_prefix[i], "_samtools_DP_chr", chr, ".txt.gz")) # import and unzip
  colnames(tmp) <- c("pos", "DP") 
  tmp <- tmp$pos[which(tmp$DP < 10)] # subset pos with DP < 10X
  list_less10X[[i]] <- tmp
  print(i)
}

dt_less10X <- as.data.table(tabulate(unlist(list_less10X))) # count occurences of each pos in chromosome
colnames(dt_less10X) <- "n_less10X"
dt_less10X$pos <- 1:nrow(dt_less10X) # tabulate counts occurences for all integers from 1 to max(vec)
dt_less10X <- dt_less10X[n_less10X != 0] # remove 0s
dt_less10X$chromosome <- chr
dt_less10X$p_less10X <- dt_less10X$n_less10X / length(dt_prefix) # calculate fraction < 10X
dt_more0.1 <- long.to.bed(dt_less10X[p_less10X > 0.1][,c("chromosome", "pos")]) # bed file for pos with fraction less than 10X > 0.1
dt_less10X <- dt_less10X[,c("pos", "n_less10X", "p_less10X")]
rm(tmp, list_less10X)

### OUTPUT
fwrite(dt_less10X, out.long.file, compress = "gzip")
fwrite(dt_more0.1, out.bed.file, compress = "gzip")



