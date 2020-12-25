rm(list = ls())
graphics.off()

library(data.table)


### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

castaneus.file <- args[1]
domesticus.file <- args[2]
musculus.file <- args[3]
functions.file <- args[4]
chr <- as.character(args[5])
out.bed.file <- args[6]

# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"

### FUNCTIONS
source(functions.file)

### IMPORT
castaneus <- fread(paste0("gunzip -cq ", castaneus.file))
domesticus <- fread(paste0("gunzip -cq ", domesticus.file))
musculus <- fread(paste0("gunzip -cq ", musculus.file))

### FORMAT
out <- rbind(castaneus, domesticus, musculus)
out <- collapse.overlap(out)

### EXPORT
fwrite(out, out.bed.file, compress = "gzip")









