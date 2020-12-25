rm(list = ls())
graphics.off()

library(data.table)


### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

classical.file <- args[1]
wild.file <- args[2]
functions.file <- args[3]
chr <- as.character(args[4])
out.bed.file <- args[5]

# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"

### FUNCTIONS
source(functions.file)

### IMPORT
classical <- fread(classical.file)
wild <- fread(wild.file)

### FORMAT
out <- rbind(classical, wild)
out <- collapse.overlap(out)

### EXPORT
fwrite(out, out.bed.file)









