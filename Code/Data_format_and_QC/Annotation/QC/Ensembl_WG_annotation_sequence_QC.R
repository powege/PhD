### SCRIPT that formats and QCs genomic sequences

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("arguments must be supplied", call.=FALSE)
} 

### SET ARGS

ref_file <- args[1]
cat_file <- args[2]
out_file <- args[3]
chr <- as.integer(args[4])

# ref_file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr19.txt"
# cat_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- ""
# chr <- 19

interest <- c("Exon - UTR", # annotations of interest
              "Exon - other",
              "Promoter", 
              "Enhancer - proximal", 
              "Enhancer - distal",
              "CTCF binding")


### FUNCTIONS

# vectorised version of substr
substr2 <- Vectorize(substr, vectorize.args = c("start", "stop"))


### IMPPORT

ref_dt <- fread(ref_file)
cat <- fread(cat_file)


### FORMAT

colnames(cat) <- c("chromosome", "start", "end", "category", "strand") # colnames
cat <- cat[chromosome == chr] # subset chromosome
cat <- cat[category %in% interest] # subset annotations of interest

if (ref_dt$POS[1] == 1 & ref_dt$POS[nrow(ref_dt)] == nrow(ref_dt)){ # if reference genome continuous
ref <- paste0(ref_dt$REF, collapse = "") # collapse referecne to string
mask <- paste0(ref_dt$REP_MASK, collapse = "") # collapse mask to string
}

cat$sequence <- unlist(substr2(x = ref, start = cat$start, stop = cat$end)) # extract sequence
cat$length <- nchar(cat$sequence) # get sequence length
cat$n_mask <- nchar(gsub("0", "", unlist(substr2(x = mask, start = cat$start, stop = cat$end))))
cat$mask_frac <- cat$n_mask / cat$length # get fraction masked
cat$n_CG <- (str_count(cat$sequence, "CG") + str_count(cat$sequence, "GC"))
cat$CG_frac <- cat$n_CG / ((cat$length/2)*2) # get fraction CG dinucleotides
cat$sequence <- NULL

### QC

removed <- data.frame() # set dt for removed seq

# filter by length
rm.id <- which(cat$length < 100 | cat$length > 10000)
if (length(rm.id) != 0){
  removed <- rbind(removed, cat[rm.id,])
  cat <- cat[-rm.id,]
}

# filter by mask fraction
rm.id <- which(cat$mask_frac > 0.9)
if (length(rm.id) != 0){
  removed <- rbind(removed, cat[rm.id,])
  cat <- cat[-rm.id,]
}


### EXPORT
fwrite(cat, out_file, col.names = F)





