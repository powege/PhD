### SCRIPT that formats Ensembl gene and multicell Regulatory Build annotation
### prioritises top ranking annotation for each nucleotide
### Outputs cols: category; chromosome; start; end

rm(list=ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
ann_file <- args[1] 
out_file <- args[2] # .csv file
chr <- args[3]
species <- args[4]
# ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr1.csv"
# out_file <- "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/Human_GENCODE_RegBuild_annotation_ranked.csv"
# chr <- as.integer(1)
# species <- "homo_sapien"


### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# GRCh38 chromosomes
H_chr <- c(1:22)
M_chr <- c(1:19)
if (species == "homo_sapien"){ 
  chromosomes <- H_chr }
if (species == "mus_musculus"){ 
  chromosomes <- M_chr }

### IMPORT 
unranked <-  fread(ann_file)

### FORMAT 
colnames(unranked) <- c("category", "chromosome", "start", "end")
unranked$category[unranked$category == "TF binding" | unranked$category == "Open chromatin"] <- "Miscellaneous"

### Set rank and output list
rank <- c(
  # "Exon - CDS", # remove for for loop 
  "Exon - UTR",  
  "Exon - other", 
  "Intron - proximal", 
  "Promoter", 
  "Enhancer - proximal", 
  "Enhancer - distal", 
  "CTCF binding",
  "Miscellaneous", 
  "Intron - distal", 
  "Unannotated"
)
out_list <- list()
ranked <- subset(unranked, unranked$category == "Exon - CDS")

### Loop through annotation categories
for(cat in rank){
  
  # subset category
  cat_sub <- subset(unranked, unranked$category == cat)
  
  # vector of POS in cat_sub, not in ranked df
  unique_pos <- setdiff(unlist(seq2(from = cat_sub$start, to = cat_sub$end)), 
                        unlist(seq2(from = ranked$start, to = ranked$end)))
  
  # format vector to dataframe
  unique_pos <- sort(unique_pos, decreasing = F)
  unique_pos <- t(sapply(split(unique_pos, findInterval(unique_pos, unique_pos[which(c(1, diff(unique_pos)) > 1)])), range))
  unique_pos <- as.data.table(unique_pos)
  colnames(unique_pos) <- c("start", "end")
  unique_pos$chromosome <- chr
  unique_pos$category <- cat
  
  # add non-duplicated pos to df
  ranked <- rbind(ranked, unique_pos) 
  print(paste0(cat, " added!"))
}

### OUTPUT
fwrite(ranked, out_file, col.names = F)


