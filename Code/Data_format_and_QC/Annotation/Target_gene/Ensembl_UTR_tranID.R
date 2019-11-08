### SCRIPT that 

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
gencode.gtf <- args[1] # .gtf file
out_file <- args[2] # .csv file
# gencode.gtf <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf"
# out_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Target_gene/Human_UTR_tranID.csv"


### IMPORT annotation raw data
gene <- fread(gencode.gtf)

### FORMAT gene annotation
colnames(gene) <- c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
# pull out transcript_biotype and transcript_id from attribute
gene$transcript_biotype <- str_match(string=gene$attribute, pattern="transcript_biotype\\s+\"([^\"]+)\"")[,2]
gene$transcript_id <- str_match(string=gene$attribute, pattern="transcript_id\\s+\"([^\"]+)\"")[,2]
# subset UTR for protein_coding biotype
utr5 <- subset(gene, gene$type == "five_prime_utr" & gene$transcript_biotype == "protein_coding")
utr3 <- subset(gene, gene$type == "three_prime_utr" & gene$transcript_biotype == "protein_coding")
# rbind 
utr <- rbind(utr5, utr3)
# add label
utr$category <- "Exon - UTR"
# subset coulmns
utr <- utr[,c("chromosome", "transcript_id", "transcript_biotype", "type", "category", "start", "end")]
# complete cases
utr <- utr[complete.cases(utr),]


### EXPORT 
fwrite(utr, out_file, col.names = F)