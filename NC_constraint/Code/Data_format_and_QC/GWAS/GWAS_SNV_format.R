### SCRIPT that formats GWAS catalog data 

rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) {
  stop("Arguments required!", call.=FALSE)
} 

# set args variables
# raw.file <- args[1] # .gff file
# formatted.file <- args[2] # .gtf file
# ann.file <- args[3]
raw.file <- "~/Dropbox/PhD/Data/GWAS/raw/gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv"
formatted.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_formatted.txt"
ann.rank.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
ann.unrank.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv"

# Import GWAS data
CV <- fread(raw.file)
CV <- CV[,c("CHR_ID", "CHR_POS", "SNPS", "DISEASE/TRAIT")]
colnames(CV) <- c("CHR", "POS", "dbSNP", "DISEASE/TRAIT")
# Subset chromosomes
CV <- CV[CV$CHR %in% c(1:22),]
# Remove duplicates
CV <- CV[!duplicated(CV),]

# import annotation data
r.ann <- fread(ann.rank.file)
ur.ann <- fread(ann.unrank.file)

# combine
tmp <- subset(r.ann, r.ann$category == "Intron" | r.ann$category == "Unannotated")
ur.ann <- subset(ur.ann, ur.ann$category != "Intron")
ur.ann <- rbind(ur.ann, tmp)
rm(tmp)

colnames(r.ann) <- c("CAT", "CHR", "START", "END")
colnames(ur.ann) <- c("CAT", "CHR", "START", "END")
CV$CHR <- as.integer(CV$CHR)
CV$POS <- as.integer(CV$POS)


# identify SNVs annotatoins
CV <- setDT(CV)[r.ann, CAT := CAT, on = .(CHR, POS >= START, POS <= END)]
# unrank_test <- setDT(CV)[ur.ann, CAT := CAT, on = .(CHR, POS >= START, POS <= END)]

# export as vcf
fwrite(CV, formatted.file, sep = "\t", col.names = T)

######

######################
### STACK OVERFLOW ###
######################

# set.seed(1)
# chr_vec <- c(sample(1:100000, 10), "12to145", "15:19", sample(1:100000, 10), "111.333", "567.1")
# int_vec <- chr_vec[c(1:10, 13:22)]
# num_vec <- chr_vec[c(1:10, 13:24)]




