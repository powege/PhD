### SCRIPT that formats ClinVar data 

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
raw.file <- "~/Dropbox/PhD/Data/ClinVar/raw/variant_summary.txt"
formatted.file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf"
ann.rank.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
ann.unrank.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv"


# Import ClinVar data
CV <- fread(raw.file)
CV <- CV[,c("Chromosome", "Start", "RS# (dbSNP)", "ReferenceAllele", "AlternateAllele",
            "Type","ClinicalSignificance", "ClinSigSimple", "Assembly")]
# Subset GRCh38
CV <- subset(CV, CV$Assembly == "GRCh38")
# Subset SNVs
CV <- subset(CV, CV$Type == "single nucleotide variant")
# Subset Pathogenic variants
# table(CV$ClinicalSignificance)
# table(CV$ClinSigSimple)
# CV <- subset(CV, CV$ClinSigSimple == 1)
CV <- subset(CV, CV$ClinicalSignificance %like% "Pathogenic" | 
                 CV$ClinicalSignificance %like% "Likely pathogenic")
# Subset chromosomes
CV <- CV[CV$Chromosome %in% c(1:22),]
# Remove duplicates
CV <- CV[!duplicated(CV),]
CV <- CV[,c("Chromosome", "Start", "RS# (dbSNP)", "ReferenceAllele", "AlternateAllele")]
colnames(CV) <- c("CHR", "POS", "dbSNP", "REF", "ALT")

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

# identify SNVs annotatoins
CV <- setDT(CV)[r.ann, CAT := CAT, on = .(CHR, POS >= START, POS <= END)]
unrank_test <- setDT(CV)[ur.ann, CAT := CAT, on = .(CHR, POS >= START, POS <= END)]

# export as vcf
fwrite(CV, formatted.file, sep = "\t", col.names = T)


