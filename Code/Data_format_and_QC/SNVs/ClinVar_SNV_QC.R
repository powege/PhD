### SCRIPT that QCs ClinVar SNVs

rm(list = ls())
graphics.off()

library(data.table)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
# if (length(args)==0) {
#   stop("Arguments required!", call.=FALSE)
# } 

# set args variables
# raw.file <- args[1] # .gff file
# formatted.file <- args[2] # .gtf file
# ann.file <- args[3]
raw.file <- "~/Dropbox/PhD/Data/ClinVar/raw/variant_summary.txt"
formatted.file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_snps_QCed.vcf"
# ann.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"

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

# export as vcf
fwrite(CV, formatted.file)


