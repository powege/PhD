### SCRIPT that formats ClinVar data for VEP

# Download ClinVar variants from:
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

rm(list = ls())
graphics.off()

library(data.table)

# Import ClinVar data
CV <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/ClinVar_variant_summary.txt")
CV <- CV[,c("Chromosome", "Start", "RS# (dbSNP)", "ReferenceAllele", "AlternateAllele",
            "Type","ClinicalSignificance", "ClinSigSimple")]
# Subset SNVs
CV <- subset(CV, CV$Type == "single nucleotide variant")
# Subset Pathogenic variants
CV <- subset(CV, CV$ClinicalSignificance == "Pathogenic" | CV$ClinicalSignificance == "Likely pathogenic")
# Subset chromosomes
CV <- CV[CV$Chromosome %in% c(1:22, "X"),]
# Remove duplicates
CV <- CV[!duplicated(CV),]
CV <- CV[,c("Chromosome", "Start", "RS# (dbSNP)", "ReferenceAllele", "AlternateAllele")]
# export as vcf
fwrite(CV, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/ClinVar_VEP_input.vcf", sep = "\t", col.names = F)


