### SCRIPT that formats GWAS catalog data 
# https://www.ebi.ac.uk/gwas/docs/file-downloads
# https://www.ebi.ac.uk/gwas/docs/file-downloads

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
gwas.file <- "~/Dropbox/PhD/Data/GWAS/raw/gwas_catalog_v1.0.2-associations_e96_r2019-09-24.tsv"
index.file <- "~/Dropbox/PhD/Data/GWAS/raw/gwas_catalog_trait-mappings_r2019-09-24.tsv"
out.ddc.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_DDC_QCed.txt"
out.other.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_nonDDC_QCed.txt"


### Import GWAS data
CV <- fread(gwas.file)
index <- fread(index.file)

### Format
CV <- CV[,c("CHR_ID", "CHR_POS", "DISEASE/TRAIT", "P-VALUE")]
colnames(CV) <- c("CHR", "POS", "DISEASE_TRAIT", "P_VALUE")
CV <- CV[CV$CHR %in% c(1:22),]
index <- index[,c("Disease trait", "Parent term")]
colnames(index) <- c("DISEASE_TRAIT", "Parent_term")
dt <- CV[index, on = "DISEASE_TRAIT", allow.cartesian=TRUE]
dt <- dt[complete.cases(dt),]
dt <- dt[dt$P_VALUE <= 5e-08]

# Subset traits with top level ontology of disease, disorder, or cancer. 
top_level_terms <- as.data.table(table(dt$Parent_term))
disease_terms <- c("Cardiovascular disease", 
                   "Neurological disorder", 
                   "Immune system disorder",
                   "Cancer",
                   "Digestive system disorder",
                   "Metabolic disorder",
                   "Other disease")
dt_ddc <- dt[dt$Parent_term %in% disease_terms,]
dt_other <- dt[!dt$Parent_term %in% disease_terms,]

# Remove duplicates
dt_ddc <- dt_ddc[,c("CHR", "POS", "Parent_term")]
dt_ddc <- dt_ddc[!duplicated(dt_ddc),]
dt_other <- dt_other[,c("CHR", "POS", "Parent_term")]
dt_other <- dt_other[!duplicated(dt_other),]

# export as vcf
fwrite(dt_ddc, out.ddc.file, sep = "\t", col.names = T)
fwrite(dt_other, out.other.file, sep = "\t", col.names = T)

######

######################
### STACK OVERFLOW ###
######################




