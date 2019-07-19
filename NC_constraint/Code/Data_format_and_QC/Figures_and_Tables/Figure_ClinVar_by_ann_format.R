rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

# H.clinvar.file <- args[1]
# H.ranked.file <- args[2]
# out.file <- args[3]
H.clinvar.file <- "~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf"
H.ranked.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
out.file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_ClinVar_by_ann_human.csv"

### IMPORT
H.ann <- fread(H.ranked.file)
H_CV <- fread(H.clinvar.file)


### FORMAT

annotation <- unique(H_CV$CAT)
total_bp <- rep(NA, length(annotation))
total_ClinVar <- rep(NA, length(annotation))
for (i in 1:length(annotation)){
sub1 <- subset(H.ann, H.ann$category == annotation[i])
total_bp[i] <- sum((sub1$end - sub1$start) + 1)
sub2 <- subset(H_CV, H_CV$CAT == annotation[i])
total_ClinVar[i] <- nrow(sub2)
}

out <- data.frame(CAT = annotation,
                  total_bp = total_bp,
                  total_ClinVar = total_ClinVar,
                  ClinVar_kb = (total_ClinVar/total_bp)*1000)
out <- subset(out, out$CAT != "")

### EXPORT
fwrite(out, out.file)
