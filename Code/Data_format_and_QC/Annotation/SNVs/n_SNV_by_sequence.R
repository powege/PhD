### SCRIPT 
# INPUT: dt of SNV sites; bed file of sequences
# OUTPUT: bed file with n SNV sites per sequence

rm(list = ls())
graphics.off()

### LIBRARYS

library(data.table)

### SET ARGS

snv_file <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_v94_controls_canPC_chr21.vcf"
bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_seqQC.csv"

### IMPORT

snv <- fread(snv_file)
bed <- fread(bed_file)

### FORMAT 

### EXPORT