rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

### Script that counts the sex-specific n. synonymous and the n. LoF variants for each transcript

### FUNCTIONS

# Function that counts the total synonymous variants
synonymous_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "synonymous_variant" |
                               names(consequence.vec) %like% "stop_retained_variant" |
                               names(consequence.vec) %like% "start_retained_variant"])
  return(out)
}

# Function that counts the total LoF variants
LoF_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "transcript_ablation" |
                               names(consequence.vec) %like% "splice_acceptor_variant" |
                               names(consequence.vec) %like% "splice_donor_variant"|
                               names(consequence.vec) %like% "stop_gained" |
                               names(consequence.vec) %like% "frameshift_variant"|
                               names(consequence.vec) %like% "stop_lost" |
                               names(consequence.vec) %like% "start_lost"|
                               names(consequence.vec) %like% "transcript_amplification"])
  return(out)
}


### IMPORT data
dt <- fread(paste0("/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/gnomAD_QCed_VEP_canPC_chr", args[1], ".txt"))
col_names <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AF", "Gene",
               "Feature", "Feature_type", "Consequence", "IMPACT", "SYMBOL", "SYMBOL_SOURCE",
               "BIOTYPE", "CANONICAL", "CCDS", "AF_male", "AF_female", "non_neuro_AF_male", 
               "non_neuro_AF_female", "controls_AF_male", "controls_AF_female", "InbreedingCoeff",
               "VQSLOD", "n_alt_alleles")
colnames(dt) <- col_names

# subset by MAF
dt <- subset(dt, dt$AF < 0.00004)
# sub <- subset(dt, dt$Feature == "ENST00000395811")

### get gene list
gene_detail <- dt[, c("CHROM","Gene","Feature","SYMBOL")]
gene_detail <- gene_detail[!duplicated(gene_detail$Feature),]
colnames(gene_detail) <- c("CHROM", "ensembl_gene_id", "ensembl_transcript_id", "HGNC_symbol")

### Count ALL synonymous SNVs
synonymous_all <- ddply(dt, "Feature", synonymous_count)
colnames(synonymous_all) <- c("ensembl_transcript_id", "n_synonymous_all")

### Count MALE synonymous SNVs
tmp <- subset(dt, dt$AF_male != 0)
synonymous_male <- ddply(tmp, "Feature", synonymous_count)
colnames(synonymous_male) <- c("ensembl_transcript_id", "n_synonymous_male")
rm(tmp)

### Count FEMALE synonymous SNVs
tmp <- subset(dt, dt$AF_female != 0)
synonymous_female <- ddply(tmp, "Feature", synonymous_count)
colnames(synonymous_female) <- c("ensembl_transcript_id", "n_synonymous_female")
rm(tmp)

### Count MALE synonymous SNVs
tmp <- subset(dt, dt$AF_male != 0)
lof_male <- ddply(tmp, "Feature", LoF_count)
colnames(lof_male) <- c("ensembl_transcript_id", "n_lof_male")
rm(tmp)

### Count FEMALE synonymous SNVs
tmp <- subset(dt, dt$AF_female != 0)
lof_female <- ddply(tmp, "Feature", LoF_count)
colnames(lof_female) <- c("ensembl_transcript_id", "n_lof_female")
rm(tmp)

### combine outputs
output <- gene_detail
output <- merge(output, synonymous_all, all = T)
output <- merge(output, synonymous_male, all = T)
output <- merge(output, synonymous_female, all = T)
output <- merge(output, lof_male, all = T)
output <- merge(output, lof_female, all = T)

# Remove all transcripts < 1 SNV
output <- subset(output, output$n_synonymous_all > 0 | 
                     output$n_synonymous_male > 0 |
                     output$n_synonymous_female > 0 |
                     output$n_lof_male > 0 |
                     output$n_lof_female > 0)

### EXPORT 
fwrite(output, "/well/lindgren/George/Data/SS_constraint/Variant_counts/SS_gnomAD_variant_counts.csv", append = T)





