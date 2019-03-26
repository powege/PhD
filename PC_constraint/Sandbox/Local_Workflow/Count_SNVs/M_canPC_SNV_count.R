# Script that counts functional and nonfunctional variants for each canonical transcript in QCed vcf

# SNV annotation categories -- 
# (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
# Non-functional:
# synonymous annotatiosn = "synonymous_variant", "stop_retained_variant","start_retained_variant"
# Functional: 
# missense annotations = "missense_variant"
# nonsnese annotations = "stop_gained", "start_lost", "stop_lost"
# spliceAD annotations = "splice_donor_variant", "splice_acceptor_variant"

rm(list=ls())
graphics.off()

library(data.table)
library(plyr)

### FUNCTIONS

# Function that counts the total synonymous SNVs
synonymous_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "synonymous_variant" |
                               names(consequence.vec) %like% "stop_retained_variant" |
                               names(consequence.vec) %like% "start_retained_variant"])
  return(out)
}

# Function that counts the total missense SNVs
missense_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "missense_variant"])
  return(out)
}

# Function that counts the total nonsense SNVs
nonsense_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "stop_gained" |
                               names(consequence.vec) %like% "start_lost" |
                               names(consequence.vec) %like% "stop_lost"])
  return(out)
}

# Function that counts the total spliceAD SNVs
spliceAD_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "splice_donor_variant" |
                               names(consequence.vec) %like% "splice_acceptor_variant"])
  return(out)
}

# Function that counts the total intron SNVs
intron_count <- function(sub){
  consequence.vec <- summary(factor(sub$Consequence))
  out <- sum(consequence.vec[names(consequence.vec) %like% "intron_variant"])
  return(out)
}

# Function that counts the non-functional SNVs, functional SNVs, and intron SNVs per Ensembl transcript
SNV_count <- function(data){
  
  ### Subset HGNC_name, Ensembl_gene_id, Ensembl_transcript_id by unique Ensembl_transcript_id 
  gene_detail <- data[, c("CHROM","Gene","Feature","SYMBOL")]
  gene_detail <- gene_detail[!duplicated(gene_detail$Feature),]
  colnames(gene_detail) <- c("CHROM", "ensembl_gene_id", "ensembl_transcript_id", "HGNC_symbol")
  
  ### Count synonymous SNVs
  synonymous_all <- ddply(data, "Feature", synonymous_count)
  colnames(synonymous_all) <- c("ensembl_transcript_id", "n_synonymous")
  
  ### Count missense SNVs 
  missense_all <- ddply(data, "Feature", missense_count)
  colnames(missense_all) <- c("ensembl_transcript_id", "n_missense")
  
  ### Count nonsense SNVs 
  nonsense_all <- ddply(data, "Feature", nonsense_count)
  colnames(nonsense_all) <- c("ensembl_transcript_id", "n_nonsense")
  
  ### Count spliceAD SNVs 
  spliceAD_all <- ddply(data, "Feature", spliceAD_count)
  colnames(spliceAD_all) <- c("ensembl_transcript_id", "n_spliceAD")
  
  ### Count spliceAD SNVs 
  intron_all <- ddply(data, "Feature", intron_count)
  colnames(intron_all) <- c("ensembl_transcript_id", "n_intron")
  
  ### combine outputs
  output <- gene_detail
  output <- merge(output, synonymous_all, all = T)
  output <- merge(output, missense_all, all = T)
  output <- merge(output, nonsense_all, all = T)
  output <- merge(output, spliceAD_all, all = T)
  output <- merge(output, intron_all, all = T)
  
  return(output)
}



### SET PATHS
input.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/"
output.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/"

### RUN MOUSE

### IMPORT 
M_data <- fread(paste0(input.path, "M_MGP_QCed_VEP_canPC_all.txt"))
# add column names
colnames(M_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                    "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS")

### Get nonfunctional, functional and intron SNV counts (any MAF) per Ensembl transcript ID
M_SNV_counts <- SNV_count(M_data)

# Remove all transcripts < 1 SNV
M_output <- subset(M_SNV_counts, M_SNV_counts$n_synonymous > 0 | 
                       M_SNV_counts$n_missense > 0 |
                       M_SNV_counts$n_nonsense > 0 |
                       M_SNV_counts$n_spliceAD > 0 |
                        M_SNV_counts$n_intron > 0)
M_removed_tran <- subset(M_SNV_counts, M_SNV_counts$n_synonymous == 0 & 
                   M_SNV_counts$n_missense == 0 &
                   M_SNV_counts$n_nonsense == 0 &
                   M_SNV_counts$n_spliceAD == 0 &
                   M_SNV_counts$n_intron == 0)

# Get total genes
M_total_ens_tran <- length(unique(M_output$ensembl_transcript_id))
M_total_ens_gene <- length(unique(M_output$ensembl_gene_id))
M_total_HGNC <- length(unique(M_output$HGNC_symbol))

# Get total SNV counts
M_total_synonymous <- sum(M_output$n_synonymous)
M_total_missense <- sum(M_output$n_missense)
M_total_nonsense <- sum(M_output$n_nonsense)
M_total_spliceAD <- sum(M_output$n_spliceAD)
M_total_intron <- sum(M_output$n_intron)

# Get meedian SNVs per gene
M_median_synonymous <- median(M_output$n_synonymous)
M_median_missense <- median(M_output$n_missense)
M_median_nonsense <- median(M_output$n_nonsense)
M_median_spliceAD <- median(M_output$n_spliceAD)
M_median_intron <- median(M_output$n_intron)

M_SNV_summary <- data.frame(n_genes = M_total_ens_gene,
                            n_synonymous = M_total_synonymous,
                            n_missense = M_total_missense,
                            n_nonsense = M_total_nonsense,
                            n_spliceAD = M_total_spliceAD,
                            n_inton = M_total_intron)

# list genes with > 0 nonsense or spliceAD SNV
M_nonsense_spliceAD_genes <- subset(M_output, M_output$n_nonsense > 0 | M_output$n_spliceAD > 0)

### OUTPUT
fwrite(M_output, paste0(output.path, "MGP_canPC_SNV_counts.txt"))
fwrite(M_SNV_summary, paste0(output.path, "MGP_canPC_SNV_totals.txt"))
# fwrite(M_nonsense_spliceAD_genes, paste0(output.path, "MGP_canPC_SNV_nonsense_spliceAD_counts.txt"))



