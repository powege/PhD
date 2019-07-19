# Script that counts synonymous, missense, and nonsense SNVs for each canonical transcript in QCed vcf

# SNV annotation categories -- 
# (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
# Non-functional:
# synonymous annotatiosn = "synonymous_variant", "stop_retained_variant","start_retained_variant"
# Functional: 
# missense annotations = "missense_variant"
# nonsnese annotations = "stop_gained", "start_lost", "stop_lost"

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### SET PATHS
M_consequence_file_path <- "~/Dropbox/PhD/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allMUSMUS_snps_QCed_VEP_v94_canPC_chr"
out_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_canPC_N_SNVs.csv"
summary_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_N_SNV_summary.csv"

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
  out <- sum(consequence.vec[names(consequence.vec) %like% "missense_variant" |
                               names(consequence.vec) %like% "protein_altering_variant"])
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

# Function that counts the non-functional SNVs, functional SNVs, and intron SNVs per Ensembl transcript
SNV_count <- function(data){
  
  ### Subset HGNC_name, Ensembl_gene_id, Ensembl_transcript_id by unique Ensembl_transcript_id 
  gene_detail <- data[, c("chromosome_name","Gene","Feature","SYMBOL")]
  gene_detail <- gene_detail[!duplicated(gene_detail$Feature),]
  colnames(gene_detail) <- c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name")
  
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
  
  ### combine outputs
  output <- gene_detail
  output <- merge(output, synonymous_all, all = T)
  output <- merge(output, missense_all, all = T)
  output <- merge(output, nonsense_all, all = T)
  output <- merge(output, spliceAD_all, all = T)
  
  return(output)
}

### IMPORT
m.con.list <- list()
for(chr in 1:19){
  m.con.list[[chr]] <- fread(paste0(M_consequence_file_path, chr, ".vcf"))
}
m.con <- do.call("rbind", m.con.list)


### FORMAT

# add column names
colnames(m.con) <- c("chromosome_name", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                      "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                      "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS")

### Get nonfunctional, functional and intron SNV counts (any MAF) per Ensembl transcript ID
M_SNV_counts <- SNV_count(m.con)

# Get total genes
M_total_ens_tran <- length(unique(M_SNV_counts$ensembl_transcript_id))
M_total_ens_gene <- length(unique(M_SNV_counts$ensembl_gene_id))
M_total_HGNC <- length(unique(M_SNV_counts$HGNC_symbol))

# Get total SNV counts
M_total_synonymous <- sum(M_SNV_counts$n_synonymous)
M_total_missense <- sum(M_SNV_counts$n_missense)
M_total_nonsense <- sum(M_SNV_counts$n_nonsense)
M_total_spliceAD <- sum(M_SNV_counts$n_spliceAD)

# Get meedian SNVs per gene
M_median_synonymous <- median(M_SNV_counts$n_synonymous)
M_median_missense <- median(M_SNV_counts$n_missense)
M_median_nonsense <- median(M_SNV_counts$n_nonsense)
M_median_spliceAD <- median(M_SNV_counts$n_spliceAD)


M_SNV_summary <- data.frame(n_genes = M_total_ens_gene,
                            n_synonymous = M_total_synonymous,
                            n_missense = M_total_missense,
                            n_nonsense = M_total_nonsense,
                            n_spliceAD = M_total_spliceAD)

# list genes with 0 nonsense or spliceAD SNV
M_nonsense_spliceAD_genes <- subset(M_SNV_counts, M_SNV_counts$n_nonsense == 0 & M_SNV_counts$n_spliceAD == 0)


### OUTPUT
fwrite(M_SNV_counts, out_file_path)
fwrite(M_SNV_summary, summary_file_path)





