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
H_consequence_file_path <- "~/Dropbox/PhD/Data/1KGP/Variants/vcf_QCed_VEP/1KGP_phase3_snps_QCed_VEP_v94_canPC_chr"
out001_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF001.csv"
out0005_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF0005.csv"
out0001_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF0001.csv"
summary_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_N_SNV_summary.csv"

# H_consequence_file_path <- "~/Dropbox/PhD/Data/gnomAD/vcf_QCed_VEP/gnomAD_v2.1.1_GRC38_snps_QCed_VEP_v94_controls_canPC_chr"
# out001_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/gnomAD_v2.1.1_QCed_canPC_N_SNVs_MAF001.csv"
# out0005_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/gnomAD_v2.1.1_QCed_canPC_N_SNVs_MAF0005.csv"
# out0001_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/gnomAD_v2.1.1_QCed_canPC_N_SNVs_MAF0001.csv"
# summary_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/gnomAD_v2.1.1_QCed_N_SNV_summary.csv"


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


### FOR LOOP CHR
H_SNV_0001_list <- list()
H_SNV_0005_list <- list()
H_SNV_001_list <- list()
for (i in 1:22){
  
  ### IMPORT 
  data <- fread(paste0(H_consequence_file_path, i, ".vcf"))
  
  # subset columns 1:19
  data <- data[,1:19]
  
  # add column names
  colnames(data) <- c("chromosome_name", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC",
                      "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                      "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS",
                      "AF")
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0001) per Ensembl transcript ID
  sub_0001 <- subset(data, data$AF >= 0.0001)
  H_SNV_0001_list[[i]] <- SNV_count(sub_0001)
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0005) per Ensembl transcript ID
  sub_0005 <- subset(data, data$AF >= 0.0005)
  H_SNV_0005_list[[i]] <- SNV_count(sub_0005)
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0005) per Ensembl transcript ID
  sub_001 <- subset(data, data$AF >= 0.001)
  H_SNV_001_list[[i]] <- SNV_count(sub_001)
  
  print(i)
}

H_SNV_001 <- do.call("rbind", H_SNV_001_list)
H_SNV_0005 <- do.call("rbind", H_SNV_0005_list)
H_SNV_0001 <- do.call("rbind", H_SNV_0001_list)

# grep "chr" from chromosome_name
H_SNV_001$chromosome_name <- gsub("chr", "", H_SNV_001$chromosome_name)
H_SNV_0005$chromosome_name <- gsub("chr", "", H_SNV_0005$chromosome_name)
H_SNV_0001$chromosome_name <- gsub("chr", "", H_SNV_0001$chromosome_name)

# subset autosomes
H_SNV_001 <- H_SNV_001[which(H_SNV_001$chromosome_name %in% c(1:22))]
H_SNV_0005 <- H_SNV_0005[which(H_SNV_0005$chromosome_name %in% c(1:22))]
H_SNV_0001 <- H_SNV_0001[which(H_SNV_0001$chromosome_name %in% c(1:22))]


### SUMMARY 
H_SNV_summary <- data.frame(
                            variable = c("n_ensembl_transcripts", 
                                         "n_ensembl_genes", 
                                         "n_external_genes", 
                                         "n_synonymous",
                                         "n_missense",
                                         "n_nonsense",
                                         "n_spliceAD",
                                         "median_synonymous",
                                         "median_missense",
                                         "median_nonsense",
                                         "median_spliceAD"),
                            value_MAF001 = c(length(unique(H_SNV_001$ensembl_transcript_id)), 
                                             length(unique(H_SNV_001$ensembl_gene_id)), 
                                             length(unique(H_SNV_001$external_gene_name)), 
                                             sum(H_SNV_001$n_synonymous),
                                             sum(H_SNV_001$n_missense),
                                             sum(H_SNV_001$n_nonsense),
                                             sum(H_SNV_001$n_spliceAD),
                                             median(H_SNV_001$n_synonymous),
                                             median(H_SNV_001$n_missense),
                                             median(H_SNV_001$n_nonsense),
                                             median(H_SNV_001$n_spliceAD)),
                            value_MAF0005 = c(length(unique(H_SNV_0005$ensembl_transcript_id)), 
                                              length(unique(H_SNV_0005$ensembl_gene_id)), 
                                                     length(unique(H_SNV_0005$external_gene_name)), 
                                                     sum(H_SNV_0005$n_synonymous),
                                                     sum(H_SNV_0005$n_missense),
                                                     sum(H_SNV_0005$n_nonsense),
                                                     sum(H_SNV_0005$n_spliceAD),
                                                     median(H_SNV_0005$n_synonymous),
                                                     median(H_SNV_0005$n_missense),
                                                     median(H_SNV_0005$n_nonsense),
                                                     median(H_SNV_0005$n_spliceAD)),
                            value_MAF0001 = c(length(unique(H_SNV_0001$ensembl_transcript_id)), 
                                              length(unique(H_SNV_0001$ensembl_gene_id)), 
                                              length(unique(H_SNV_0001$external_gene_name)), 
                                              sum(H_SNV_0001$n_synonymous),
                                              sum(H_SNV_0001$n_missense),
                                              sum(H_SNV_0001$n_nonsense),
                                              sum(H_SNV_0001$n_spliceAD),
                                              median(H_SNV_0001$n_synonymous),
                                              median(H_SNV_0001$n_missense),
                                              median(H_SNV_0001$n_nonsense),
                                              median(H_SNV_0001$n_spliceAD))
                            )


# list genes with 0 nonsense or spliceAD SNV
H_nonsense_spliceAD_genes_001 <- subset(H_SNV_001, H_SNV_001$n_nonsense == 0 & H_SNV_001$n_spliceAD == 0)
H_nonsense_spliceAD_genes_0005 <- subset(H_SNV_0005, H_SNV_0005$n_nonsense == 0 & H_SNV_0005$n_spliceAD == 0)
H_nonsense_spliceAD_genes_0001 <- subset(H_SNV_0001, H_SNV_0001$n_nonsense == 0 & H_SNV_0001$n_spliceAD == 0)

### OUTPUT
fwrite(H_SNV_001, out001_file_path)
fwrite(H_SNV_0005, out0005_file_path)
fwrite(H_SNV_0001, out0001_file_path)
fwrite(H_SNV_summary, summary_file_path)





