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

### HUMAN

### FOR LOOP CHR
M_SNV_0001_list <- list()
M_SNV_0005_list <- list()
M_SNV_001_list <- list()
chr <- c(1:22, "X")
for (i in 1:length(chr)){

  ### IMPORT 
  file_path <- paste0(input.path, "H_1000GP_QCed_VEP_all_canPC_chr", chr[i], ".txt")
  data <- fread(file_path)
  
  # subset columns
  data <- data[,1:19]
  
  # add column names
  colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC",
                      "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                      "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS",
                      "AF")
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0001) per Ensembl transcript ID
  sub_0001 <- subset(data, data$AF >= 0.0001)
  temp <- SNV_count(sub_0001)
  colnames(temp) <- c(colnames(temp)[1:4], paste0(colnames(temp)[5:9], "_0.0001"))
  M_SNV_0001_list[[i]] <- temp
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0005) per Ensembl transcript ID
  sub_0005 <- subset(data, data$AF >= 0.0005)
  temp <- SNV_count(sub_0005)
  colnames(temp) <- c(colnames(temp)[1:4], paste0(colnames(temp)[5:9], "_0.0005"))
  M_SNV_0005_list[[i]] <- temp
  
  ### Get nonfunctional, functional and intron SNV counts (MAF >= 0.0005) per Ensembl transcript ID
  sub_001 <- subset(data, data$AF >= 0.001)
  temp <- SNV_count(sub_001)
  colnames(temp) <- c(colnames(temp)[1:4], paste0(colnames(temp)[5:9], "_0.001"))
  M_SNV_001_list[[i]] <- temp
  
  print(i)
}

# rbind dataframes in list
M_SNV_0001 <- do.call("rbind", M_SNV_0001_list)
M_SNV_0005 <- do.call("rbind", M_SNV_0005_list)
M_SNV_001 <- do.call("rbind", M_SNV_001_list)

# merge dataframes
SNV_counts <- merge(M_SNV_0001, M_SNV_0005, all = T)
SNV_counts <- merge(SNV_counts, M_SNV_001, all = T)


# Remove all transcripts < 1 SNV
output <- subset(SNV_counts, SNV_counts$n_synonymous_0.0001 > 0 | 
                     SNV_counts$n_missense_0.0005 > 0 |
                     SNV_counts$n_nonsense_0.0005 > 0 |
                     SNV_counts$n_spliceAD_0.0005 > 0 |
                     SNV_counts$n_intron_0.0005 > 0)
removed_tran <- subset(SNV_counts, SNV_counts$n_synonymous_0.0001 == 0 & 
                           SNV_counts$n_missense_0.0005 == 0 &
                           SNV_counts$n_nonsense_0.0005 == 0 &
                           SNV_counts$n_spliceAD_0.0005 == 0 &
                           SNV_counts$n_intron_0.0005 == 0)

# list genes with > 0 nonsense or spliceAD SNV
# nonsense_spliceAD_genes <- subset(output, output$n_nonsense_0.0005 > 1 | output$n_spliceAD_0.0005 > 1)

### OUTPUT
fwrite(output, paste0(output.path, "1000GP_canPC_SNV_counts.txt"))
# fwrite(M_nonsense_spliceAD_genes, paste0(output.path, "MGP_canPC_SNV_nonsense_spliceAD_counts.txt"))



