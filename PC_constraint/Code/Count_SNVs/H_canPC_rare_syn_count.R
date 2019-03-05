# Script that counts synonymous variants for each canonical transcript in QCed vcf

rm(list = ls())
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

# Function that counts the non-functional SNVs, functional SNVs, and intron SNVs per Ensembl transcript
SNV_count <- function(data){
  
  ### Subset HGNC_name, Ensembl_gene_id, Ensembl_transcript_id by unique Ensembl_transcript_id 
  gene_detail <- data[, c("CHROM","Gene","Feature","SYMBOL")]
  gene_detail <- gene_detail[!duplicated(gene_detail$Feature),]
  colnames(gene_detail) <- c("CHROM", "ensembl_gene_id", "ensembl_transcript_id", "HGNC_symbol")
  
  ### Count synonymous SNVs
  synonymous_all <- ddply(data, "Feature", synonymous_count)
  colnames(synonymous_all) <- c("ensembl_transcript_id", "n_synonymous")
  
  ### combine outputs
  output <- gene_detail
  output <- merge(output, synonymous_all, all = T)
  
  return(output)
}


### SET PATHS
input.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/"
output.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/"

### HUMAN

### FOR LOOP CHR
H_SNV_0002_list <- list()
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
  sub_0002 <- subset(data, data$AF <= 0.0002)
  temp <- SNV_count(sub_0002)
  H_SNV_0002_list[[i]] <- temp

  print(i)
}

# rbind dataframes in list
SNV_counts <- do.call("rbind", H_SNV_0002_list)

# Remove all transcripts < 1 SNV
output <- subset(SNV_counts, SNV_counts$n_synonymous > 0)
removed_tran <- subset(SNV_counts, SNV_counts$n_synonymous == 0)

# list genes with > 0 nonsense or spliceAD SNV
# nonsense_spliceAD_genes <- subset(output, output$n_nonsense_0.0005 > 1 | output$n_spliceAD_0.0005 > 1)

### OUTPUT
fwrite(output, paste0(output.path, "1000GP_canPC_rare_syn_SNV_counts.txt"))
# fwrite(M_nonsense_spliceAD_genes, paste0(output.path, "MGP_canPC_SNV_nonsense_spliceAD_counts.txt"))

########
rm(list = ls())
graphics.off()

counts <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_rare_syn_SNV_counts.txt")
k3 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")
k5 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu.csv")

df <- data.table(ensembl_transcript_id = k3$ensembl_transcript_id,
                 cds_length = k3$cds_length,
                 k3_p_syn = k3$p_syn,
                 k5_p_syn = k5$p_syn)
df <- df[counts, on = "ensembl_transcript_id"]

plot(df$k3_p_syn, df$k5_p_syn)
cor.test(df$k3_p_syn, df$k5_p_syn)

cds_mod <- lm(df$n_synonymous ~ df$cds_length)
summary(cds_mod)

k3_mod <- lm(df$n_synonymous ~ df$k3_p_syn)
summary(k3_mod)

k5_mod <- lm(df$n_synonymous ~ df$k5_p_syn)
summary(k5_mod)

mod <- lm(df$n_synonymous ~ df$k3_p_syn + df$k5_p_syn)
summary(mod)
