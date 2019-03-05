# Script that counts the total functional and nonfunctional variants at differnet MAF thresholds 

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

# Function that returns the total SNVs (all transcripts) above a given MAF threshold
alakazam <- function(data, MAF_threshold){
  
  # Subset data by MAF cutoff
  data_maf <- subset(data, data$AF >= MAF_threshold)
  
  ### Count unique Ensembl transcripts
  transcript_total <- length(unique(data_maf$Feature))
  
  ### Count synonymous SNVs 
  synonymous_all <- ddply(data_maf, "Feature", synonymous_count)
  synonymous_total <- sum(synonymous_all$V1, na.rm = T)
  
  ### Count missense SNVs
  missense_all <- ddply(data_maf, "Feature", missense_count)
  missense_total <- sum(missense_all$V1, na.rm = T)
  
  ### Count nonsense SNVs
  nonsense_all <- ddply(data_maf, "Feature", nonsense_count)
  nonsense_total <- sum(nonsense_all$V1, na.rm = T)
  
  ### Count spliceAD SNVs
  spliceAD_all <- ddply(data_maf, "Feature", spliceAD_count)
  spliceAD_total <- sum(spliceAD_all$V1, na.rm = T)
  
  ### Count intron SNVs
  intron_all <- ddply(data_maf, "Feature", intron_count)
  intron_total <- sum(intron_all$V1, na.rm = T)
  
  ### combine outputs
  output <- c(MAF_threshold, transcript_total, synonymous_total, missense_total, nonsense_total, spliceAD_total, intron_total)
  names(output) <- c("MAF_threshold", "n_genes", "n_synonymous", "n_missense", "n_nonsense", "n_spliceAD", "n_intron")
  
  return(output)
}

### SET PATHS
data.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/"
out.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/"

### FOR LOOP CHR
MAF_counts_list <- list()
chr <- c(1:22, "X")
for (i in 1:length(chr)){
  
### IMPORT 
file_path <- paste0(data.path, "H_1000GP_QCed_VEP_all_canPC_chr", chr[i], ".txt")
data <- fread(file_path)

# subset columns
data <- data[,1:19]

# add column names
colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC",
                    "Gene", "Feature", "Feature_type", "Consequence", "IMPACT", 
                    "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS",
                    "AF")

### Calculate total SNVs at different MAF thresholds
# no MAF threshold
MAF_counts <- alakazam(data, 0)
MAF_counts[1] <- 0
maf_seq <- seq(from = 0.0001, to = 0.001, by = 0.0001)
for(j in maf_seq){
  temp <- alakazam(data, j)
  MAF_counts <- rbind(MAF_counts, temp)
}
colnames(MAF_counts) <- c("MAF_threshold", "n_genes", "n_synonymous", "n_missense", "n_nonsense", "n_spliceAD", "n_intron")
rownames(MAF_counts) <- NULL

MAF_counts_list[[i]] <- MAF_counts
}

# Sum all matrices in list:
output <- Reduce('+', MAF_counts_list)
output[,1] <- c(0, maf_seq)
output <- as.data.table(output)

### OUTPUT
out_file_path <- paste0(out.path, "1000GP_canPC_SNV_totals.txt")
fwrite(x = output, file = out_file_path)


#####
