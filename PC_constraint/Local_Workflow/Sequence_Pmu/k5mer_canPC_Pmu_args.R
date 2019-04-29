#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args: 
# 1. mutation probabiliy table file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table")
# 2. mutation coding table file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/AA_mutation_table.csv")
# 3. transcript sequence file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_QCed.csv")
# 4. chromosome (eg 21)
# 5. output file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu_chr21.csv")

# SCRIPT that takes a list of CDS and returns the relative probabilities of synonymous, 
# missense, and nonsense (stop-gained) mutations, given 5-mer 
# substitution probabilities.

# rm(list=ls())
# graphics.off()

library(data.table)

### Test if there is five arguments: if not, return an error
if (length(args)!=5) {
  stop("Five arguments must be supplied", call.=FALSE)
} 


### FUNCTIONS

### FUNCTION that splits string into vector of 5-mers
## INPUT: character string
## OUTPUT: vector of 5-mer character strings
five_split <- function(sequence){
  
  vec <- strsplit(sequence, "")[[1]]
  out <- rep(NA, nchar(sequence)-4)
  
  for (i in 1:(nchar(sequence)-4)){
    out[i] <- paste0(vec[i], vec[i+1], vec[i+2], vec[i+3], vec[i+4])
  }
  
  return(out)
}

### FUNCTION that splits string into list of codons with adjasent two bases (ie base, base, codon, base, base)
## INPUT: character string (multiple of 3, including start and stop codons)
## OUTPUT: list of codons with adjasent bases (character string)
seven_split <- function(sequence){
  
  vec <- strsplit(sequence, "")[[1]]
  vec <- vec[2:(length(vec)-1)]
  
  out <- list()
  
  for (i in 1:(length(vec)-3)){
    if(i %% 3 == 0){
      out[[i/3]] <- paste(vec[(i-2):(i+4)], collapse = "")
    }
  }
  
  return(out)
}

### FUNCTION that calculates the sequence specific probability of synonymous, 
### missense, and nonsense mutations, based on a 5-mer mutation rate table
## INPUT: CT = Mutation table with all possible codon to codon point mutations 
##            "Codon_from","Codon_to","AA_from","AA_to","Mutation_type"
##        RT = Mutation rate table 
## sequence = genetic sequence string (char)
## OUTPUT: Named vector with probability of each mutation type (P_syn, P_mis, P_non)
Moaning_Myrtle <- function(sevenbase, CT, RT){
  
  codon <- substring(sevenbase, 3, 5)
  fivelets <- five_split(sevenbase)
  
  # table of potential base changes and the coding consequence
  # Codon_from; Codon_to, AA_from; AA_to
  df1 <- subset(CT, CT$Codon_from == codon)
  df1$seq <- sevenbase
  
  # Identify first, second, and third bases in codon
  B1 <- strsplit(codon, "")[[1]][1]
  B2 <- strsplit(codon, "")[[1]][2]
  B3 <- strsplit(codon, "")[[1]][3]
  
  # create a dataframe of potential codon changes for each 5-mer substitution. 
  df2 <- data.frame()
  for ( i in 1:length(fivelets)){
    sub <- subset(RT, RT$k5_from == fivelets[i])
    sub$codon_from <- codon
    if (i == 1){
      sub$codon_to <- substring(sub$k5_to, 3) # take last three bases of k5_to
    }
    if (i == 2){
      sub$codon_to <- substring(sub$k5_to, 2, 4) # take middle three bases of k5_to
    }
    if (i == 3){
      sub$codon_to <- substring(sub$k5_to, 1, 3) # take first three bases of k5_to 
    }
    df2 <- rbind(df2, sub)
  }
  
  # combine codon coding consequences with codon mutation probabilities
  df1 <- df1[order(df1$Codon_to),]
  df2 <- df2[order(df2$codon_to),]
  df <- cbind(df1, df2)
  
  # sum consequence probabilities
  for (i in 1:nrow(df)){
    syn <- sum(df$k5_mu_rr[df$Mutation_type == "syn"])
    mis <- sum(df$k5_mu_rr[df$Mutation_type == "mis"])
    non <- sum(df$k5_mu_rr[df$Mutation_type == "non"])
  }
  
  output <- c(syn, mis, non)
  names(output) <- c("p_syn", "p_mis", "p_non")
  
  return(output)
}


alakazam <- function(gene.seq, CT, RT){
  
  # split sequence into codons +- 1bp
  codons <- seven_split(gene.seq)
  
  # apply Moaning Myrtle to all codons +- 1bp
  p_list <- lapply(codons, Moaning_Myrtle, CT=CT, RT=RT)
  
  # sum all probabilites in list for each mutation type
  out <- Reduce(`+`, p_list)
  
  return(out)
}


### RUN

# mutation probabiliy table
RT <- fread(args[1])

# mutation coding table
CT <- fread(args[2])

# transcript sequence
seq <- fread(args[3])

# subset sequences by chromosome
seq_sub <- subset(seq, seq$chromosome_name == args[4])

input <- as.list(seq_sub$coding)
output <- lapply(input, alakazam, CT = CT, RT = RT)
df_out <- cbind(data.frame(chromosome_name = seq_sub$chromosome_name,
                           external_gene_name = seq_sub$external_gene_name,
                           ensembl_gene_id = seq_sub$ensembl_gene_id,
                           ensembl_transcript_id = seq_sub$ensembl_transcript_id,
                           cds_length = seq_sub$cds_length),
                as.data.frame(do.call(rbind, output)))

print(paste("CHROMOSOME", args[4], "DONE!", sep = " "))

fwrite(df_out, args[5])


# test
# args <- c("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table",
#           "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/AA_mutation_table.csv",
#           "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_QCed.csv",
#           21,
#           "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu_chr21.csv")
