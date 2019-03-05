# SCRIPT that takes a list of CDS and returns the relative probabilities of synonymous, 
# missense, and nonsense (stop-gained) mutations, given 3-mer 
# substitution probabilities.

rm(list=ls())
graphics.off()

library(data.table)

### FUNCTIONS

### FUNCTION that splits string into vector of trinucleotides
## INPUT: character string
## OUTPUT: vector of trinucleotide character strings
three_split <- function(sequence){
  
  vec <- strsplit(sequence, "")[[1]]
  out <- rep(NA, nchar(sequence)-2)
  
  for (i in 1:(nchar(sequence)-2)){
    out[i] <- paste0(vec[i], vec[i+1], vec[i+2])
  }
  
  return(out)
}

### FUNCTION that splits string into list of codons with adjasent bases (ie base, codon, base)
## INPUT: character string (multiple of 3, including start and stop codons)
## OUTPUT: list of codons with adjasent bases (character string)
five_split <- function(sequence){
  
  vec <- strsplit(sequence, "")[[1]]
  vec <- vec[3:(length(vec)-2)]
  
  out <- list()
  
  for (i in 1:length(vec)){
    if(i %% 3 == 0){
      out[[i/3]] <- paste(vec[(i-2):(i+2)], collapse = "")
    }
  }
  
  return(out)
}

### FUNCTION that calculates the sequence specific probability of synonymous, 
### missense, and nonsense mutations, based on a trinucleotide mutation rate table
## INPUT: CT = Mutation table with all possible codon to codon point mutations 
##            "Codon_from","Codon_to","AA_from","AA_to","Mutation_type"
##        RT = Mutation rate table 
## sequence = genetic sequence string (char)
## OUTPUT: Named vector with probability of each mutation type (P_syn, P_mis, P_non)
Moaning_Myrtle <- function(fivebase, CT, RT){
  
  codon <- substring(fivebase, 2, 4)
  triplets <- three_split(fivebase)
  
  # table of potential base changes and the coding consequence
  # Codon_from; Codon_to, AA_from; AA_to
  df1 <- subset(CT, CT$Codon_from == codon)
  df1$seq <- fivebase
  
  # Identify first, second, and third bases in codon
  B1 <- strsplit(codon, "")[[1]][1]
  B2 <- strsplit(codon, "")[[1]][2]
  B3 <- strsplit(codon, "")[[1]][3]
  
  # create a dataframe of sequence specific mutation probabilities
  df2 <- data.frame()
  for ( i in 1:length(triplets)){
    sub <- subset(RT, RT$k3_from == triplets[i])
    sub$k3_from <- codon
    if (i == 1){
      sub$k3_to <- substring(sub$k3_to, 2)
      sub$k3_to <- paste0(sub$k3_to, B3)
    }
    if (i == 3){
      sub$k3_to <- substr(sub$k3_to, 1, nchar(as.character(sub$k3_to))-1)
      sub$k3_to <- paste0(B1, sub$k3_to)
    }
    df2 <- rbind(df2, sub)
  }
  
  df1 <- df1[order(df1$Codon_to),]
  df2 <- df2[order(df2$k3_to),]
  df <- cbind(df1, df2)
  
  for (i in 1:nrow(df)){
    syn <- sum(df$k3_mu_rr[df$Mutation_type == "syn"])
    mis <- sum(df$k3_mu_rr[df$Mutation_type == "mis"])
    non <- sum(df$k3_mu_rr[df$Mutation_type == "non"])
  }
  
  output <- c(syn, mis, non)
  names(output) <- c("p_syn", "p_mis", "p_non")
  
  return(output)
}


alakazam <- function(gene.seq, CT, RT){
  
  # split sequence into codons +- 1bp
  codons <- five_split(gene.seq)
  
  # apply Moaning Myrtle to all codons +- 1bp
  p_list <- lapply(codons, Moaning_Myrtle, CT=CT, RT=RT)
  
  # sum all probabilites in list for each mutation type
  out <- Reduce(`+`, p_list)
  
  return(out)
}


### IMPORT DATA

# mutation probabiliy table
H_RT <- fread("/well/lindgren/George/Workflows/PC_constraint/Data/Mu_rates/H_3mer_mu_rate.table")
M_RT <- fread("/well/lindgren/George/Workflows/PC_constraint/Data/Mu_rates/M_3mer_mu_rate.table")

# mutation coding table
CT <- fread("/well/lindgren/George/Workflows/PC_constraint/Data/Mu_rates/AA_mutation_table.csv")

# transcript sequence
H_seq <- fread("/well/lindgren/George/Workflows/PC_constraint/Data/Ensembl/H_canPC_SEQ_QCed.csv")
M_seq <- fread("/well/lindgren/George/Workflows/PC_constraint/Data/Ensembl/M_canPC_SEQ_QCed.csv")

### LOOP BY CHR
H_CHR <- c(1:22, "X")
H_output_list <- list()
for (chr in 1:length(H_CHR)){
  print(paste("HUMAN CHROMOSOME", H_CHR[chr], "BEING PROCESSED", sep = " "))
  H_seq_sub <- subset(H_seq, H_seq$chromosome_name == H_CHR[chr])
  H_input <- as.list(H_seq_sub$coding)
  H_output <- lapply(H_input, alakazam, CT = CT, RT = H_RT) 
  H_df <- cbind(data.frame(chromosome_name = H_seq_sub$chromosome_name,
                           external_gene_name = H_seq_sub$external_gene_name,
                           ensembl_gene_id = H_seq_sub$ensembl_gene_id,
                           ensembl_transcript_id = H_seq_sub$ensembl_transcript_id,
                           cds_length = H_seq_sub$cds_length),
                as.data.frame(do.call(rbind, H_output)))
  H_output_list[[chr]] <- H_df
  print(paste("HUMAN CHROMOSOME", H_CHR[chr], "DONE!", sep = " "))
}
H_all_out <- do.call("rbind", H_output_list)
fwrite(H_all_out, "/well/lindgren/George/Workflows/PC_constraint/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")

M_CHR <- c(1:19, "X")
M_output_list <- list()
for (chr in 1:length(M_CHR)){
  print(paste("MOUSE CHROMOSOME", M_CHR[chr], "BEING PROCESSED", sep = " "))
  M_seq_sub <- subset(M_seq, M_seq$chromosome_name == M_CHR[chr])
  M_input <- as.list(M_seq_sub$coding)
  M_output <- lapply(M_input, alakazam, CT = CT, RT = M_RT)
  M_df <- cbind(data.frame(chromosome_name = M_seq_sub$chromosome_name,
                           external_gene_name = M_seq_sub$external_gene_name,
                           ensembl_gene_id = M_seq_sub$ensembl_gene_id,
                           ensembl_transcript_id = M_seq_sub$ensembl_transcript_id,
                           cds_length = M_seq_sub$cds_length),
                as.data.frame(do.call(rbind, M_output)))
  M_output_list[[chr]] <- M_df
  print(paste("MOUSE CHROMOSOME", M_CHR[chr], "DONE!", sep = " "))
}
M_all_out <- do.call("rbind", M_output_list)
fwrite(M_all_out, "/well/lindgren/George/Workflows/PC_constraint/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")

