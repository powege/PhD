### SCRIPT that totals the potential synonymous, missense and nonnsense substitutions for each transcript

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("Args required", call.=FALSE)
# }

### SET VARS
mu_table_file_path <- args[1]
seq_file_path <- args[2]
out_file_path <- args[3]
species <- args[4]
# mu_table_file_path <- "~/Dropbox/PhD/Data/Mu_rates/AA_mutation_table.csv"
# seq_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCpass.csv"
# out_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv"
# species <- "mouse"

### FUNCTIONS

### FUNCTION that splits string into vector of trinucleotides
## INPUT: character string
## OUTPUT: vector of trinucleotide character strings
three_split <- function(sequence){
  
  # method 1
  # system.time({
  # out1 <- substring(sequence, seq(1, nchar(sequence), 3), seq(3, nchar(sequence), 3))
  # })
  
  # method 2
  # system.time({
  vec <- strsplit(sequence, "")[[1]]
  out2 <- paste0(vec[c(TRUE, FALSE, FALSE)], vec[c(FALSE, TRUE, FALSE)], vec[c(FALSE, FALSE, TRUE)])
  # })
  
  return(out2)
}

### FUNCTION that calculates the probability of synonymous, 
### missense, and nonsense mutations, based on a codon mutation rate table
## INPUT: CT = Mutation table with all possible codon to codon point mutations 
##            "Codon_from","Codon_to","AA_from","AA_to","Mutation_type"
## sequence = genetic sequence string (char)
## OUTPUT: Named vector with probability of each mutation type (P_syn, P_mis, P_non)
# codon <- codons[[3]]
Moaning_Myrtle <- function(codon, CT){
  
  # table of potential base changes and the coding consequence
  # Codon_from; Codon_to, AA_from; AA_to
  df <- subset(CT, CT$Codon_from == codon)
  syn <- nrow(df[df$Mutation_type == "syn",])
  mis <- nrow(df[df$Mutation_type == "mis",])
  non <- nrow(df[df$Mutation_type == "non",])
  
  output <- c(syn, mis, non)
  names(output) <- c("p_syn", "p_mis", "p_non")
  
  return(output)
}

# gene.seq <- H_input[[1]]
alakazam <- function(gene.seq, CT){
  
  # split sequence into codons
  codons <- three_split(gene.seq)
  
  # cut first and last codons
  start_codon <- codons[1]
  stop_codon <- codons[length(codons)]
  codons <- codons[2:(length(codons)-1)]
  
  # apply Moaning Myrtle to codons (excluding first and last)
  p_list <- lapply(codons, Moaning_Myrtle, CT=CT)
  
  # calculate P stop-lost for last codons
  # add P start-lost for first codon
  
  # sum all probabilites in list for each mutation type
  out <- Reduce(`+`, p_list)
  
  return(out)
}


### IMPORT

# mutation coding table
CT <- fread(mu_table_file_path)

# transcript sequence
A_seq <- fread(seq_file_path)

# SET CHR
if (species == "human"){ CHR <- 1:22}
if (species == "mouse"){ CHR <- 1:19}

### LOOP BY CHR

output_list <- list()
for (chr in 1:length(CHR)){
  print(paste("CHROMOSOME", CHR[chr], "BEING PROCESSED", sep = " "))
  seq_sub <- subset(A_seq, A_seq$chromosome_name == CHR[chr])
  input <- as.list(seq_sub$coding)
  output <- lapply(input, alakazam, CT = CT) 
  df <- cbind(data.frame(chromosome_name = seq_sub$chromosome_name,
                           external_gene_name = seq_sub$external_gene_name,
                           ensembl_gene_id = seq_sub$ensembl_gene_id,
                           ensembl_transcript_id = seq_sub$ensembl_transcript_id,
                           cds_length = seq_sub$cds_length),
                as.data.frame(do.call(rbind, output)))
  output_list[[chr]] <- df
  print(paste("CHROMOSOME", CHR[chr], "DONE!", sep = " "))
}
all_out <- do.call("rbind", output_list)


### OUTPUT 
fwrite(all_out, out_file_path)


