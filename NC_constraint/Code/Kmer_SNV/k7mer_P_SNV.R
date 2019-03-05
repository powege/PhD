# SCRIPT that returns the relative probabilities of SNV at each POS in a sequence, 

rm(list=ls())
graphics.off()

library(data.table)

### FUNCTIONS

### FUNCTION that returnes the probability of SNV for each k7mer in a sequence. 
## INPUT: vector of bases
## OUTPUT: vector of 7-mer probabilities of SNV
slender_worm <- function(seq, P_SNV){
  
  kmer <- rep(NA, length(seq)-6)
  
  for (i in 1:(length(seq)-6)){
    kmer[i] <- paste0(seq[i], seq[i+1], seq[i+2], seq[i+3], seq[i+4], seq[i+5], seq[i+6])
  }
  
  # identify indices in P_SNV
  ind <- match(kmer, P_SNV$k7_from, nomatch = NA)
  pSNV <- c(rep(NA, 3), P_SNV$p_any_snp_given_k7[ind], rep(NA, 3))

  return(pSNV)
}

### IMPORT DATA

# mutation probabiliy table
P_SNV <- fread("~/Dropbox/BitBucket_repos/phd/NC_constraint/Data/SNV_rates/MGP_7mer_SNV_rate.table")

# sequence (try reading in only REF column)
# seq <- c(sample(c("A","C","G","T"),500000,rep=TRUE,prob=c(0.4,0.1,0.1,0.4)), NA, sample(c("A","C","G","T"),500000,rep=TRUE,prob=c(0.4,0.1,0.1,0.4)))
seq <- fread("~/Dropbox/BitBucket_repos/phd/NC_constraint/Data/Alignment/M_chr19_formatted.txt", select = "REF")

### FORMAT
P_SNV <- P_SNV[,c("k7_from", "p_any_snp_given_k7")]
P_SNV <- unique(P_SNV)
seq <- seq$REF

### RUN
system.time({ 
  out <- slender_worm(seq, P_SNV)
  })

### OUTPUT 