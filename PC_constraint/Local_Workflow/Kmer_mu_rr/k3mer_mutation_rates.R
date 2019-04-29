# SCRIPT that estimates the relative probabilities of substitution for the middle base in 3-mers.

rm(list = ls())
graphics.off()

library(data.table)
library(dplyr)

### FUNCTIONS

### function that returns all possible kmers and middle base substitutions
kmer_changes <- function(kmer){
  
  tmp <- do.call(CJ, replicate(n = kmer, expr = c("A", "C", "G", "T"), FALSE))
  from <- do.call(paste, c(tmp, sep=""))
  tmp$V2 <- "A"
  to.A <- do.call(paste, c(tmp, sep=""))
  tmp$V2 <- "C"
  to.C <- do.call(paste, c(tmp, sep=""))
  tmp$V2 <- "G"
  to.G <- do.call(paste, c(tmp, sep=""))
  tmp$V2 <- "T"
  to.T <- do.call(paste, c(tmp, sep=""))
  from <- rep(from, 4)
  to <- c(to.A, to.C, to.G, to.T)
  out <- data.table(from, to)
  out <- out[-which(out$from == out$to),]
  out <- out[order(out$from),]
  colnames(out) <- paste0("k", kmer, "_", colnames(out))
  
  return(out)
}

### FUNCTION that counts k3mers in sequence
count_3mer <- function(anc.SEQ, k3_all){
  
  b1 <- anc.SEQ[1:(length(anc.SEQ)-2)]
  b2 <- anc.SEQ[2:(length(anc.SEQ)-1)]
  b3 <- anc.SEQ[3:(length(anc.SEQ))]
  tri.total.temp <- paste(b1, b2, b3, sep = "")
  tri.total.temp <- tri.total.temp[!grepl("-", tri.total.temp)]
  tri.total.temp <- data.frame(table(tri.total.temp))
  colnames(tri.total.temp) <- c(colnames(k3_all), "count")
  tri.total <- merge(k3_all, tri.total.temp, all = T)
  
  return(tri.total)
}


### FUNCTION that counts k3mer subsetitutions from ancestral to mutant sequence
count_3mer_sub <- function(anc.SEQ, mu.SEQ, k3_mu_all){
  
  # identify SNV POS
  SNV.POS <- which(anc.SEQ != mu.SEQ)
  
  # remove SNV at fist or last base
  rm.id <- c(1,length(anc.SEQ))
  SNV.POS <- SNV.POS[!SNV.POS %in% rm.id]
  
  if (length(SNV.POS) != 0){
  
  # identify ancestral kmer
  b1 <- anc.SEQ[SNV.POS - 1]
  b2 <- anc.SEQ[SNV.POS]
  b3 <- anc.SEQ[SNV.POS + 1]
  tri.from <- paste(b1, b2, b3, sep = "")
  
  # identify mutant kmer
  b1 <- anc.SEQ[SNV.POS - 1]
  b2 <- mu.SEQ[SNV.POS]
  b3 <- anc.SEQ[SNV.POS + 1]
  tri.to <- paste(b1, b2, b3, sep = "")
  
  # remove all indicies in tri.from and tri.to that contain chrs other than ACGT or nchar != 3
  # tri.to <- tri.to[!grepl("-", tri.to)]
  
  mu.count.temp <- data.frame(tri.from, tri.to)
  mu.count.temp <- data.frame(table(mu.count.temp))
  mu.count.temp <- mu.count.temp[mu.count.temp$Freq != 0, ]
  colnames(mu.count.temp) <- c("k3_from", "k3_to", "count")
  mu.count <- merge(k3_mu_all, mu.count.temp, all = T)
  
  return(mu.count)
  } else {
    k3_mu_all$count <- 0
    return(k3_mu_all)
  }
}

### FUNCTION that calculates 3mer relative substitution rates
k3_mu_rr <- function(chr, in_file_path){
  
# Initiate output lists
k3_mu_count_list <- list()
k3_total_list <- list()

# Identify all possible 3-mer substitutions
k3_mu_all <- kmer_changes(3)
k3_all <- data.table(k3_from = unique(k3_mu_all$k3_from))

# count by chromosome
for(i in 1:length(chr)){

  ### IMPORT
  in_seq_path <- paste0(in_file_path, chr[i], ".txt")
  df.seq <- fread(in_seq_path)
  # df.snp <- fread("/Volumes/Untitled/PC_constraint/Paper/Data/Ancestral_mutant/M_ancestral_mutant_SNP_chr19.txt")
  colnames(df.seq) <- c("POS", "ANCESTRAL", "MUTANT")
  # colnames(df.snp) <- c("POS", "ANCESTRAL", "MUTANT")
  
  ### Split into vectors of sequential POS
  ANCESTRAL <- split(df.seq$ANCESTRAL, cumsum(c(TRUE, diff(df.seq$POS)!=1)))
  MUTANT <- split(df.seq$MUTANT, cumsum(c(TRUE, diff(df.seq$POS)!=1)))
  
  ### Remove vectors with bases < 3
  ANCESTRAL <- ANCESTRAL[lapply(ANCESTRAL, length) >= 3]
  MUTANT <- MUTANT[lapply(MUTANT, length) >= 3]
  
  ### Count 3-mers in ancestral sequence
  kmer3_totals <- lapply(ANCESTRAL, count_3mer, k3_all = k3_all)
  # sum all counts in list
  k3_total_list[[i]] <- rowSums(sapply(kmer3_totals, `[[`, 2), na.rm = TRUE)
  
  ### Count the 3-mer substitutions with mutant sequence
  kmer3_sub_totals <- list()
  for (j in 1:length(ANCESTRAL)){
    kmer3_sub_totals[[j]] <- count_3mer_sub(anc.SEQ = ANCESTRAL[[j]], 
                                            mu.SEQ = MUTANT[[j]],
                                            k3_mu_all = k3_mu_all)
    # print(j)
    }
  # sum all counts in list
  k3_mu_count_list[[i]] <- rowSums(sapply(kmer3_sub_totals, `[[`, 3), na.rm = TRUE)
  
  print(paste("chr", chr[i], "done!"))
}

# Sum total kmer across chromosomes
k3_total <- do.call(cbind, k3_total_list)
k3_total <- apply(k3_total, 1, sum, na.rm = T)
k3_total <- data.frame(k3_from = k3_all$k3_from,
                       k3_from_N = k3_total)

# Sum total kmer changes across chromosomes
k3_mu_total <- do.call(cbind, k3_mu_count_list)
k3_mu_total <- apply(k3_mu_total, 1, sum, na.rm = T)
k3_mu_total <- data.frame(k3_from = k3_mu_all$k3_from,
                          k3_to = k3_mu_all$k3_to,
                          k3_mu_N = k3_mu_total)

### Calculate relative mu rates

k3_out <- merge(k3_mu_total, k3_total, all = T)

k3_out$p_any_snp_given_k3 <- NA
for (ii in 1:nrow(k3_out)){
  k3_out$p_any_snp_given_k3[ii] <- sum(k3_out$k3_mu_N[k3_out$k3_from == k3_out$k3_from[ii]])/k3_out$k3_from_N[ii]
}

k3_out$p_snp_given_k3 <- NA
for (ii in 1:nrow(k3_out)){
  k3_out$p_snp_given_k3[ii] <- k3_out$k3_mu_N[ii]/sum(k3_out$k3_mu_N[k3_out$k3_from == k3_out$k3_from[ii]])
}

k3_out$k3_mu_rr <- NA
for (ii in 1:nrow(k3_out)){
  k3_out$k3_mu_rr[ii] <- k3_out$k3_mu_N[ii]/k3_out$k3_from_N[ii]
}
print(paste("CHROMOSOME", chr[i], "DONE!", sep = " "))
return(k3_out)
}


### RUN MOUSE

# chr <- c(1:19, "X")
# in_file_path   <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ancestral_mutant/M_ancestral_mutant_SEQ_chr"
# k3_out <- k3_mu_rr(chr, in_file_path)
# write.csv(k3_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_3mer_mu_rate.table", row.names = F)


### RUN HUMAN

chr <- c(1:22, "X")
in_file_path <- "/Volumes/HarEx/PC_constraint/Paper/Data/Ancestral_mutant/H_ancestral_mutant_SEQ_chr"
k3_out <- k3_mu_rr(chr, in_file_path)
write.csv(k3_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_3mer_mu_rate.table", row.names = F)




