# SCRIPT that estimates the relative probabilities of substitution for the middle base in 5-mers.

rm(list = ls())
graphics.off()

library(data.table)
library(dplyr)

### FUNCTIONS

### function that returns all possible kmers and middle base substitutions
k5mer_changes <- function(kmer){
  
  tmp <- do.call(CJ, replicate(n = kmer, expr = c("A", "C", "G", "T"), FALSE))
  from <- do.call(paste, c(tmp, sep=""))
  tmp$V3 <- "A"
  to.A <- do.call(paste, c(tmp, sep=""))
  tmp$V3 <- "C"
  to.C <- do.call(paste, c(tmp, sep=""))
  tmp$V3 <- "G"
  to.G <- do.call(paste, c(tmp, sep=""))
  tmp$V3 <- "T"
  to.T <- do.call(paste, c(tmp, sep=""))
  from <- rep(from, 4)
  to <- c(to.A, to.C, to.G, to.T)
  out <- data.table(from, to)
  out <- out[-which(out$from == out$to),]
  out <- out[order(out$from),]
  colnames(out) <- paste0("k", kmer, "_", colnames(out))
  
  return(out)
}

### FUNCTION that counts k5mers in sequence
count_5mer <- function(anc.SEQ, k5_all){
  
  b1 <- anc.SEQ[1:(length(anc.SEQ)-4)]
  b2 <- anc.SEQ[2:(length(anc.SEQ)-3)]
  b3 <- anc.SEQ[3:(length(anc.SEQ)-2)]
  b4 <- anc.SEQ[4:(length(anc.SEQ)-1)]
  b5 <- anc.SEQ[5:(length(anc.SEQ))]
  tri.total.temp <- paste(b1, b2, b3, b4, b5, sep = "")
  tri.total.temp <- tri.total.temp[!grepl("-", tri.total.temp)]
  tri.total.temp <- data.frame(table(tri.total.temp))
  colnames(tri.total.temp) <- c(colnames(k5_all), "count")
  tri.total <- merge(k5_all, tri.total.temp, all = T)
  
  return(tri.total)
}


### FUNCTION that counts k5mer subsetitutions from ancestral to mutant sequence
count_5mer_sub <- function(anc.SEQ, mu.SEQ, k5_mu_all){
  
  # identify SNV POS
  SNV.POS <- which(anc.SEQ != mu.SEQ)
  
  # remove SNV at fist or last base
  rm.id <- c((1:2), ((length(anc.SEQ)-1):length(anc.SEQ)))
  SNV.POS <- SNV.POS[!SNV.POS %in% rm.id]
  
  if (length(SNV.POS) != 0){
    
    # # identify ancestral kmer
    # b1 <- anc.SEQ[SNV.POS - 2]
    # b2 <- anc.SEQ[SNV.POS - 1]
    # b3 <- anc.SEQ[SNV.POS]
    # b4 <- anc.SEQ[SNV.POS + 1]
    # b5 <- anc.SEQ[SNV.POS + 2]
    # tri.from <- paste(b1, b2, b3, b4, b5, sep = "")
    # 
    # # identify mutant kmer
    # b1 <- anc.SEQ[SNV.POS - 2]
    # b2 <- anc.SEQ[SNV.POS - 1]
    # b3 <- mu.SEQ[SNV.POS]
    # b4 <- anc.SEQ[SNV.POS + 1]
    # b5 <- anc.SEQ[SNV.POS + 2]
    # tri.to <- paste(b1, b2, b3, b4, b5, sep = "")
    
    # identify flanking bases
    b1 <- anc.SEQ[SNV.POS - 2]
    b2 <- anc.SEQ[SNV.POS - 1]
    b4 <- anc.SEQ[SNV.POS + 1]
    b5 <- anc.SEQ[SNV.POS + 2]
    
    b.anc <- anc.SEQ[SNV.POS]
    b.mu <- mu.SEQ[SNV.POS]
    
    tri.from <- paste(b1, b2, b.anc, b4, b5, sep = "")
    tri.to <- paste(b1, b2, b.mu, b4, b5, sep = "")
    
    # remove all indicies in tri.from and tri.to that contain chrs other than ACGT or nchar != 5
    # tri.to <- tri.to[!grepl("-", tri.to)]
    
    mu.count.temp <- data.frame(tri.from, tri.to)
    mu.count.temp <- data.frame(table(mu.count.temp))
    mu.count.temp <- mu.count.temp[mu.count.temp$Freq != 0, ]
    colnames(mu.count.temp) <- c("k5_from", "k5_to", "count")
    mu.count <- merge(k5_mu_all, mu.count.temp, all = T)
    
    return(mu.count)
  } else {
    k5_mu_all$count <- 0
    return(k5_mu_all)
  }
}

### FUNCTION that calculates 3mer relative substitution rates
k5_mu_rr <- function(chr, in_file_path){
  
  # Initiate output lists
  k5_mu_count_list <- list()
  k5_total_list <- list()
  
  # Identify all possible 3-mer substitutions
  k5_mu_all <- k5mer_changes(5)
  k5_all <- data.table(k5_from = unique(k5_mu_all$k5_from))
  
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
    
    ### Remove vectors with bases < 5
    ANCESTRAL <- ANCESTRAL[lapply(ANCESTRAL, length) >= 5]
    MUTANT <- MUTANT[lapply(MUTANT, length) >= 5]
    
    ### Count 5-mers in ancestral sequence
    kmer5_totals <- lapply(ANCESTRAL, count_5mer, k5_all = k5_all)
    # sum all counts in list
    k5_total_list[[i]] <- rowSums(sapply(kmer5_totals, `[[`, 2), na.rm = TRUE)
    
    ### Count the 3-mer substitutions with mutant sequence
    kmer5_sub_totals <- list()
    for (j in 1:length(ANCESTRAL)){
      kmer5_sub_totals[[j]] <- count_5mer_sub(anc.SEQ = ANCESTRAL[[j]], 
                                              mu.SEQ = MUTANT[[j]],
                                              k5_mu_all = k5_mu_all)
      # print(j)
    }
    # sum all counts in list
    k5_mu_count_list[[i]] <- rowSums(sapply(kmer5_sub_totals, `[[`, 3), na.rm = TRUE)
    
    print(paste("chr", chr[i], "done!"))
  }
  
  # Sum total kmer across chromosomes
  k5_total <- do.call(cbind, k5_total_list)
  k5_total <- apply(k5_total, 1, sum, na.rm = T)
  k5_total <- data.frame(k5_from = k5_all$k5_from,
                         k5_from_N = k5_total)
  
  # Sum total kmer changes across chromosomes
  k5_mu_total <- do.call(cbind, k5_mu_count_list)
  k5_mu_total <- apply(k5_mu_total, 1, sum, na.rm = T)
  k5_mu_total <- data.frame(k5_from = k5_mu_all$k5_from,
                            k5_to = k5_mu_all$k5_to,
                            k5_mu_N = k5_mu_total)
  
  ### Calculate relative mu rates
  
  k5_out <- merge(k5_mu_total, k5_total, all = T)
  
  k5_out$p_any_snp_given_k5 <- NA
  for (ii in 1:nrow(k5_out)){
    k5_out$p_any_snp_given_k5[ii] <- sum(k5_out$k5_mu_N[k5_out$k5_from == k5_out$k5_from[ii]])/k5_out$k5_from_N[ii]
  }
  
  k5_out$p_snp_given_k5 <- NA
  for (ii in 1:nrow(k5_out)){
    k5_out$p_snp_given_k5[ii] <- k5_out$k5_mu_N[ii]/sum(k5_out$k5_mu_N[k5_out$k5_from == k5_out$k5_from[ii]])
  }
  
  k5_out$k5_mu_rr <- NA
  for (ii in 1:nrow(k5_out)){
    k5_out$k5_mu_rr[ii] <- k5_out$k5_mu_N[ii]/k5_out$k5_from_N[ii]
  }
  print(paste("CHROMOSOME", chr[i], "DONE!", sep = " "))
  return(k5_out)
}


### RUN 

# MOUSE
chr <- c(1:19, "X")
in_file_path   <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ancestral_mutant/M_ancestral_mutant_SEQ_chr"
k5_out <- k5_mu_rr(chr, in_file_path)
write.csv(k5_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_5mer_mu_rate.table", row.names = F)

# HUMAN
chr <- c(1:22, "X")
in_file_path   <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ancestral_mutant/H_ancestral_mutant_SEQ_chr"
k5_out <- k5_mu_rr(chr, in_file_path)
write.csv(k5_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table", row.names = F)

