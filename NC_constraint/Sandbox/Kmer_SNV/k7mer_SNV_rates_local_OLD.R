# SCRIPT that estimates the relative probabilities of substitution for the middle base in 7-mers.

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)

### FUNCTIONS

### function that returns all possible kmers and middle base substitutions
k7mer_changes <- function(kmer){
  
  tmp <- do.call(CJ, replicate(n = kmer, expr = c("A", "C", "G", "T"), FALSE))
  from <- do.call(paste, c(tmp, sep=""))
  tmp$V4 <- "A"
  to.A <- do.call(paste, c(tmp, sep=""))
  tmp$V4 <- "C"
  to.C <- do.call(paste, c(tmp, sep=""))
  tmp$V4 <- "G"
  to.G <- do.call(paste, c(tmp, sep=""))
  tmp$V4 <- "T"
  to.T <- do.call(paste, c(tmp, sep=""))
  from <- rep(from, 4)
  to <- c(to.A, to.C, to.G, to.T)
  out <- data.table(from, to)
  out <- out[-which(out$from == out$to),]
  out <- out[order(out$from),]
  colnames(out) <- paste0("k", kmer, "_", colnames(out))
  
  return(out)
}

### FUNCTION that counts k7mers in sequence
count_7mer <- function(anc.SEQ, k7_all){
  
  b1 <- anc.SEQ[1:(length(anc.SEQ)-6)]
  b2 <- anc.SEQ[2:(length(anc.SEQ)-5)]
  b3 <- anc.SEQ[3:(length(anc.SEQ)-4)]
  b4 <- anc.SEQ[4:(length(anc.SEQ)-3)]
  b5 <- anc.SEQ[5:(length(anc.SEQ)-2)]
  b6 <- anc.SEQ[6:(length(anc.SEQ)-1)]
  b7 <- anc.SEQ[7:(length(anc.SEQ))]
  tri.total.temp <- paste0(b1, b2, b3, b4, b5, b6, b7)
  tri.total.temp <- tri.total.temp[!grepl("N", tri.total.temp)]
  tri.total.temp <- data.frame(table(tri.total.temp))
  colnames(tri.total.temp) <- c(colnames(k7_all), "count")
  tri.total <- merge(k7_all, tri.total.temp, all = T)
  tri.total$count[is.na(tri.total$count)] <- 0
  
  return(tri.total)
}


### FUNCTION that counts k7mer subsetitutions from ancestral to mutant sequence
count_7mer_sub <- function(anc.SEQ, mu.SEQ, k7_mu_all){
  
  # identify SNV POS
  SNV.POS <- which(anc.SEQ != mu.SEQ)
  
  # remove SNV at fist or last base
  rm.id <- c((1:3), ((length(anc.SEQ)-2):length(anc.SEQ)))
  SNV.POS <- SNV.POS[!SNV.POS %in% rm.id]
  
  if (length(SNV.POS) != 0){
    
    # identify ancestral kmer
    b1 <- anc.SEQ[SNV.POS - 3]
    b2 <- anc.SEQ[SNV.POS - 2]
    b3 <- anc.SEQ[SNV.POS - 1]
    b4 <- anc.SEQ[SNV.POS]
    b5 <- anc.SEQ[SNV.POS + 1]
    b6 <- anc.SEQ[SNV.POS + 2]
    b7 <- anc.SEQ[SNV.POS + 3]
    tri.from <- paste0(b1, b2, b3, b4, b5, b6, b7)
    
    # identify mutant kmer
    b1 <- anc.SEQ[SNV.POS - 3]
    b2 <- anc.SEQ[SNV.POS - 2]
    b3 <- anc.SEQ[SNV.POS - 1]
    b4 <- mu.SEQ[SNV.POS]
    b5 <- anc.SEQ[SNV.POS + 1]
    b6 <- anc.SEQ[SNV.POS + 2]
    b7 <- anc.SEQ[SNV.POS + 3]
    tri.to <- paste0(b1, b2, b3, b4, b5, b6, b7)
    
    # remove all indicies in tri.from and tri.to that contain chrs other than ACGT or nchar != 7
    tri.to <- tri.to[!grepl("N", tri.to)]
    tri.from <- tri.from[!grepl("N", tri.from)]
    
    mu.count.temp <- data.table(tri.from, tri.to)
    mu.count.temp <- as.data.table(table(mu.count.temp))
    colnames(mu.count.temp) <- c("k7_from", "k7_to", "count")
    mu.count.temp <- mu.count.temp[mu.count.temp$count != 0, ]
    mu.count <- mu.count.temp[k7_mu_all, on = c("k7_from", "k7_to")]
    mu.count$count[is.na(mu.count$count)] <- 0
    
    return(mu.count)
  } else {
    k7_mu_all$count <- 0
    return(k7_mu_all)
  }
}

### FUNCTION that calculates 7mer relative substitution rates
k7_mu_rr <- function(ref_files, vcf_files){
  
  # Initiate output lists
  k7_mu_count_list <- list()
  k7_total_list <- list()
  
  # Identify all possible 7-mer substitutions
  k7_mu_all <- k7mer_changes(7)
  k7_all <- data.table(k7_from = unique(k7_mu_all$k7_from))
  
  
  
  # count by chromosome
  for(i in 1:length(ref_files)){
    
    ### IMPORT reference genome and SNV POS
    ref <- fread(ref_files[i])
    vcf <- fread(vcf_files[i])
  
  ### FORMAT
  # merge
  vcf <- vcf[,c(2,4,5)]
  colnames(vcf) <- c("POS", "REF", "ALT")
  ref <- vcf[ref, on = c("POS", "REF")]
  rm(vcf)
  
  # generate ALT sequence
  ind <- which(ref$ALT == "A" | ref$ALT == "T" | ref$ALT == "C" | ref$ALT == "G")
  tmp <- ref$REF
  tmp[ind] <- ref$ALT[ind]
  ref$ALT <- tmp
  rm(tmp)
  
  #####################################
  # cut masked regions
  # ref <- subset(ref, ref$REP_MASK == 0)
  #####################################
  
  ### Split into vectors of sequential POS
  ANCESTRAL <- split(ref$REF, cumsum(c(TRUE, diff(ref$POS)!=1)))
  MUTANT <- split(ref$ALT, cumsum(c(TRUE, diff(ref$POS)!=1)))
  rm(ref)
  
  ### Remove vectors with bases < 7
  ANCESTRAL <- ANCESTRAL[lapply(ANCESTRAL, length) >= 7]
  MUTANT <- MUTANT[lapply(MUTANT, length) >= 7]
  
  ### Count 7-mers in ancestral sequence
  kmer7_totals <- lapply(ANCESTRAL, count_7mer, k7_all = k7_all)
  # sum all counts in list
  k7_total_list[[i]] <- rowSums(sapply(kmer7_totals, `[[`, 2), na.rm = TRUE)
  rm(kmer7_totals)
  
  ### Count the 7-mer substitutions with mutant sequence
  kmer7_sub_totals <- list()
  for (j in 1:length(ANCESTRAL)){
    kmer7_sub_totals[[j]] <- count_7mer_sub(anc.SEQ = ANCESTRAL[[j]],
                                            mu.SEQ = MUTANT[[j]],
                                            k7_mu_all = k7_mu_all)
    # print(j)
  }
  
  # sum all counts in list
  k7_mu_count_list[[i]] <- rowSums(sapply(kmer7_sub_totals, `[[`, 3), na.rm = TRUE)
  
  print(paste(i, "done!"))
  }
  
  # Sum total kmer across chromosomes
  k7_total <- do.call(cbind, k7_total_list)
  k7_total <- apply(k7_total, 1, sum, na.rm = T)
  k7_total <- data.frame(k7_from = k7_all$k7_from,
                         k7_from_N = k7_total)
  
  # Sum total kmer changes across chromosomes
  k7_mu_total <- do.call(cbind, k7_mu_count_list)
  k7_mu_total <- apply(k7_mu_total, 1, sum, na.rm = T)
  k7_mu_total <- data.frame(k7_from = k7_mu_all$k7_from,
                            k7_to = k7_mu_all$k7_to,
                            k7_mu_N = k7_mu_total)
  
  ### Calculate relative mu rates
  
  k7_out <- merge(k7_mu_total, k7_total, all = T)
  
  p_any_snp_given_k7 <- function(sub){
    sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
  }
  tmp <- ddply(k7_out, "k7_from", p_any_snp_given_k7)
  colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
  k7_out <- merge(k7_out, tmp, all = T)
  
  k7_out$p_snp_given_k7 <- k7_out$k7_mu_N/k7_out$k7_from_N
  
  k7_out$k7_mu_rr <- k7_out$k7_mu_N/k7_out$k7_from_N
  
  return(k7_out)
}


### RUN

# MOUSE
ref_files <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr19.txt"
vcf_files <- '~/Dropbox/PhD/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr19.vcf'
# # k7_out <- k7_mu_rr(ref_files = ref_files, vcf_files = vcf_files)
# fwrite(k7_out, "~/Dropbox/MGP_7mer_SNV_rate.table")

### TEST



#### TO DO
# edit script to:
# get kmer for each POS 
# cut masked POS
# count kmers in unmasked POS 
# count kmer changes in unmasked POS

