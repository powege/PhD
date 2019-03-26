# SCRIPT that estimates the relative probabilities of substitution for the middle base in 7-mers.
# ARGS: REF file, vcf file, output file

#############
### Error: cannot allocate vector of size 1.5 Gb
#############

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Three argumants must be supplied: REF file path, vcf file path, out file path", call.=FALSE)
} 

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
k7_mu_rr <- function(ref_file, vcf_file){

  # Identify all possible 7-mer substitutions
  k7_mu_all <- k7mer_changes(7)
  k7_all <- data.table(k7_from = unique(k7_mu_all$k7_from))

  # count by chromosome

    ### IMPORT reference genome and SNV POS
    ref <- fread(ref_file)
    vcf <- fread(vcf_file)

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
    k7_total <- rowSums(sapply(kmer7_totals, `[[`, 2), na.rm = TRUE)
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
    k7_mu_count_list <- rowSums(sapply(kmer7_sub_totals, `[[`, 3), na.rm = TRUE)
    rm(kmer7_sub_totals)

    # format output
    k7_total <- data.table(k7_from = k7_all$k7_from,
                           k7_from_N = k7_total)
    k7_mu_total <- data.frame(k7_from = k7_mu_all$k7_from,
                              k7_to = k7_mu_all$k7_to,
                              k7_mu_N = k7_mu_count_list)
    return(k7_mu_total)
}


### RUN

out <- k7_mu_rr(ref_file = args[1], vcf_file = args[2])
fwrite(out, args[3])

### TEST

# out <- k7_mu_rr(ref_file = '~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr19.txt',
#                 vcf_file = '~/Dropbox/PhD/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr19.vcf')
# fwrite(out, args[3])
