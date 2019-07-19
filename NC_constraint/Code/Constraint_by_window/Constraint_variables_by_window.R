# SCRIPT that creates data table for each chromosome

rm(list = ls())
graphics.off()

library(data.table)
# library(evobiR) ##### NOT ON CLUSTER!!!
library(zoo)
library(MASS)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args
species <- args[1]
chromo <- as.integer(args[2])
window.size <- as.integer(args[3])
window.shift <- as.integer(args[4])
ref.file <- args[5]
vcf.file <- args[6]
mask.file <- args[7]
alignment.file <- args[8]
mu_table.file <- args[9]
out.file <- args[10]
# chromo <- 19
# ref.file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr19.txt"
# vcf.file <- "/well/lindgren/George/Data/1KGP/Variants/vcf_QCed_VEP/1000GP_phase3_QCed_VEP_v94_allPASS_chr19.txt"
# mask.file <- "/well/lindgren/George/Data/MGP/Masks/MGP_Mask_chr19.txt"
# mu_table.file <- "/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_StrictMask_chr19.txt"
# species <- "human"
# window.size <- 950
# window.shift <- 50
# alignment.file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt"


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

# FUNCTION that calculates the number of CpG dinucleotides in a sequence
n_CpG <- function(seq){
  string <- paste0(seq, collapse = "")
  CpG_count <- str_count(string, "CG")
  # to calculate proportion:
  # total.dinuc <- nchar(string)/2
  # CG_prop <- CG_count/total.dinuc
  return(CpG_count)
}

# FUNCTION that calculates the weighted sum of a window
# the middle numbers are weighted 1
# all other numbers are weighted lit a linear function between 0 and 1 (ie 0 outside the window)
sum_weight <- function(sequence, window.size, window.shift){

  outer.length <- (window.size - window.shift )/2
  
  weighted <- sum((sequence[1:outer.length]) * seq(from = 0, to = 1, length.out = outer.length)) +
                   sum((sequence[(outer.length + 1):(outer.length + window.shift)] * 1)) +
                   sum((sequence[(outer.length + window.shift + 1):(window.size)] * seq(from = 1, to = 0, length.out = outer.length))
  )
  return(weighted)
}

# FUNCTION that calculates the weighted number of CpG dinucleotides in a sequence
n_CpG_weight <- function(sequence, window.size, window.shift){
  
  CG.int <- unlist(str_locate_all(paste0(sequence, collapse = ""), "CG"))
  binary <- rep(0, length(sequence))
  binary[CG.int] <- 1
  
  outer.length <- (window.size - window.shift )/2
  
  weighted <- sum((binary[1:outer.length]) * seq(from = 0, to = 1, length.out = outer.length)) +
    sum((binary[(outer.length + 1):(outer.length + window.shift)] * 1)) +
    sum((binary[(outer.length + window.shift + 1):(window.size)] * seq(from = 1, to = 0, length.out = outer.length))
  )
  return(weighted)
}

# FUNCTION that sums the middle numbers in a window
n_central <- function(sequence, window.size, window.shift){
  outer.length <- (window.size - window.shift )/2
  central <- sum(sequence[(outer.length + 1):(outer.length + window.shift)])
  return(central)
}

# FUNCTION that vectorises seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### IMPORT

# reference genome and soft mask
ref <- fread(ref.file, fill = T)
# SNVs
vcf <- fread(vcf.file, fill = T)
# mask
mask <- fread(mask.file, fill = T)
# SNV rate table
P_SNV <- fread(mu_table.file, fill = T)
# alignment
align <- fread(alignment.file, fill = T)
  
  ### FORMAT
  
  vcf <- vcf[,c(2,4,5)]
  colnames(vcf) <- c("POS", "REF", "ALT")
  
  # get min and max POS for refenrecne
  min.POS <-  min(which(ref$REF != "N"))
  max.POS <-   max(which(ref$REF != "N"))
  
  # subset REF by min and max POS 
  dt <- subset(ref, ref$POS >= min.POS & ref$POS <= max.POS)
  rm(ref)
  
  # merge REF with ALT
  dt <- vcf[dt, on = c("POS", "REF")]

  # identify POS with SNV
  dt$SNV <- 0
  dt$SNV[which(dt$POS %in% vcf$POS)] <- 1
  rm(vcf)
  
  # merge mask
  colnames(mask) <- c("POS", "MASK")
  dt <- mask[dt, on = c("POS")]
  
  # identify POS that are aligned
  align <- subset(align, align$V1 == chromo)
  dt$aligned <- 0
  dt$aligned[which(dt$POS %in% unlist(seq2(from = align$V2, to = align$V3)))] <- 1
  rm(align)
  
  # get unique p_any_snp_given_kmer
  P_SNV <- P_SNV[,c("k7_from", "p_any_snp_given_k7")]
  P_SNV <- unique(P_SNV)
  
  # get probability of SNV based on kmer. 
  dt$P_SNV <- slender_worm(dt$REF, P_SNV)

  # calculate variables by sliding window using zoo
  x.out <- data.table(
    n_SNV = rollapply(dt$SNV, width = window.size, by = window.shift, FUN = sum, align = "left"),
    n_SNV_weighted = rollapply(dt$SNV, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
    
    p_SNV_given_kmers = rollapply(dt$P_SNV, width = window.size, by = window.shift, FUN = sum, align = "left"),
    p_SNV_given_kmers_weighted = rollapply(dt$P_SNV, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left"),
    
    n_CpG = rollapply(dt$REF, width = window.size, by = window.shift, FUN = n_CpG, align = "left"),
    n_CpG_weighted = rollapply(dt$REF, width = window.size, by = window.shift, FUN = n_CpG_weight, window.size = window.size, window.shift = window.shift, align = "left"),
    
    n_mask = rollapply(dt$MASK, width = window.size, by = window.shift, FUN = sum, align = "left"),
    n_mask_central = rollapply(dt$MASK, width = window.size, by = window.shift, FUN = n_central, window.size = window.size, window.shift = window.shift, align = "left"),
    n_RepMask = rollapply(dt$REP_MASK, width = window.size, by = window.shift, FUN = sum, align = "left"),
    
    n_aligned_central = rollapply(dt$aligned, width = window.size, by = window.shift, FUN = n_central, window.size = window.size, window.shift = window.shift, align = "left")
    # n_aligned = rollapply(dt$aligned, width = window.size, by = window.shift, FUN = sum, align = "left"),
    # n_aligned_weighted = rollapply(dt$aligned, width = window.size, by = window.shift, FUN = sum_weight, window.size = window.size, window.shift = window.shift, align = "left")
  )

  # Add chromosome, POS from, POS to for each window
  outer.length <- (window.size - window.shift )/2
  x.out$CHR <- rep(chromo, nrow(x.out))
  x.out$POS_from <- seq(from = (dt$POS[1] + outer.length),
                       to = (dt$POS[1] + (outer.length - 1)) + (window.shift * nrow(x.out)),
                       by = window.shift)
  x.out$POS_to <- seq(from <- (dt$POS[1] + (outer.length + window.shift - 1)),
                     to = (dt$POS[1] + (outer.length + window.shift - 1)) + (window.shift * (nrow(x.out)-1)),
                     by = window.shift)
  
  
  print(paste0(chromo, " done!"))

  
  ### EXPORT  
fwrite(x.out, out.file)


#####

# calculate variables by window 500 bp sliding 100 bp using evobiR
# x5.1 <- data.table(
#   n_SNV = SlidingWindow("sum", dt$SNV, 500, 100),
#   p_SNV_given_kmers = SlidingWindow("sum", dt$P_SNV, 500, 100),
#   # p.CpG <- SlidingWindow("CG_worm", dt$REF, 500, 100),
#   Read_depth = SlidingWindow("sum", dt$RD, 500, 100),
#   Repeats = SlidingWindow("sum", dt$REP_MASK, 500,100)
# )

# # calculate variables by window 500 bp sliding 100 bp using zoo
# x5.1 <- data.table(
#   n_SNV = rollapply(dt$SNV, width = 500, by = 100, FUN = sum, align = "left"),
#   p_SNV_given_kmers = rollapply(dt$P_SNV, width = 500, by = 100, FUN = sum, align = "left"),
#   Read_depth = rollapply(dt$RD, width = 500, by = 100, FUN = sum, align = "left"),
#   Repeats = rollapply(dt$REP_MASK, width = 500, by = 100, FUN = sum, align = "left")
# )
# 
# x5.1$n_SNV <- (x5.1$n_SNV/500)
# x5.1$p_SNV_given_kmers <- (x5.1$p_SNV_given_kmers/500)
# x5.1$Read_depth <- (x5.1$Read_depth/500)
# x5.1$Repeats <- (x5.1$Repeats/500)
# 
# # Add chromosome, POS from, POS to for each window
# x5.1$CHR <- rep(args[1], nrow(x5.1))
# x5.1$POS_from <- seq(from = (dt$POS[1] + 200),
#                      to = (dt$POS[1] + 199) + (100*nrow(x5.1)),
#                      by = 100)
# x5.1$POS_to <- seq(from <- (dt$POS[1] + 299),
#                    to = (dt$POS[1] + 299) + (100*(nrow(x5.1)-1)),
#                    by = 100)
# 
# # calculate variables by window 700 bp sliding 100 bp using evobiR
# # x7.1 <- data.table(
# #   n_SNV = SlidingWindow("sum", dt$SNV, 700, 100),
# #   p_SNV_given_kmers = SlidingWindow("sum", dt$P_SNV, 700, 100),
# #   # p.CpG <- SlidingWindow("CG_worm", dt$REF, 700, 100),
# #   Read_depth = SlidingWindow("sum", dt$RD, 700, 100),
# #   Repeats = SlidingWindow("sum", dt$REP_MASK, 700,100)
# # )
# 
# # calculate variables by window 700 bp sliding 100 bp using zoo
# x7.1 <- data.table(
#   n_SNV = rollapply(dt$SNV, width = 700, by = 100, FUN = sum, align = "left"),
#   p_SNV_given_kmers = rollapply(dt$P_SNV, width = 700, by = 100, FUN = sum, align = "left"),
#   Read_depth = rollapply(dt$RD, width = 700, by = 100, FUN = sum, align = "left"),
#   Repeats = rollapply(dt$REP_MASK, width = 700, by = 100, FUN = sum, align = "left")
# )
# 
# x7.1$n_SNV <- (x7.1$n_SNV/700)
# x7.1$p_SNV_given_kmers <- (x7.1$p_SNV_given_kmers/700)
# x7.1$Read_depth <- (x7.1$Read_depth/700)
# x7.1$Repeats <- (x7.1$Repeats/700)
# 
# # Add chromosome, POS from, POS to for each window
# x7.1$CHR <- rep(args[1], nrow(x7.1))
# x7.1$POS_from <- seq(from = (dt$POS[1] + 300),
#                       to = (dt$POS[1] + 299) + (100*nrow(x7.1)),
#                       by = 100)
# x7.1$POS_to <- seq(from <- (dt$POS[1] + 399),
#                     to = (dt$POS[1] + 399) + (100*(nrow(x7.1)-1)),
#                     by = 100)
# 
# # calculate variables by window 900 bp sliding 100 bp using evobiR
# # x9.1 <- data.table(
# #   n_SNV = SlidingWindow("sum", dt$SNV, 900, 100),
# #   p_SNV_given_kmers = SlidingWindow("sum", dt$P_SNV, 900, 100),
# #   # p.CpG <- SlidingWindow("CG_worm", dt$REF, 900, 100),
# #   Read_depth = SlidingWindow("sum", dt$RD, 900, 100),
# #   Repeats = SlidingWindow("sum", dt$REP_MASK, 900,100)
# # )
# 
# # calculate variables by window 900 bp sliding 100 bp using zoo
# x9.1 <- data.table(
#   n_SNV = rollapply(dt$SNV, width = 900, by = 100, FUN = sum, align = "left"),
#   p_SNV_given_kmers = rollapply(dt$P_SNV, width = 900, by = 100, FUN = sum, align = "left"),
#   Read_depth = rollapply(dt$RD, width = 900, by = 100, FUN = sum, align = "left"),
#   Repeats = rollapply(dt$REP_MASK, width = 900, by = 100, FUN = sum, align = "left")
# )
# 
# x9.1$n_SNV <- (x9.1$n_SNV/900)
# x9.1$p_SNV_given_kmers <- (x9.1$p_SNV_given_kmers/900)
# x9.1$Read_depth <- (x9.1$Read_depth/900)
# x9.1$Repeats <- (x9.1$Repeats/900)
# 
# # Add chromosome, POS from, POS to for each window
# x9.1$CHR <- rep(args[1], nrow(x9.1))
# x9.1$POS_from <- seq(from = (dt$POS[1] + 400),
#                    to = (dt$POS[1] + 399) + (100*nrow(x9.1)),
#                    by = 100)
# x9.1$POS_to <- seq(from <- (dt$POS[1] + 499),
#                  to = (dt$POS[1] + 499) + (100*(nrow(x9.1)-1)),
#                  by = 100)
# 
# # calculate variables by window 1100 bp sliding 100 bp using evobiR
# # x11.1 <- data.table(
# #   n_SNV = SlidingWindow("sum", dt$SNV, 1100, 100),
# #   p_SNV_given_kmers = SlidingWindow("sum", dt$P_SNV, 1100, 100),
# #   # p.CpG <- SlidingWindow("CG_worm", dt$REF, 1100, 100),
# #   Read_depth = SlidingWindow("sum", dt$RD, 1100, 100),
# #   Repeats = SlidingWindow("sum", dt$REP_MASK, 1100,100)
# # )
# 
# # calculate variables by window 1100 bp sliding 100 bp using zoo
# x11.1 <- data.table(
#   n_SNV = rollapply(dt$SNV, width = 1100, by = 100, FUN = sum, align = "left"),
#   p_SNV_given_kmers = rollapply(dt$P_SNV, width = 1100, by = 100, FUN = sum, align = "left"),
#   Read_depth = rollapply(dt$RD, width = 1100, by = 100, FUN = sum, align = "left"),
#   Repeats = rollapply(dt$REP_MASK, width = 1100, by = 100, FUN = sum, align = "left")
# )
# 
# x11.1$n_SNV <- (x11.1$n_SNV/1100)
# x11.1$p_SNV_given_kmers <- (x11.1$p_SNV_given_kmers/1100)
# x11.1$Read_depth <- (x11.1$Read_depth/1100)
# x11.1$Repeats <- (x11.1$Repeats/1100)
# 
# # Add chromosome, POS from, POS to for each window
# x11.1$CHR <- rep(args[1], nrow(x11.1))
# x11.1$POS_from <- seq(from = (dt$POS[1] + 500),
#                      to = (dt$POS[1] + 499) + (100*nrow(x11.1)),
#                      by = 100)
# x11.1$POS_to <- seq(from <- (dt$POS[1] + 599),
#                    to = (dt$POS[1] + 599) + (100*(nrow(x11.1)-1)),
#                    by = 100)

# fwrite(x5.1, paste0("/well/lindgren/George/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", chromo, "_500_100.csv"))
# fwrite(x7.1, paste0("/well/lindgren/George/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", args[1], "_700_100.csv"))
# fwrite(x9.1, paste0("/well/lindgren/George/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", args[1], "_900_100.csv"))
# fwrite(x11.1, paste0("/well/lindgren/George/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", args[1], "_1100_100.csv"))

