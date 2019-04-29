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
  stop("At least one argument must be supplied", call.=FALSE)
} 

# set args
chromo <- args[1]
ref.file <- args[2]
vcf.file <- args[3]
mask.file <- args[4]
mu_table.file <- args[5]
species <- args[6]
# chromo <- 19
# ref.file <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr19.txt"
# vcf.file <- "/well/lindgren/George/Data/1KGP/Variants/vcf_QCed_VEP/1000GP_phase3_QCed_VEP_v94_allPASS_chr19.txt" 
# mask.file <- "/well/lindgren/George/Data/MGP/Masks/MGP_Mask_chr19.txt"
# mu_table.file <- "/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_StrictMask_chr19.txt" 
# species <- "human"


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

# FUNCTION that calculates the weighted sum of a vector of length 650. 
# the middle 50 numbers are weighted 1
# all other numbers are weighted 0.1
sum_650_50_weight <- function(seq){
  central <- sum(seq[c(301:350)]) * 1
  outer <- sum(seq[c(1:300,351:650)]) * 0.1
  weighted <- (central + outer)
  return(weighted)
}

# FUNCTION that calculates the weighted number of CpG dinucleotides in a sequence (650 long).
# the middle 52 numbers are weighted 1 
# all other numbers are weighted 0.1
n_CpG_650_50_weight <- function(seq){
  central_CG_count <- str_count(paste0(seq[c(300:351)], collapse = ""), "CG") * 1
  outer_CG_count1 <- str_count(paste0(seq[c(1:300)], collapse = ""), "CG") * 0.1
  outer_CG_count2 <- str_count(paste0(seq[c(351:650)], collapse = ""), "CG") * 0.1
  weighted <- (central_CG_count + outer_CG_count1 + outer_CG_count2)
  return(weighted)
}

# FUNCTION that sums the middle 50 numbers in a vector length 650
n_mask_650_50_central <- function(seq){
  central <- sum(seq[c(301:350)])
  return(central)
}


### IMPORT

# reference genome
ref <- fread(ref.file, fill = T)
# SNVs
vcf <- fread(vcf.file, fill = T)
# mask
mask <- fread(mask.file, fill = T)
# SNV rate table
P_SNV <- fread(mu_table.file, fill = T)

  
  ### FORMAT
  
  vcf <- vcf[,c(2,4,5)]
  colnames(vcf) <- c("POS", "REF", "ALT")
  
  # get min and max POS for refenrecne
  min.POS <-  min(which(ref$REF != "N"))
  max.POS <-   max(which(ref$REF != "N"))
  
  # subset REF by min and max POS 
  dt <- subset(ref, ref$POS >= min.POS & ref$POS <= max.POS)
  dt$REP_MASK <- NULL
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
  
  # get unique p_any_snp_given_kmer
  P_SNV <- P_SNV[,c("k7_from", "p_any_snp_given_k7")]
  P_SNV <- unique(P_SNV)
  
  # get probability of SNV based on kmer. 
  dt$P_SNV <- slender_worm(dt$REF, P_SNV)

  # calculate variables by window 650 bp sliding 50 bp using zoo
  x65.05 <- data.table(
    n_SNV = rollapply(dt$SNV, width = 650, by = 50, FUN = sum, align = "left"),
    n_SNV_weighted = rollapply(dt$SNV, width = 650, by = 50, FUN = sum_650_50_weight, align = "left"),
    
    p_SNV_given_kmers = rollapply(dt$P_SNV, width = 650, by = 50, FUN = sum, align = "left"),
    p_SNV_given_kmers_weighted = rollapply(dt$P_SNV, width = 650, by = 50, FUN = sum_650_50_weight, align = "left"),
    
    n_CpG = rollapply(dt$REF, width = 650, by = 50, FUN = n_CpG, align = "left"),
    n_CpG_weighted = rollapply(dt$REF, width = 650, by = 50, FUN = n_CpG_650_50_weight, align = "left"),
    
    n_mask = rollapply(dt$MASK, width = 650, by = 50, FUN = sum, align = "left"),
    n_mask_central = rollapply(dt$MASK, width = 650, by = 50, FUN = n_mask_650_50_central, align = "left")
  )

  # Add chromosome, POS from, POS to for each window
  x65.05$CHR <- rep(args[1], nrow(x65.05))
  x65.05$POS_from <- seq(from = (dt$POS[1] + 300),
                       to = (dt$POS[1] + 299) + (50*nrow(x65.05)),
                       by = 50)
  x65.05$POS_to <- seq(from <- (dt$POS[1] + 349),
                     to = (dt$POS[1] + 349) + (50*(nrow(x65.05)-1)),
                     by = 50)
  
  
  print(paste0(chromo, " done!"))

  
  ### EXPORT  
fwrite(x65.05, paste0("/well/lindgren/George/Data/NC_constraint/Constraint/", species, "_constraint_variables_by_window_chr", chromo, "_650_50.csv"))


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

