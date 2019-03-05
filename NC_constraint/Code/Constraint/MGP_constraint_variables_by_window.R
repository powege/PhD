# script that creates data table for each chromosome
# Header: POS; REF; ALT; SNV; P_SNV; Read Depth; Human REF; 

rm(list = ls())
graphics.off()

library(data.table)
library(evobiR)
library(MASS)
library(stringr)

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

# # FUNCTION that calculates the proportion of CpG dinucleotides in a sequence
# CG_worm <- function(seq){
#   string <- paste0(seq, collapse = "")
#   CG_count <- str_count(string, "CG")
#   total.dinuc <- nchar(string)/2
#   CG_prop <- CG_count/total.dinuc
#   return(CG_prop)
# }


# set output list
output <- list()

# set chr
chr <- c(1:19, "X")

for (i in 1:length(chr)){
  
  ### IMPORT
  ref <- fread(paste0("/Volumes/HarEx/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_v94_chr", chr[i], ".txt"))
  vcf <- fread(paste0("~/Dropbox/PhD/Data/MGP/vcf_QCed_VEP/M_MGP_QCed_VEP_all_chr", chr[i], ".txt"))
  # alignment <- fread(paste0("/Volumes/HarEx/PC_constraint/Paper/Data/Ensembl_alignments/Mouse_Muscaroli_alignment_chr", chr[i], ".txt"))
  rd <- fread(paste0("/Volumes/HarEx/Data/MGP/Read_depth/MGP_low_depth_POS_chr", chr[i], ".txt"))
  P_SNV <- fread("~/Dropbox/BitBucket_repos/phd/NC_constraint/Data/SNV_rates/MGP_7mer_SNV_rate.table")
  
  ### FORMAT
  
  # get min and max POS for refenrecne
  vcf <- vcf[,c(2,4,5)]
  colnames(vcf) <- c("POS", "REF", "ALT")
  min.POS <-  min(which(ref$REF != "N"))
  max.POS <-   max(which(ref$REF != "N"))
  
  # subset REF by min and max POS 
  dt <- subset(ref, ref$POS >= min.POS & ref$POS <= max.POS)
  rm(ref)
  
  # merge REF with ALT
  dt <- vcf[dt, on = c("POS", "REF")]
  
  # identify POS wit SNV
  dt$SNV <- 0
  dt$SNV[which(dt$POS %in% vcf$POS)] <- 1
  rm(vcf)
  
  # get POS with low coverage
  rd <- rd$V2
  dt$RD <- 1
  dt$RD[dt$POS %in% rd] <- 0
  
  # get unique p_any_snp_given_kmer
  P_SNV <- P_SNV[,c("k7_from", "p_any_snp_given_k7")]
  P_SNV <- unique(P_SNV)
  
  # get probability of SNV based on kmer. 
  dt$P_SNV <- slender_worm(dt$REF, P_SNV)

  # calculate variables by window
  x1 <- data.table(
                  n_SNV = SlidingWindow("sum", dt$SNV, 900, 100),
                  p_SNV_given_kmers = SlidingWindow("sum", dt$P_SNV, 900, 100),
                  # p.CpG <- SlidingWindow("CG_worm", dt$REF, 900, 100),
                  Read_depth = SlidingWindow("sum", dt$RD, 900, 100)
                  )

  x1$Read_depth <- (x1$Read_depth/900)
  
  # Add chromosome, POS from, POS to for each window
  x1$CHR <- rep(chr[i], nrow(x1))
  x1$POS_from <- seq(from = (dt$POS[1] + 400),
                  to = (dt$POS[1] + 399) + (100*nrow(x1)),
                  by = 100)
  x1$POS_to <- seq(from <- (dt$POS[1] + 499),
                to = (dt$POS[1] + 499) + (100*(nrow(x1)-1)),
                by = 100)

  # subset windows with >= 50% POS have RD >= 10x for all strains. 
  output[[i]] <- subset(x1, x1$Read_depth >= 0.5)

  print(paste0(chr[i], " done!"))
  
}

out <- do.call("rbind", output)
fwrite(out, "~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_variables_by_window.csv")


#####


# dt <- fread(paste0("~/Dropbox/BitBucket_repos/phd/NC_constraint/Data/Formatted/M_POS_REF_ALT_SNV_pSNV_RD_ALIGN_chr", chr, ".txt"))

# get_POS_401 <- function(POS){
#   return(POS[401])
# }
# id.2 <- SlidingWindow("get_POS_401", dt$POS, 900, 100)
# x <- cbind(id, x)
# 
# hist(x1$RD)
# x.sub <- subset(x1, x1$RD >= 450)
# hist(x.sub$RD)
# 
# mod <- lm(x.sub$n.SNV ~ x.sub$p.Kmer + x.sub$RD)
# summary(mod)
# constraint <- studres(mod)
# 
# sub.int <- sample(nrow(x.sub), 100000, replace = F)
# sub <- x.sub[sub.int,]
# plot(sub$n.SNV ~ sub$p.Kmer)
# summary(lm((sub$n.SNV ~ sub$p.Kmer)))
# plot(sub$n.SNV ~ sub$RD)
# plot(sub$n.SNV~sub$p.CpG)


