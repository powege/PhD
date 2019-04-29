# SCRIPT that estimates the relative probabilities of substitution for the middle base in 7-mers.
# ASSUMES NO MISSING POS in REFERENCE
# ASSUMES MISSING BASES FILLED WITH N
# ARGS: REF file, vcf file, output file, chromosome, all/unmasked

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  stop("Five argumants must be supplied: REF file path, vcf file path, out file path, chromosome, all/unmasked", call.=FALSE)
} 

# set args
ref_file <- args[1]
vcf_file <- args[2]
out_file <- args[3]
chr <- args[4]
mask <- args[5]

#################
### FUNCTIONS ###    
#################

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

##############
### SCRIPT ###    
##############

    ### IMPORT reference genome and SNV POS
    ref <- fread(ref_file, fill = T)
    vcf <- fread(vcf_file, fill = T)
    
    ### merge
    vcf <- vcf[,c(2,4,5)]
    colnames(vcf) <- c("POS", "REF", "ALT")
    ref <- vcf[ref, on = c("POS", "REF")]
    rm(vcf)
    
    ### calculate ancestral 7mer for each POS
    b1 <- ref$REF[1:(length(ref$REF)-6)]
    b2 <- ref$REF[2:(length(ref$REF)-5)]
    b3 <- ref$REF[3:(length(ref$REF)-4)]
    b4 <- ref$REF[4:(length(ref$REF)-3)]
    b5 <- ref$REF[5:(length(ref$REF)-2)]
    b6 <- ref$REF[6:(length(ref$REF)-1)]
    b7 <- ref$REF[7:(length(ref$REF))]
    kmer <- paste0(b1, b2, b3, b4, b5, b6, b7)
    kmer[grep("N", kmer)] <- NA
    ref$ANCESTRAL_KMER <- c(rep(NA, 3), kmer, rep(NA, 3)) # ensure output equal length to input
    rm(b1, b2, b3, b4, b5, b6, b7, kmer)
    
    ### calculate mutant 7mer for each POS
    SNV.POS <- which(ref$REF != ref$ALT) # identify SNV POS
    b1 <- ref$REF[1:(length(ref$REF)-6)]
    b2 <- ref$REF[2:(length(ref$REF)-5)]
    b3 <- ref$REF[3:(length(ref$REF)-4)]
    b4 <- ref$ALT[4:(length(ref$REF)-3)]
    b5 <- ref$REF[5:(length(ref$REF)-2)]
    b6 <- ref$REF[6:(length(ref$REF)-1)]
    b7 <- ref$REF[7:(length(ref$REF))]
    kmer <- paste0(b1, b2, b3, b4, b5, b6, b7)
    kmer[grep("N", kmer)] <- NA
    ref$MUTANT_KMER <- c(rep(NA, 3), kmer, rep(NA, 3)) # ensure output equal length to input
    rm(b1, b2, b3, b4, b5, b6, b7, kmer)
    
    if(mask == "unmasked"){
    ### cut masked regions
    ref <- subset(ref, ref$REP_MASK == 0)
    }

    ### count ancestral kmers
    k7_mu_all <- k7mer_changes(7) # Identify all possible 7-mer substitutions
    k7_all <- data.table(k7_from = unique(k7_mu_all$k7_from))
    kmer_from <- as.data.table(table(ref$ANCESTRAL_KMER))
    colnames(kmer_from) <- c(colnames(k7_all), "k7_from_N")
    kmer_from <- kmer_from[k7_all, on = "k7_from"]
    kmer_from$k7_from_N[is.na(kmer_from$k7_from_N)] <- 0
    
    ### count kmer changes
    SNV.POS <- which(ref$REF != ref$ALT) # identify SNV POS
    mu.count <- data.table(kmer.from = ref$ANCESTRAL_KMER[SNV.POS],
                                kmer.to = ref$MUTANT_KMER[SNV.POS])
    mu.count <- as.data.table(table(mu.count))
    colnames(mu.count) <- c("k7_from", "k7_to", "k7_mu_N")
    mu.count <- mu.count[mu.count$k7_mu_N != 0, ]
    mu.count <- mu.count[k7_mu_all, on = c("k7_from", "k7_to")]
    mu.count$k7_mu_N[is.na(mu.count$k7_mu_N)] <- 0
    
    ### OUTPUT
    k7_out <- merge(kmer_from, mu.count, all = T)
    k7_out$Chromosome <- chr # add chomosome
    fwrite(k7_out, out_file)
    
    
   