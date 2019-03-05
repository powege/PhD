# SCRIPT that converts Multiple Allignment Format (.MAF) to text file with each base and genomic coordinates

rm(list = ls())
graphics.off()

library(stringr)
library(dplyr)
library(utils)
library(data.table)
library(zoo)

### FUNCTIONS

# Function that retuns reference sequence (B6J) from .MAF input
turnip1 <- function(ref1){
  
  # regular expression!
  start <- str_extract(ref1, "(?<=\\s)\\d+(?=\\s)")
  refseq <- str_extract(ref1, "[\\w\\-]+(?=\\s*)$")
  refseq <- strsplit(refseq, "")[[1]]
  
  return(refseq)
}

# Function that retuns reference sequence start coordinate from .MAF input
turnip2 <- function(ref1){
  
  # regular expression!
  start <- str_extract(ref1, "(?<=\\s)\\d+(?=\\s)")
  refseq <- str_extract(ref1, "[\\w\\-]+(?=\\s*)$")
  refseq <- strsplit(refseq, "")[[1]]
  
  refbp <- rep(NA_integer_, length(refseq))
  refbp[refseq != "-"] <- seq(from = start, by = 1, length.out = sum(refseq != "-"))
  
  return(refbp)
}

# Function that retuns alligned sequence (Mus caroli) from .MAF input
turnip3 <- function(car1){
  
  # regular expression!
  start <- str_extract(car1, "(?<=\\s)\\d+(?=\\s)")
  carseq <- str_extract(car1, "[\\w\\-]+(?=\\s*)$")
  carseq <- strsplit(carseq, "")[[1]]
  carseq2 <- replace(carseq, carseq=="-", NA)
  
  return(carseq)
}

# function that counts the nember of .maf files for each chromosome
file_count <- function(allignment.path, dir, chr){
  files <- list.files(paste0(allignment.path, dir))
  files.chr <- list()
  for (i in chr){
    files.chr[[i]] <- length(grep(paste0("chr", i, "_"), files))
  }
  files.chr <- unlist(files.chr)
  return(files.chr)
}

# Function that converts .maf files to to tab delimited file with each base and genomic coordinates
# output cols: (POS; REF; SPECIES)
MAF_convert <- function(chr, files.chr, allignment.path, dir, in.file.prefix, out.file.prefix){
  
  for (i in chr){
    
    for (j in 1:files.chr[i]){
      
      # read
      x <- read.table(gzfile(paste0(allignment.path, dir, in.file.prefix, i, "_", j, ".maf.gz")),
                      sep = "\t")
      
      # format
      x <- as.character(x$V1)
      
      ind.ref <- seq(2, length(x), by=3)
      ref <- x[ind.ref]
      ind.car <- seq(3, length(x), by=3)
      car <- x[ind.car]
      
      ref.bp <- unlist(lapply(ref, turnip2))
      ref.seq <- unlist(lapply(ref, turnip1))
      car.seq <- unlist(lapply(car, turnip3))
      
      ref.seq <- toupper(ref.seq)
      car.seq <- toupper(car.seq)
      
      bases <- c("A", "T", "C", "G")
      car.seq[!car.seq %in% bases] <- NA
      ref.seq[!ref.seq %in% bases] <- NA
      
      df <- data.frame(ref.bp, ref.seq, car.seq)
      df <- df[!(is.na(df$ref.bp)),]
      
      write.table(df, paste0(allignment.path, out.file.prefix, i, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      
      print(i)
      print(j)
    }
  }
}


### SET VARS (chr, files.chr, allignment.path, dir, in.file.prefix, out.file.prefix)

# path
allignment.path <- "../../Data/Ensembl_allignments/"

# directories
M_dir <- "mus_musculus_GRCm38_vs_mus_caroli_CAROLI_EIJ_v1_1_lastz_net/"
H_dir <- "hsap_grch38.v.ptro_pan_tro_3.0.lastz_net/"

# file prefix
M_in_prefix <- "mus_musculus_GRCm38_vs_mus_caroli_CAROLI_EIJ_v1_1_lastz_net.chr"
M_out_prefix <- "M_caroli_allignment_chr"
H_in_prefix <- "homo_sapiens_GRCh38_vs_pan_troglodytes_Pan_tro_3_0_lastz_net.chr"
H_out_prefix <- "H_chimpanzee_allignment_chr"

# chromosomes
M_chr <- c(1:19, "X")
H_chr <- c(1:22, "X")

### Count files
M_files.chr <- file_count(allignment.path, M_dir, M_chr)
H_files.chr <- file_count(allignment.path, H_dir, H_chr)

### Convert

# MAF_convert(chr = M_chr, 
#             files.chr = M_files.chr, 
#             allignment.path = allignment.path, 
#             dir = M_dir, 
#             in.file.prefix = M_in_prefix, 
#             out.file.prefix = M_out_prefix)

# MAF_convert(chr = H_chr, 
#             files.chr = H_files.chr, 
#             allignment.path = allignment.path, 
#             dir = H_dir, 
#             in.file.prefix = H_in_prefix, 
#             out.file.prefix = H_out_prefix)

