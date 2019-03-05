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
allignment.path <- "/Volumes/HarEx/PC_constraint/Paper/Data/Ensembl_alignments/"

# directories
dir <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/"

# file prefix
in.file.prefix <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr"
out.file.prefix <- "Mouse_human_alignment_chr"

# chromosomes
chr <- c(1:22, "X")

### Count files
files.chr <- file_count(allignment.path, dir, chr)

### Convert

MAF_convert(chr = chr,
            files.chr = files.chr,
            allignment.path = allignment.path,
            dir = dir,
            in.file.prefix = in.file.prefix,
            out.file.prefix = out.file.prefix)


