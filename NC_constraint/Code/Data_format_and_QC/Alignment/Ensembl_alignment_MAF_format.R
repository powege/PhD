# SCRIPT that converts Multiple Allignment Format (.MAF) to text file with each base and genomic coordinates

rm(list = ls())
graphics.off()

library(stringr)
library(plyr)
library(dplyr)
library(utils)
library(data.table)
library(zoo)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

# set args variables
# (chr, files.chr, input.path, dir, in.file.prefix, out.file.prefix)
input.path <- args[1]
output.path <- args[2]
dir <- args[3]
in_prefix <- args[4]
out_prefix <- args[5]
# input.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Raw/"
# output.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/"
# dir <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/"
# in_prefix <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr"
# out_prefix <- "H_HtoM_alignment_long_chr"


### FUNCTIONS

# Function that returns sequence chromosome as string from .maf input
turnip_CHR <- function(ref1){
 
  vars <- sub("^[^.]*\\.", "", ref1)
  vars <- strsplit(vars, "\\s+")
  chr <- vars[[1]][1]
  refseq <- vars[[1]][6]
  refseq <- strsplit(refseq, "")[[1]]
  chr_vec <- rep(NA, length(refseq))
  chr_vec[refseq != "-"] <- chr
  
  return(chr_vec)
}


# Function that retuns sequence coordinates as string from .maf input
turnip_POS <- function(ref1){
  
  vars <- sub("^[^.]*\\.", "", ref1)
  vars <- strsplit(vars, "\\s+")
  
  # identify strand direction
  strand <- vars[[1]][4]
  # strand <- str_extract(ref1, "(?<=\\s)[+-](?=\\s)")
  
  # if + strand:
  if (strand == "+"){
    start <- as.integer(vars[[1]][2]) + 1
    # start <- str_extract(ref1, "(?<=\\s)\\d+(?=\\s)") 
    # start <- as.integer(start) + 1
    refseq <- vars[[1]][6]
    # refseq <- str_extract(ref1, "[\\w\\-]+(?=\\s*)$")
    refseq <- strsplit(refseq, "")[[1]]
    refbp <- rep(NA_integer_, length(refseq))
    refbp[refseq != "-"] <- seq(from = start, by = 1, length.out = sum(refseq != "-"))
  }
  
  # if - strand: 
  if (strand == "-"){
    chr_length <- as.integer(vars[[1]][5])
    minus <- as.integer(vars[[1]][2])
    start <- chr_length - minus
    refseq <- vars[[1]][6]
    refseq <- strsplit(refseq, "")[[1]]
    refbp <- rep(NA_integer_, length(refseq))
    refbp[refseq != "-"] <- seq(from = start, by = -1, length.out = sum(refseq != "-"))
  }
  
  return(refbp)
}

# Function that retuns sequence bases as string from .maf input
turnip_REF <- function(ref1){
  
  refseq <- str_extract(ref1, "[\\w\\-]+(?=\\s*)$")
  refseq <- strsplit(refseq, "")[[1]]
  refseq <- replace(refseq, refseq == "-", NA)
  
  return(refseq)
}


# Function that counts the nember of .maf files for each chromosome
file_count <- function(input.path, dir, chr){
  files <- list.files(paste0(input.path, dir))
  files.chr <- list()
  for (i in chr){
    files.chr[[i]] <- length(grep(paste0("chr", i, "_"), files))
  }
  files.chr <- unlist(files.chr)
  return(files.chr)
}

### FUNCTION that converts dataframe with cols c("CHR", "POS") from long format to short (c("CHR", "start", "end"))
# sub <- subset(m_df, m_df$spec2.chr == "1")
long_to_short <- function(sub){
  tmp_pos <- sort(sub[,2], decreasing = F)
  tmp_pos <- as.data.table(t(sapply(split(tmp_pos, findInterval(tmp_pos, tmp_pos[which(c(1, diff(tmp_pos)) > 1)])), range)))
  colnames(tmp_pos) <- c("start", "end")
  # tmp_pos$CHR <- sub[1,1]
  return(tmp_pos)
}


# Function that converts .maf files to to tab delimited file with each base and genomic coordinates
# output cols: (POS; REF; SPECIES)
MAF_convert <- function(chr, files.chr, input.path, output.path, dir, in.file.prefix, out.file.prefix){
  
  for (i in chr){
    
    for (j in 1:files.chr[i]){
      
      # read
      x <- read.table(gzfile(paste0(input.path, dir, in.file.prefix, i, "_", j, ".maf.gz")),
                      sep = "\t")
      
      # format
      x <- as.character(x$V1)
      
      ind.spec1 <- seq(2, length(x), by=3)
      spec1 <- x[ind.spec1]
      ind.spec2 <- seq(3, length(x), by=3)
      spec2 <- x[ind.spec2]
      
      spec1.pos <- unlist(lapply(spec1, turnip_POS))
      spec2.pos <- unlist(lapply(spec2, turnip_POS))
      
      spec1.chr <- unlist(lapply(spec1, turnip_CHR))
      spec2.chr <- unlist(lapply(spec2, turnip_CHR))
      
      spec1.ref <- unlist(lapply(spec1, turnip_REF))
      spec2.ref <- unlist(lapply(spec2, turnip_REF))
      spec1.ref <- toupper(spec1.ref)
      spec2.ref <- toupper(spec2.ref)
      
      bases <- c("A", "T", "C", "G")
      spec1.ref[!spec1.ref %in% bases] <- NA
      spec2.ref[!spec2.ref %in% bases] <- NA
      
      df <- data.frame(spec1.chr, spec1.pos, spec1.ref, spec2.chr, spec2.pos, spec2.ref)
      df <- df[complete.cases(df),]
      write.table(df, paste0(output.path, out.file.prefix, i, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)

      h_df_long <- df[,c("spec1.chr", "spec1.pos")]
      m_df_long <- df[,c("spec2.chr", "spec2.pos")]
      colnames(h_df_long) <- c("CHR", "POS")
      colnames(m_df_long) <- c("CHR", "POS")
      
      m_df_short <- ddply(m_df_long, "CHR", long_to_short)
      h_df_short <- ddply(h_df_long, "CHR", long_to_short)
      
      # short one file per species 
      write.table(h_df_short, paste0(output.path, "H_HtoM_alignment_short.txt"), sep = "\t", col.names = F, row.names = F, append = T)
      write.table(m_df_short, paste0(output.path, "M_HtoM_alignment_short.txt"), sep = "\t", col.names = F, row.names = F, append = T)
      
      # short by chr
      # write.table(h_df_short, paste0(input.path, "HtoM_alignment_H_POS_chr", i, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # for (kk in c(1:19, "X")){
      #   sub_out <- subset(m_df_short, m_df_short == kk)
      #   write.table(sub_out, paste0(input.path, "HtoM_alignment_M_POS_chr", kk, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # }
      print(i)
      print(j)
    }
  }
}

# set chr
H_chr <- c(1:22, "X")

### Count files
files.chr <- file_count(input.path, dir, H_chr)

### Convert

MAF_convert(chr = H_chr,
            files.chr = files.chr,
            input.path = input.path,
            output.path = output.path,
            dir = dir,
            in.file.prefix = in_prefix,
            out.file.prefix = out_prefix)


##########

# #### CHECK ALIGNS WITH 38 REF
# m.ref <- fread("/Volumes/HarEx/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_v94_chr10.txt")
# test_df <- subset(df, df$spec2.chr == "10")
# test_df <- test_df[,c("spec2.pos", "spec2.ref")]
# colnames(test_df) <- c("POS", "REF2")
# test <- m.ref[test_df, on="POS"]
# no_match <- test[which(test$REF != test$REF2),]
# match <- test[which(test$REF == test$REF2),]



