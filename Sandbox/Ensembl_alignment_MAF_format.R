# SCRIPT that converts Multiple Allignment Format (.MAF) to text file with each base and genomic coordinates

rm(list = ls())
graphics.off()

library(stringr)
library(plyr)
library(dplyr)
library(utils)
library(data.table)
library(zoo)
library(reshape2)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

# set args variables
input.path <- args[1]
output.path <- args[2]
dir <- args[3]
in_prefix <- args[4]
out_prefix_long <- args[5]
out_prefix_short <- args[6]
chr_id <- as.integer(args[7])
# input.path <- "/well/lindgren/George/Data/Ensembl/Alignment/Raw/"
# output.path <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/"
# dir <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/"
# in_prefix <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr"
# out_prefix_long <- "H_HtoM_alignment_long_chr"
# out_prefix_short <- "HM_alignment_short_chr"
# chr_id <- 1
input.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Raw/"
output.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/"
dir <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/"
in_prefix <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr"
out_prefix_long <- "H_HtoM_alignment_long_chr"
out_prefix_short <- "HM_alignment_short_chr"
chr_id <- 1

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
MAF_convert <- function(chr, files.chr, input.path, output.path, dir, in_prefix, out_prefix_long, out_prefix_short){
  
    # j=3
    for (j in 1:files.chr){
      
      # read
      x <- read.table(gzfile(paste0(input.path, dir, in_prefix, chr, "_", j, ".maf.gz")),
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
      
      spec1.ID <- c()
      for (id in 1:length(spec1)){
        vars <- sub("^[^.]*\\.", "", spec1[id])
        vars <- strsplit(vars, "\\s+")
        refseq <- vars[[1]][6]

        tmp <- regmatches(refseq, gregexpr("[^-]+|-+", refseq))[[1]]
        tmp_len <- unlist(lapply(tmp, nchar))

        ID_vec <- list()
        for (idd in 1:length(tmp_len)){
          ID_vec[[idd]] <- rep(paste0(id, ".", idd), tmp_len[idd])
        }
        ID_vec <- unlist(ID_vec)
        spec1.ID <- c(spec1.ID, ID_vec)
      }

      spec2.ID <- c()
      for (id in 1:length(spec2)){
        vars <- sub("^[^.]*\\.", "", spec2[id])
        vars <- strsplit(vars, "\\s+")
        refseq <- vars[[1]][6]

        tmp <- regmatches(refseq, gregexpr("[^-]+|-+", refseq))[[1]]
        tmp_len <- unlist(lapply(tmp, nchar))

        ID_vec <- list()
        for (idd in 1:length(tmp_len)){
          ID_vec[[idd]] <- rep(paste0(id, ".", idd), tmp_len[idd])
        }
        ID_vec <- unlist(ID_vec)
        spec2.ID <- c(spec2.ID, ID_vec)
      }
      
      align.ID <- paste0(spec1.ID, spec2.ID)
      
      spec1.ref <- toupper(spec1.ref)
      spec2.ref <- toupper(spec2.ref)
      bases <- c("A", "T", "C", "G")
      spec1.ref[!spec1.ref %in% bases] <- NA
      spec2.ref[!spec2.ref %in% bases] <- NA
      
      df <- data.frame(spec1.chr, spec1.pos, spec1.ref, spec2.chr, spec2.pos, spec2.ref, align.ID)
      df <- df[complete.cases(df),]
      
      # convert to short
      n_ID <- length(unique(df$align.ID))
      S1_CHR <- rep(NA, n_ID)
      S2_CHR <- rep(NA, n_ID)
      S1_start <- rep(NA, n_ID)
      S2_start <- rep(NA, n_ID)
      S1_end <- rep(NA, n_ID)
      S2_end <- rep(NA, n_ID)
      for (k in 1:n_ID){
        sub <- subset(df, df$align.ID == unique(df$align.ID)[k])
        S1_CHR[k] <- as.character(sub$spec1.chr[1])
        S2_CHR[k] <- as.character(sub$spec2.chr[1])
        S1_start[k] <- sub$spec1.pos[1]
        S2_start[k] <- sub$spec2.pos[1]
        S1_end[k] <- sub$spec1.pos[nrow(sub)]
        S2_end[k] <- sub$spec2.pos[nrow(sub)]
        # print(k)
      }
      df_short <- data.frame(S1_CHR, S1_start, S1_end,
                             S2_CHR, S2_start, S2_end)
      # identical(abs(df_short$S1_end - df_short$S1_start), abs(df_short$S2_end - df_short$S2_start))

      ### INCREASE PERFORMANCE WITH DATA.TABLE!!! 
      # dt_long <- as.data.table(df)
      # dt_long <- dt_long[,c("spec1.chr", "spec1.pos", "spec2.chr", "spec2.pos")]
      # test <- dt_long[, .(
      #   S1_CHR=spec1.chr[1L], S1_start=spec1.pos[1L], S1_end=spec1.pos[.N],
      #   S2_CHR=spec2.chr[1L], S2_start=spec2.pos[1L], S2_end=spec2.pos[.N]),
      #   by=rleid(spec1.chr, spec2.chr,
      #            c(0L, cumsum(diff(spec1.pos) > 1L)),
      #            c(0L, cumsum(diff(spec2.pos) > 1L)))][, -1L]
      # identical(abs(test$S1_end - test$S1_start), abs(test$S2_end - test$S2_start))
      # tmp <- rbind(df_short, test)
      # tmp <- tmp[!duplicated(tmp),]
      # tmp <- subset(tmp, tmp$S1_start==1026176)
      
      df <- df[,c("spec1.chr", "spec1.pos", "spec1.ref", "spec2.chr", "spec2.pos", "spec2.ref")]
      
      # # write tables
      write.table(df_short, paste0(output.path, out_prefix_short, chr, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # write.table(df, paste0(output.path, out_prefix_long, chr, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)

      # h_df_long <- df[,c("spec1.chr", "spec1.pos")]
      # m_df_long <- df[,c("spec2.chr", "spec2.pos")]
      # colnames(h_df_long) <- c("CHR", "POS")
      # colnames(m_df_long) <- c("CHR", "POS")
      # 
      # m_df_short <- ddply(m_df_long, "CHR", long_to_short)
      # h_df_short <- ddply(h_df_long, "CHR", long_to_short)
      # 
      # # short one file per species 
      # write.table(h_df_short, paste0(output.path, "H_HtoM_alignment_short.txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # write.table(m_df_short, paste0(output.path, "M_HtoM_alignment_short.txt"), sep = "\t", col.names = F, row.names = F, append = T)
      
      # short by chr
      # write.table(h_df_short, paste0(input.path, "HtoM_alignment_H_POS_chr", i, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # for (kk in c(1:19, "X")){
      #   sub_out <- subset(m_df_short, m_df_short == kk)
      #   write.table(sub_out, paste0(input.path, "HtoM_alignment_M_POS_chr", kk, ".txt"), sep = "\t", col.names = F, row.names = F, append = T)
      # }
      print(j)
    }
  }

# set chr
H_chr <- c(1:22, "X")[chr_id]

### Count files
files.chr <- file_count(input.path, dir, H_chr)


### Convert

MAF_convert(chr = H_chr,
            files.chr = files.chr,
            input.path = input.path,
            output.path = output.path,
            dir = dir,
            in_prefix = in_prefix,
            out_prefix_long = out_prefix_long,
            out_prefix_short = out_prefix_short)


##########

# #### CHECK ALIGNS WITH 38 REF
# m.ref <- fread("/Volumes/HarEx/Data/Ensembl/Reference/Formatted/Mouse_REF_Ensembl_v94_chr10.txt")
# test_df <- subset(df, df$spec2.chr == "10")
# test_df <- test_df[,c("spec2.pos", "spec2.ref")]
# colnames(test_df) <- c("POS", "REF2")
# test <- m.ref[test_df, on="POS"]
# no_match <- test[which(test$REF != test$REF2),]
# match <- test[which(test$REF == test$REF2),]


# ### STACK OVERFLOW
# 
# dt_long <- data.table(LABEL_A = c(rep("A", 20), rep("A", 15), rep ("A", 10), rep ("A", 15), rep ("A", 10)),
#                       SEQ_A = c(11:25, 31:35, 61:75, 76:85, 86:100, 110:119),
#                       LABEL_B= c(rep("C", 20), rep("D", 15), rep("F", 10), rep("G",15), rep("D", 10)),
#                       SEQ_B = c(1:20, 25:11, 16:25, 15:1, 1:5, 8:12))
# 
# dt_short <- data.table(LABEL_A = c("A", "A", "A", "A", "A", "A", "A"),
#                        Start_A = c(11, 31, 61, 76, 86, 110, 115),
#                        End_A = c(25, 35, 75, 85, 100, 114, 119),
#                        LABEL_B= c("C", "C", "D", "F", "G", "D", "D"),
#                        Start_B = c(1, 16, 25, 16, 15, 1, 8),
#                        End_B = c(15, 20, 11, 25, 1, 5, 12))
# identical(abs(dt_short$End_A - dt_short$Start_A), abs(dt_short$End_B - dt_short$Start_B))
# 
# test <- dt_long[, .(
#   LABEL_A=LABEL_A[1L], Start_A=SEQ_A[1L], End_A=SEQ_A[.N],
#   LABEL_B=LABEL_B[1L], Start_B=SEQ_B[1L], End_B=SEQ_B[.N]),
#   by=rleid(LABEL_A, LABEL_B,
#            c(0L, cumsum(diff(SEQ_A) > 1L)),
#            c(0L, cumsum(diff(SEQ_B) > 1L)))][, -1L]
# 
# identical(abs(test$End_A - test$Start_A), abs(test$End_B - test$Start_B))
# 
# 
# dt_long[, ID := cumsum(c(0, +(LABEL_A != shift(LABEL_A) | LABEL_B != shift(LABEL_B))[-1])) + 1]
# 
# 
# 
# dt_long_ID <- data.table(LABEL_A = c(rep("A", 20), rep("A", 15), rep ("A", 10), rep ("A", 15), rep ("A", 10)),
#                          SEQ_A = c(11:30, 61:75, 76:85, 86:100, 110:119),
#                          LABEL_B= c(rep("C", 20), rep("D", 15), rep("F", 10), rep("G",15), rep("D", 10)),
#                          SEQ_B = c(1:20, 25:11, 16:25, 15:1, 1:5, 8:12),
#                          ID = c(rep(1, 20), rep(2, 15), rep(3, 10), rep(4, 15), rep(5, 5), rep(6, 5) ))
# 
# 
# 
# 
# identical(abs(dt_short$End_A - dt_short$Start_A), abs(dt_short$End_B - dt_short$Start_B))


# dt_long[, ID := .GRP, by = .(rleid(LABEL_A), rleid(LABEL_B))]
# dt_long[, ID := c(0, cumsum(diff(LABEL_A) != 0 | diff(LABEL_B) != 0)) + 1]
# 
# 
# 
# tmp_pos <- sort(dt_long$SEQ_B, decreasing = F)
# tmp_pos <- as.data.table(t(sapply(split(tmp_pos, findInterval(tmp_pos, tmp_pos[which(c(1, diff(tmp_pos)) > 1)])), range)))
# 
# 
# dcast(dt_long[, .SD[c(1, .N)], .(LABEL_A, LABEL_B)], 
#       LABEL_A + LABEL_B ~ c("Start", "End")[rowid(LABEL_A, LABEL_B)], 
#       value.var = c("SEQ_A", "SEQ_B"))
# 
# 
# 
# 
# sub <- subset(dt_long, dt_long$LABEL_B == 3)
# ballbag <- function(sub){
#   tmp_pos <- sort(sub$SEQ_B, decreasing = F)
#   tmp_pos <- as.data.table(t(sapply(split(tmp_pos, findInterval(tmp_pos, tmp_pos[which(c(1, diff(tmp_pos)) > 1)])), range)))
# }
# 
# 
# dcast(dt_long[, .SD[c(1, .N)], .(LABEL_A, LABEL_B)], 
#       LABEL_A + LABEL_B ~ c("Start", "End")[rowid(LABEL_A, LABEL_B)], 
#       value.var = c("SEQ_A", "SEQ_B"))
# 
# test <- dt_long[, .(Start_A = first(SEQ_A), End_A = last(SEQ_A), Start_B = first(SEQ_B), End_B = last(SEQ_B)), by = .(LABEL_A, LABEL_B)][]
# 
# 
# 
# 
# 
# 
# 
# df <- as.data.table(df)
# test <- df[, .(spec1.start = first(spec1.pos), spec1.end = last(spec1.pos), spec2.start = first(spec2.pos), spec2.end = last(spec2.pos)), 
#            by = .(spec1.chr, spec2.chr)][]
# test$len_S1 <- test$spec1.end - test$spec1.start
# test$len_S2 <- test$spec2.end - test$spec2.start
# 
# 
# string <- "AATTGGCGCTAG---AT-TTACG----"
# 
# string1 <- "AATTGGCGCTAG"
# string2 <- "---"
# string3 <- "AT"
# string4 <- "-"
# string5 <- "TTACG"
# string6 <- "----"
# 
# strsplit(string, "[-]+")
# res <- regmatches(string, gregexpr("[^-]+|-+", string))





