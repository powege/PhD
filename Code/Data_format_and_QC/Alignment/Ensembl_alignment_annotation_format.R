# SCRIPT that converts Multiple Allignment Format (.MAF) to text file 

rm(list = ls())
graphics.off()

library(stringr)
library(stringi)
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

### set args variables
input.path <- args[1]
output.path <- args[2]
dir <- args[3]
in_prefix <- args[4]
out_long_prefix <- args[5]
out_short_prefix <- args[6]
h_ann_file <- args[7]
m_ann_file <- args[8]
chr_id <- as.integer(args[9])

# input.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Raw/"
# output.path <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/"
# dir <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/"
# in_prefix <- "homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr"
# out_long_prefix <- "HM_alignment_annotation_long_Hchr"
# out_short_prefix <- "HM_alignment_annotation_short_Hchr"
# h_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# m_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv"
# chr_id <- 1

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

# Function that converts long to short by ind
# sub <- subset(dt_long, dt_long$ID == "1427119318314353")
long_to_short <- function(sub){
  data.table(
    CHR_A = sub$spec1.chr[1],
    POS_A_START = sub$spec1.pos[1],
    POS_A_END = sub$spec1.pos[nrow(sub)],
    ANN_A =  sub$spec1.ann[1],
    CHR_B = sub$spec2.chr[1],
    POS_B_START = sub$spec2.pos[1],
    POS_B_END = sub$spec2.pos[nrow(sub)],
    ANN_B =  sub$spec2.ann[1]
  )
}

# FUNCTION that sorts the characters in a string, or vercor of strings
sort_cat3 <- function(strings){
  apply(stri_extract_all_regex(strings, "\\p{L}", simplify = TRUE), 1, function(i){
    stri_join(stri_sort(i), collapse = "")
  })
}

# Function that converts .maf files to to tab delimited file with each base and genomic coordinates
# output cols: (POS; REF; SPECIES)
MAF_convert <- function(chr, files.chr, input.path, output.path, dir, in_prefix, out_long_prefix, out_short_prefix){
  
  # j=30
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
    
    spec1.ref <- toupper(spec1.ref) # remove softmask
    spec2.ref <- toupper(spec2.ref)
    bases <- c("A", "T", "C", "G")
    spec1.ref[!spec1.ref %in% bases] <- NA
    spec2.ref[!spec2.ref %in% bases] <- NA
    
    df <- data.frame(spec1.chr, spec1.pos, spec1.ref, spec2.chr, spec2.pos, spec2.ref)
    df <- subset(df, df$spec2.chr %in% c(1:19, "X")) # subset automomes and X (remove patches)
    df <- df[complete.cases(df),]
    
    # data.table to get annotations by POS
    h_dt1 <- data.table(
      CHR = as.character(df$spec1.chr),
      INT = df$spec1.pos)
    h_dt2 <- data.table(
      CHR = as.character(h_ann$chromosome),
      INT_START = h_ann$start,
      INT_END = h_ann$end,
      CAT = as.character(h_ann$category)) # category as.character
    # str(h_dt1)
    # str(h_dt2)
    # which(h_dt2$INT_START > h_dt2$INT_END)
    h_dt1 <- h_dt2[h_dt1, on=.(CHR=CHR, INT_START<=INT, INT_END>=INT), allow.cartesian = T]
    if (identical(h_dt1$INT_END, h_dt1$INT_START)){
      h_dt1 <- h_dt1[,c("CHR", "INT_END", "CAT")]
      colnames(h_dt1) <- c("CHR", "INT", "CAT")
      h_dt1 <- h_dt1[!duplicated(h_dt1),]
    }
    h_dt1 <- dcast(h_dt1, CHR + INT ~ CAT, value.var = "CAT")
    h_dt1[is.na(h_dt1)] <- ""
    if (ncol(h_dt1) > 3){
      h_dt1 <- data.table(
        spec1.chr = h_dt1$CHR,
        spec1.pos = h_dt1$INT,
        spec1.ann = apply(h_dt1[,3:ncol(h_dt1)], 1, paste0, collapse = ""))
    } else {
      h_dt1 <- as.data.table(h_dt1)
      colnames(h_dt1) <- c("spec1.chr", "spec1.pos", "spec1.ann")
    }
    # h_dt1$spec1.ann <- sort_cat3(h_dt1$spec1.ann) # sort strings
      # table(h_dt1$spec1.ann)
    
    m_dt1 <- data.table(
      CHR = as.character(df$spec2.chr),
      INT = df$spec2.pos)
    m_dt2 <- data.table(
      CHR = as.character(m_ann$chromosome),
      INT_START = m_ann$start,
      INT_END = m_ann$end,
      CAT = as.character(m_ann$category))
    # str(m_dt1)
    # str(m_dt2)
    # which(m_dt2$INT_START > m_dt2$INT_END)
    m_dt1 <- m_dt2[m_dt1, on=.(CHR=CHR, INT_START<=INT, INT_END>=INT), allow.cartesian = T]
    if (identical(m_dt1$INT_END, m_dt1$INT_START)){
      m_dt1 <- m_dt1[,c("CHR", "INT_END", "CAT")]
      colnames(m_dt1) <- c("CHR", "INT", "CAT")
      m_dt1 <- m_dt1[!duplicated(m_dt1),]
    }
    m_dt1 <- dcast(m_dt1, CHR + INT ~ CAT, value.var = "CAT")
    m_dt1[is.na(m_dt1)] <- ""
    if (ncol(m_dt1) > 3){
      m_dt1 <- data.table(
        spec2.chr = m_dt1$CHR,
        spec2.pos = m_dt1$INT,
        spec2.ann = apply(m_dt1[,3:ncol(m_dt1)], 1, paste0, collapse = ""))
    } else {
      m_dt1 <- as.data.table(m_dt1)
      colnames(m_dt1) <- c("spec2.chr", "spec2.pos", "spec2.ann")
    }
    # m_dt1$spec2.ann <- sort_cat(m_dt1$spec2.ann) # sort strings
    # table(m_dt1$spec2.ann)
    
    # m_dt1 may be smaller then h_dt1 if h_POS align to multiple m_POS
    # merge m_dt1 and h_dt1 with df
    dt_long <- as.data.table(df)
    dt_long <- h_dt1[dt_long, on = c("spec1.chr", "spec1.pos")]
    dt_long <- m_dt1[dt_long, on = c("spec2.chr", "spec2.pos")]
    
    dt_long <- dt_long[,c("spec1.chr", "spec1.pos", "spec1.ref", "spec1.ann",
                        "spec2.chr", "spec2.pos", "spec2.ref", "spec2.ann")]
    
    # convert long to short format
    INT_A_split <- cumsum(c(TRUE, abs(diff(dt_long$spec1.pos))!=1))
    INT_B_split <- cumsum(c(TRUE, abs(diff(dt_long$spec2.pos))!=1))
    CAT1_A_split <- cumsum(c(1, diff(as.numeric(as.factor(dt_long$spec1.chr))) != 0))
    CAT1_B_split <- cumsum(c(1, diff(as.numeric(as.factor(dt_long$spec2.chr))) != 0))
    CAT2_A_split <- cumsum(c(1, diff(as.numeric(as.factor(dt_long$spec1.ann))) != 0))
    CAT2_B_split <- cumsum(c(1, diff(as.numeric(as.factor(dt_long$spec2.ann))) != 0))
    dt_long_ind <- data.table(INT_A_split,
                              INT_B_split,
                              CAT1_A_split,
                              CAT1_B_split,
                              CAT2_A_split,
                              CAT2_B_split)
    dt_long$ID <- apply(dt_long_ind, 1, paste0, collapse = "") # unique ID for changes
    dt_long$ID <- as.factor(dt_long$ID)
    dt_short <- ddply(dt_long, "ID", long_to_short)
    dt_long$ID <- NULL
    dt_short$ID <- NULL
    
    # write tables
    fwrite(dt_long, paste0(output.path, out_long_prefix, chr, ".csv"), col.names = F, row.names = F, append = T)
    fwrite(dt_short, paste0(output.path, out_short_prefix, chr, ".csv"), col.names = F, row.names = F, append = T)
    
    print(j)
  }
}


### SCRIPT

### IMPORT

# annotation
h_ann <- fread(h_ann_file)
m_ann <- fread(m_ann_file)

### FORMAT

# colnames 
colnames(h_ann) <- c("category", "chromosome", "start", "end")
colnames(m_ann) <- c("category", "chromosome", "start", "end")

# categories
h_ann$category <- as.character(h_ann$category)
m_ann$category <- as.character(m_ann$category)
h_ann$category[h_ann$category == "TF binding" | h_ann$category == "Open chromatin"] <- "Miscellaneous"
m_ann$category[m_ann$category == "TF binding" | m_ann$category == "Open chromatin"] <- "Miscellaneous"

# number code annotations (plyr)
h_ann$category <- mapvalues(h_ann$category, from=c("Exon - CDS",
                                                   "Exon - UTR",
                                                   "Exon - other",
                                                   "Promoter",
                                                   "Enhancer - proximal",
                                                   "Enhancer - distal",
                                                   "CTCF binding",
                                                   "Miscellaneous",
                                                   "Intron - proximal",
                                                   "Intron - distal",
                                                   "Unannotated"), 
                   to=c("A","B","C","D","E","F","G","H","I","J","K"))
m_ann$category <- mapvalues(m_ann$category, from=c("Exon - CDS",
                                                   "Exon - UTR",
                                                   "Exon - other",
                                                   "Promoter",
                                                   "Enhancer - proximal",
                                                   "Enhancer - distal",
                                                   "CTCF binding",
                                                   "Miscellaneous",
                                                   "Intron - proximal",
                                                   "Intron - distal",
                                                   "Unannotated"), 
                            to=c("A","B","C","D","E","F","G","H","I","J","K"))

# ensure start <= end 
tmp1 <- subset(h_ann, h_ann$end < h_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(h_ann, h_ann$end >= h_ann$start)
h_ann <- rbind(tmp1, tmp2)
tmp1 <- subset(m_ann, m_ann$end < m_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(m_ann, m_ann$end >= m_ann$start)
m_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)

# set chr
chr <- c(1:22, "X")[chr_id]

# count files
files.chr <- file_count(input.path, dir, chr)


### Convert

MAF_convert(chr = chr,
            files.chr = files.chr,
            input.path = input.path,
            output.path = output.path,
            dir = dir,
            in_prefix = in_prefix,
            out_long_prefix = out_long_prefix,
            out_short_prefix = out_short_prefix)


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
# 
# ### CREATE 
# set.seed(11)
# dt1 <- data.table(
#                   CHR = c(rep("1", 5000), rep("1", 5000)),
#                   INT = sample(1:110000, 10000, replace = F))
# dt2 <- data.table(
#                   CHR = c(rep("1", 500), rep("1", 500)),
#                   INT_START = sample(1:100000, 1000, replace = F),
#                   INT_END = sample(1:100000, 1000, replace = F),
#                   CAT = sample(c("A", "B", "C", "D", "E"), 1000, replace = T))
# 
# ### FORMAT
# # ensure start <= end 
# tmp1 <- subset(dt2, dt2$INT_END < dt2$INT_START)
# colnames(tmp1) <- c("CHR", "INT_END", "INT_START", "CAT")
# tmp2 <- subset(dt2, dt2$INT_END >= dt2$INT_START)
# dt2 <- rbind(tmp1, tmp2)
# rm(tmp1, tmp2)
# 
# dt1 <- dt1[order(INT),]
# dt2 <- dt2[order(CAT, INT_START),]
# 
#
# dt1$CAT <- NA
# for_test <- rep(NA, nrow(dt1))
# for (i in 1:nrow(dt1)){
#   ann.ind <- which(dt2$INT_START <= dt1$INT[i] & dt2$INT_END >= dt1$INT[i])
#   # dt1$CAT[i] <- paste0(unique(dt2$CAT[ann.ind]), collapse = "")
#   for_test[i] <- paste0(unique(dt2$CAT[ann.ind]), collapse = "")
#   print(i)
# }

# dt1 <- dt2[dt1, on=.(CHR=CHR, INT_START<=INT, INT_END>=INT), allow.cartesian = T]
# if (identical(dt1$INT_END, dt1$INT_START)){
#   dt1 <- dt1[,c("CHR", "INT_END", "CAT")]
#   colnames(dt1) <- c("CHR", "INT", "CAT")
#   dt1 <- dt1[!duplicated(dt1),]
# }
# dt1 <- dcast(dt1, INT + CHR ~ CAT, value.var = "CAT")
# dt1[is.na(dt1)] <- ""
# dt1 <- data.table(INT = dt1[,1],
#                   CHR = dt1[,2],
#                     CAT = apply(dt1[,3:ncol(dt1)], 1, paste0, collapse = ""))
# 
# table(for_test)
# table(dt1$CAT)
# 
# 
# dt1 <- data.table(
#   INT = df$spec1.pos)
# 
# dt2 <- subset(h_ann, h_ann$chromosome == "1")
# dt2 <- data.table(
#   INT_START = dt2$start,
#   INT_END = dt2$end,
#   CAT = dt2$category)
# 
# str(dt1)
# str(dt2)
# which(dt2$INT_START > dt2$INT_END)
# 
# dt1 <- dt2[dt1, on=.(INT_START<=INT, INT_END>=INT), allow.cartesian = T]
# if (identical(dt1$INT_END, dt1$INT_START)){
#   dt1 <- dt1[,c("INT_END", "CAT")]
#   colnames(dt1) <- c("INT", "CAT")
#   dt1 <- dt1[!duplicated(dt1),]
# }
# dt1 <- dcast(dt1, INT ~ CAT, value.var = "CAT")
# dt1[is.na(dt1)] <- ""
# dt1 <- data.table(INT = dt1[,1],
#                   CAT = apply(dt1[,2:ncol(dt1)], 1, paste0, collapse = ""))
# # table(dt1$CAT)




# ### for loop
# # sort by m_chr
# df <- df[order(spec2.chr),]
# df <- df[complete.cases(df),]
# # subset df
# h_df <- df[,1:3]
# # subset h_ann by h_chr
# h_sub <- subset(h_ann, h_ann$chromosome == h_df$spec1.chr[1])
# # annotation for each position
# h_df$spec1.ann <- NA
# for (i in 1:nrow(h_df)){
#   ann.ind <- which(h_sub$start <= h_df$spec1.pos[i] & h_sub$end >= h_df$spec1.pos[i])
#   h_df$spec1.ann[i] <- paste0(unique(sort(h_sub$category[ann.ind])), collapse = "")
#   print(i)
# }
# # table(h_df$spec1.ann)
# 
# # subset df 
# m_df <- df[,4:6]
# m_df_ann <- data.frame()
# for (chr_i in unique(m_df$spec2.chr)){
#   # subset m_df by chr
#   m_chr <- subset(m_df, m_df$spec2.chr == chr_i)
#   # subset m_ann by chr
#   m_sub <- subset(m_ann, m_ann$chromosome == m_chr$spec2.chr[1])
#   # annotation for each position
#   m_chr$spec2.ann <- NA
#   for (i in 1:nrow(m_chr)){
#     ann.ind <- which(m_sub$start <= m_chr$spec2.pos[i] & m_sub$end >= m_chr$spec2.pos[i])
#     m_chr$spec2.ann[i] <- paste0(unique(sort(m_sub$category[ann.ind])), collapse = "")
#     print(i)
#   }
#   m_df_ann <- rbind(m_df_ann, m_chr)
# }
# # table(h_df$spec1.ann)
# # table(m_df_ann$spec2.ann)



# INT_A = c(1:10, 16:30, 41:70, 72:116)
# CAT1_A = rep("X", 100)
# CAT2_A = c(rep("A", 10), rep("C", 5), rep("B", 10), rep("A", 20), rep("B", 55))
# INT_B = c(10:1, 20:6, 41:115)
# CAT1_B = c(rep("Y", 10), rep("X", 15), rep("Y", 75))
# CAT2_B = c(rep("A", 10), rep("C", 15), rep("B", 5), rep("A", 25), rep("B", 45))
# 
# dt_long <- data.table(INT_A,
#                       INT_B,
#                       CAT1_A,
#                       CAT1_B,
#                       CAT2_A,
#                       CAT2_B)
# 
# 
# INT_A_split <- cumsum(c(TRUE, abs(diff(INT_A))!=1))
# INT_B_split <- cumsum(c(TRUE, abs(diff(INT_B))!=1))
# CAT1_A_split <- cumsum(c(1, diff(as.numeric(as.factor(CAT1_A))) != 0))
# CAT1_B_split <- cumsum(c(1, diff(as.numeric(as.factor(CAT1_B))) != 0))
# CAT2_A_split <- cumsum(c(1, diff(as.numeric(as.factor(CAT2_A))) != 0))
# CAT2_B_split <- cumsum(c(1, diff(as.numeric(as.factor(CAT2_B))) != 0))
# 
# dt_long_ind <- data.table(INT_A_split,
#                       INT_B_split,
#                       CAT1_A_split,
#                       CAT1_B_split,
#                       CAT2_A_split,
#                       CAT2_B_split)
# dt_long$ID <- apply(dt_long_ind, 1, paste0, collapse = "")
# dt_long$ID <- as.factor(dt_long$ID)
# 
# # sub <- subset(dt_long, dt_long$ID == "331344331344")
# long_to_short <- function(sub){
#   data.table(
#     INT_A_START = sub$INT_A[1],
#     INT_A_END = sub$INT_A[nrow(sub)],
#     CAT1_A = sub$CAT1_A[1],
#     CAT2_A =  sub$CAT2_A[1],
#     INT_B_START = sub$INT_B[1],
#     INT_B_END = sub$INT_B[nrow(sub)],
#     CAT1_B = sub$CAT1_B[1],
#     CAT2_B =  sub$CAT2_B[1]
#   )
# }
# 
# out <- ddply(dt_long, "ID", long_to_short)
# out$ID <- NULL


