### Script that identifies all bp not annotated in GenCode or the Ensembl Regulatory Build. 

rm(list=ls())
graphics.off()

library(data.table)
library(dplyr)
library(purrr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)!=3) {
  stop("Exactly three arguments must be supplied", call.=FALSE)
} 

# set args variables
annotation.csv <- args[1] 
out_file <- args[2] 
species <- args[3] 
# annotation.csv <- "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/Human_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/Human_GENCODE_RegBuild_unannotated.csv"
# species <- "homo_sapien"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### IMPORT
# annotation
annotation <- fread(annotation.csv)

# GRCh38 chromosome lengths
H_chr_length <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                  145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                  101991189,  90338345,  83257441,  80373285,  58617616,  64444167,  46709983,
                  50818468)
M_chr_length <- c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 
                  129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 
                  104043685, 98207768, 94987271, 90702639, 61431566)
if (species == "homo_sapien"){ chr_length <- H_chr_length }
if (species == "mus_musculus"){ chr_length <- M_chr_length }


### FORMAT unannotated POS

out_list <- list()
for(chr in 1:length(chr_length)){
  
  # subset annotation by chr
  tmp <- subset(annotation, annotation$chromosome == chr)
  tmp <- tmp[,c("start", "end")]
  colnames(tmp) <- c("from", "to")
  
  # vector of POS not in tmp
  unique_pos <- setdiff(x = 1:chr_length[chr], 
                        y = unlist(seq2(from = tmp$from, to = tmp$to)))
  # unique_pos <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
  
  # format vector to dataframe
  unique_pos <- sort(unique_pos, decreasing = F)
  unique_pos <- t(sapply(split(unique_pos, findInterval(unique_pos, unique_pos[which(c(1, diff(unique_pos)) > 1)])), range))
  unique_pos <- as.data.table(unique_pos)
  colnames(unique_pos) <- c("start", "end")
  unique_pos$chromosome <- chr
  unique_pos$category <- "Unannotated"
  
  out_list[[chr]] <- unique_pos
  
  # see: https://stackoverflow.com/questions/55711684/how-to-identify-all-sequential-numbers-not-covered-by-to-and-from-positions/55712517#55712517
  #   miss <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
  # i <- 
  #   miss %>% 
  #   diff %>% 
  #   `>`(1) %>% 
  #   rev %>%
  #   cumsum %>% 
  #   rev 
  # out <- map_df(split(miss, c(i, 0)), ~list(from = head(.x, 1), to = tail(.x, 1))) %>% 
  #   dplyr::arrange(from)
  # 
  # colnames(out) <- c("start", "end")
  # out$chromosome <- chr
  # out$category <- "Unannotated"
  # out_list[[chr]] <- out
}

unannotated <- do.call("rbind", out_list)
unannotated <- unannotated[,c("category", "chromosome", "start", "end")]


### EXPORT
fwrite(unannotated, out_file)


########

# df1 <- data.frame(from = c(7, 22, 35, 21, 50),
#                   to = c(13, 29, 43, 31, 60))
# df2 <- data.frame(from = c(1, 14, 32, 44, 61),
#                   to = c(6, 20, 34, 49, 100))
# 
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# seq <- c(min(df1$from):max(df1$to))
# df1_POS <- unlist(seq2(from = df1$from, to = df1$to))
# df1_POS <- unique(df1_POS)
# system.time({ df2_POS <- seq[!seq %in% df1_POS] })
# system.time({ df2_POS <- setdiff(seq, df1_POS) })
# # split vector into list of sequences #### ASK ON STACK OVERFLOW! 
# # calculate differnece between each element in vector
# df2_POS <- df2_POS[order(df2_POS)]
# dif <- c(NA, df2_POS[2:(length(df2_POS))] - df2_POS[1:(length(df2_POS)-1)])
# # identify indicies where b == a+1 
# a <- 1
# z <- 100
# to_ind <- which(dif != 1) -1
# to <- c(min(df1$from)-1, df2_POS[to_ind], max(df1$from)-1, z)
# from_ind <- which(dif != 1)
# from <- c(a, min(df1$to)+1, df2_POS[from_ind], max(df1$to)+1)
# out <- data.frame(from = from,
#                   to = to)

### Get REF lengths
## using R
## Human
# chr <- c(1:22)
# length <- rep(NA, length(chr))
# for (i in 1:length(chr)){
#   length[i] <- nrow(fread(paste0("/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr", i, ".txt"), select = 1L))
# }
## using bash
## Human
# for CHR in {1..22}
# do
# wc -l /well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr"$CHR".txt
# done
## Mouse
# for CHR in {1..19}
# do
# wc -l /well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$CHR".txt
# done          


