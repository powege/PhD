rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)!=4) {
  stop("Exactly four arguments must be supplied", call.=FALSE)
} 

# set args variables
annotated_file <- args[1] 
unannotated_file <- args[2] 
out_file <- args[3] 
species <- args[4]
# annotated_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv"
# unannotated_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_unannotated.csv"

### Set species chromosomes
if (species == "mouse"){ species_chr <- c(1:19)}
if (species == "human"){ species_chr <- c(1:22)}


### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT
annotated <- fread(annotated_file)
unannotated <- fread(unannotated_file)

### FORMAT 
unranked <- rbind(annotated, unannotated)
rm(annotated, unannotated)
unranked$start <- as.integer(unranked$start)
unranked$end <- as.integer(unranked$end)

### get total bases in each cat
# unranked$length <- abs(unranked$end - unranked$start)
# cats <- unique(unranked$category)
# total <- rep(NA, length(cats))
# for(cat in 1:length(cats)){
#   sub <- subset(unranked, unranked$category == cats[cat])
#   total[cat] <- sum(sub$length)
# }
# totals <- data.frame(cat = cats, total = total)
# totals <- totals[order(totals$total),]

### LOOP through chromosomes

# Set rank and output list
rank <- c(
          # "Exon - CDS", # remove for for loop 
          "Exon - UTR", 
          "Exon - non-coding", 
          "Promoter", 
          "Enhancer",
          "Open chromatin",
          "TF binding",
          "Promoter flanking", 
          "Intron",
          "Unannotated"
          )
out_list <- list()

for(chr in species_chr){
  
  unranked_chr <- subset(unranked, unranked$chromosome == chr)
  ranked_chr <- subset(unranked_chr, unranked_chr$category == "Exon - CDS")

    # loop through annotation categories
    for(cat in rank){
  
      # subset category
      cat_sub <- subset(unranked_chr, unranked_chr$category == cat)

      # vector of POS in cat_sub, not in ranked df
      unique_pos <- setdiff(unlist(seq2(from = cat_sub$start, to = cat_sub$end)), 
                            unlist(seq2(from = ranked_chr$start, to = ranked_chr$end)))
  
      # format vector to dataframe
      unique_pos <- sort(unique_pos, decreasing = F)
      unique_pos <- t(sapply(split(unique_pos, findInterval(unique_pos, unique_pos[which(c(1, diff(unique_pos)) > 1)])), range))
      unique_pos <- as.data.table(unique_pos)
      colnames(unique_pos) <- c("start", "end")
      unique_pos$chromosome <- chr
      unique_pos$category <- cat

      # add non-duplicated pos to df
      ranked_chr <- rbind(ranked_chr, unique_pos) 
      print(paste0(cat, " added!"))
}

  out_list[[chr]] <- ranked_chr
  print(paste0(chr, " done!"))

}

# merge to one dataframe
ranked_out <- do.call("rbind", out_list)

### EXPORT
fwrite(ranked_out, out_file)




#####

### Stack Overflow example code

# df1 <- data.frame(from = c(1, 18, 53, 65),
#                   to = c(10, 23, 47, 70))
# 
# df2 <- data.frame(from = c(7, 29, 35, 21, 50),
#                   to = c(13, 22, 43, 31, 60))
# 
# # all integers in df2 not in df1
# df3 <- data_frame(from = c(11, 24, 35, 54),
#                   to = c(13, 31, 43, 60))

