### Script that formats fraction of genomic annotations by percentile

rm(list=ls())
graphics.off()

library(data.table)

### ARGS 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)==0) {
  stop("Args required", call.=FALSE)
} 

perc <- as.integer(args[1])
annotation.file <- args[2]
ranked.annotation.file <- args[3]
constraint.file <- args[4]
out.file <- args[5]
species <- args[6]
# perc <- 1
# annotation.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv"
# ranked.annotation.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
# constraint.file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_human_650_50.csv"
# out.file <- "~/Dropbox/"
# species <- "human"

### FUNCTIONS

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT

# annotation data
ann <- fread(annotation.file)
rank.ann <- fread(ranked.annotation.file)
# constraint data
con <- fread(constraint.file)


### FORMAT 

# Subset constraint percentile
con <- subset(con, con$Constraint_percentile_CpG == perc)

# Combine annotation files
ann <- rbind(subset(ann, ann$category != "Intron"), 
              subset(rank.ann, rank.ann$category == "Intron" | rank.ann$category == "Unannotated"))
rm(rank.ann)

# Set output list
out_list <- list()

# Set chromosome number
if (species == "human"){ n_chr <- 22}
if (species == "mouse"){ n_chr <- 19}

# For each annotation:
for (a in 1:length(unique(ann$category))){

 # Subset annotation
  ann.sub <- subset(ann, ann$category == unique(ann$category)[a])

 # Set vectors
  percentile <- rep(perc, n_chr)
  annotation <- rep(unique(ann$category)[a], n_chr)
  chromosome <- 1:n_chr
  total_percentile_POS <- rep(NA, n_chr)
  total_annotation_POS <- rep(NA, n_chr) 
  

  # For each chromosome:
  for (chr in chromosome){

    # Subset annotation chromosome
    ann.chr.sub <- subset(ann.sub, ann.sub$chromosome == chr)
    # Subset constraint percentile chromsome
    con.chr.sub <- subset(con, con$CHR == chr)
    # Get total percentile chromosome POS
    total_percentile_POS[chr] <- sum((con.chr.sub$POS_to + 1) - con.chr.sub$POS_from)
    # Vectorise all percentile chromosome POS
    per_chr_pos <- unlist(seq2(from = con.chr.sub$POS_from, to = con.chr.sub$POS_to))
    # Vectorise all annotation chromosome POS
    ann_chr_pos <- unlist(seq2(from = ann.chr.sub$start, to = ann.chr.sub$end))
    # Identify total annotation POS in percentile
    total_annotation_POS[chr] <- total_percentile_POS[chr] - length(setdiff(per_chr_pos, ann_chr_pos))
    print(chr)
  }

  out_list[[a]] <- data.frame(percentile = c(percentile, percentile[1]), 
                       annotation = c(annotation, annotation[1]),
                       chromosome = c(chromosome, "total"),
                       total_percentile_POS = c(total_percentile_POS, sum(total_percentile_POS, na.rm = T)),
                       total_annotation_POS = c(total_annotation_POS, sum(total_annotation_POS, na.rm = T))
                       )

  print(unique(ann$category)[a])
}


### OUTPUT

# Append output: Percentile; Annotation; Chromosome; Total_Annotation; Total_POS
output <- do.call("rbind", out_list)
fwrite(output, out.file, append = T)

