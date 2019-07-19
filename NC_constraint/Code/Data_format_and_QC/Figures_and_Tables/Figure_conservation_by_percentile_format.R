### SCRIPT that gets the fraction of gerp and human-mouse alignment for each constraint percentile. 

rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

perc <- as.integer(args[1])
gerp.file <- args[2]
align.file <- args[3]
constraint.file <- args[4]
out.file <- args[5]
species <- args[6]

# perc <- 1
# gerp.file <- "~/Dropbox/PhD/Data/Ensembl/GERP/gerp_constrained_elements.homo_sapiens.bed"
# align.file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt"
# constraint.file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_human_950_50_MAF001.csv"
# # out.file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/"
# species <- "human"


### FUNCTIONS

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT

# conservation data
align <- fread(align.file)
gerp <- fread(gerp.file)
# constraint data
con <- fread(constraint.file)


### FORMAT 

# Set chromosome number
if (species == "human"){ n_chr <- 22}
if (species == "mouse"){ n_chr <- 19}

# Subset autosomes
align <- subset(align, align$V1 %in% as.character(1:n_chr))
gerp <- subset(gerp, gerp$V1 %in% as.character(1:n_chr))

# Subset constraint percentile
con <- subset(con, con$Constraint_percentile_CpG == perc)

# Set output list
out_list <- list()


  # Set vectors
  percentile <- rep(perc, n_chr)
  chromosome <- 1:n_chr
  total_percentile_POS <- rep(NA, n_chr)
  total_GERP_POS <- rep(NA, n_chr) 
  total_aligned_POS <- rep(NA, n_chr) 
  
  
  # For each chromosome:
  for (chr in chromosome){
    
    # Subset gerp chromosome
    gerp.chr.sub <- subset(gerp, gerp$V1 == chr)
    # Subset alignment chromosome
    align.chr.sub <- subset(align, align$V1 == chr)
    # Subset constraint percentile chromsome
    con.chr.sub <- subset(con, con$CHR == chr)
    # Get total percentile chromosome POS
    total_percentile_POS[chr] <- sum((con.chr.sub$POS_to + 1) - con.chr.sub$POS_from)
    # Vectorise all percentile chromosome POS
    per_chr_pos <- unlist(seq2(from = con.chr.sub$POS_from, to = con.chr.sub$POS_to))
    # Vectorise all gerp chromosome POS
    gerp_chr_pos <- unlist(seq2(from = gerp.chr.sub$V2, to = gerp.chr.sub$V3))
    # Vectorise all aligned chromosome POS
    align_chr_pos <- unlist(seq2(from = align.chr.sub$V2, to = align.chr.sub$V3))
    # Identify total gerp POS in percentile
    total_GERP_POS[chr] <- total_percentile_POS[chr] - length(setdiff(per_chr_pos, gerp_chr_pos))
    # Identify total aligned POS in percentile
    total_aligned_POS[chr] <- total_percentile_POS[chr] - length(setdiff(per_chr_pos, align_chr_pos))
    print(chr)
  }
  
  output <- data.frame(percentile = c(percentile, percentile[1]), 
                    chromosome = c(chromosome, "total"),
                    total_percentile_POS = c(total_percentile_POS, sum(total_percentile_POS, na.rm = T)),
                    total_GERP_POS = c(total_GERP_POS, sum(total_GERP_POS, na.rm = T)),
                    total_aligned_POS = c(total_aligned_POS, sum(total_aligned_POS, na.rm = T))
  )
  


### OUTPUT
fwrite(output, out.file, append = T)






####################
# species <- args[1]
# gerp.file <- args[2]
# alogn.file <- args[3]
# out.file <- args[4]
# species <- args[3]
# annotation <- args[4]
# bed.file <- "~/Dropbox/PhD/Data/Ensembl/GERP/gerp_constrained_elements.homo_sapiens.bed"
# out.file <- "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_GERP_by_percentile_650_50.csv"
# species <- "human"
# 
# 
# ### FUNCTIONS
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# 
# ### IMPORT
# gerp <- fread(bed.file, skip = 1)
# 
# ### FORMAT
# colnames(gerp) <- c("chromosome", "start", "end") 
# 
# # get vector of all GERP POS_ID 
# out <- list()
# if (species == "human") { chr <- c(1:22) }
# if (species == "mouse") { chr <- c(1:19) }
# for(i in chr){
#   sub <- subset(gerp, gerp$chromosome == chr[i])
#   all.POS <- unlist(seq2(from = sub$start, to = sub$end))
#   all.CHR <- rep(chr[i], length(all.POS))
#   out[[i]] <- paste0(all.CHR, "_", all.POS)
#   print(i)
# }
# out <- unique(unlist(out))
# 
# # get proportion of annotation POS in each percentile
# prop <- rep(NA, 100)
# percentile <- c(1:100)
# for (i in percentile){
#   # p <- fread("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/human_POS_ID_percentile_1_650_50.csv", header = F)
#   p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", i, "_650_50.csv"), header = F)
#   # prop[i] <- table(p$V1 %in% out)[2]/length(p$V1)
#   prop[i] <- sum(p$V1 %chin% out)/length(p$V1)
#     print(i)
# }
# 
# df_out <- data.frame(Percentile = percentile, 
#                      Fraction = prop,
#                      Annotation = rep(annotation, 100))
# 
# ### EXPORT
# fwrite(df_out, out.file, col.names = F)