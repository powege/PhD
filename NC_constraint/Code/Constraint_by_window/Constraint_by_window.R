rm(list = ls())
graphics.off()

library(data.table)
library(MASS)
library(plyr)
# library(dplyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied", call.=FALSE)
# } 

# set args
# contraint_variables_path <- args[1]
# contraint_variables_window_size <- args[2]
# output_file <- args[3]
# species <- args[4]
contraint_variables_path <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/mouse_constraint_variables_by_window_chr"
contraint_variables_window_size <- "_650_50.csv"
output_file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window"
species <- "mouse"

if (species == "mouse") { chr_vec <- c(1:19)}
if (species == "human") { chr_vec <- c(1:22)}


### IMPORT

# all chr
dt.list <- list()
for (chr in chr_vec){
  dt.list[[chr]] <- fread(paste0(contraint_variables_path, chr, contraint_variables_window_size))
}
dt <- rbind.fill(dt.list)
rm(dt.list)


### QC

dt <- dt[complete.cases(dt),]
len_dt_all <- nrow(dt)
if (species == "human"){
dt <- subset(
            dt, 
            dt$n_mask <= (650*0.8)
            &
            dt$n_mask_central <= (50*0.2)
             )
}
if (species == "mouse"){
dt <- subset(
            dt, 
            dt$n_mask <= (650*0.8)
            &
            dt$n_mask_central <= (50*0.2)
            )
}

### CONSTRAINT

# Calculate constraint with kmer
mod_kmer <- lm(
  dt$n_SNV_weighted ~ dt$p_SNV_given_kmers_weighted
          )
summary(mod_kmer)
dt$Constraint_score_kmer <- studres(mod_kmer)


# Calculate constraint with CpG
mod_CpG <- lm(
  dt$n_SNV_weighted ~ dt$n_CpG_weighted
)
summary(mod_CpG)
dt$Constraint_score_CpG <- studres(mod_CpG)

# percentile rank
percentile_rank <- function(x) ceiling((rank(x)/length(x))*100) 
dt$Constraint_percentile_kmer <- percentile_rank(dt$Constraint_score_kmer)
dt$Constraint_percentile_CpG <- percentile_rank(dt$Constraint_score_CpG)

### OUTPUT
# dt <- dt[, c("CHR", "POS_from", "POS_to", "Constraint_score", "Constraint_percentile")]
fwrite(dt, paste0(output_file, species, "_", contraint_variables_window_size))

#####

### TEST

# # subset top percentiles
# p1 <- subset(dt, dt$Constraint_percentile_CpG == 1)
# p2 <- subset(dt, dt$Constraint_percentile_CpG == 2)
# p3 <- subset(dt, dt$Constraint_percentile_CpG == 3)
# p98 <- subset(dt, dt$Constraint_percentile_CpG == 98)
# p99 <- subset(dt, dt$Constraint_percentile_CpG == 99)
# p100 <- subset(dt, dt$Constraint_percentile_CpG == 100)
# 
# hist(dt$n_SNV, breaks = 700)
# hist(dt$n_SNV_weighted, breaks = 700)
# 
# # plot variables for top percentiles
# hist(p1$n_SNV_weighted, breaks = 700)
# hist(p2$n_SNV_weighted, breaks = 700)
# hist(p3$n_SNV_weighted, breaks = 700)
# hist(p98$n_SNV_weighted, breaks = 700)
# hist(p99$n_SNV_weighted, breaks = 700)
# hist(p100$n_SNV_weighted, breaks = 700)
# 
# hist(p1$p_SNV_given_kmers, breaks = 700)
# hist(p2$p_SNV_given_kmers, breaks = 700)
# hist(p3$p_SNV_given_kmers, breaks = 700)
# hist(p98$p_SNV_given_kmers, breaks = 700)
# hist(p99$p_SNV_given_kmers, breaks = 700)
# hist(p100$p_SNV_given_kmers, breaks = 700)
# 
# hist(p1$n_mask, breaks = 700)
# hist(p2$n_mask, breaks = 700)
# hist(p3$n_mask, breaks = 700)
# hist(p98$n_mask, breaks = 700)
# hist(p99$n_mask, breaks = 700)
# hist(p100$n_mask, breaks = 700)



# # get POS_ID for percentiles
# get_ID_all <- function(sub){
#   all.POS <- as.character(rep(sub$POS_from, each = 100) + 0:99)
#   all.CHR <- rep(sub$CHR, each = 100)
#   all.ID <- paste0(all.CHR, "_", all.POS)
#   return(all.ID)
# }
# 
# chromo <- 2
# ID.list <- list()
# for(i in c(1:5, 96:100)){
#   sub <- subset(dt, dt$CHR == chromo)
#   sub <- subset(sub, sub$Constraint_percentile_CpG == i)
#   ID.list[[i]] <- as.data.table(get_ID_all(sub))
#   print(i)
# }
# 
# # get POS_ID for CDS
# gene.pos <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv")
# table(gene.pos$category)
# gene.pos <- subset(gene.pos, gene.pos$category == "Exon - CDS")
# sub <- subset(gene.pos, gene.pos$chromosome == chromo)
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# all.POS <- unlist(seq2(from = sub$start, to = sub$end))
# all.CHR <- rep(chromo, length(all.POS))
# CDS_out <- paste0(all.CHR, "_", all.POS)
# 
# # get fraction of CDS in each percentile
# prop <- rep(NA, 100)
# for (i in c(1:5, 96:100)){
#   p <- ID.list[[i]]
#   prop[i] <- table(p$V1 %in% CDS_out)[2]/length(p$V1)
#   print(i)
# }
# 
# # plot
# percentile <- c(1:100)
# df.plot <- data.frame(Percentile = percentile, CDS_fraction = prop)
# plot(df.plot$CDS_fraction~df.plot$Percentile)
