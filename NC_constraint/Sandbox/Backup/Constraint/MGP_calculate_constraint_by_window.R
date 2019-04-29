rm(list = ls())
graphics.off()

library(data.table)
library(MASS)
library(plyr)
# library(dplyr)

# set variables
contraint_variables_path <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/mouse_constraint_variables_by_window_chr"
contraint_variables_window_size <- "_650_50.csv"
output_file_path <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_by_window"

### import

# chr
chrom = 1
dt <- fread(paste0(contraint_variables_path, chrom, contraint_variables_window_size))

# all chr
dt.list <- list()
for (chr in 1:19){
  dt.list[[chr]] <- fread(paste0(contraint_variables_path, chr, contraint_variables_window_size))
}
dt <- rbind.fill(dt.list)
rm(dt.list)

### QC

len_dt_all <- nrow(dt)
dt <- dt[complete.cases(dt),]
# subset windows with n_mask >= 0.5
dt <- subset(
            dt, dt$n_mask >= 0.5
             )

# Calculate constraint with kmer
mod_kmer <- lm(
  dt$n_SNV ~ dt$p_SNV_given_kmers
          )
summary(mod_kmer)
dt$Constraint_score_kmer <- studres(mod_kmer)


# Calculate constraint with CpG
mod_CpG <- lm(
  dt$n_SNV ~ dt$n_CpG
)
summary(mod_CpG)
dt$Constraint_score_CpG <- studres(mod_CpG)

# percentile rank
percentile_rank <- function(x) ceiling((rank(x)/length(x))*100) 
dt$Constraint_percentile_kmer <- percentile_rank(dt$Constraint_score_kmer)
dt$Constraint_percentile_CpG <- percentile_rank(dt$Constraint_score_CpG)

# Output
# dt <- dt[, c("CHR", "POS_from", "POS_to", "Constraint_score", "Constraint_percentile")]
fwrite(dt, paste0(output_file_path, contraint_variables_window_size))

#####

### TEST 

# # subset top percentiles
# p1 <- subset(dt, dt$Constraint_percentile == 1)
# p2 <- subset(dt, dt$Constraint_percentile == 2)
# p3 <- subset(dt, dt$Constraint_percentile == 3)
# p98 <- subset(dt, dt$Constraint_percentile == 98)
# p99 <- subset(dt, dt$Constraint_percentile == 99)
# p100 <- subset(dt, dt$Constraint_percentile == 100)
# 
# # plot variables for top percentiles
# hist(p1$n_SNV, breaks = 700)
# hist(p2$n_SNV, breaks = 700)
# hist(p3$n_SNV, breaks = 700)
# hist(p98$n_SNV, breaks = 700)
# hist(p99$n_SNV, breaks = 700)
# hist(p100$n_SNV, breaks = 700)
# 
# hist(p1$p_SNV_given_kmers, breaks = 700)
# hist(p2$p_SNV_given_kmers, breaks = 700)
# hist(p3$p_SNV_given_kmers, breaks = 700)
# hist(p98$p_SNV_given_kmers, breaks = 700)
# hist(p99$p_SNV_given_kmers, breaks = 700)
# hist(p100$p_SNV_given_kmers, breaks = 700)
# 
# hist(p1$Read_depth, breaks = 700)
# hist(p2$Read_depth, breaks = 700)
# hist(p3$Read_depth, breaks = 700)
# hist(p98$Read_depth, breaks = 700)
# hist(p99$Read_depth, breaks = 700)
# hist(p100$Read_depth, breaks = 700)
# 
# hist(p1$Repeats, breaks = 700)
# hist(p2$Repeats, breaks = 700)
# hist(p3$Repeats, breaks = 700)
# hist(p1$Repeats, breaks = 700)
# hist(p2$Repeats, breaks = 700)
# hist(p3$Repeats, breaks = 700)

# # get POS_ID for percentiles
# get_ID_all <- function(sub){
#   all.POS <- as.character(rep(sub$POS_from, each = 100) + 0:99)
#   all.CHR <- rep(sub$CHR, each = 100)
#   all.ID <- paste0(all.CHR, "_", all.POS)
#   return(all.ID)
# }
# ID.list <- list()
# for(i in c(1:5, 96:100)){
#   sub <- subset(dt, dt$Constraint_percentile == i)
#   ID.list[[i]] <- as.data.table(get_ID_all(sub))
#   print(i)
# }
# 
# # get POS_ID for CDS
# gene.pos <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/M_exon_POS.csv")
# sub <- subset(gene.pos, gene.pos$chromosome_name == chrom)
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# all.POS <- unlist(seq2(from = sub$exon_chrom_start, to = sub$exon_chrom_end))
# all.CHR <- rep(chrom, length(all.POS))
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
