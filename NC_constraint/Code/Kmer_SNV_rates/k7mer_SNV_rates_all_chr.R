# SCRIPT that combines chromosome kmer counts, and calculates autosome wide 
# 7mer specific probabilities of SNV 

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Three argumants must be supplied: input file path, N chromosomes, output file path", call.=FALSE)
} 

# set args
in_file <- args[1]
n_chr <- args[2]
out_file <- args[3]
# in_file <- "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_chr"
# n_chr <- 19
# out_file <- "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table"


### IMPORT chromosome counts into list
dt_list <- list()
for (chr in 1:n_chr){
  dt_list[[chr]] <- fread(paste0(in_file, chr, ".table"))
  }

# sum rows across list
k7_from <- dt_list[[1]]$k7_from
k7_from_N <- rowSums(sapply(dt_list, `[[`, 2), na.rm = TRUE)
k7_to <- dt_list[[1]]$k7_to
k7_mu_N <- rowSums(sapply(dt_list, `[[`, 4), na.rm = TRUE)

# output datatable
out <- data.table(k7_from = k7_from,
                  k7_from_N = k7_from_N,
                  k7_to = k7_to,
                  k7_mu_N = k7_mu_N)

# calculate p_any_snp_given_k7 and k7_mu_rr
p_any_snp_given_k7 <- function(sub){
  sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
}
tmp <- ddply(out, "k7_from", p_any_snp_given_k7)
colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
out <- merge(out, tmp, all = T)
rm(tmp)
out$p_snp_given_k7 <- out$k7_mu_N/out$k7_from_N
out$k7_mu_rr <- out$k7_mu_N/out$k7_from_N

### EXPORT
fwrite(out, out_file)


############

# ### PLOT
# 
# plot(H_um_out$k7_mu_rr, H_out$k7_mu_rr)
# cor.test(H_um_out$k7_mu_rr, H_out$k7_mu_rr)
# 
# plot(M_um_out$k7_mu_rr, M_out$k7_mu_rr)
# cor.test(M_um_out$k7_mu_rr, M_out$k7_mu_rr)
# 
# plot(M_um_out$k7_mu_rr, H_um_out$k7_mu_rr)
# cor.test(M_um_out$k7_mu_rr, H_um_out$k7_mu_rr)
# 
# plot(M_out$k7_mu_rr, H_out$k7_mu_rr)
# cor.test(M_out$k7_mu_rr, H_out$k7_mu_rr)
# 
# 
# plot_dt <- data.table(label = paste(H_out$k7_from, "to", H_out$k7_to, sep = " "),
#                       H_umREF_k7_mu = H_um_out$k7_mu_rr,
#                       H_REF_k7_mu = H_out$k7_mu_rr,
#                       M_umREF_k7_mu = M_um_out$k7_mu_rr,
#                       M_REF_k7_mu = M_out$k7_mu_rr
#                       )
# library("PerformanceAnalytics")
# chart.Correlation(plot_dt[,c(2:5)], histogram=TRUE, pch=12, method = "spearman")


