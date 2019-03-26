# SCRIPT that combines chromosome kmer counts, and calculates autosome wide 
# 7mer specific probabilities of SNV 

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### ALL MOUSE ###

### IMPORT chromosome counts into list
M_dt_list <- list()
for (chr in 1:19){
  M_dt_list[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_chr", chr, ".table"))
  }

# sum rows across list
k7_from <- M_dt_list[[1]]$k7_from
k7_from_N <- rowSums(sapply(M_dt_list, `[[`, 2), na.rm = TRUE)
k7_to <- M_dt_list[[1]]$k7_to
k7_mu_N <- rowSums(sapply(M_dt_list, `[[`, 4), na.rm = TRUE)

# output datatable
M_out <- data.table(k7_from = k7_from,
                  k7_from_N = k7_from_N,
                  k7_to = k7_to,
                  k7_mu_N = k7_mu_N)

# calculate p_any_snp_given_k7 and k7_mu_rr
p_any_snp_given_k7 <- function(sub){
  sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
}
tmp <- ddply(M_out, "k7_from", p_any_snp_given_k7)
colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
M_out <- merge(M_out, tmp, all = T)
rm(tmp)
M_out$p_snp_given_k7 <- M_out$k7_mu_N/M_out$k7_from_N
M_out$k7_mu_rr <- M_out$k7_mu_N/M_out$k7_from_N

### EXPORT
fwrite(M_out, "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table")

### UNMASKED MOUSE ###

### IMPORT chromosome counts into list
M_um_dt_list <- list()
for (chr in 1:19){
  M_um_dt_list[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_unmasked_chr", chr, ".table"))
}

# sum rows across list
k7_from <- M_um_dt_list[[1]]$k7_from
k7_from_N <- rowSums(sapply(M_um_dt_list, `[[`, 2), na.rm = TRUE)
k7_to <- M_um_dt_list[[1]]$k7_to
k7_mu_N <- rowSums(sapply(M_um_dt_list, `[[`, 4), na.rm = TRUE)

# output datatable
M_um_out <- data.table(k7_from = k7_from,
                  k7_from_N = k7_from_N,
                  k7_to = k7_to,
                  k7_mu_N = k7_mu_N)

# calculate p_any_snp_given_k7 and k7_mu_rr
p_any_snp_given_k7 <- function(sub){
  sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
}
tmp <- ddply(M_um_out, "k7_from", p_any_snp_given_k7)
colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
M_um_out <- merge(M_um_out, tmp, all = T)
rm(tmp)
M_um_out$p_snp_given_k7 <- M_um_out$k7_mu_N/M_um_out$k7_from_N
M_um_out$k7_mu_rr <- M_um_out$k7_mu_N/M_um_out$k7_from_N

### EXPORT
fwrite(M_um_out, "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_unmasked.table")

### PLOT
plot(M_um_out$k7_mu_rr, out$k7_mu_rr)
cor.test(M_um_out$k7_mu_rr, out$k7_mu_rr)

### ALL HUMAN ###

### IMPORT chromosome counts into list
H_dt_list <- list()
for (chr in 1:19){
  H_dt_list[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_chr", chr, ".table"))
}

# sum rows across list
k7_from <- H_dt_list[[1]]$k7_from
k7_from_N <- rowSums(sapply(H_dt_list, `[[`, 2), na.rm = TRUE)
k7_to <- H_dt_list[[1]]$k7_to
k7_mu_N <- rowSums(sapply(H_dt_list, `[[`, 4), na.rm = TRUE)

# output datatable
H_out <- data.table(k7_from = k7_from,
                    k7_from_N = k7_from_N,
                    k7_to = k7_to,
                    k7_mu_N = k7_mu_N)

# calculate p_any_snp_given_k7 and k7_mu_rr
p_any_snp_given_k7 <- function(sub){
  sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
}
tmp <- ddply(H_out, "k7_from", p_any_snp_given_k7)
colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
H_out <- merge(H_out, tmp, all = T)
rm(tmp)
H_out$p_snp_given_k7 <- H_out$k7_mu_N/H_out$k7_from_N
H_out$k7_mu_rr <- H_out$k7_mu_N/H_out$k7_from_N

### EXPORT
fwrite(H_out, "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table")

### UNMASKED MOUSE ###

### IMPORT chromosome counts into list
H_um_dt_list <- list()
for (chr in 1:19){
  H_um_dt_list[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_unmasked_chr", chr, ".table"))
}

# sum rows across list
k7_from <- H_um_dt_list[[1]]$k7_from
k7_from_N <- rowSums(sapply(H_um_dt_list, `[[`, 2), na.rm = TRUE)
k7_to <- H_um_dt_list[[1]]$k7_to
k7_mu_N <- rowSums(sapply(H_um_dt_list, `[[`, 4), na.rm = TRUE)

# output datatable
H_um_out <- data.table(k7_from = k7_from,
                       k7_from_N = k7_from_N,
                       k7_to = k7_to,
                       k7_mu_N = k7_mu_N)

# calculate p_any_snp_given_k7 and k7_mu_rr
p_any_snp_given_k7 <- function(sub){
  sub$p_any_snp_given_k7 <- sum(sub$k7_mu_N)/sub$k7_from_N[1]
}
tmp <- ddply(H_um_out, "k7_from", p_any_snp_given_k7)
colnames(tmp) <- c("k7_from", "p_any_snp_given_k7")
H_um_out <- merge(H_um_out, tmp, all = T)
rm(tmp)
H_um_out$p_snp_given_k7 <- H_um_out$k7_mu_N/H_um_out$k7_from_N
H_um_out$k7_mu_rr <- H_um_out$k7_mu_N/H_um_out$k7_from_N

### EXPORT
fwrite(H_um_out, "~/Dropbox/PhD/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_unmasked.table")

### PLOT
plot(H_um_out$k7_mu_rr, out$k7_mu_rr)
cor.test(H_um_out$k7_mu_rr, out$k7_mu_rr)

