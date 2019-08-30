# SCRIPT that calculates constraint scores

rm(list=ls())
graphics.off()

library(data.table)
library(MASS)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

############
### SET VARS 
############

N_SNV_file <- args[1]
P_SNV_file <- args[2]
N_CpG_file <- args[3]
rd_file <- args[4]
out_file <- args[5]

N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF001.csv"
P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv"
N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv"
rd_file <- "~/Dropbox/PhD/Data/gnomAD/coverage/Ensembl_v94_human_canPC_gnomad_mask.csv"
out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF001_v2.csv"

M_N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_canPC_N_SNVs.csv"
M_P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv"
M_N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv"
M_rd_file <- "~/Dropbox/PhD/Data/MGP/coverage/Ensembl_v94_mouse_canPC_MGP_mask.csv"
M_out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv"


##########
### HUMAN
#########

### IMPORT DATA 
###############

# Human
H_N_SNV <- fread(H_N_SNV_file)
H_P_SNV <- fread(H_P_SNV_file)
H_N_CpG <- fread(H_N_CpG_file)
H_rd <- fread(H_rd_file)

### FORMAT
##########

# merge
H_df <- H_N_SNV[H_P_SNV, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
H_df <- H_df[H_N_CpG, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
H_df <- H_rd[H_df, on = c("chromosome_name", "external_gene_name", "ensembl_gene_id",      
                          "ensembl_transcript_id", "cds_length")]

# remove na's
rm.ind <- which(is.na(H_df$n_synonymous) | is.na(H_df$n_missense) | is.na(H_df$n_nonsense) | 
                  is.na(H_df$p_syn) | is.na(H_df$p_mis) | is.na(H_df$p_non))
if (length(rm.ind) != 0) { H_df <- H_df[-rm.ind,] }

# remove cds > 15000
rm.id <- which(H_df$cds_length>15000)
if (length(rm.id) != 0){ H_df <- H_df[-rm.id,] }

# remove cds < 300
rm.id <- which(H_df$cds_length<300)
if (length(rm.id) != 0){ H_df <- H_df[-rm.id,] }

# remove cds > 10% masked
rm.id <- which(H_df$cds_fraction_masked > 0.1)
if (length(rm.id) != 0){  H_df <- H_df[-rm.id,] }

# check ratio of synonymous to total 
H_syn_ratio_total <- sum(H_df$n_synonymous, na.rm = T)/(sum(H_df$n_synonymous, na.rm = T) +
                                                          sum(H_df$n_missense, na.rm = T) +
                                                          sum(H_df$n_nonsense, na.rm = T))


# calculate "n_functional"
H_df$n_functional <- H_df$n_missense + H_df$n_nonsense

# calculate "p_functional"
H_df$p_functional <- H_df$p_mis + H_df$p_non


### CALCULATE Functional Z-score
################################

# fit model for synonymous variants
H_mod <- lm(n_synonymous ~ p_syn + n_CpG,
            data = H_df)
# H_mod <- lm(n_synonymous ~ cds_length,
#           data = H_df)
# H_mod <- lm(n_synonymous ~ p_syn,
#           data = H_df)
# plot(H_mod)
# summary(H_mod)

# model coeff
H_coef <- summary(H_mod)$coef
H_r2 <- summary(H_mod)$adj.r.squared
H_mse <- mean(H_mod$residuals^2) # mean squared error
# 95% confidence interval
H_error_95 <- qnorm(0.975)*sd(H_mod$residuals)/sqrt(length(H_mod$residuals))


### PREDICT
H_df$exp_functional <- predict(H_mod, data.frame(p_syn = H_df$p_functional,
                                                 n_CpG = H_df$n_CpG))
H_df$exp_functional_95CI_upper <- H_df$exp_functional + H_error_95
H_df$exp_functional_95CI_lower <- H_df$exp_functional - H_error_95

### funZ SCORES
H_dif <- H_df$exp_functional - H_df$n_functional
H_dif_upper <- H_df$exp_functional_95CI_upper - H_df$n_functional
H_dif_lower <- H_df$exp_functional_95CI_lower - H_df$n_functional

H_df$fun_Z <- (H_dif - mean(H_dif, na.rm = T))/sd(H_dif, na.rm = T)
H_df$fun_Z_95CI_upper <- (H_dif_upper - mean(H_dif, na.rm = T))/sd(H_dif, na.rm = T)
H_df$fun_Z_95CI_lower <- (H_dif_lower - mean(H_dif, na.rm = T))/sd(H_dif, na.rm = T)


### CALCULATE OE ratio
######################

H_df$OE_ratio <- H_df$n_functional/H_df$exp_functional


### OUTPUT
##########

fwrite(H_df, H_out_file)


#########
### MOUSE
#########

### IMPORT DATA
###############

M_N_SNV <- fread(M_N_SNV_file)
M_P_SNV <- fread(M_P_SNV_file)
M_N_CpG <- fread(M_N_CpG_file)
M_rd <- fread(M_rd_file)

### FORMAT
##########

# merge
M_df <- M_N_SNV[M_P_SNV, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
M_df <- M_df[M_N_CpG, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
M_df <- M_rd[M_df, on = c("chromosome_name", "external_gene_name", "ensembl_gene_id",
                          "ensembl_transcript_id", "cds_length")]
# remove na's
rm.ind <- which(is.na(M_df$n_synonymous) | is.na(M_df$n_missense) | is.na(M_df$n_nonsense) |
                  is.na(M_df$p_syn) | is.na(M_df$p_mis) | is.na(M_df$p_non))
if (length(rm.ind) != 0) { M_df <- M_df[-rm.ind,] }

# remove cds > 15000
rm.id <- which(M_df$cds_length>15000)
if (length(rm.id) != 0){ M_df <- M_df[-rm.id,] }

# remove cds < 300
rm.id <- which(M_df$cds_length<300)
if (length(rm.id) != 0){ M_df <- M_df[-rm.id,] }

# remove cds > 10% masked
rm.id <- which(M_df$cds_fraction_masked > 0.1)
if (length(rm.id) != 0){  M_df <- M_df[-rm.id,] }

# check ratio of synonymous to total
M_syn_ratio_total <- sum(M_df$n_synonymous, na.rm = T)/(sum(M_df$n_synonymous, na.rm = T) +
                                                          sum(M_df$n_missense, na.rm = T) +
                                                          sum(M_df$n_nonsense, na.rm = T))
# divide n_synonymou by 2
M_df$n_synonymous <- as.numeric(M_df$n_synonymous)/2

# calculate "n_functional"
M_df$n_functional <- M_df$n_missense + M_df$n_nonsense

# calculate "p_functional"
M_df$p_functional <- M_df$p_mis + M_df$p_non


### CALCULATE Functional Z-score
################################

# fit model for synonymous variants
M_mod <- lm(n_synonymous ~ p_syn + n_CpG,
            data = M_df)
# M_mod <- lm(n_synonymous ~ p_syn,
#           data = M_df)
# M_mod <- lm(n_synonymous ~ cds_length,
#           data = M_df)
# plot(M_mod)
# summary(M_mod)

# model coeff
M_coef <- summary(M_mod)$coef
M_r2 <- summary(M_mod)$adj.r.squared
M_mse <- mean(M_mod$residuals^2) # mean squared error
# 95% confidence interval
M_error_95 <- qnorm(0.975)*sd(M_mod$residuals)/sqrt(length(M_mod$residuals))


### PREDICT
M_df$exp_functional <- predict(M_mod, data.frame(p_syn = M_df$p_functional,
                                                 n_CpG = M_df$n_CpG))
M_df$exp_functional_95CI_upper <- M_df$exp_functional + M_error_95
M_df$exp_functional_95CI_lower <- M_df$exp_functional - M_error_95

### funZ SCORES
M_dif <- M_df$exp_functional - M_df$n_functional
M_dif_upper <- M_df$exp_functional_95CI_upper - M_df$n_functional
M_dif_lower <- M_df$exp_functional_95CI_lower - M_df$n_functional

M_df$fun_Z <- (M_dif - mean(M_dif, na.rm = T))/sd(M_dif, na.rm = T)
M_df$fun_Z_95CI_upper <- (M_dif_upper - mean(M_dif, na.rm = T))/sd(M_dif, na.rm = T)
M_df$fun_Z_95CI_lower <- (M_dif_lower - mean(M_dif, na.rm = T))/sd(M_dif, na.rm = T)

length(unique(M_df$ensembl_gene_id))


### CALCULATE OE ratio
######################

M_df$OE_ratio <- M_df$n_functional/M_df$exp_functional


### OUTPUT
fwrite(M_df, M_out_file)


##########

### CALCULATE Functional Z-score adjusted for CDS length
########################################################

# calculate cds_length adjusted variables 
H_df$n_synonymous_adj <- H_df$n_synonymous / H_df$cds_length
H_df$n_functional_adj <- H_df$n_functional / H_df$cds_length
H_df$p_synonymous_adj <- H_df$p_syn / H_df$cds_length
H_df$p_functional_adj <- H_df$p_functional / H_df$cds_length
H_df$n_CpG_adj <- H_df$n_CpG / H_df$cds_length

# fit model for synonymous variants
H_mod_adj <- lm(n_synonymous_adj ~ p_synonymous_adj + n_CpG_adj,
                data = H_df)
# plot(H_mod_adj)
# summary(H_mod_adj)

# model coeff
H_coef_adj <- summary(H_mod_adj)$coef
H_r2_adj <- summary(H_mod_adj)$adj.r.squared
H_mse_adj <- mean(H_mod_adj$residuals^2) # mean squared error
# 95% confidence interval
H_error_95_adj <- qnorm(0.975)*sd(H_mod_adj$residuals)/sqrt(length(H_mod_adj$residuals))


### PREDICT
H_df$exp_functional_adj <- predict(H_mod_adj, data.frame(p_synonymous_adj = H_df$p_functional_adj,
                                                         n_CpG_adj = H_df$n_CpG_adj))

### funZ SCORES
H_dif_adj <- H_df$exp_functional_adj - H_df$n_functional_adj
H_df$fun_Z_adj <- (H_dif_adj - mean(H_dif_adj, na.rm = T))/sd(H_dif_adj, na.rm = T)


### CALCULATE Functional Z-score adjusted for CDS length
########################################################

# calculate cds_length adjusted variables 
M_df$n_synonymous_adj <- M_df$n_synonymous / M_df$cds_length
M_df$n_functional_adj <- M_df$n_functional / M_df$cds_length
M_df$p_synonymous_adj <- M_df$p_syn / M_df$cds_length
M_df$p_functional_adj <- M_df$p_functional / M_df$cds_length
M_df$n_CpG_adj <- M_df$n_CpG / M_df$cds_length

# fit model for synonymous variants
M_mod_adj <- lm(n_synonymous_adj ~ p_synonymous_adj + n_CpG_adj,
                data = M_df)
# plot(M_mod_adj)
# summary(M_mod_adj)

# model coeff
M_coef_adj <- summary(M_mod_adj)$coef
M_r2_adj <- summary(M_mod_adj)$adj.r.squared
M_mse_adj <- mean(M_mod_adj$residuals^2) # mean squared error
# 95% confidence interval
M_error_95_adj <- qnorm(0.975)*sd(M_mod_adj$residuals)/sqrt(length(M_mod_adj$residuals))


### PREDICT
M_df$exp_functional_adj <- predict(M_mod_adj, data.frame(p_synonymous_adj = M_df$p_functional_adj,
                                                         n_CpG_adj = M_df$n_CpG_adj))

### funZ SCORES
M_dif_adj <- M_df$exp_functional_adj - M_df$n_functional_adj
M_df$fun_Z_adj <- (M_dif_adj - mean(M_dif_adj, na.rm = T))/sd(M_dif_adj, na.rm = T)




# old <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_v1.csv")
# colnames(old)
# old <- old[,c("ensembl_transcript_id", "fun_Z_0.0001")]
# tmp <- old[M_df, on="ensembl_transcript_id"]
# plot(tmp$fun_Z, tmp$fun_Z_0.0001)
# cor.test(tmp$fun_Z, tmp$fun_Z_0.0001)
# 
# plot_tmp <- subset(tmp, tmp$fun_Z < 10 & tmp$fun_Z > -10)
# plot_tmp <- subset(plot_tmp, plot_tmp$fun_Z_0.0001 < 10 & plot_tmp$fun_Z_0.0001 > -10)
# plot(plot_tmp$fun_Z, plot_tmp$fun_Z_0.0001)
# cor.test(plot_tmp$fun_Z, plot_tmp$fun_Z_0.0001)
# 
# 
# ### FUNCTION for plotting Pearson's correlation on ggplot
# corr_eqn <- function(x,y, digits = 2) {
#   corr_coef <- round(cor(x, y), digits = digits)
#   paste("italic(r) == ", corr_coef)
# }
# pearsons <- corr_eqn(plot_tmp$fun_Z, plot_tmp$fun_Z_0.0001, digits = 3)
# library(ggplot2)
# Fig2 <- ggplot(plot_tmp, aes(x = fun_Z, y = fun_Z_0.0001)) +
#   geom_point() +
#   geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
#   # annotate("text", x = -4, y = 8, label = lm_equation, colour="black", size = 5, parse=TRUE) +
#   annotate("text", x = -7, y = 8, label = pearsons, colour="black", size = 5, parse=TRUE) +
#   # scale_color_manual(breaks = c("Most constrained", "Least constrained"), values=c('gray80', 'red', 'blue')) +
#   xlab("funZ (Mus musculus strains)") +
#   ylab('funZ (all strains)') +
#   xlim(-10, 10) +
#   ylim(-10, 10) +
#   theme_bw() +
#   theme(legend.position="top",
#         legend.title=element_blank(),
#         legend.text=element_text(size=14),
#         text = element_text(size = 14),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.background=element_blank(),
#         plot.margin=unit(c(1,1,1,1),"cm"))
# Fig2
# ggsave("~/Dropbox/PhD/Data/Figures/FunZ_spretus_musculus_cor.jpg", Fig2, height = 6, width = 6)

# library(gmodels)
# test <- ci(H_mod$residuals)


### test kmer r2
# h_kmer <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")
# # h_kmer <- subset(h_kmer, h_kmer$chromosome_name != "X")
# # h_kmer$chromosome_name <- as.integer(h_kmer$chromosome_name)
# h_kmer <- h_kmer[,c("external_gene_name", "p_syn")]
# colnames(h_kmer) <- c("external_gene_name", "p_syn_kmer")
# test <- merge(H_df, h_kmer)
# test_mod1 <- lm(n_synonymous ~ p_syn + n_CpG,
#             data = test)
# summary(test_mod1)
# test_mod2 <- lm(n_synonymous ~ p_syn_kmer,
#                 data = test)
# summary(test_mod2)
# 
# m_kmer <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")
# # m_kmer <- subset(m_kmer, m_kmer$chromosome_name != "X")
# # m_kmer$chromosome_name <- as.integer(m_kmer$chromosome_name)
# m_kmer <- m_kmer[,c("external_gene_name", "p_syn")]
# colnames(m_kmer) <- c("external_gene_name", "p_syn_kmer")
# test <- merge(M_df, m_kmer)
# test_mod1 <- lm(n_synonymous ~ p_syn + n_CpG,
#                 data = test)
# summary(test_mod1)
# test_mod2 <- lm(n_synonymous ~ p_syn_kmer,
#                 data = test)
# summary(test_mod2)
