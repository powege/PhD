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
SNV_source <- args[6]

# N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF001.csv"
# P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv"
# N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv"
# rd_file <- "~/Dropbox/PhD/Data/1KGP/Masks/Formatted/Ensembl_v94_human_canPC_1KGP_mask.csv"
# out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF001_v2.csv"
# SNV_source <- "1KGP"
# 
# N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_canPC_N_SNVs.csv"
# P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv"
# N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv"
# rd_file <- "~/Dropbox/PhD/Data/MGP/coverage/Ensembl_v94_mouse_canPC_MGP_mask.csv"
# out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv"
# SNV_source <- "MGP"

################
### IMPORT DATA 
###############

N_SNV <- fread(N_SNV_file)
P_SNV <- fread(P_SNV_file)
N_CpG <- fread(N_CpG_file)
rd <- fread(rd_file)

##################
### FORMAT and QC
################

# merge
df <- N_SNV[P_SNV, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
df <- df[N_CpG, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
df <- rd[df, on = c("chromosome_name", "external_gene_name", "ensembl_gene_id",      
                           "ensembl_transcript_id", "cds_length")]

# remove na's
rm.ind <- which(is.na(df$n_synonymous) | is.na(df$n_missense) | is.na(df$n_nonsense) | 
  is.na(df$p_syn) | is.na(df$p_mis) | is.na(df$p_non))
if (length(rm.ind) != 0) { df <- df[-rm.ind,] }

# remove cds > 15000
rm.id <- which(df$cds_length>15000)
if (length(rm.id) != 0){ df <- df[-rm.id,] }

# remove cds < 300
rm.id <- which(df$cds_length<300)
if (length(rm.id) != 0){ df <- df[-rm.id,] }

# remove cds > 10% masked
rm.id <- which(df$cds_fraction_masked > 0.1)
if (length(rm.id) != 0){  df <- df[-rm.id,] }

# check ratio of synonymous to total 
syn_ratio_total <- sum(df$n_synonymous, na.rm = T)/(sum(df$n_synonymous, na.rm = T) +
                                                                      sum(df$n_missense, na.rm = T) +
                                                                      sum(df$n_nonsense, na.rm = T))
# divide n_syn by 2
if (SNV_source == "MGP"){ df$n_synonymous <- as.numeric(df$n_synonymous)/2 }

# calculate "n_functional"
df$n_functional <- df$n_missense + df$n_nonsense

# calculate "p_functional"
df$p_functional <- df$p_mis + df$p_non

################################
### CALCULATE Functional Z-score
################################

# fit model for synonymous variants
mod <- lm(n_synonymous ~ p_syn + n_CpG,
          data = df)
# mod <- lm(n_synonymous ~ cds_length,
#           data = df)
# mod <- lm(n_synonymous ~ p_syn,
#           data = df)
# plot(mod)
# summary(mod)

# model coeff
coef <- summary(mod)$coef
r2 <- summary(mod)$adj.r.squared
mse <- mean(mod$residuals^2) # mean squared error
# 95% confidence interval
error_95 <- qnorm(0.975)*sd(mod$residuals)/sqrt(length(mod$residuals))

### PREDICT
df$exp_functional <- predict(mod, data.frame(p_syn = df$p_functional,
                                               n_CpG = df$n_CpG))
# df$exp_functional_95CI_upper <- df$exp_functional + error_95
# df$exp_functional_95CI_lower <- df$exp_functional - error_95

### funZ SCORES
dif <- df$exp_functional - df$n_functional
# dif_upper <- df$exp_functional_95CI_upper - df$n_functional
# dif_lower <- df$exp_functional_95CI_lower - df$n_functional

df$fun_Z <- (dif - mean(dif, na.rm = T))/sd(dif, na.rm = T)
# df$fun_Z_95CI_upper <- (dif_upper - mean(dif, na.rm = T))/sd(dif, na.rm = T)
# df$fun_Z_95CI_lower <- (dif_lower - mean(dif, na.rm = T))/sd(dif, na.rm = T)

######################
### CALCULATE OE ratio
######################

df$OE_ratio <- df$n_functional/df$exp_functional

###########
### OUTPUT
##########

fwrite(df, out_file)


# ########################################################
# ### CALCULATE Functional Z-score adjusted for CDS length
# ########################################################
# 
# # calculate cds_length adjusted variables 
# df$n_synonymous_adj <- df$n_synonymous / df$cds_length
# df$n_functional_adj <- df$n_functional / df$cds_length
# df$p_synonymous_adj <- df$p_syn / df$cds_length
# df$p_functional_adj <- df$p_functional / df$cds_length
# df$n_CpG_adj <- df$n_CpG / df$cds_length
# 
# # fit model for synonymous variants
# mod_adj <- lm(n_synonymous_adj ~ p_synonymous_adj + n_CpG_adj,
#                 data = df)
# # plot(mod_adj)
# # summary(mod_adj)
# 
# # model coeff
# coef_adj <- summary(mod_adj)$coef
# r2_adj <- summary(mod_adj)$adj.r.squared
# mse_adj <- mean(mod_adj$residuals^2) # mean squared error
# # 95% confidence interval
# error_95_adj <- qnorm(0.975)*sd(mod_adj$residuals)/sqrt(length(mod_adj$residuals))
# 
# 
# ### PREDICT
# df$exp_functional_adj <- predict(mod_adj, data.frame(p_synonymous_adj = df$p_functional_adj,
#                                                          n_CpG_adj = df$n_CpG_adj))
# 
# ### funZ SCORES
# dif_adj <- df$exp_functional_adj - df$n_functional_adj
# df$fun_Z_adj <- (dif_adj - mean(dif_adj, na.rm = T))/sd(dif_adj, na.rm = T)
# 
# ###########
# ### OUTPUT
# ##########
# 
# fwrite(df, out_file)


##################################################
### CALCULATE Functional Z-score MSE upper limmit
#################################################


### funZ SCORES MSE upper limit 
dif_MSE <- (df$exp_functional + mse) - df$n_functional
df$fun_Z_MSE <- (dif_MSE - mean(dif_MSE, na.rm = T))/sd(dif_MSE, na.rm = T)

### OUTPUT
fwrite(df, out_file)
