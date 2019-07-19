# SCRIPT that calculates constraint scores

rm(list=ls())
graphics.off()

library(data.table)
library(MASS)

### SET VARS 

H_N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF001.csv"
H_P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv"
H_N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv"
H_out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_MAF001_v2.csv"

M_N_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_canPC_N_SNVs.csv"
M_P_SNV_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv"
M_N_CpG_file <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv"
M_out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv"



### HUMAN

### IMPORT DATA 

# Human
H_N_SNV <- fread(H_N_SNV_file)
H_P_SNV <- fread(H_P_SNV_file)
H_N_CpG <- fread(H_N_CpG_file)

### FORMAT

# merge
H_df <- H_N_SNV[H_P_SNV, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
H_df <- H_df[H_N_CpG, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]

# remove na's
rm.ind <- which(is.na(H_df$n_synonymous) | is.na(H_df$n_missense) | is.na(H_df$n_nonsense) | 
  is.na(H_df$p_syn) | is.na(H_df$p_mis) | is.na(H_df$p_non))
if (length(rm.ind) != 0) { H_df <- H_df[-rm.ind,] }

# check ratio of synonymous to total 
H_syn_ratio_total <- sum(H_df$n_synonymous, na.rm = T)/(sum(H_df$n_synonymous, na.rm = T) +
                                                                      sum(H_df$n_missense, na.rm = T) +
                                                                      sum(H_df$n_nonsense, na.rm = T))


# calculate "n_functional"
H_df$n_functional <- H_df$n_missense + H_df$n_nonsense

# calculate "p_functional"
H_df$p_functional <- H_df$p_mis + H_df$p_non


### CALCULATE Functional Z-score

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


### OUTPUT
fwrite(H_df, H_out_file)


### MOUSE

### IMPORT DATA
M_N_SNV <- fread(M_N_SNV_file)
M_P_SNV <- fread(M_P_SNV_file)
M_N_CpG <- fread(M_N_CpG_file)

### FORMAT

# merge
M_df <- M_N_SNV[M_P_SNV, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]
M_df <- M_df[M_N_CpG, on = c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id")]

# remove na's
rm.ind <- which(is.na(M_df$n_synonymous) | is.na(M_df$n_missense) | is.na(M_df$n_nonsense) |
                  is.na(M_df$p_syn) | is.na(M_df$p_mis) | is.na(M_df$p_non))
if (length(rm.ind) != 0) { M_df <- M_df[-rm.ind,] }

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


### OUTPUT
fwrite(M_df, M_out_file)

##########
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
