rm(list = ls())
graphics.off()

counts <- fread("~/Dropbox/PhD/Data/SS_constraint/Variant_counts/SS_gnomAD_variant_counts.csv")

p_mu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")

colnames(counts) <- c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",   
                      "external_gene_name", "n_synonymous_all", "n_synonymous_male",  
                      "n_synonymous_female", "n_lof_male", "n_lof_female")
counts$chromosome_name <- as.character(counts$chromosome_name)
dt <- p_mu[counts, on = c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",   
                         "external_gene_name")]


### QC

dt$n_synonymous_mean <- (dt$n_synonymous_female + dt$n_synonymous_male)/2

# remove na's
rm.ind <- which(is.na(dt$n_synonymous_mean) | is.na(dt$p_syn) | is.na(dt$p_non) | 
                  is.na(dt$n_lof_male) | is.na(dt$n_lof_female))
if (length(rm.ind) != 0) { df_temp <- dt[-rm.ind,] }

# remove genes n_total > 300
df_temp <- subset(df_temp, df_temp$n_synonymous_mean <= 400)

# remove genes with length > 15000
df_temp <- subset(df_temp, df_temp$cds_length <= 15000)


mod <- lm(n_synonymous_mean ~ p_syn, data = df_temp)
# plot(mod)
# summary(mod)

# r2 <- summary(mod)$adj.r.squared
mse <- mean(mod$residuals^2)

### PREDICT
df_temp$exp_syn <- predict(mod,  data.frame(p_syn = df_temp$p_syn))
df_temp$exp_lof <- predict(mod,  data.frame(p_syn = df_temp$p_non))


### Z scores
dif <- df_temp$exp_syn - df_temp$n_synonymous_male
mean_dif <- mean(dif, na.rm = T)
sd_dif <- sd(dif, na.rm = T)
df_temp$syn_male_Z <- (dif - mean_dif)/sd_dif

dif <- df_temp$exp_syn - df_temp$n_synonymous_female
mean_dif <- mean(dif, na.rm = T)
sd_dif <- sd(dif, na.rm = T)
df_temp$syn_female_Z <- (dif - mean_dif)/sd_dif

dif <- df_temp$exp_lof - df_temp$n_lof_male
mean_dif <- mean(dif, na.rm = T)
sd_dif <- sd(dif, na.rm = T)
df_temp$lof_male_Z <- (dif - mean_dif)/sd_dif

dif <- df_temp$exp_lof - df_temp$n_lof_female
mean_dif <- mean(dif, na.rm = T)
sd_dif <- sd(dif, na.rm = T)
df_temp$lof_female_Z <- (dif - mean_dif)/sd_dif

plot(df_temp$syn_female_Z, df_temp$syn_male_Z)
plot(df_temp$lof_female_Z, df_temp$lof_male_Z)

percentile <- ecdf(df_temp$lof_male_Z)
df_temp$lof_male_Z_p <- percentile(df_temp$lof_male_Z)
percentile <- ecdf(df_temp$lof_female_Z)
df_temp$lof_female_Z_p <- percentile(df_temp$lof_female_Z)

df_temp$lof_dif <- df_temp$lof_male_Z_p - df_temp$lof_female_Z_p
hist(df_temp$lof_dif)
summary(df_temp$lof_dif)
x <- subset(df_temp, df_temp$lof_dif > 0.65 | df_temp$lof_dif < -0.65)


impc <- c("Mpo", "Zmynd8", "Snrnp200", "Gpr139", "Ankrd1", "Spag4", "Dact2", "Gpa33",
"Ofcc1", "Clec9a", "Sacs", "Prdm15")
impc <- toupper(impc)
y <- df_temp[which(df_temp$external_gene_name %in% impc),]



# try X chromosome
