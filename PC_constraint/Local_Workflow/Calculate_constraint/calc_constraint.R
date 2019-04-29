# SCRIPT that estimates constraint z-scores and RVIS for genes. 

rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)
library(MASS)
source("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Code/Calculate_constraint/function_funZ.R")
source("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Code/Calculate_constraint/function_misZ.R")
source("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Code/Calculate_constraint/function_RVIS.R")


### IMPORT DATA 

# Human
H_SNV_counts <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_SNV_counts.txt")
H_Pmu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")
H_intron <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_intron_length.csv")
H_regional <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_regional_mu_rates.csv")

# Mouse
M_SNV_counts <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/MGP_canPC_SNV_counts.txt")
M_Pmu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")
M_intron <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/M_intron_length.csv")
M_regional <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_regional_mu_rates.csv")


### FORMAT

# Human
H_SNV_counts <- H_SNV_counts[,c("ensembl_transcript_id", "n_synonymous_0.0001","n_missense_0.0001",   
                                "n_nonsense_0.0001", "n_spliceAD_0.0001", "n_intron_0.0001",      
                                "n_synonymous_0.0005", "n_missense_0.0005", "n_nonsense_0.0005",    
                                "n_spliceAD_0.0005", "n_intron_0.0005", "n_synonymous_0.001",   
                                "n_missense_0.001", "n_nonsense_0.001", "n_spliceAD_0.001",     
                                "n_intron_0.001")]
H_regional <- H_regional[,c("ensembl_gene_id", "r_mu")]
H_df <- H_SNV_counts[H_Pmu, on = "ensembl_transcript_id"]
H_df <- H_intron[H_df, on = "ensembl_transcript_id"]
H_df <- H_regional[H_df, on = "ensembl_gene_id"]

# Mouse 
M_SNV_counts <- M_SNV_counts[,c("ensembl_transcript_id", "n_synonymous","n_missense",   
                                "n_nonsense", "n_spliceAD", "n_intron")]
colnames(M_SNV_counts) <- c("ensembl_transcript_id", 
                            paste0(c("n_synonymous", "n_missense", "n_nonsense", "n_spliceAD", "n_intron"), "_0.0001"))
M_regional <- M_regional[,c("ensembl_gene_id", "r_mu")]
M_df <- M_SNV_counts[M_Pmu, on = "ensembl_transcript_id"]
M_df <- M_intron[M_df, on = "ensembl_transcript_id"]
M_df <- M_regional[M_df, on = "ensembl_gene_id"]

### Check ratio of synonymous to total 
H_syn_0005_ratio_median <- median(H_df$n_synonymous_0.0005, na.rm = T)/(median(H_df$n_synonymous_0.0005, na.rm = T) +
                                                  median(H_df$n_missense_0.0005, na.rm = T) +
                                                  median(H_df$n_nonsense_0.0005, na.rm = T))
M_syn_ratio_median <- median(M_df$n_synonymous_0.0001, na.rm = T)/(median(M_df$n_synonymous_0.0001, na.rm = T) +
                                                             median(M_df$n_missense_0.0001, na.rm = T) +
                                                             median(M_df$n_nonsense_0.0001, na.rm = T))
H_syn_0005_ratio_total <- sum(H_df$n_synonymous_0.0005, na.rm = T)/(sum(H_df$n_synonymous_0.0005, na.rm = T) +
                                                                        sum(H_df$n_missense_0.0005, na.rm = T) +
                                                                        sum(H_df$n_nonsense_0.0005, na.rm = T))
M_syn_ratio_total <- sum(M_df$n_synonymous_0.0001, na.rm = T)/(sum(M_df$n_synonymous_0.0001, na.rm = T) +
                                                                   sum(M_df$n_missense_0.0001, na.rm = T) +
                                                                   sum(M_df$n_nonsense_0.0001, na.rm = T))

# adjust Mouse n_syn (n/2) 
M_df_RVIS <- M_df
M_df_Z <- M_df
M_df_Z$n_synonymous_0.0001 <- (M_df_Z$n_synonymous_0.0001/2)


### Calculate funZ 

# Human
H_funZ_0.001 <- fun_Z(H_df, "0.001",  "0.001")
summary(H_funZ_0.001[[2]])
sqrt(H_funZ_0.001[[3]])
H_df_funZ_0.001 <- H_funZ_0.001[[1]]

H_funZ_0.0005 <- fun_Z(H_df, "0.0005", "0.0005")
summary(H_funZ_0.0005[[2]])
sqrt(H_funZ_0.0005[[3]])
H_df_funZ_0.0005 <- H_funZ_0.0005[[1]]

H_funZ_0.0001 <- fun_Z(H_df, "0.0001", "0.0001")
summary(H_funZ_0.0001[[2]])
sqrt(H_funZ_0.0001[[3]])
H_df_funZ_0.0001 <- H_funZ_0.0001[[1]]

# Mouse
M_funZ <- fun_Z(M_df_Z, "0.0001", "0.0001")
summary(M_funZ[[2]])
sqrt(M_funZ[[3]])
M_df_funZ <- M_funZ[[1]]


### Calculate misZ 

# Human
# H_misZ_0.001 <- mis_Z(H_df, "0.001")
# summary(H_misZ_0.001[[2]])
# sqrt(H_misZ_0.001[[3]])
# H_df_misZ_0.001 <- H_misZ_0.001[[1]]

H_misZ_0.0005 <- mis_Z(H_df, "0.0005")
summary(H_misZ_0.0005[[2]])
sqrt(H_misZ_0.0005[[3]])
H_df_misZ_0.0005 <- H_misZ_0.0005[[1]]

# H_misZ_0.0001 <- mis_Z(H_df, "0.0001")
# summary(H_misZ_0.0001[[2]])
# sqrt(H_misZ_0.0001[[3]])
# H_df_misZ_0.0001 <- H_misZ_0.0001[[1]]

# Mouse
M_misZ <- mis_Z(M_df_Z, "0.0001")
summary(M_misZ[[2]])
sqrt(M_misZ[[3]])
M_df_misZ <- M_misZ[[1]]


### Calculate RVIS

# Human
H_RVIS_0.001 <- RVIS(H_df, fun_maf = "0.001", all_maf = "0.0001")
summary(H_RVIS_0.001[[2]])
H_df_RVIS_0.001 <- H_RVIS_0.001[[1]]

# H_RVIS_0.0005 <- RVIS(H_df, fun_maf = "0.0005", all_maf = "0.0001")
# summary(H_RVIS_0.0005[[2]])
# H_df_RVIS_0.0005 <- H_RVIS_0.0005[[1]]

# Mouse
M_RVIS <- RVIS(M_df_RVIS, "0.0001", "0.0001")
summary(M_RVIS[[2]])
M_df_RVIS <- M_RVIS[[1]]


### Save datasets

# # Human RVIS
fwrite(H_df_RVIS_0.001, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_RVIS_0.001.csv")
# Human funZ
fwrite(H_df_funZ_0.001, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
fwrite(H_df_funZ_0.0005, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.0005.csv")
fwrite(H_df_funZ_0.0001, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.0001.csv")
# Human misZ
fwrite(H_df_misZ_0.0005, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_misZ_0.0005.csv")


# Mouse RVIS
fwrite(M_df_RVIS, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_RVIS.csv")
# Mouse funZ
fwrite(M_df_funZ, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")
# Mouse misZ
fwrite(M_df_misZ, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_misZ.csv")


### COMBINE BY OTHOLOGS

# subset columns

H_out <- H_df_funZ_0.001[,c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                                    "fun_Z_0.001", "fun_Z_percentile_0.001")]
M_out <- M_df_funZ[,c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                                       "fun_Z_0.0001", "fun_Z_percentile_0.0001")]

orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/HM_Ensembl_orthologs.csv")
colnames(M_out) <- paste0("M_", colnames(M_out))
colnames(H_out) <- paste0("H_", colnames(H_out))
orth_out <- orths[M_out, on = c("M_ensembl_transcript_id", "M_ensembl_gene_id")]
orth_out <- orth_out[H_out, on = c("H_ensembl_transcript_id", "H_ensembl_gene_id")]
orth_out <- orth_out[,c("H_external_gene_name", "M_external_gene_name",
                        "H_ensembl_transcript_id", "M_ensembl_transcript_id",
                        "M_fun_Z_0.0001", "M_fun_Z_percentile_0.0001", 
                        # "M_RVIS_0.0001", "M_RVIS_percentile_0.0001", 
                        # "M_mis_Z_0.0001", "M_mis_Z_percentile_0.0001", 
                        "H_fun_Z_0.001", "H_fun_Z_percentile_0.001", 
                        # "H_RVIS_0.001", "H_RVIS_percentile_0.001",
                        # "H_mis_Z_0.0005", "H_mis_Z_percentile_0.0005",
                        "orthology_type")]
# orth_out <- orth_out[!is.na(orth_out$H_external_gene_name),]
summary(orth_out)
orth_out <- orth_out[complete.cases(orth_out),]
fwrite(orth_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")

summary(orth_out)
length(unique(orth_out$H_external_gene_name))
length(unique(orth_out$H_ensembl_transcript_id))
length(unique(orth_out$M_external_gene_name))
length(unique(orth_out$M_ensembl_transcript_id))
table(orth_out$orthology_type)
o2o <- subset(orth_out, orth_out$orthology_type == "ortholog_one2one")
length(unique(o2o$H_external_gene_name))
length(unique(o2o$H_ensembl_transcript_id))
length(unique(o2o$M_external_gene_name))
length(unique(o2o$M_ensembl_transcript_id))
o2m <- subset(orth_out, orth_out$orthology_type == "ortholog_one2many")
length(unique(o2m$H_external_gene_name))
length(unique(o2m$H_ensembl_transcript_id))
length(unique(o2m$M_external_gene_name))
length(unique(o2m$M_ensembl_transcript_id))
m2m <- subset(orth_out, orth_out$orthology_type == "ortholog_many2many")
length(unique(m2m$H_external_gene_name))
length(unique(m2m$H_ensembl_transcript_id))
length(unique(m2m$M_external_gene_name))
length(unique(m2m$M_ensembl_transcript_id))


##########

### COMBINE BY OTHOLOGS

# # subset columns
# H_df_RVIS_0.001 <- H_df_RVIS_0.001[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                                       "RVIS_0.001", "RVIS_percentile_0.001")]
# H_df_funZ_0.001 <- H_df_funZ_0.001[,c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                                       "fun_Z_0.001", "fun_Z_percentile_0.001")]
# H_df_misZ_0.0005 <- H_df_misZ_0.0005[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                                         "mis_Z_0.0005", "mis_Z_percentile_0.0005")]
# M_df_RVIS <- M_df_RVIS[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                           "RVIS_0.0001", "RVIS_percentile_0.0001")]
# M_df_funZ <- M_df_funZ[,c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                           "fun_Z_0.0001", "fun_Z_percentile_0.0001")]
# M_df_misZ <- M_df_misZ[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
#                           "mis_Z_0.0001", "mis_Z_percentile_0.0001")]
# H_out <- merge(H_df_RVIS_0.001, H_df_funZ_0.001, all = T)
# H_out <- merge(H_out, H_df_misZ_0.0005, all = T)
# 
# M_out <- merge(M_df_RVIS, M_df_funZ, all = T)
# M_out <- merge(M_out, M_df_misZ, all = TRUE)
# 
# orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/HM_Ensembl_orthologs.csv")
# colnames(M_out) <- paste0("M_", colnames(M_out))
# colnames(H_out) <- paste0("H_", colnames(H_out))
# orth_out <- orths[M_out, on = c("M_ensembl_transcript_id", "M_ensembl_gene_id")]
# orth_out <- orth_out[H_out, on = c("H_ensembl_transcript_id", "H_ensembl_gene_id")]
# orth_out <- orth_out[,c("H_external_gene_name", "M_external_gene_name",
#                         "H_ensembl_transcript_id", "M_ensembl_transcript_id",
#                         "M_fun_Z_0.0001", "M_fun_Z_percentile_0.0001", 
#                         # "M_RVIS_0.0001", "M_RVIS_percentile_0.0001", 
#                         # "M_mis_Z_0.0001", "M_mis_Z_percentile_0.0001", 
#                         "H_fun_Z_0.001", "H_fun_Z_percentile_0.001", 
#                         # "H_RVIS_0.001", "H_RVIS_percentile_0.001",
#                         # "H_mis_Z_0.0005", "H_mis_Z_percentile_0.0005",
#                         "orthology_type")]
# # orth_out <- orth_out[!is.na(orth_out$H_external_gene_name),]
# orth_out <- orth_out[complete.cases(orth_out),]
# fwrite(orth_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")


