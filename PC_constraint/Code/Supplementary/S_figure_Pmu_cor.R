rm(list = ls())
graphics.off()

library(gridExtra)
library(data.table)
library(plyr)

### IMPORT
H <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")
M <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")


### FORMAT

H_df <- H[,c("ensembl_transcript_id","cds_length","p_syn","p_mis","p_non")]
M_df <- M[,c("ensembl_transcript_id","cds_length","p_syn","p_mis","p_non")]

# orthologous genes
orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
orths <- orths[,c("H_ensembl_transcript_id", "M_ensembl_transcript_id")]

colnames(H_df) <- paste0("H_", colnames(H_df))
colnames(M_df) <- paste0("M_", colnames(M_df))
df <- H_df[orths, on = "H_ensembl_transcript_id"]
df <- M_df[df, on = "M_ensembl_transcript_id"]

plot(df$M_p_syn, df$H_p_syn)
cor.test(df$M_p_syn, df$H_p_syn)
plot(df$M_p_mis, df$H_p_mis)
cor.test(df$M_p_mis, df$H_p_mis)
plot(df$M_p_non, df$H_p_non)
cor.test(df$M_p_non, df$H_p_non)
plot(df$M_cds_length, df$H_cds_length)
cor.test(df$M_cds_length, df$H_cds_length)
