rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)


### IMPORT DATA 

# Human
H_SNV_counts <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_SNV_counts.txt")

# Mouse
M_SNV_counts <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/MGP_canPC_SNV_counts.txt")

### FORMAT

# Human
H_df <- H_SNV_counts[,c("ensembl_transcript_id", "n_synonymous_0.001",   
                                "n_missense_0.001", "n_nonsense_0.001", "n_spliceAD_0.001",     
                                "n_intron_0.001")]
# H_df <- H_SNV_counts[,c("ensembl_transcript_id", "n_synonymous_0.0001","n_missense_0.0001",   
#                         "n_nonsense_0.0001", "n_spliceAD_0.0001", "n_intron_0.0001",      
#                         "n_synonymous_0.0005", "n_missense_0.0005", "n_nonsense_0.0005",    
#                         "n_spliceAD_0.0005", "n_intron_0.0005", "n_synonymous_0.001",   
#                         "n_missense_0.001", "n_nonsense_0.001", "n_spliceAD_0.001",     
#                         "n_intron_0.001")]

# Mouse 
M_df <- M_SNV_counts[,c("ensembl_transcript_id", "n_synonymous","n_missense",   
                                "n_nonsense", "n_spliceAD", "n_intron")]

# adjust Mouse n_syn (n/2) 
M_df$n_synonymous <- (M_df$n_synonymous/2)

# orthologous genes
orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
orths <- orths[,c("H_ensembl_transcript_id", "M_ensembl_transcript_id")]

colnames(H_df) <- paste0("H_", colnames(H_df))
colnames(M_df) <- paste0("M_", colnames(M_df))
df <- H_df[orths, on = "H_ensembl_transcript_id"]
df <- M_df[df, on = "M_ensembl_transcript_id"]

plot(df$M_n_synonymous, df$H_n_synonymous_0.001)
cor.test(df$M_n_synonymous, df$H_n_synonymous_0.001, method = "spearman")
plot(df$M_n_missense, df$H_n_missense_0.001)
cor.test(df$M_n_missense, df$H_n_missense_0.001, method = "spearman")
plot(df$M_n_nonsense, df$H_n_nonsense_0.001)
cor.test(df$M_n_nonsense, df$H_n_nonsense_0.001, method = "spearman")

