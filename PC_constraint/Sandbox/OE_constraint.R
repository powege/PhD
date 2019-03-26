rm(list = ls())
graphics.off()

library(data.table)

### Script that calculates the OE ratio as a measure of canstraint for human and mouse genes. 
# OE ratio = sequence specific probability of functional (missense and nonsense) SNV / observed number of functional SNVs

m_pSNV <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")
h_pSNV <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")

m_SNVs <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/MGP_canPC_SNV_counts.txt")
h_SNVs <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_SNV_counts.txt")


colnames(m_SNVs) <- c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                      "external_gene_name", "n_synonymous", "n_missense",
                      "n_nonsense", "n_spliceAD", "n_intron")
m_df <- m_pSNV[m_SNVs, on = c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                              "external_gene_name")]
m_df$n_functional <- m_df$n_missense + m_df$n_nonsense
m_df$p_functional <- m_df$p_mis + m_df$p_non
m_df <- m_df[,c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                "external_gene_name", "p_functional", "n_functional")]
m_df$OE <- m_df$n_functional/m_df$p_functional

h_SNVs <- h_SNVs[,c("ensembl_transcript_id","CHROM","ensembl_gene_id","HGNC_symbol",
                    "n_missense_0.001", "n_nonsense_0.001")]
colnames(h_SNVs) <- c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                      "external_gene_name", "n_missense",
                      "n_nonsense")
h_df <- h_pSNV[h_SNVs, on = c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                              "external_gene_name")]
h_df$n_functional <- h_df$n_missense + h_df$n_nonsense
h_df$p_functional <- h_df$p_mis + h_df$p_non
h_df <- h_df[,c("ensembl_transcript_id", "chromosome_name", "ensembl_gene_id",
                "external_gene_name", "p_functional", "n_functional")]
h_df$OE <- h_df$n_functional/h_df$p_functional
hist(h_df$OE)

orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/HM_Ensembl_orthologs.csv")
orths <- orths[, c("H_external_gene_name", "M_external_gene_name", "orthology_type")]
colnames(h_df) <- paste0("H_", colnames(h_df))  
colnames(m_df) <- paste0("M_", colnames(m_df))  
df <- orths[h_df, on = "H_external_gene_name"]
df <- df[m_df, on = "M_external_gene_name"]
df <- df[complete.cases(df),]
plot(df$M_OE, df$H_OE)
cor.test(df$M_OE, df$H_OE)
