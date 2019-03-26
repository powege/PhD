rm(list = ls())
graphics.off()

H_df_funZ_0.001 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
H_df_funZ_0.0005 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.0005.csv")
H_df_funZ_0.0001 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.0001.csv")


a <- H_df_funZ_0.001[,c("ensembl_transcript_id", "fun_Z_0.001", "fun_Z_percentile_0.001")]
b <- H_df_funZ_0.0005[,c("ensembl_transcript_id", "fun_Z_0.0005", "fun_Z_percentile_0.0005")]
c <- H_df_funZ_0.0001[,c("ensembl_transcript_id", "fun_Z_0.0001", "fun_Z_percentile_0.0001")]
df <- merge(a, b)
df <- merge(df, c)

cor.test(df$fun_Z_0.001, df$fun_Z_0.0005)
cor.test(df$fun_Z_0.001, df$fun_Z_0.0001)
cor.test(df$fun_Z_0.0005, df$fun_Z_0.0001)

cor.test(df$fun_Z_percentile_0.001, df$fun_Z_percentile_0.0005)
cor.test(df$fun_Z_percentile_0.001, df$fun_Z_percentile_0.0001)
cor.test(df$fun_Z_percentile_0.0005, df$fun_Z_percentile_0.0001)
