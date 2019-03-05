rm(list = ls())
graphics.off()

library(data.table)
library(pROC)
library(MASS)

### IMPORT 

funZ <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
RVIS_esp <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/RVIS_Petrovski2013_ESP.csv")
Lek <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")

### FORMAT

# colnames
RVIS_esp <- RVIS_esp[,c(1,2)]
colnames(RVIS_esp) <- c("external_gene_name", "RVIS_ESP")
Lek <- Lek[,c("gene", "mis_z", "lof_z", "pLI")]
colnames(Lek) <- c("external_gene_name", "mis_z_Lek", "lof_z_Lek", "pLI_Lek")

funZ <- funZ[,c("external_gene_name", "fun_Z_0.001")]

df <- merge(Lek, RVIS_esp)
df <- merge(df, funZ)
df <- df[!is.na(df$external_gene_name),]
df <- df[df$external_gene_name != "",]

# percentiles
percentile <- ecdf(df$RVIS_ESP)
df$RVIS_ESP_percentile <- percentile(df$RVIS_ESP)
percentile <- ecdf(df$mis_z_Lek)
df$mis_z_Lek_percentile <- percentile(df$mis_z_Lek)
percentile <- ecdf(df$pLI_Lek)
df$pLI_Lek_percentile <- percentile(df$pLI_Lek)
percentile <- ecdf(df$fun_Z_0.001)
df$fun_Z_0.001_percentile <- percentile(df$fun_Z_0.001)

cor.test(df$fun_Z_0.001_percentile, df$RVIS_ESP_percentile)
cor.test(df$fun_Z_0.001_percentile, df$mis_z_Lek_percentile)
cor.test(df$fun_Z_0.001_percentile, df$pLI_Lek)

cor.test(df$fun_Z_0.001, df$RVIS_ESP, method = "spearman")
cor.test(df$fun_Z_0.001, df$mis_z_Lek, method = "spearman")
cor.test(df$fun_Z_0.001, df$pLI_Lek, method = "spearman")

cor.test(df$RVIS_ESP, df$mis_z_Lek, method = "spearman")
cor.test(df$RVIS_ESP, df$pLI_Lek, method = "spearman")

cor.test(df$mis_z_Lek, df$pLI_Lek, method = "spearman")

