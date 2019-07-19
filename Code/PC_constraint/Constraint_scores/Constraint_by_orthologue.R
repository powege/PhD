# SCRIPT that formats orthologue output

rm(list=ls())
graphics.off()

library(data.table)

### SET VARS 

H_con_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_MAF001_v2.csv"
M_con_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv"
orth_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv"
out_file <- "~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv"

### IMPORT
H_out <- fread(H_con_file)
M_out <- fread(M_con_file)
orths <- fread(orth_file)

### FORMAT
H_out <- H_out[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                  "fun_Z", "fun_Z_95CI_upper", "fun_Z_95CI_lower")]
M_out <- M_out[,c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                  "fun_Z", "fun_Z_95CI_upper", "fun_Z_95CI_lower")]
colnames(M_out) <- paste0("M_", colnames(M_out))
colnames(H_out) <- paste0("H_", colnames(H_out))
orth_out <- orths[M_out, on = c("M_ensembl_transcript_id", "M_ensembl_gene_id", "M_external_gene_name")]
orth_out <- orth_out[H_out, on = c("H_ensembl_transcript_id", "H_ensembl_gene_id", "H_external_gene_name")]
orth_out <- orth_out[,c("H_external_gene_name", "M_external_gene_name",
                        "H_ensembl_gene_id", "M_ensembl_gene_id",
                        "H_ensembl_transcript_id", "M_ensembl_transcript_id",
                        "M_fun_Z", "H_fun_Z", 
                        "M_fun_Z_95CI_lower", "H_fun_Z_95CI_lower",
                        "M_fun_Z_95CI_upper", "H_fun_Z_95CI_upper",
                        "orthology_type")]
# orth_out <- orth_out[!is.na(orth_out$H_external_gene_name),]
summary(orth_out)
orth_out <- orth_out[complete.cases(orth_out),]

### OUTPUT
fwrite(orth_out, out_file)


######

# cor.test(orth_out$M_fun_Z, orth_out$H_fun_Z)
# plot(orth_out$M_fun_Z, orth_out$H_fun_Z)
