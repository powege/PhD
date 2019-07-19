rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

### SET PATHS
H_in_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_human_canPC_seq_QCpass.csv"
M_in_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCpass.csv"
H_out_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv"
M_out_file_path <- "~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv"

### IMPORT
m.seq <- fread(M_in_file_path)
h.seq <- fread(H_in_file_path)

# get CpG proportion for each annotation
m.seq$n_CpG <- sapply(1:nrow(m.seq), function(i) {
  str_count(m.seq$coding[i], "CG")
})
h.seq$n_CpG <- sapply(1:nrow(h.seq), function(i) {
  str_count(h.seq$coding[i], "CG")
})

m.out <- m.seq[,c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id","n_CpG")]
h.out <- h.seq[,c("chromosome_name","external_gene_name","ensembl_gene_id","ensembl_transcript_id","n_CpG")]

### OUTPUT
fwrite(m.out, M_out_file_path)
fwrite(h.out, H_out_file_path)

