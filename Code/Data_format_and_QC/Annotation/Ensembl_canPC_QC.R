rm(list=ls())
graphics.off()

library(data.table)


#### SET PATHS
H_seq_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_seq.csv"
M_seq_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_seq.csv"
H_rd_file <- "~/Dropbox/PhD/Data/gnomAD/coverage/Ensembl_v94_human_canPC_gnomad_mask.csv"
M_rd_file <- "~/Dropbox/PhD/Data/MGP/coverage/Ensembl_v94_mouse_canPC_MGP_mask.csv"
M_rem_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCfail.csv"
H_rem_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_human_canPC_seq_QCfail.csv"
M_QCed_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_QCpass.csv"
H_QCed_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_human_canPC_QCpass.csv"


### IMPORT 
H_seq <-  fread(H_seq_file)
M_seq <-  fread(M_seq_file)
H_rd <- fread(H_rd_file)
M_rd <- fread(M_rd_file)

### FORMAT

# merge datasets
M_dt <- M_rd[M_seq, on = c("chromosome_name", "external_gene_name", "ensembl_gene_id",
                           "ensembl_transcript_id", "cds_length")]
H_dt <- H_rd[H_seq, on = c("chromosome_name", "external_gene_name", "ensembl_gene_id",      
                           "ensembl_transcript_id", "cds_length")]

### QC 

M_removed <- data.frame()
H_removed <- data.frame()

# # remove cds > 10% masked
# rm.id <- which(M_dt$cds_fraction_masked > 0.1)
# if (length(rm.id) != 0){
#   M_removed <- rbind(M_removed, M_dt[rm.id,])
#   M_dt <- M_dt[-rm.id,]
# }
# rm.id <- which(H_dt$cds_fraction_masked > 0.1)
# if (length(rm.id) != 0){
#   H_removed <- rbind(H_removed, H_dt[rm.id,])
#   H_dt <- H_dt[-rm.id,]
# }

# remove seq not divisible by 3
ncod <- M_dt$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  M_removed <- rbind(M_removed, M_dt[rm.id,])
  M_dt <- M_dt[-rm.id,]
}
ncod <- H_dt$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  H_removed <- rbind(H_removed, H_dt[rm.id,])
  H_dt <- H_dt[-rm.id,]
}

# remove seq with nchar < 9
rm.id <- which(M_dt$cds_length<9)
if (length(rm.id) != 0){
  M_removed <- rbind(M_removed, M_dt[rm.id,])
  M_dt <- M_dt[-rm.id,]
}
rm.id <- which(H_dt$cds_length<9)
if (length(rm.id) != 0){
  H_removed <- rbind(H_removed, H_dt[rm.id,])
  H_dt <- H_dt[-rm.id,]
}

# # remove seq with nchar > 15000
# rm.id <- which(M_dt$cds_length>15000)
# if (length(rm.id) != 0){
#   M_removed <- rbind(M_removed, M_dt[rm.id,])
#   M_dt <- M_dt[-rm.id,]
# }
# rm.id <- which(H_dt$cds_length>15000)
# if (length(rm.id) != 0){
#   H_removed <- rbind(H_removed, H_dt[rm.id,])
#   H_dt <- H_dt[-rm.id,]
# }
# 
# # remove seq with nchar < 300
# rm.id <- which(M_dt$cds_length<300)
# if (length(rm.id) != 0){
#   M_removed <- rbind(M_removed, M_dt[rm.id,])
#   M_dt <- M_dt[-rm.id,]
# }
# rm.id <- which(H_dt$cds_length<300)
# if (length(rm.id) != 0){
#   H_removed <- rbind(H_removed, H_dt[rm.id,])
#   H_dt <- H_dt[-rm.id,]
# }

# remove seq with characters other than ATCG 
rm.id <- grep('[^ATGC]', M_dt$coding)
if (length(rm.id) != 0){
  M_removed <- rbind(M_removed, M_dt[rm.id,])
  M_dt <- M_dt[-rm.id,]
}
rm.id <- grep('[^ATGC]', H_dt$coding)
if (length(rm.id) != 0){
  H_removed <- rbind(H_removed, H_dt[rm.id,])
  H_dt <- H_dt[-rm.id,]
}

# remove seq that do not start with start codon
c1 <- substring(M_dt$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  M_removed <- rbind(M_removed, M_dt[rm.id,])
  M_dt <- M_dt[-rm.id,]
}
c1 <- substring(H_dt$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  H_removed <- rbind(H_removed, H_dt[rm.id,])
  H_dt <- H_dt[-rm.id,]
}

# remove seq that do not end with stop codon
cn <- substring(M_dt$coding, M_dt$cds_length-2, M_dt$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  M_removed <- rbind(M_removed, M_dt[rm.id,])
  M_dt <- M_dt[-rm.id,]
}
cn <- substring(H_dt$coding, H_dt$cds_length-2, H_dt$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  H_removed <- rbind(H_removed, H_dt[rm.id,])
  H_dt <- H_dt[-rm.id,]
}

### REORDER
M_dt <- M_dt[order(M_dt$chromosome_name)]
H_dt <- H_dt[order(H_dt$chromosome_name)]
M_removed <- M_removed[order(M_removed$chromosome_name)]
H_removed <- H_removed[order(H_removed$chromosome_name)]

M_dt <- M_dt[complete.cases(M_dt),]
H_dt <- H_dt[complete.cases(H_dt),]
M_removed <- M_removed[complete.cases(M_removed),]
H_removed <- H_removed[complete.cases(H_removed),]


### OUTPUT 
fwrite(M_removed, M_rem_file)
fwrite(M_dt, M_QCed_file)
fwrite(H_removed, H_rem_file)
fwrite(H_dt, H_QCed_file)
