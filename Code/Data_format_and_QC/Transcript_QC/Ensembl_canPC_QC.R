rm(list=ls())
graphics.off()

library(data.table)


#### SET PATHS
H_in_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_seq.csv"
M_in_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_seq.csv"
M_rem_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCfail.csv"
H_rem_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_human_canPC_seq_QCfail.csv"
M_seq_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCpass.csv"
H_seq_file <- "~/Dropbox/PhD/Data/Ensembl/BioMart/QCed/Ensembl_v94_human_canPC_seq_QCpass.csv"


### IMPORT 
H_seq <-  fread(H_in_file_path)
M_seq <-  fread(M_in_file_path)


### QC 

M_removed <- data.frame()
H_removed <- data.frame()

# remove seq not divisible by 3
ncod <- M_seq$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
ncod <- H_seq$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

# remove seq with nchar < 9
rm.id <- which(M_seq$cds_length<9)
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
rm.id <- which(H_seq$cds_length<9)
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

# remove seq with nchar > 15000
rm.id <- which(M_seq$cds_length>15000)
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
rm.id <- which(H_seq$cds_length>15000)
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

# remove seq with characters other than ATCG 
rm.id <- grep('[^ATGC]', M_seq$coding)
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
rm.id <- grep('[^ATGC]', H_seq$coding)
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

# remove seq that do not start with start codon
c1 <- substring(M_seq$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
c1 <- substring(H_seq$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

# remove seq that do not end with stop codon
cn <- substring(M_seq$coding, M_seq$cds_length-2, M_seq$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  M_seq <- M_seq[-rm.id,]
  M_removed <- rbind(M_removed, M_seq[rm.id,])
}
cn <- substring(H_seq$coding, H_seq$cds_length-2, H_seq$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  H_seq <- H_seq[-rm.id,]
  H_removed <- rbind(H_removed, H_seq[rm.id,])
}

### REORDER
M_seq <- M_seq[order(M_seq$chromosome_name)]
H_seq <- H_seq[order(H_seq$chromosome_name)]
M_removed <- M_removed[order(M_removed$chromosome_name)]
H_removed <- H_removed[order(H_removed$chromosome_name)]

M_seq <- M_seq[complete.cases(M_seq),]
H_seq <- H_seq[complete.cases(H_seq),]
M_removed <- M_removed[complete.cases(M_removed),]
H_removed <- H_removed[complete.cases(H_removed),]


### OUTPUT 
fwrite(M_removed, M_rem_file)
fwrite(M_seq, M_seq_file)
fwrite(H_removed, H_rem_file)
fwrite(H_seq, H_seq_file)
