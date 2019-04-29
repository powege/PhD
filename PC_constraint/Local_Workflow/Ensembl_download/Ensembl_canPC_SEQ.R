# Script that pulls and QCs Ensembl transctipts

rm(list = ls())
graphics.off()

# library(pegas)
# library(vcfR)
# library(seqinr)
# require(plyr)
library(biomaRt)
library(data.table)
# library(reshape2)


### SET RELATIVE PATHS and FILES
M_tran.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/MGP_canPC_SNV_counts.txt"
M_seq.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/M_canPC_SEQ_QCed.csv"
M_rem.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/M_canPC_SEQ_removed.csv"

H_tran.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_SNV_counts.txt"
H_seq.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_QCed.csv"
H_rem.file <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_removed.csv"


### IMPORT canonical transcript IDs
M_canT <- fread(M_tran.file)
H_canT <- fread(H_tran.file)

### Pull Ensembl transcritpt CDS

M_ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

# x <- listFilters(ensembl)
# attributePages(ensembl)
# x <- listAttributes(ensembl, page = "sequences")
M_seq <- getBM(attributes=c('chromosome_name',
                            'external_gene_name',
                            'ensembl_gene_id',
                            'ensembl_transcript_id',
                            'cds_length',
                            'coding'),
             filters = c('chromosome_name',
                         'biotype',
                         'ensembl_transcript_id'),
             values = list(as.character(c(1:19, "X")),
                           "protein_coding",
                           M_canT$ensembl_transcript_id),
             mart = M_ensembl)

H_ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
H_seq <- getBM(attributes=c('chromosome_name',
                            'external_gene_name',
                            'ensembl_gene_id',
                            'ensembl_transcript_id',
                            'cds_length',
                            'coding'),
               filters = c('chromosome_name',
                           'biotype',
                           'ensembl_transcript_id'),
               values = list(as.character(c(1:22, "X")),
                             "protein_coding",
                             H_canT$ensembl_transcript_id),
               mart = H_ensembl)


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


### OUTPUT 
fwrite(M_removed, M_rem.file)
fwrite(M_seq, M_seq.file)
fwrite(H_removed, H_rem.file)
fwrite(H_seq, H_seq.file)

