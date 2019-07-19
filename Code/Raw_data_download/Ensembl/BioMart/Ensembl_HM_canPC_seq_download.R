# Script that pulls and QCs Ensembl transctipts

rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
orth_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv"
H_out_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_seq.csv"
M_out_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_seq.csv"

### IMPORT orthologues
orths <- fread(orth_file_path)

### List canPC
h_canPC <- unique(orths$H_ensembl_transcript_id)
m_canPC <- unique(orths$M_ensembl_transcript_id)

### Pull Ensembl transcritpt CDS

## Human

# listEnsemblArchives()
# listMarts(host = 'Oct2018.archive.ensembl.org')
h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
# filters <- listFilters(h_ensembl94)
# attributePages(h_ensembl94)
# sequences <- listAttributes(h_ensembl94, page = "sequences")
H_seq <- getBM(attributes=c('chromosome_name',
                            'external_gene_name',
                            'ensembl_gene_id',
                            'ensembl_transcript_id',
                            'cds_length',
                            'coding'),
               filters = c('chromosome_name',
                           'biotype',
                           'ensembl_transcript_id'),
               values = list(as.character(c(1:22)),
                             "protein_coding",
                             h_canPC),
               mart = h_ensembl94)

## Mouse

m_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='mmusculus_gene_ensembl')
M_seq <- getBM(attributes=c('chromosome_name',
                            'external_gene_name',
                            'ensembl_gene_id',
                            'ensembl_transcript_id',
                            'cds_length',
                            'coding'),
               filters = c('chromosome_name',
                           'biotype',
                           'ensembl_transcript_id'),
               values = list(as.character(c(1:19)),
                             "protein_coding",
                             m_canPC),
               mart = m_ensembl94)



### OUTPUT 
fwrite(M_seq, M_out_file_path)
fwrite(H_seq, H_out_file_path)


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


#######

# # Get list of canonical transcripts
# m_vcf <- list()
# for (i in 1:19){
# m_vcf[[i]] <- fread(paste0("~/Dropbox/PhD/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allMUSMUS_snps_QCed_VEP_v94_canPC_chr", i, ".vcf"))
# }
# m_vcf <- do.call("rbind", m_vcf)
# m_canPC <- unique(m_vcf$V9)
# 
# h_vcf <- list()
# for (i in 1:22){
#   h_vcf[[i]] <- fread(paste0("~/Dropbox/PhD/Data/1KGP/Variants/vcf_QCed_VEP/1KGP_phase3_snps_QCed_VEP_v94_canPC_chr", i, ".vcf"))
# }
# h_vcf <- do.call("rbind", h_vcf)
# h_canPC <- unique(h_vcf$V10)
