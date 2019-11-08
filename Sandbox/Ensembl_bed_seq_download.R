# Script that pulls sequences given bed  file

rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
bed_file_path <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
species <- "human"
seq_file_path <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_seq.csv"

### IMPORT 
bed <- fread(bed_file_path)



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
