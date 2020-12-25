rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
in.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_m.mus_orthologues.csv"
out.h.sap.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_canPC_seq.csv"

### IMPORT orthologues
orths <- fread(in.orth.file)

### List canPC
h_canPC <- unique(orths$H_ensembl_transcript_id)

### SET MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')

### Pull Ensembl transcritpt CDS

## Human

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
                             h_canPC),
               mart = h_ensembl94)


### OUTPUT 
fwrite(H_seq, out.h.sap.file)

