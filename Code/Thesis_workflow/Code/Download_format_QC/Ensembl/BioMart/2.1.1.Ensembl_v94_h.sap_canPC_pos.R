rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
in.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_m.mus_orthologues.csv"
out.h.sap.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_canPC_pos.csv"

### IMPORT orthologues
orths <- fread(inn.orth.file)

### List canPC
h_canPC_orth <- unique(orths$H_ensembl_transcript_id)

### SET MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')


### Pull Ensembl gene IDs, transcript IDs (v94)
# https://www.ensembl.org/Help/Glossary

h_orth <- getBM(attributes=c('chromosome_name',
                             'external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_exon_id',
                             'rank',
                             'transcription_start_site',
                             'genomic_coding_start',
                             'genomic_coding_end',
                             'cds_length',
                             'strand'),
                filters = c('chromosome_name',
                            'biotype',
                            'ensembl_transcript_id'),
                values = list(as.character(c(1:22, "X")),
                              "protein_coding",
                              h_canPC_orth),
                mart = h_ensembl94)


### OUTPUT
fwrite(h_orth, H_out_file_path)


