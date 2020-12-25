rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
in.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_m.mus_orthologues.csv"
out.m.mus.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_pos.csv"

### IMPORT orthologues
orths <- fread(inn.orth.file)

### List canPC
m_canPC_orth <- unique(orths$M_ensembl_transcript_id)

### SET MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

m_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='mmusculus_gene_ensembl')

### Pull Ensembl gene IDs, transcript IDs (v94)
# https://www.ensembl.org/Help/Glossary

m_orth <- getBM(attributes=c('chromosome_name',
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
                              m_canPC_orth),
                mart = m_ensembl94)

### OUTPUT
fwrite(m_orth, M_out_file_path)


