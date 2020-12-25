rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS

out.h.file.path <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_TSS.csv"


### SET MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')


### Pull Ensembl gene IDs, transcript IDs (v94)
# https://www.ensembl.org/Help/Glossary

h_dt <- getBM(attributes=c('chromosome_name',
                             'external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'transcription_start_site',
                             'transcript_biotype',
                             'strand'),
                filters = c('chromosome_name',
                            'transcript_biotype'),
                values = list(as.character(c(1:22)),
                              "protein_coding"),
                mart = h_ensembl94)


### OUTPUT
fwrite(h_dt, out.h.file.path)

