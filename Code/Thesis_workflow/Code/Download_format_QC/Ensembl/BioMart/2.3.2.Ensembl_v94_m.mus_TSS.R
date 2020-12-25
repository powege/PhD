rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS

out.m.file.path <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_TSS.csv"


### SET MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

m_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='mmusculus_gene_ensembl')


### Pull Ensembl gene IDs, transcript IDs (v94)
# https://www.ensembl.org/Help/Glossary

m_dt <- getBM(attributes=c('chromosome_name',
                           'external_gene_name',
                           'ensembl_gene_id',
                           'ensembl_transcript_id',
                           'transcription_start_site',
                           'transcript_biotype',
                           'strand'),
              filters = c('chromosome_name',
                          'transcript_biotype'),
              values = list(as.character(c(1:19)),
                            "protein_coding"),
              mart = m_ensembl94)

### OUTPUT
fwrite(m_dt, out.m.file.path)

