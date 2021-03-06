rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

#### SET PATHS
H_out_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_TSS.csv"
M_out_file_path <- "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_TSS.csv"

### Pull Ensembl gene IDs, transcript IDs (v94)
# https://www.ensembl.org/Help/Glossary

## Human

# listEnsemblArchives()
# listMarts(host = 'Oct2018.archive.ensembl.org')
h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='hsapiens_gene_ensembl')
# filters <- listFilters(h_ensembl94)
# attributePages(h_ensembl94)
# attributes <- listAttributes(h_ensembl94, page = "structure")
# attributes2 <- listAttributes(h_ensembl94, page = "feature_page")
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


## Mouse

m_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
                       biomart='ENSEMBL_MART_ENSEMBL', 
                       dataset='mmusculus_gene_ensembl')
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
fwrite(h_dt, H_out_file_path)
fwrite(m_dt, M_out_file_path)



#################

# ## Human
# 
# # listEnsemblArchives()
# # listMarts(host = 'Oct2018.archive.ensembl.org')
# h_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
#                        biomart='ENSEMBL_MART_ENSEMBL', 
#                        dataset='hsapiens_gene_ensembl')
# # filters <- listFilters(h_ensembl94)
# # attributePages(h_ensembl94)
# # attributes <- listAttributes(h_ensembl94, page = "structure")
# # attributes2 <- listAttributes(h_ensembl94, page = "feature_page")
# h.orth1 <- getBM(attributes=c('chromosome_name',
#                               'external_gene_name',
#                              'ensembl_gene_id',
#                              'ensembl_transcript_id',
#                              'ccds'),
#                 filters = c('chromosome_name',
#                             'biotype'),
#                 values = list(as.character(c(1:22)),
#                               "protein_coding"),
#                 mart = h_ensembl94)
# h.orth2 <- getBM(attributes=c('external_gene_name',
#                               'ensembl_gene_id',
#                               'ensembl_transcript_id',
#                               'ensembl_exon_id',
#                               'rank',
#                               'genomic_coding_start',
#                               'genomic_coding_end',
#                               'cds_length'),
#                  filters = c('chromosome_name',
#                              'biotype'),
#                  values = list(as.character(c(1:22)),
#                                "protein_coding"),
#                  mart = h_ensembl94)
# 
# ## Mouse
# 
# m_ensembl94 <- useMart(host='Oct2018.archive.ensembl.org', 
#                        biomart='ENSEMBL_MART_ENSEMBL', 
#                        dataset='mmusculus_gene_ensembl')
# m.orth <- getBM(attributes=c('external_gene_name',
#                              'ensembl_gene_id',
#                              'ensembl_transcript_id',
#                              'ensembl_peptide_id',
#                              'hsapiens_homolog_ensembl_gene',
#                              'hsapiens_homolog_associated_gene_name',
#                              'hsapiens_homolog_canonical_transcript_protein',
#                              'hsapiens_homolog_ensembl_peptide',
#                              'hsapiens_homolog_orthology_type',
#                              'hsapiens_homolog_perc_id',
#                              'hsapiens_homolog_perc_id_r1',
#                              'hsapiens_homolog_goc_score',
#                              'hsapiens_homolog_wga_coverage',
#                              'hsapiens_homolog_dn',
#                              'hsapiens_homolog_ds',
#                              'hsapiens_homolog_orthology_confidence'),
#                 filters = c('biotype'),
#                 values = list("protein_coding"),
#                 mart = m_ensembl94)
# 
# ### FORMAT
# 
# # Subset all genes with CCDS
# h.orth1 <- subset(h.orth1, h.orth1$ccds != "")
# 
# # Merge datasets
# h_out <- merge(h.orth1, h.orth2)
# 
# # Subset canonnical transcripts (ie longest cds for each gene)
# h_out <- merge(aggregate(cds_length ~ external_gene_name, h_out, max), h_out, all = F)
# 
# 
# ### OUTPUT
# # write.csv(out, "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv", row.names = F)


