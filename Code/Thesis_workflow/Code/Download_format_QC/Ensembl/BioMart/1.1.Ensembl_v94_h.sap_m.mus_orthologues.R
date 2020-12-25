rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

### SET VARS 
out.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_h.sap_m.mus_orthologues.csv"

### SER MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

## Human
h_ensembl94 <- useMart(
  host='http://oct2018.archive.ensembl.org', 
  biomart='ENSEMBL_MART_ENSEMBL', 
  dataset='hsapiens_gene_ensembl'
  )

# Mouse
m_ensembl94 <- useMart(
  host='Oct2018.archive.ensembl.org', 
  biomart='ENSEMBL_MART_ENSEMBL', 
  dataset='mmusculus_gene_ensembl'
  )

### PULL DATA 

# filters <- listFilters(h_ensembl94)
# attributePages(h_ensembl94)
# attributes <- listAttributes(h_ensembl94, page = "homologs")
# attributes <- listAttributes(h_ensembl94, page = "feature_page")

# Human
h.orth <- getBM(attributes=c('external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_peptide_id',
                             'mmusculus_homolog_ensembl_gene',
                             'mmusculus_homolog_associated_gene_name',
                             'mmusculus_homolog_canonical_transcript_protein',
                             'mmusculus_homolog_ensembl_peptide',
                             'mmusculus_homolog_orthology_type',
                             'mmusculus_homolog_perc_id',
                             'mmusculus_homolog_perc_id_r1',
                             'mmusculus_homolog_goc_score',
                             'mmusculus_homolog_wga_coverage',
                             'mmusculus_homolog_dn',
                             'mmusculus_homolog_ds',
                             'mmusculus_homolog_orthology_confidence'),
                filters = c('biotype'),
                values = list("protein_coding"),
                mart = h_ensembl94)

## Mouse
m.orth <- getBM(attributes=c('external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_peptide_id',
                             'hsapiens_homolog_ensembl_gene',
                             'hsapiens_homolog_associated_gene_name',
                             'hsapiens_homolog_canonical_transcript_protein',
                             'hsapiens_homolog_ensembl_peptide',
                             'hsapiens_homolog_orthology_type',
                             'hsapiens_homolog_perc_id',
                             'hsapiens_homolog_perc_id_r1',
                             'hsapiens_homolog_goc_score',
                             'hsapiens_homolog_wga_coverage',
                             'hsapiens_homolog_dn',
                             'hsapiens_homolog_ds',
                             'hsapiens_homolog_orthology_confidence'),
                filters = c('biotype'),
                values = list("protein_coding"),
                mart = m_ensembl94)


### FORMAT

# mmusculus_homolog_canonical_transcript_protein is the canonical human peptide ID
# remove rows with no canonical peptide ID
h.orth <- subset(h.orth, h.orth$mmusculus_homolog_ensembl_gene != "")
# h.orth <- h.orth[which(h.orth$ensembl_peptide_id %in%
#                        m.orth$hsapiens_homolog_canonical_transcript_protein),]

m.orth <- subset(m.orth, m.orth$hsapiens_homolog_ensembl_gene != "")
# m.orth <- m.orth[which(m.orth$ensembl_peptide_id %in% 
#                          h.orth$mmusculus_homolog_canonical_transcript_protein),]

colnames(h.orth) <- c("H_external_gene_name",                            
                      "H_ensembl_gene_id",                               
                      "H_ensembl_transcript_id",                         
                      "H_ensembl_peptide_id",                            
                      "M_ensembl_gene_id",                
                      "M_external_gene_name",        
                      "mmusculus_homolog_canonical_transcript_protein",
                      "M_ensembl_peptide_id",             
                      "orthology_type",              
                      "Maa_match_Haa",                     
                      "Haa_match_Maa",                  
                      "goc_score",                   
                      "wga_coverage",                
                      "dn",                          
                      "ds",                          
                      "orthology_confidence")  

h.orth <- h.orth[c("H_external_gene_name",                            
                   "H_ensembl_gene_id",                               
                   "H_ensembl_transcript_id",                         
                   "H_ensembl_peptide_id",                            
                   "M_ensembl_gene_id",                
                   "M_external_gene_name",        
                   "M_ensembl_peptide_id",             
                   "orthology_type",              
                   "Maa_match_Haa",                     
                   "Haa_match_Maa",                  
                   "goc_score",                   
                   "wga_coverage",                
                   "dn",                          
                   "ds",                          
                   "orthology_confidence")]

colnames(m.orth) <- c("M_external_gene_name",                            
                      "M_ensembl_gene_id",                               
                      "M_ensembl_transcript_id",                         
                      "M_ensembl_peptide_id",                            
                      "H_ensembl_gene_id",                
                      "H_external_gene_name",        
                      "hsapiens_homolog_canonical_transcript_protein",
                      "H_ensembl_peptide_id",             
                      "orthology_type",              
                      "Haa_match_Maa",                     
                      "Maa_match_Haa",                  
                      "goc_score",                   
                      "wga_coverage",                
                      "dn",                          
                      "ds",                          
                      "orthology_confidence")  

m.orth <- m.orth[c("M_external_gene_name",                            
                   "M_ensembl_gene_id",                               
                   "M_ensembl_transcript_id",                         
                   "M_ensembl_peptide_id",                            
                   "H_ensembl_gene_id",                
                   "H_external_gene_name",        
                   "H_ensembl_peptide_id",             
                   "orthology_type",              
                   "Haa_match_Maa",                     
                   "Maa_match_Haa",                  
                   "goc_score",                   
                   "wga_coverage",                
                   "dn",                          
                   "ds",                          
                   "orthology_confidence")]

out <- merge(h.orth, m.orth, all = T)
out <- subset(out, !is.na(out$H_ensembl_transcript_id) & !is.na(out$M_ensembl_transcript_id))


### OUTPUT
fwrite(out, out.orth.file)


