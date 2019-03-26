# Script that pulls a dataframe of human and mouse orthologs from Ensembl 

rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

### Pull ennsembl gene IDs, transcript IDs

ensembl = useMart("ensembl")
# listDatasets(ensembl)
h.ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# x <- listFilters(h.ensembl)
# attributePages(h.ensembl)
# x <- listAttributes(h.ensembl, page = "homologs")
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
                mart = h.ensembl)

ensembl = useMart("ensembl")
# x <- listDatasets(ensembl)
m.ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
# x <- listFilters(m.ensembl)
# attributePages(m.ensembl)
# x <- listAttributes(m.ensembl, page = "homologs")
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
                mart = m.ensembl)


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

write.csv(out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/HM_Ensembl_orthologs.csv", row.names = F)

# length(unique(out$H_external_gene_name))
# length(unique(out$M_external_gene_name))
