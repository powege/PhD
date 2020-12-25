rm(list = ls())
graphics.off()

library(biomaRt)
library(data.table)

### SET VARS 
out.orth.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_r.nor_orthologues.csv"

### SER MARTS

# listEnsemblArchives()
# listMarts(host='http://oct2018.archive.ensembl.org')

# Rat
r_ensembl94 <- useMart(
  host='http://oct2018.archive.ensembl.org', 
  biomart='ENSEMBL_MART_ENSEMBL', 
  dataset='rnorvegicus_gene_ensembl'
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

# Rat
r.orth <- getBM(attributes=c('external_gene_name',
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
                mart = r_ensembl94)

# Mouse
m.orth <- getBM(attributes=c('external_gene_name',
                             'ensembl_gene_id',
                             'ensembl_transcript_id',
                             'ensembl_peptide_id',
                             'rnorvegicus_homolog_ensembl_gene',
                             'rnorvegicus_homolog_associated_gene_name',
                             'rnorvegicus_homolog_canonical_transcript_protein',
                             'rnorvegicus_homolog_ensembl_peptide',
                             'rnorvegicus_homolog_orthology_type',
                             'rnorvegicus_homolog_perc_id',
                             'rnorvegicus_homolog_perc_id_r1',
                             'rnorvegicus_homolog_goc_score',
                             'rnorvegicus_homolog_wga_coverage',
                             'rnorvegicus_homolog_dn',
                             'rnorvegicus_homolog_ds',
                             'rnorvegicus_homolog_orthology_confidence'),
                filters = c('biotype'),
                values = list("protein_coding"),
                mart = m_ensembl94)

### FORMAT

# mmusculus_homolog_canonical_transcript_protein is the canonical rat peptide ID
# remove rows with no canonical peptide ID
r.orth <- subset(r.orth, r.orth$mmusculus_homolog_ensembl_gene != "")
# r.orth <- r.orth[which(r.orth$ensembl_peptide_id %in%
#                        m.orth$rnorvegicus_homolog_canonical_transcript_protein),]

m.orth <- subset(m.orth, m.orth$rnorvegicus_homolog_ensembl_gene != "")
# m.orth <- m.orth[which(m.orth$ensembl_peptide_id %in% 
#                          r.orth$mmusculus_homolog_canonical_transcript_protein),]

colnames(r.orth) <- c("R_external_gene_name",                            
                      "R_ensembl_gene_id",                               
                      "R_ensembl_transcript_id",                         
                      "R_ensembl_peptide_id",                            
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

r.orth <- r.orth[c("R_external_gene_name",                            
                   "R_ensembl_gene_id",                               
                   "R_ensembl_transcript_id",                         
                   "R_ensembl_peptide_id",                            
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
                      "R_ensembl_gene_id",                
                      "R_external_gene_name",        
                      "rnorvegicus_homolog_canonical_transcript_protein",
                      "R_ensembl_peptide_id",             
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
                   "R_ensembl_gene_id",                
                   "R_external_gene_name",        
                   "R_ensembl_peptide_id",             
                   "orthology_type",              
                   "Haa_match_Maa",                     
                   "Maa_match_Haa",                  
                   "goc_score",                   
                   "wga_coverage",                
                   "dn",                          
                   "ds",                          
                   "orthology_confidence")]

out <- merge(r.orth, m.orth, all = T)
out <- subset(out, !is.na(out$R_ensembl_transcript_id) & !is.na(out$M_ensembl_transcript_id))


### OUTPUT
fwrite(out, out.orth.file)


