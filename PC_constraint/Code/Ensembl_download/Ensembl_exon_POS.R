### Script that pulls the start and end coordinates for each exon in human and mice from Ensembl 

rm(list=ls())
graphics.off()

library(data.table)
library(dplyr)
library(biomaRt)

### PULL EXON COORDINATES FROM ENSEMBL

get_coords <- function(species, chr){
  ensembl <- useEnsembl(biomart="ensembl", dataset=species)
  gene.pos <- getBM(attributes=c('chromosome_name',
                                 'external_gene_name',
                                 'ensembl_gene_id',
                                 'ensembl_transcript_id',
                                 'ensembl_exon_id',
                                 'rank',
                                 'exon_chrom_start',
                                 'exon_chrom_end'),
                      filters = c('chromosome_name',
                                  'biotype'),
                      values = list(as.character(chr),
                                    "protein_coding"),
                      mart = ensembl)
  return(gene.pos)
}

M_exon.pos <- get_coords(species = "mmusculus_gene_ensembl", chr = c(1:19, "X"))
H_exon.pos <- get_coords(species = "hsapiens_gene_ensembl", chr = c(1:22, "X"))

write.csv(M_exon.pos, "../../Data/Ensembl/M_exon_POS.csv", row.names = F)
write.csv(H_exon.pos, "../../Data/Ensembl/H_exon_POS.csv", row.names = F)





