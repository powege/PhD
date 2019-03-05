### Script that pulls the start and end coordinates for each gene in human and mice from Ensembl 

rm(list=ls())
graphics.off()

library(data.table)
library(biomaRt)

### PULL GENE COORDINATES FROM ENSEMBL

get_coords <- function(species, chr){
  ensembl <- useEnsembl(biomart="ensembl", dataset=species)
  gene.pos <- getBM(attributes=c('chromosome_name',
                                 'external_gene_name',
                                 'ensembl_gene_id',
                                 'start_position',
                                 'end_position'),
                    filters = c('chromosome_name',
                                'biotype'),
                    values = list(as.character(chr),
                                  "protein_coding"),
                    mart = ensembl)
  return(gene.pos)
}

M_gene.pos <- get_coords(species = "mmusculus_gene_ensembl", chr = c(1:19, "X"))
H_gene.pos <- get_coords(species = "hsapiens_gene_ensembl", chr = c(1:22, "X"))

fwrite(M_gene.pos, "../../Data/Ensembl/M_gene_POS.csv")
fwrite(H_gene.pos, "../../Data/Ensembl/H_gene_POS.csv")

