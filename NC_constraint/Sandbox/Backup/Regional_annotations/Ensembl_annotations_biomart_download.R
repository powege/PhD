rm(list=ls())
graphics.off()

### Script that pulls from Ensembl BiomaRt:
# Exon CDS coordinates
# Exon ncRNA coordinates
# UTR coordinates

library(biomaRt)
library(data.table)

### PULL EXON COORDINATES FROM ENSEMBL

exon_coords <- function(species, chr){
  ensembl <- useEnsembl(biomart="ensembl",   # GRCh38 is default
                        dataset=species)
  exon.pos <- getBM(attributes=c('chromosome_name',
                                 'external_gene_name',
                                 'ensembl_gene_id',
                                 'ensembl_transcript_id',
                                 'ensembl_exon_id',
                                 'rank',
                                 'exon_chrom_start',
                                 'exon_chrom_end',
                                 'gene_biotype'),
                    filters = c('chromosome_name'),
                    values = list(as.character(chr)),
                    mart = ensembl)
  return(exon.pos)
}

UTR_coords <- function(species, chr){
  ensembl <- useEnsembl(biomart="ensembl", 
                        dataset=species)
  UTR.pos <- getBM(attributes=c('chromosome_name',
                                 'external_gene_name',
                                 'ensembl_gene_id',
                                 'ensembl_transcript_id',
                                 '5_utr_start',
                                 '5_utr_end',
                                 '3_utr_start',
                                 '3_utr_end'),
                    filters = c('chromosome_name',
                                'biotype'),
                    values = list(as.character(chr),
                                  "protein_coding"),
                    mart = ensembl)
  return(UTR.pos)
}

# M_exon.pos <- exon_coords(species = "mmusculus_gene_ensembl", chr = c(1:19, "X"))
H_exon.pos <- exon_coords(species = "hsapiens_gene_ensembl", chr = c(1:22))

# M_UTR.pos <- UTR_coords(species = "mmusculus_gene_ensembl", chr = c(1:19, "X"))
H_UTR.pos <- UTR_coords(species = "hsapiens_gene_ensembl", chr = c(1:22))

# M_UTR.pos <- M_UTR.pos[!is.na(M_UTR.pos$`5_utr_start`) | !is.na(M_UTR.pos$`3_utr_start`),]
H_UTR.pos <- H_UTR.pos[!is.na(H_UTR.pos$`5_utr_start`) | !is.na(H_UTR.pos$`3_utr_start`),]

# fwrite(M_exon.pos, "~/Dropbox/PhD/Data/Ensembl/BioMart/M_exon_POS.csv")
fwrite(H_exon.pos, "~/Dropbox/PhD/Data/Ensembl/BioMart/H_exon_POS.csv")

# fwrite(M_UTR.pos, "~/Dropbox/PhD/Data/Ensembl/BioMart/M_UTR_POS.csv")
fwrite(H_UTR.pos, "~/Dropbox/PhD/Data/Ensembl/BioMart/H_UTR_POS.csv")

#####

