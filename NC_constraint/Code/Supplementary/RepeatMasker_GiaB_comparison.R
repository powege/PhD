# Script that compares overlap between repeatmasked regions in Ensembl and UCSC

rm(list = ls())
graphics.off()

library(data.table)

ensembl <- fread("~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr19.txt")
ucsc <- fread("~/Dropbox/PhD/Data/UCSC/M_rmsk.txt")

ensembl_POS <- ensembl$POS[ensembl$REP_MASK == 1]

ucsc <- subset(ucsc, ucsc$V6 == "chr19")
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
ucsc_POS <- unlist(seq2(from = ucsc$V7, to = ucsc$V8))
ucsc_POS <- unique(ucsc_POS)

table(ensembl_POS %in% ucsc_POS)
table(ucsc_POS %in% ensembl_POS)
table(ucsc$V12)

giab <- fread("~/Dropbox/PhD/Data/GiaB/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed")
H_ucsc <- fread("~/Dropbox/PhD/Data/UCSC/H_rmsk.txt")

H_ucsc <- subset(H_ucsc, H_ucsc$V6 == "chr19")
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
H_ucsc_POS <- unlist(seq2(from = H_ucsc$V7, to = H_ucsc$V8))
H_ucsc_POS <- unique(H_ucsc_POS)

giab <- subset(giab, giab$V1 == "chr19")
giab_POS <- unlist(seq2(from = giab$V2, to = giab$V3))
table(giab_POS %in% H_ucsc_POS)
