rm(list = ls())
graphics.off()

library(data.table)
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

ann.rank <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv")
ann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv")

# intron
int.rank <- ann.rank[ann.rank$category == "Intron" & ann.rank$chromosome == 1,]
int <- ann[ann$category == "Intron" & ann$chromosome == 1,]

pos.int.rank <- unlist(seq2(from = int.rank$start, to = int.rank$end))
pos.int.rank.unq <- unique(pos.int.rank)

pos.int <- unlist(seq2(from = int$start, to = int$end))
pos.int.unq <- unique(pos.int)

dt_sub <- ann.rank[ann.rank$category == "Intron",]
out <- list()
chr <- c(1:22) 
for(i in chr){
  sub <- subset(dt_sub, dt_sub$chromosome == chr[i])
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(chr[i], length(all.POS))
  out[[i]] <- paste0(all.CHR, "_", all.POS)
  print(i)
}
out <- unique(unlist(out))





# unannotated
un.rank <- ann.rank[ann.rank$category == "Unannotated" & ann.rank$chromosome == 1,]
un <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_unannotated.csv")
un <- un[un$chromosome == 1,]

pos.un.rank <- unlist(seq2(from = un.rank$start, to = un.rank$end))
pos.un.rank.unq <- unique(pos.un.rank)

pos.un <- unlist(seq2(from = un$start, to = un$end))
pos.un.unq <- unique(pos.un)


