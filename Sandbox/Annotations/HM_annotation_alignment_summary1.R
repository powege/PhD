rm(list=ls())
graphics.off()

library(data.table)
library(plyr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)<1) {
  stop("More that one argument must be supplied", call.=FALSE)
} 

### set args variables
ann_align_prefix  <- args[1]
H_ann_prefix <- args[2]
out_file <- args[3]

ann_align_prefix <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_short_Hchr"
H_ann_prefix <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr"
# out_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/test.csv"

# ann_align_prefix  <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_short_Hchr"
# H_ann_prefix <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr"
# out_file <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_summary.csv"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

dt_out <- data.table()
h_chr <- c(1:22, "X")
for (chr in 1:length(h_chr)){

### IMPORT
dt_align <- fread(paste0(ann_align_prefix, h_chr[chr], ".csv"))
dt_Hann <- fread(paste0(H_ann_prefix, h_chr[chr], ".csv"))

### FORMAT

# colnames 
colnames(dt_Hann) <- c("category", "chromosome", "start", "end")
colnames(dt_align) <- c("H_chromosome", "H_start", "H_end", "H_category", "M_chromosome", "M_start", "M_end", "M_category")
# group misc 
dt_Hann$category[dt_Hann$category == "TF binding" | dt_Hann$category == "Open chromatin"] <- "Miscellaneous"
# number code annotations (plyr)
dt_Hann$category <- mapvalues(dt_Hann$category, from=c("Exon - CDS",
                                                   "Exon - UTR",
                                                   "Exon - other",
                                                   "Promoter",
                                                   "Enhancer - proximal",
                                                   "Enhancer - distal",
                                                   "CTCF binding",
                                                   "Miscellaneous",
                                                   "Intron - proximal",
                                                   "Intron - distal",
                                                   "Unannotated"), 
                            to=c("A","B","C","D","E","F","G","H","I","J","K"))
# sequence length
dt_Hann$len <- abs(dt_Hann$end-dt_Hann$start)+1
dt_align$len <- abs(dt_align$H_start-dt_align$H_end)+1

annotation <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")
H_chromosome <- rep(h_chr[chr], length(annotation)*length(annotation))
H_annotation <- list()
M_annotation <- list()
H_ann_total <- list()
H_ann_align_total <- list()
HM_ann_align_total <- list()

for (i in 1:length(annotation)){
  
  sub <- subset(dt_Hann, dt_Hann$category == annotation[i])
  H_ann_total[[i]] <- rep(length(unique(unlist(seq2(from = sub$start, to = sub$end)))), length(annotation))
  
  sub1 <- subset(dt_align, dt_align$H_category %like% annotation[i])
  H_ann_align_total[[i]] <- rep(length(unique(unlist(seq2(from = sub1$H_start, to = sub1$H_end)))), length(annotation))
  
  tmp_vec <- rep(NA, length(annotation))
  
  for (j in 1:length(annotation)){
    
    sub2 <- subset(sub1, sub1$M_category %like% annotation[j])
    tmp_vec[j] <- length(unique(unlist(seq2(from = sub2$H_start, to = sub2$H_end))))
    
    # print(c(i,j))
  }
  
  HM_ann_align_total[[i]] <- tmp_vec
  H_annotation[[i]] <- rep(annotation[i], length(annotation))
  M_annotation[[i]] <- annotation
}  
  

dt_tmp <- data.table(H_chromosome = H_chromosome,
                     H_annotation = unlist(H_annotation),
                     H_ann_total = unlist(H_ann_total),
                     H_ann_align_total = unlist(H_ann_align_total),
                     M_annotation = unlist(M_annotation),
                     HM_ann_align_total = unlist(HM_ann_align_total)
                     )

dt_tmp$H_annotation <- mapvalues(dt_tmp$H_annotation, from=c("A","B","C","D","E","F","G","H","I","J","K"),
                              to=c("Exon - CDS",
                                   "Exon - UTR",
                                   "Exon - other",
                                   "Promoter",
                                   "Enhancer - proximal",
                                   "Enhancer - distal",
                                   "CTCF binding",
                                   "Miscellaneous",
                                   "Intron - proximal",
                                   "Intron - distal",
                                   "Unannotated"))
dt_tmp$M_annotation <- mapvalues(dt_tmp$M_annotation, from=c("A","B","C","D","E","F","G","H","I","J","K"),
                                 to=c("Exon - CDS",
                                      "Exon - UTR",
                                      "Exon - other",
                                      "Promoter",
                                      "Enhancer - proximal",
                                      "Enhancer - distal",
                                      "CTCF binding",
                                      "Miscellaneous",
                                      "Intron - proximal",
                                      "Intron - distal",
                                      "Unannotated"))


dt_out <-  rbind(dt_out, dt_tmp)

print(h_chr[chr])
}

dt_out$ind <- paste0(dt_out$H_annotation, dt_out$M_annotation)

# sub <- subset(dt_out, dt_out$ind == "PromoterExon - CDS")
format_fun <- function(sub){
  data.table(H_annotation = sub$H_annotation[1],
             H_ann_total = sum(sub$H_ann_total),
             H_ann_align_total = sum(sub$H_ann_align_total),
             M_annotation = sub$M_annotation[1],
             HM_ann_align_total = sum(sub$HM_ann_align_total)
  )
}

dt_plot <- ddply(dt_out, "ind", format_fun)
dt_plot$ind <- NULL
dt_plot$H_ann_align_frac <- dt_plot$H_ann_align_total/dt_plot$H_ann_total
dt_plot$HM_ann_align_frac <- dt_plot$HM_ann_align_total/dt_plot$H_ann_align_total 

### EXPORT
fwrite(dt_plot, out_file)


#####
