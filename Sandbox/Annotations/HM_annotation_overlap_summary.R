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
ann_prefix <- args[1]
out_file <- args[2]
species <- args[3]

# ann_prefix <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr"
# # out_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/test.csv"
# species <- "human"


### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


if (species == "human"){ chr <- c(1:22)}
if (species == "mouse"){ chr <- c(1:19)}
dt_out_list <- list()
for (k in 1:length(chr)){

### IMPORT
dt <- fread(paste0(ann_prefix, chr[k], ".csv"))

### FORMAT

colnames(dt) <- c("category", "chromosome", "start", "end")

# number code annotations (plyr)
dt$category <- mapvalues(dt$category, from=c("Exon - CDS",
                                                       "Exon - UTR",
                                                       "Exon - other",
                                                       "Promoter",
                                                       "Enhancer - proximal",
                                                       "Enhancer - distil",
                                                       "CTCF binding",
                                                       "Open chromatin",
                                                       "TF binding",
                                                       "Intron - proximal",
                                                       "Intron - distil",
                                                       "Unannotated"), 
                              to=c("A","B","C","D","E","F","G","H","I","J", "K", "L"))
dt$length <- abs(dt$end-dt$start)+1


annotation <- c("A","B","C","D","E","F","G","H","I","J", "K", "L")
chromosome <- rep(chr[k], (length(annotation)*length(annotation)))
annotation_primary <- rep(annotation, each = length(annotation))
annotation_secondary <- rep(annotation, times = length(annotation))
primary_total_bp <- list()
secondary_total_bp <- list()

for (i in 1:length(annotation)){
  primary_sub <- subset(dt, dt$category == annotation[i])
  primary_sub_bp <- unique(unlist(seq2(from = primary_sub$start, to = primary_sub$end)))
  primary_total_bp[[i]] <- rep(length(primary_sub_bp), length(annotation))
  
  overlap_vec <- rep(NA, length(annotation))
  for (j in 1:length(annotation)){
    secondary_sub <- subset(dt, dt$category == annotation[j])
    secondary_sub_bp <- unique(unlist(seq2(from = secondary_sub$start, to = secondary_sub$end)))
    overlap_bp <- primary_sub_bp[which(primary_sub_bp %in% secondary_sub_bp)]
    overlap_vec[j] <- length(overlap_bp)
    print(j)
  }
  secondary_total_bp[[i]] <- overlap_vec
  print(i)
}
  
dt_out_list[[k]] <- data.table(chromosome = chromosome,
                     annotation_primary = annotation_primary,
                     annotation_secondary = annotation_secondary,
                     primary_total_bp = unlist(primary_total_bp),
                     secondary_total_bp = unlist(secondary_total_bp))
print(k)
}

dt_out <- do.call("rbind", dt_out_list)
dt_out$annotation_primary <- mapvalues(dt_out$annotation_primary, from=c("A","B","C","D","E","F","G","H","I","J", "K", "L"),
                                 to=c("Exon - CDS",
                                      "Exon - UTR",
                                      "Exon - other",
                                      "Promoter",
                                      "Enhancer - proximal",
                                      "Enhancer - distil",
                                      "CTCF binding",
                                      "Open chromatin",
                                      "TF binding",
                                      "Intron - proximal",
                                      "Intron - distil",
                                      "Unannotated"))
dt_out$annotation_secondary <- mapvalues(dt_out$annotation_secondary, from=c("A","B","C","D","E","F","G","H","I","J", "K", "L"),
                                       to=c("Exon - CDS",
                                            "Exon - UTR",
                                            "Exon - other",
                                            "Promoter",
                                            "Enhancer - proximal",
                                            "Enhancer - distil",
                                            "CTCF binding",
                                            "Open chromatin",
                                            "TF binding",
                                            "Intron - proximal",
                                            "Intron - distil",
                                            "Unannotated"))

### EXPORT
fwrite(dt_out, out_file)


#####

