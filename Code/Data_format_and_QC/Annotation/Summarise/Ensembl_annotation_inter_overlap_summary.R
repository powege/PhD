### SCRIPT that calculates the fraction of between annotaton overlap. 

rm(list = ls())
graphics.off()

library(data.table)

### SET ARGS
# in_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
# out_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_annotation_interoverlap_summary.csv"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

overlap <- function(in_file, out_file, species){

### IMPORT
ann_all <- fread(in_file)

### FORMAT
colnames(ann_all) <- c("category", "chromosome", "start", "end")
# table(ann_all$chromosome)
annotation <- unique(ann_all$category)[1:10] # exclude Intron - distal and Unannotated

### FOR LOOP

n_ann_A_LIST <- list()
n_ann_B_LIST <- list()
n_overlap_LIST <- list()

for (chr in 1:22){ # For each chromosome
  
  ann_chr <- subset(ann_all, ann_all$chromosome == chr) # subset by chromosome
  
  n_ann_A_tmp <- list() # set output lists
  n_ann_B_tmp <- list()
  n_overlap_tmp <- list()
  
  for(i in 1:length(annotation)){ # for each annotation
    
    sub_A <- subset(ann_chr, ann_chr$category == annotation[i]) # subset by annotation
    tmp_vec_A <- unique(as.vector(unlist(seq2(from = sub_A$start, to = sub_A$end)))) # vector of unique POS
    n_ann_A_tmp[[i]] <- rep(length(unique(tmp_vec_A)), length(annotation)) # n unique POS
    
    n_ann_B_vec <- rep(NA, length(annotation)) # set output vectors
    n_overlap_vec <- rep(NA, length(annotation))
    
    for (j in 1:length(annotation)){ # for each annotation
      
      sub_B <- subset(ann_chr, ann_chr$category == annotation[j]) # subset by annotation
      tmp_vec_B <- unique(as.vector(unlist(seq2(from = sub_B$start, to = sub_B$end)))) # vector of unique POS
      n_ann_B_vec[j] <- length(tmp_vec_B) # n unique POS
      n_overlap_vec[j] <- length(intersect(x = tmp_vec_A, y = tmp_vec_B)) # n overlapping POS
      
    }
  
    n_ann_B_tmp[[i]] <- n_ann_B_vec
    n_overlap_tmp[[i]] <- n_overlap_vec
    
    print(annotation[i])
  }
  
  n_ann_A_LIST[[chr]] <- unlist(n_ann_A_tmp)
  n_ann_B_LIST[[chr]] <- unlist(n_ann_B_tmp)
  n_overlap_LIST[[chr]] <- unlist(n_overlap_tmp)
  
  print(chr)
}

out_dt <- data.table(cat_A = rep(annotation, each = length(annotation)),
                     cat_B = rep(annotation, times = length(annotation)),
                     n_cat_A = apply(do.call("cbind", n_ann_A_LIST), 1, sum),
                     n_cat_B = apply(do.call("cbind", n_ann_B_LIST), 1, sum),
                     n_overlap = apply(do.call("cbind", n_overlap_LIST), 1, sum))
                     
out_dt$fraction <- out_dt$n_overlap/out_dt$n_cat_A
out_dt$species <- species

### EXPORT
fwrite(out_dt, out_file)
}

### SCRIPT 

overlap(in_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv",
        out_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_annotation_interoverlap_summary.csv",
        species = "Human")

overlap(in_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv",
        out_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_annotation_interoverlap_summary.csv",
        species = "Mouse")

