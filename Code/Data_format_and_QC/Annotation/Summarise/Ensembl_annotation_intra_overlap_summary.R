### SCRIPT that calculates the fraction of within annotaton overlap. 

rm(list = ls())
graphics.off()

library(data.table)

### SET ARGS
in_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
out_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_annotation_intraoverlap_summary.csv"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### IMPORT
ann_all <- fread(in_file)

### FORMAT
colnames(ann_all) <- c("category", "chromosome", "start", "end")
# table(ann_all$chromosome)
annotation <- unique(ann_all$category)

### FOR LOOP

# output lists
n_pos_list <- list()
n_unique_list <- list()

for (chr in 1:22){
  
  ann <- subset(ann_all, ann_all$chromosome == chr)
  
  all_pos <- rep(NA, (length(annotation)-2))
  unique_pos <- rep(NA, (length(annotation)-2))
  
  for(i in 1:(length(annotation)-2)){ # exclude Intron - distal and Unannotated
    sub <- subset(ann, ann$category == annotation[i])
    tmp_vec <- unlist(seq2(from = sub$start, to = sub$end))
    all_pos[i] <- length(tmp_vec)
    unique_pos[i] <- length(unique(tmp_vec))
    print(annotation[i])
  }
  
  n_pos_list[[chr]] <- data.table(n_pos = all_pos)
  n_unique_list[[chr]] <- data.table(n_unique = unique_pos)

  print(chr)
}

out_dt <- data.table(category = annotation[1:(length(annotation)-2)],
                              n_pos = apply(do.call("cbind", n_pos_list), 1, sum),
                              n_unique = apply(do.call("cbind", n_unique_list), 1, sum))
out_dt$fraction <- out_dt$n_unique/out_dt$n_pos

### EXPORT
fwrite(out_dt, out_file)
