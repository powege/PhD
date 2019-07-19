rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

annotation.file <- args[1]
out.file <- args[2]
species <- args[3]

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# IMPORT
dt <- fread(annotation.file)
dt_sub <- subset(dt, dt$category == "Intron")

if (species == "human") { 

  chr <- c(1:22) 
  for(i in chr){
    sub <- subset(dt_sub, dt_sub$chromosome == i)
    all.POS <- unlist(seq2(from = sub$start, to = sub$end))
    all.CHR <- rep(i, length(all.POS))
    out <- paste0(all.CHR, "_", all.POS)
    out <- as.data.table(unique(out)) # ensure unique
    fwrite(
      out, 
      "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_intron_POS_ID_long.csv",
      append = T,
      col.names = F
    )
    print(i)
  }
  rm(out)
  
  # import
  out <- fread("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_intron_POS_ID_long.csv", header = F)
  
  # get proportion of annotation POS in each percentile
  prop <- rep(NA, 100)
  percentile <- c(1:100)
  for (i in percentile){
    p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", i, "_650_50.csv"), header = F)
    prop[i] <- sum(p$V1 %chin% out)/length(p$V1)
    print(i)
  }
  
  df_out <- data.frame(Percentile = percentile, 
                       Fraction = prop,
                       Annotation = rep(annotation, 100))
  fwrite(df_out, out.file, append = T, col.names = F)
}

if (species == "mouse") { 
  
  chr <- c(1:19) 
  for(i in chr){
    sub <- subset(dt_sub, dt_sub$chromosome == i)
    all.POS <- unlist(seq2(from = sub$start, to = sub$end))
    all.CHR <- rep(i, length(all.POS))
    out <- paste0(all.CHR, "_", all.POS)
    out <- as.data.table(unique(out)) # ensure unique
    fwrite(
      out, 
      "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/mouse_intron_POS_ID_long.csv",
      append = T,
      col.names = F
    )
    print(i)
  }
  rm(out)
  
  # import
  out <- fread("/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/mouse_intron_POS_ID_long.csv", header = F)
  
  # get proportion of annotation POS in each percentile
  prop <- rep(NA, 100)
  percentile <- c(1:100)
  for (i in percentile){
    p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", i, "_650_50.csv"), header = F)
    prop[i] <- sum(p$V1 %chin% out)/length(p$V1)
    print(i)
  }
  
  df_out <- data.frame(Percentile = percentile, 
                       Fraction = prop,
                       Annotation = rep(annotation, 100))
  fwrite(df_out, out.file, append = T, col.names = F)
}
