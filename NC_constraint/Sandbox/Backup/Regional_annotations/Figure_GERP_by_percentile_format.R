rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

gerp.file <- args[1]
out.file <- args[2]
species <- args[3]
# gerp.file <- "~/Dropbox/PhD/Data/Ensembl/GERP/gerp_constrained_elements.homo_sapiens.bed"
# out.file <- "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_GERP_by_percentile_650_50.csv"
# species <- "human"


### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

### IMPORT
gerp <- fread(gerp.file, skip = 1)

### FORMAT
colnames(gerp) <- c("chromosome", "start", "end") 

# get vector of all GERP POS_ID 
out <- list()
if (species == "human") { chr <- c(1:22) }
if (species == "mouse") { chr <- c(1:19) }
for(i in chr){
  sub <- subset(gerp, gerp$chromosome == chr[i])
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(chr[i], length(all.POS))
  out[[i]] <- paste0(all.CHR, "_", all.POS)
  print(i)
}
out <- unique(unlist(out))

# get proportion of annotation POS in each percentile
prop <- rep(NA, 100)
percentile <- c(1:100)
for (i in percentile){
  # p <- fread("~/Dropbox/PhD/Data/NC_constraint/Percentile_POS_ID/human_POS_ID_percentile_1_650_50.csv", header = F)
  p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", i, "_650_50.csv"), header = F)
  prop[i] <- table(p$V1 %in% out)[2]/length(p$V1)
    print(i)
}

df_out <- data.frame(Percentile = percentile, 
                     Fraction = prop,
                     Annotation = rep("GERP", 100))

### EXPORT
fwrite(df_out, out.file, col.names = F)